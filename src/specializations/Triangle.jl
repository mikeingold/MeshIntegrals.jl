################################################################################
#                    Specialized Methods for Triangle
#
# Why Specialized?
#   The Triangle geometry is a surface simplex whose parametric function in
#   Meshes.jl uses barycentric coordinates on a domain {u,v} with coordinates
#   that are non-negative and bound by the surface $u + v ≤ 1$. This requires a
#   multi-step domain transformation whose derivation is detailed in the package
#   documentation.
################################################################################

function integral(
        f,
        triangle::Meshes.Triangle,
        rule::GaussLegendre;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    _guarantee_analytical(Meshes.Triangle, diff_method)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = _gausslegendre(FP, rule.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformations:
    #   xᵢ [-1,1] ↦ R [0,1]
    #   xⱼ [-1,1] ↦ φ [0,π/2]
    uR(xᵢ) = (1 // 2) * (xᵢ + 1)
    uφ(xⱼ) = (T(π) / 4) * (xⱼ + 1)

    # Integrate the Barycentric triangle by transforming it into polar coordinates
    #   with a modified radius
    #     R = r ( sinφ + cosφ )
    #   s.t. integration bounds become rectangular
    #     R ∈ [0, 1] and φ ∈ [0, π/2]
    function integrand(((wᵢ, wⱼ), (xᵢ, xⱼ)))
        R = uR(xᵢ)
        φ = uφ(xⱼ)
        a, b = sincos(φ)
        u = R * (1 - a / (a + b))
        v = R * (1 - b / (a + b))
        return wᵢ * wⱼ * f(triangle(u, v)) * R / (a + b)^2
    end

    # Calculate 2D Gauss-Legendre integral over modified-polar-Barycentric coordinates
    # Apply a linear domain-correction factor
    return FP(π / 4) * Meshes.area(triangle) .* sum(integrand, zip(wws, xxs))
end

function integral(
        f,
        triangle::Meshes.Triangle,
        rule::GaussKronrod;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    _guarantee_analytical(Meshes.Triangle, diff_method)

    # Integrate the Barycentric triangle in (u,v)-space: (0,0), (0,1), (1,0)
    #   i.e. \int_{0}^{1} \int_{0}^{1-u} f(u,v) dv du
    ∫u(u) = QuadGK.quadgk(v -> f(triangle(u, v)), zero(FP), FP(1 - u); rule.kwargs...)[1]
    ∫ = QuadGK.quadgk(∫u, zero(FP), one(FP); rule.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2 * Meshes.area(triangle) .* ∫
end

function integral(
        f,
        triangle::Meshes.Triangle,
        rule::HAdaptiveCubature;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    _guarantee_analytical(Meshes.Triangle, diff_method)

    # Integrate the Barycentric triangle by transforming it into polar coordinates
    #   with a modified radius
    #     R = r ( sinφ + cosφ )
    #   s.t. integration bounds become rectangular
    #     R ∈ [0, 1] and φ ∈ [0, π/2]
    function integrand(Rφ)
        R, φ = Rφ
        a, b = sincos(φ)
        u = R * (1 - a / (a + b))
        v = R * (1 - b / (a + b))
        return f(triangle(u, v)) * R / (a + b)^2
    end

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = _zeros(T, 2)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(ts) = Unitful.ustrip.(integrandunits, integrand(ts))
    # Integrate only the unitless values
    ∫ = HCubature.hcubature(uintegrand, _zeros(FP, 2), (T(1), T(π / 2)), rule.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2 * Meshes.area(triangle) * integrandunits .* ∫
end

################################################################################
#                               jacobian
################################################################################

_has_analytical(::Type{T}) where {T <: Meshes.Triangle} = true

################################################################################
#                              Parametric
################################################################################

function _parametric(triangle::Meshes.Triangle, t1, t2)
    if (t1 < 0 || t1 > 1) || (t2 < 0 || t2 > 1)
        msg = "triangle(t1, t2) is not defined for (t1, t2) outside [0, 1]^2."
        throw(DomainError((t1, t2), msg))
    end

    # Form a line segment between sides
    a = triangle(0, t2)
    b = triangle(t2, 0)
    segment = Meshes.Segment(a, b)

    return segment(t1)
end
