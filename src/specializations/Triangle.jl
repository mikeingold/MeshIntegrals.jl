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

"""
    integral(f, triangle::Triangle, ::GaussLegendre; FP=Float64)

Like [`integral`](@ref) but integrates over the surface of a `triangle`
by transforming the triangle into a polar-barycentric coordinate system and
using a Gauss-Legendre quadrature rule along each barycentric dimension of the
triangle.
"""
function integral(
        f::F,
        triangle::Meshes.Ngon{3},
        rule::GaussLegendre;
        diff_method::DM = FiniteDifference(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = _gausslegendre(FP, rule.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformations:
    #   xᵢ [-1,1] ↦ R [0,1]
    #   xⱼ [-1,1] ↦ φ [0,π/2]
    uR(xᵢ) = T(1 // 2) * (xᵢ + 1)
    uφ(xⱼ) = T(π / 4) * (xⱼ + 1)

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

"""
    integral(f, triangle::Triangle, ::GaussKronrod; FP=Float64)

Like [`integral`](@ref) but integrates over the surface of a `triangle` using nested
Gauss-Kronrod quadrature rules along each barycentric dimension of the triangle.
"""
function integral(
        f::F,
        triangle::Meshes.Ngon{3},
        rule::GaussKronrod;
        diff_method::DM = FiniteDifference(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    # Integrate the Barycentric triangle in (u,v)-space: (0,0), (0,1), (1,0)
    #   i.e. \int_{0}^{1} \int_{0}^{1-u} f(u,v) dv du
    ∫u(u) = QuadGK.quadgk(v -> f(triangle(u, v)), zero(FP), FP(1 - u); rule.kwargs...)[1]
    ∫ = QuadGK.quadgk(∫u, zero(FP), one(FP); rule.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2 * Meshes.area(triangle) .* ∫
end

"""
    integral(f, triangle::Triangle, ::GaussKronrod; FP=Float64)

Like [`integral`](@ref) but integrates over the surface of a `triangle` by
transforming the triangle into a polar-barycentric coordinate system and using
an h-adaptive cubature rule.
"""
function integral(
        f::F,
        triangle::Meshes.Ngon{3},
        rule::HAdaptiveCubature;
        diff_method::DM = FiniteDifference(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
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
    ∫ = HCubature.hcubature(integrand, zeros(FP, 2), FP[1, π / 2], rule.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2 * Meshes.area(triangle) .* ∫
end
