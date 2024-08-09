################################################################################
#                          Generalized 2D Methods
################################################################################

function _integral_2d(
    f,
    geometry2d,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]²
    xs, ws = _gausslegendre(FP, settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformation: x [-1,1] ↦ t [0,1]
    t(x) = FP(1/2) * x + FP(1/2)
    point(xi, xj) = geometry2d(t(xi), t(xj))

    # Integrate f over the geometry
    integrand(((wi,wj), (xi,xj))) = wi * wj * f(point(xi,xj)) * differential(geometry2d, t.([xi, xj]))
    return FP(1/4) .* sum(integrand, zip(wws,xxs))
end

function _integral_2d(
    f,
    geometry2d,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    integrand(u,v) = f(geometry2d(u,v)) * differential(geometry2d, [u,v])
    ∫₁(v) = QuadGK.quadgk(u -> integrand(u,v), FP(0), FP(1); settings.kwargs...)[1]
    return QuadGK.quadgk(v -> ∫₁(v), FP(0), FP(1); settings.kwargs...)[1]
end

function _integral_2d(
    f,
    geometry2d,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    integrand(uv) = f(geometry2d(uv...)) * differential(geometry2d, uv)

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    integrandunits = Unitful.units(integrand(FP[0.5,0.5]))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, FP[0,0], FP[1,1]; settings.kwargs...)[1]

    # Reapply units
    return value * integrandunits
end


################################################################################
#                  Specialized Methods for CylinderSurface
################################################################################

function integral(
    f::F,
    cyl::Meshes.CylinderSurface,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    error("Integrating a CylinderSurface with GaussLegendre not supported.")
    # TODO Planned to support in the future
    # Waiting for resolution on whether CylinderSurface includes the terminating disks
    # on its surface by definition, and whether there will be parametric function to
    # generate those.
end

function integral(
    f::F,
    cyl::Meshes.CylinderSurface,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Integrate the rounded sides of the cylinder's surface
    # \int ( \int f(r̄) dz ) dφ
    function sides_inner∫(φ)
        sidelength = norm(cyl(φ,FP(1)) - cyl(φ,FP(0)))
        return sidelength * QuadGK.quadgk(z -> f(cyl(φ,z)), FP(0), FP(1); settings.kwargs...)[1]
    end
    sides = (FP(2π) * cyl.radius) .* QuadGK.quadgk(φ -> sides_inner∫(φ), FP(0), FP(1); settings.kwargs...)[1]

    # Integrate the top and bottom disks
    # \int ( \int r f(r̄) dr ) dφ
    function disk_inner∫(φ,plane,z)
        # Parameterize the top surface of the cylinder
        rimedge = cyl(φ,FP(z))
        centerpoint = plane.p
        r̄ = rimedge - centerpoint
        radius = norm(r̄)
        point(r) = centerpoint + (r / radius) * r̄

        return radius^2 * QuadGK.quadgk(r -> r * f(point(r)), FP(0), FP(1); settings.kwargs...)[1]
    end
    top    = FP(2π) .* QuadGK.quadgk(φ -> disk_inner∫(φ,cyl.top,1), FP(0), FP(1); settings.kwargs...)[1]
    bottom = FP(2π) .* QuadGK.quadgk(φ -> disk_inner∫(φ,cyl.bot,0), FP(0), FP(1); settings.kwargs...)[1]

    return sides + top + bottom
end

function integral(
    f::F,
    cyl::Meshes.CylinderSurface,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    error("Integrating a CylinderSurface with HAdaptiveCubature not supported.")
    # TODO Planned to support in the future
end

################################################################################
#                      Specialized Methods for Plane
################################################################################

function integral(
    f::F,
    plane::Meshes.Plane,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]²
    xs, ws = _gausslegendre(FP, settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Normalize the Plane's orthogonal vectors
    plane = Plane(plane.p, Meshes.unormalize(plane.u), Meshes.unormalize(plane.v))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f over the Plane
    integrand(((wi,wj), (xi,xj))) = wi * wj * f(plane(t(xi), t(xj))) * t′(xi) * t′(xj)
    return sum(integrand, zip(wws,xxs))
end

function integral(
    f::F,
    plane::Meshes.Plane,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Plane's orthogonal vectors
    plane = Plane(plane.p, Meshes.unormalize(plane.u), Meshes.unormalize(plane.v))

    # Integrate f over the Plane
    inner∫(v) = QuadGK.quadgk(u -> f(plane(u,v)), FP(-Inf), FP(Inf); settings.kwargs...)[1]
    return QuadGK.quadgk(v -> inner∫(v), FP(-Inf), FP(Inf); settings.kwargs...)[1]
end

function integral(
    f::F,
    plane::Meshes.Plane,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Plane's orthogonal vectors
    plane = Plane(plane.p, Meshes.unormalize(plane.u), Meshes.unormalize(plane.v))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f over the Plane
    integrand(x::AbstractVector) = f(plane(t(x[1]), t(x[2]))) * t′(x[1]) * t′(x[2])
    return HCubature.hcubature(integrand, FP[-1, -1], FP[1, 1]; settings.kwargs...)[1]
end

################################################################################
#                    Specialized Methods for Triangle
################################################################################

"""
    integral(f, triangle::Meshes.Triangle, ::GaussLegendre)

Like [`integral`](@ref) but integrates over the surface of a `triangle`
by transforming the triangle into a polar-barycentric coordinate system and
using a Gauss-Legendre quadrature rule along each barycentric dimension of the
triangle.
"""
function integral(
    f::F,
    triangle::Meshes.Ngon{3},
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = _gausslegendre(FP, settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformations:
    #   xᵢ [-1,1] ↦ R [0,1]
    #   xⱼ [-1,1] ↦ φ [0,π/2]
    uR(xᵢ) = T(1/2) * (xᵢ + 1)
    uφ(xⱼ) = T(π/4) * (xⱼ + 1)

    # Integrate the Barycentric triangle by transforming it into polar coordinates
    #   with a modified radius
    #     R = r ( sinφ + cosφ )
    #   s.t. integration bounds become rectangular
    #     R ∈ [0, 1] and φ ∈ [0, π/2]
    function integrand(((wᵢ,wⱼ), (xᵢ,xⱼ)))
        R = uR(xᵢ)
        φ = uφ(xⱼ)
        a,b = sincos(φ)
        u = R * (1 - a / (a + b))
        v = R * (1 - b / (a + b))
        return wᵢ * wⱼ * f(triangle(u, v)) * R / (a + b)^2
    end

    # Calculate 2D Gauss-Legendre integral over modified-polar-Barycentric coordinates
    # Apply a linear domain-correction factor
    return FP(π/4) * area(triangle) .* sum(integrand, zip(wws,xxs))
end

"""
    integral(f, triangle::Meshes.Triangle, ::GaussKronrod)

Like [`integral`](@ref) but integrates over the surface of a `triangle` using nested
Gauss-Kronrod quadrature rules along each barycentric dimension of the triangle.
"""
function integral(
    f::F,
    triangle::Meshes.Ngon{3},
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Integrate the Barycentric triangle in (u,v)-space: (0,0), (0,1), (1,0)
    #   i.e. \int_{0}^{1} \int_{0}^{1-u} f(u,v) dv du
    inner∫(u) = QuadGK.quadgk(v -> f(triangle(u,v)), FP(0), FP(1-u); settings.kwargs...)[1]
    outer∫ = QuadGK.quadgk(inner∫, FP(0), FP(1); settings.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2 * area(triangle) .* outer∫
end

"""
    integral(f, triangle::Meshes.Triangle, ::GaussKronrod)

Like [`integral`](@ref) but integrates over the surface of a `triangle` by
transforming the triangle into a polar-barycentric coordinate system and using
an h-adaptive cubature rule.
"""
function integral(
    f::F,
    triangle::Meshes.Ngon{3},
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Integrate the Barycentric triangle by transforming it into polar coordinates
    #   with a modified radius
    #     R = r ( sinφ + cosφ )
    #   s.t. integration bounds become rectangular
    #     R ∈ [0, 1] and φ ∈ [0, π/2]
    function integrand(Rφ)
        R,φ = Rφ
        a,b = sincos(φ)
        u = R * (1 - a / (a + b))
        v = R * (1 - b / (a + b))
        return f(triangle(u, v)) * R / (a + b)^2
    end
    intval = HCubature.hcubature(integrand, FP[0, 0], FP[1, π/2], settings.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2 * area(triangle) .* intval
end
