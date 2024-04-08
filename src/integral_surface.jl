################################################################################
#                          Generalized 2D Methods
################################################################################

function _integral_2d(
    f,
    geometry2d::G,
    settings::GaussLegendre
) where {Dim, T, G<:Meshes.Geometry{Dim,T}}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]²
    xs, ws = _gausslegendre(T, settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformation: x [-1,1] ↦ t [0,1]
    t(x) = T(1/2) * x + T(1/2)
    point(xi, xj) = geometry2d(t(xi), t(xj))

    # Integrate f over the geometry
    integrand(((wi,wj), (xi,xj))) = wi * wj * f(point(xi,xj)) * differential(geometry2d, [u(xi), v(xj)])

    return T(1/4) .* sum(integrand, zip(wws,xxs))
end

function _integral_2d(
    f,
    geometry2d::G,
    settings::GaussKronrod
) where {Dim, T, G<:Meshes.Geometry{Dim,T}}
    integrand(u,v) = f(geometry2d(u,v)) * differential(geometry2d, [u,v])
    ∫₁(v) = QuadGK.quadgk(u -> integrand(u,v), T(0), T(1); settings.kwargs...)[1]
    return QuadGK.quadgk(v -> ∫₁(v), T(0), T(1); settings.kwargs...)[1]
end

function _integral_2d(
    f,
    geometry2d::G,
    settings::HAdaptiveCubature
) where {Dim, T, G<:Meshes.Geometry{Dim,T}}
    integrand(uv) = f(geometry2d(uv...)) * differential(geometry2d, uv)
    return HCubature.hcubature(integrand, T[0,0], T[1,1]; settings.kwargs...)[1]
end


################################################################################
#                  Specialized Methods for CylinderSurface
################################################################################

function integral(
    f::F,
    cyl::Meshes.CylinderSurface{T},
    settings::GaussLegendre
) where {F<:Function, T}
    error("Integrating a CylinderSurface{T} with GaussLegendre not supported.")
    # TODO Planned to support in the future
    # Waiting for resolution on whether CylinderSurface includes the terminating disks
    # on its surface by definition, and whether there will be parametric function to
    # generate those.
end

function integral(
    f::F,
    cyl::Meshes.CylinderSurface{T},
    settings::GaussKronrod
) where {F<:Function, T}
    # Validate the provided integrand function
    # A CylinderSurface is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    # Integrate the rounded sides of the cylinder's surface
    # \int ( \int f(r̄) dz ) dφ
    function sides_innerintegral(φ)
        sidelength = norm(cyl(φ,T(1)) - cyl(φ,T(0)))
        return sidelength * quadgk(z -> f(cyl(φ,z)), T(0), T(1); settings.kwargs...)[1]
    end
    sides = (T(2π) * cyl.radius) .* quadgk(φ -> sides_innerintegral(φ), T(0), T(1); settings.kwargs...)[1]

    # Integrate the top and bottom disks
    # \int ( \int r f(r̄) dr ) dφ
    function disk_innerintegral(φ,plane,z)
        # Parameterize the top surface of the cylinder
        rimedge = cyl(φ,T(z))
        centerpoint = plane.p
        r̄ = rimedge - centerpoint
        radius = norm(r̄)
        point(r) = centerpoint + (r / radius) * r̄

        return radius^2 * quadgk(r -> r * f(point(r)), T(0), T(1); settings.kwargs...)[1]
    end
    top    = T(2π) .* quadgk(φ -> disk_innerintegral(φ,cyl.top,1), T(0), T(1); settings.kwargs...)[1]
    bottom = T(2π) .* quadgk(φ -> disk_innerintegral(φ,cyl.bot,0), T(0), T(1); settings.kwargs...)[1]

    return sides + top + bottom
end

function integral(
    f::F,
    cyl::Meshes.CylinderSurface{T},
    settings::HAdaptiveCubature
) where {F<:Function, T}
    error("Integrating a CylinderSurface{T} with HAdaptiveCubature not supported.")
    # TODO Planned to support in the future
end

################################################################################
#                      Specialized Methods for Plane
################################################################################

function integral(
    f::F,
    plane::Meshes.Plane{T},
    settings::GaussLegendre
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Plane is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]²
    xs, ws = _gausslegendre(T, settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Change of variables: s,t [-Inf,Inf] ↦ x,y [-1,1]
    s(x) = x / (T(1) - x^2)
    t(y) = y / (T(1) - y^2)
    Δ(u) = (T(1) + u^2) / (T(1) - u^2)^2

    integrand(((wi,wj), (xi,xj))) = wi * wj * f(plane(s(xi), t(xj))) * Δ(xi) * Δ(xj)
    return sum(integrand, zip(wws,xxs))
end

function integral(
    f::F,
    plane::Meshes.Plane{T},
    settings::GaussKronrod
) where {F<:Function, T}
    integrand(u,v) = f(plane(u, v))
    innerintegral(v) = QuadGK.quadgk(u -> integrand(u,v), T(-Inf), T(Inf); settings.kwargs...)[1]
    return QuadGK.quadgk(v -> innerintegral(v), T(-Inf), T(Inf); settings.kwargs...)[1]
end

function integral(
    f::F,
    plane::Meshes.Plane{T},
    settings::HAdaptiveCubature
) where {F<:Function, T}
    # Change of variables: s,t [-Inf,Inf] ↦ x,y [-1,1]
    s(x) = x / (T(1) - x^2)
    t(y) = y / (T(1) - y^2)
    Δ(u) = (T(1) + u^2) / (T(1) - u^2)^2

    integrand(x, y) = f(plane(s(x), t(y))) * Δ(x) * Δ(y)
    integrand(xy) = integrand(xy[1], xy[2])
    return HCubature.hcubature(integrand, T[-1, -1], T[1, 1]; settings.kwargs...)[1]
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
    triangle::Meshes.Ngon{3,Dim,T},
    settings::GaussLegendre
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = _gausslegendre(T, settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformations:
    #   xᵢ [-1,1] ↦ R [0,1]
    #   xⱼ [-1,1] ↦ φ [0,π/2]
    uR(xᵢ) = T(1/2) * (xᵢ + T(1))
    uφ(xⱼ) = T(π/4) * (xⱼ + T(1))

    # Integrate the Barycentric triangle by transforming it into polar coordinates
    #   with a modified radius
    #     R = r ( sinφ + cosφ )
    #   s.t. integration bounds become rectangular
    #     R ∈ [0, 1] and φ ∈ [0, π/2]
    function integrand(((wᵢ,wⱼ), (xᵢ,xⱼ)))
        R = uR(xᵢ)
        φ = uφ(xⱼ)
        a,b = sincos(φ)
        u = R * (T(1) - a/(a+b))
        v = R * (T(1) - b/(a+b))
        return wᵢ * wⱼ * f(triangle(u,v)) * R / (a+b)^2
    end

    # Calculate 2D Gauss-Legendre integral over modified-polar-Barycentric coordinates
    # Apply a linear domain-correction factor
    return T(π/4) * area(triangle) .* sum(integrand, zip(wws,xxs))
end

"""
    integral(f, triangle::Meshes.Triangle, ::GaussKronrod)

Like [`integral`](@ref) but integrates over the surface of a `triangle` using nested
Gauss-Kronrod quadrature rules along each barycentric dimension of the triangle.
"""
function integral(
    f::F,
    triangle::Meshes.Ngon{3,Dim,T},
    settings::GaussKronrod
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Integrate the Barycentric triangle in (u,v)-space: (0,0), (0,1), (1,0)
    #   i.e. \int_{0}^{1} \int_{0}^{1-u} f(u,v) dv du
    innerintegral(u) = QuadGK.quadgk(v -> f(triangle(u,v)), 0, 1-u; settings.kwargs...)[1]
    outerintegral = QuadGK.quadgk(innerintegral, 0, 1; settings.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return T(2) * area(triangle) .* outerintegral
end

"""
    integral(f, triangle::Meshes.Triangle, ::GaussKronrod)

Like [`integral`](@ref) but integrates over the surface of a `triangle` by
transforming the triangle into a polar-barycentric coordinate system and using
an h-adaptive cubature rule.
"""
function integral(
    f::F,
    triangle::Meshes.Ngon{3,Dim,T},
    settings::HAdaptiveCubature
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Integrate the Barycentric triangle by transforming it into polar coordinates
    #   with a modified radius
    #     R = r ( sinφ + cosφ )
    #   s.t. integration bounds become rectangular
    #     R ∈ [0, 1] and φ ∈ [0, π/2]
    function integrand(Rφ)
        R,φ = Rφ
        a,b = sincos(φ)
        u = R * (T(1) - a/(a+b))
        v = R * (T(1) - b/(a+b))
        return f(triangle(u,v)) * R / (a+b)^2
    end
    intval = hcubature(integrand, T[0,0], T[1,π/2], settings.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return T(2) * area(triangle) .* intval
end
