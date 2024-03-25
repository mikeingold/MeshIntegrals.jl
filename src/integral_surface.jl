################################################################################
#                          Generalized 2D Methods
################################################################################

function _integral_2d(
    f,
    geometry2d,
    settings::GaussLegendre
)
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]²
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformation: x [-1,1] ↦ u,v [0,1]
    u(x) = 0.5x + 0.5
    v(x) = 0.5x + 0.5
    point(xi,xj) = geometry2d(u(xi), v(xj))

    function paramfactor(xi, xj)
        uv = [u(xi), v(xj)]
        J = jacobian(geometry2d, uv)
        return norm(J[1] × J[2])
    end

    integrand(((wi,wj), (xi,xj))) = wi * wj * f(point(xi,xj)) * paramfactor(xi,xj)

    return 0.25 .* sum(integrand, zip(wws,xxs))
end

function _integral_2d(
    f,
    geometry2d,
    settings::GaussKronrod
)
    function paramfactor(uv)
        J = jacobian(geometry2d, uv)
        return norm(J[1] × J[2])
    end

    integrand(u,v) = f(geometry2d(u,v)) * paramfactor([u,v])
    innerintegral(v) = QuadGK.quadgk(u -> integrand(u,v), 0, 1; settings.kwargs...)[1]
    return QuadGK.quadgk(v -> innerintegral(v), 0, 1; settings.kwargs...)[1]
end

function _integral_2d(
    f,
    geometry2d,
    settings::HAdaptiveCubature
)
    function paramfactor(uv)
        J = jacobian(geometry2d, uv)
        return norm(J[1] × J[2])
    end

    integrand(uv) = paramfactor(uv) * f(geometry2d(uv...))
    return hcubature(integrand, [0,0], [1,1]; settings.kwargs...)[1]
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
        sidelength = norm(cyl(φ,1) - cyl(φ,0))
        return sidelength * quadgk(z -> f(cyl(φ,z)), 0, 1; settings.kwargs...)[1]
    end
    sides = (2π * cyl.radius) .* quadgk(φ -> sides_innerintegral(φ), 0, 1; settings.kwargs...)[1]

    # Integrate the top and bottom disks
    # \int ( \int r f(r̄) dr ) dφ
    function disk_innerintegral(φ,plane,z)
        # Parameterize the top surface of the cylinder
        rimedge = cyl(φ,z)
        centerpoint = plane.p
        r̄ = rimedge - centerpoint
        radius = norm(r̄)
        point(r) = centerpoint + (r / radius) * r̄

        return radius^2 * quadgk(r -> r * f(point(r)), 0, 1; settings.kwargs...)[1]
    end
    top    = 2π .* quadgk(φ -> disk_innerintegral(φ,cyl.top,1), 0, 1; settings.kwargs...)[1]
    bottom = 2π .* quadgk(φ -> disk_innerintegral(φ,cyl.bot,0), 0, 1; settings.kwargs...)[1]

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
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Change of variables: s,t [-Inf,Inf] ↦ x,y [-1,1]
    s(x) = x / (1 - x^2)
    t(y) = y / (1 - y^2)
    Δ(u) = (1 + u^2) / (1 - u^2)^2

    integrand(((wi,wj), (xi,xj))) = wi * wj * f(plane(s(xi), t(xj))) * Δ(xi) * Δ(xj)
    return sum(integrand, zip(wws,xxs))
end

function integral(
    f::F,
    plane::Meshes.Plane{T},
    settings::GaussKronrod
) where {F<:Function, T}
    integrand(u,v) = f(plane(u, v))
    innerintegral(v) = QuadGK.quadgk(u -> integrand(u,v), -Inf, Inf; settings.kwargs...)[1]
    return QuadGK.quadgk(v -> innerintegral(v), -Inf, Inf; settings.kwargs...)[1]
end

function integral(
    f::F,
    plane::Meshes.Plane{T},
    settings::HAdaptiveCubature
) where {F<:Function, T}
    # Change of variables: s,t [-Inf,Inf] ↦ x,y [-1,1]
    s(x) = x / (1 - x^2)
    t(y) = y / (1 - y^2)
    Δ(u) = (1 + u^2) / (1 - u^2)^2

    integrand(x, y) = f(plane(s(x), t(y))) * Δ(x) * Δ(y)
    integrand(xy) = integrand(xy[1], xy[2])
    return HCubature.hcubature(integrand, [-1, -1], [1, 1]; settings.kwargs...)[1]
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
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformations:
    #   xᵢ [-1,1] ↦ R [0,1]
    #   xⱼ [-1,1] ↦ φ [0,π/2]
    uR(xᵢ) = 0.5 * (xᵢ + 1)
    uφ(xⱼ) = (π/4) * (xⱼ + 1)

    # Integrate the Barycentric triangle by transforming it into polar coordinates
    #   with a modified radius
    #     R = r ( sinφ + cosφ )
    #   s.t. integration bounds become rectangular
    #     R ∈ [0, 1] and φ ∈ [0, π/2]
    function integrand(((wᵢ,wⱼ), (xᵢ,xⱼ)))
        R = uR(xᵢ)
        φ = uφ(xⱼ)
        a,b = sincos(φ)
        u = R * (1 - a/(a+b))
        v = R * (1 - b/(a+b))
        return wᵢ * wⱼ * f(triangle(u,v)) * R / (a+b)^2
    end

    # Calculate 2D Gauss-Legendre integral over modified-polar-Barycentric coordinates
    # Apply a linear domain-correction factor
    return (π/4) * area(triangle) .* sum(integrand, zip(wws,xxs))
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
    return 2.0 * area(triangle) .* outerintegral
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
        u = R * (1 - a/(a+b))
        v = R * (1 - b/(a+b))
        return f(triangle(u,v)) * R / (a+b)^2
    end
    intval = hcubature(integrand, [0,0], [1,π/2], settings.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2.0 * area(triangle) .* intval
end
