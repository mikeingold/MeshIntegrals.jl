################################################################################
#                              Common Methods
################################################################################

function surfaceintegral(
    f::F,
    circle::Meshes.Circle,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Circle is parametrically 1D, whereas Disk is parametrically 2D
    # Circle and Disk are both definitionally embedded in 3D-space
    # Convert to a Disk, integrate its surface
    return surfaceintegral(f, Disk(circle.plane, circle.radius), settings)
end

################################################################################
#                               Gauss-Legendre
################################################################################

function surfaceintegral(
    f,
    box::Meshes.Box{2,T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformation: u,v [-1,1] ↦ s,t [0,1]
    s(u) = 0.5u + 0.5
    t(v) = 0.5v + 0.5
    point(xi,xj) = box(s(xi), t(xj))

    # Calculate weight-node product
    g(((wi,wj), (xi,xj))) = wi * wj * f(point(xi,xj))

    # Calculate 2D Gauss-Legendre integral of f over parametric coordinates [-1,1]^2
    # Apply a linear domain-correction factor [-1,1]^2 ↦ area(box)
    return 0.25 * area(box) .* sum(g, zip(wws,xxs))
end

function surfaceintegral(
    f,
    disk::Meshes.Disk{T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    # A Disk is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformation: u,v [-1,1] ↦ s,t [0,1]
    s(u) = 0.5u + 0.5
    t(v) = 0.5v + 0.5
    point(xi,xj) = disk(s(xi), t(xj))

    # Calculate weight-node product with curvilinear correction
    g(((wi,wj), (xρ,xϕ))) = wi * wj * f(point(xρ,xϕ)) * disk.radius * s(xρ)

    # Calculate 2D Gauss-Legendre integral of f over parametric coordinates [-1,1]²
    # Apply curvilinear domain-correction factor [-1,1]² ↦ [0,1]² ↦ [0,ρ]x[0,2π]
    return (0.25 * area(disk)) .* sum(g, zip(wws,xxs))
end

function surfaceintegral(
    f,
    disk::Meshes.Sphere{2,T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    # A Sphere{2,T} is simply a circle in 2D-space
    _validate_integrand(f,2,T)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformation: u,v [-1,1] ↦ s,t [0,1]
    s(u) = 0.5u + 0.5
    t(v) = 0.5v + 0.5
    point(xi,xj) = disk(s(xi), t(xj))

    # Calculate weight-node product with curvilinear correction
    g(((wi,wj), (xρ,xϕ))) = wi * wj * f(point(xρ,xϕ)) * disk.radius * s(xρ)

    # Calculate 2D Gauss-Legendre integral of f over parametric coordinates [-1,1]²
    # Apply curvilinear domain-correction factor [-1,1]² ↦ [0,1]² ↦ [0,ρ]x[0,2π]
    return (0.25 * area(disk)) .* sum(g, zip(wws,xxs))
end

"""
    surfaceintegral(f, triangle::Meshes.Triangle, ::GaussLegendre)

Like [`surfaceintegral`](@ref) but integrates over the surface of a `triangle`
using a Gauss-Legendre quadrature rule along each barycentric dimension of the
triangle.
"""
function surfaceintegral(
    f,
    triangle::Meshes.Ngon{3,Dim,T},
    settings::GaussLegendre
) where {Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Domain transformation: u,v [-1,1] ↦ s,t [0,1]
    s(u) = 0.5u + 0.5
    t(v) = 0.5v + 0.5
    point(xi,xj) = triangle(s(xi), t(xj))

    # Determine output type of f at a Point inside the triangle
    # Define an applicable zero value
    fzero = zero(f(point(-0.5,-0.5)))

    # Calculate weight-node product
    function g(((wi,wj), (xi,xj)))
        if 0 < (s(xi) + t(xj)) < 1
            # Valid coordinate (inside triangle)
            return wi * wj * f(point(xi,xj))
        else
            # Invalid coordinate (outside triangle)
            return fzero
        end
    end

    # Calculate 2D Gauss-Legendre integral of f over Barycentric coordinates [-1,1]^2
    # Apply a linear domain-correction factor [-1,1]^2 ↦ area(triangle)
    return 0.5 * area(triangle) .* sum(g, zip(wws,xxs))
end

################################################################################
#                               Gauss-Kronrod
################################################################################

function surfaceintegral(
    f,
    box::Meshes.Box{2,T},
    settings::GaussKronrod
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    # Integrate the box in parametric (u,v)-space
    innerintegral(u) = QuadGK.quadgk(v -> f(box(u,v)), 0, 1; settings.kwargs...)[1]
    outerintegral = QuadGK.quadgk(innerintegral, 0, 1; settings.kwargs...)[1]

    # Apply a linear domain-correction factor 1 ↦ area(box)
    return area(box) .* outerintegral
end

function surfaceintegral(
    f,
    disk::Meshes.Disk{T},
    settings::GaussKronrod
) where {T}
    # Validate the provided integrand function
    # A Disk is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    # Integrate the box in parametric (ρ,ϕ)-space
    innerintegral(ϕ) = QuadGK.quadgk(ρ -> f(disk(ρ,ϕ)) * disk.radius * ρ, 0, 1; settings.kwargs...)[1]
    outerintegral = QuadGK.quadgk(innerintegral, 0, 1; settings.kwargs...)[1]

    # Apply a linear domain-correction factor [0,1]² ↦ [0,ρ]x[0,2π]
    return (2π * disk.radius) .* outerintegral
end

function surfaceintegral(
    f,
    disk::Meshes.Sphere{2,T},
    settings::GaussKronrod
) where {T}
    # Validate the provided integrand function
    # A Sphere{2,T} is simply a circle in 2D-space
    _validate_integrand(f,2,T)

    # Integrate the box in parametric (ρ,ϕ)-space
    innerintegral(ϕ) = QuadGK.quadgk(ρ -> f(disk(ρ,ϕ)) * disk.radius * ρ, 0, 1; settings.kwargs...)[1]
    outerintegral = QuadGK.quadgk(innerintegral, 0, 1; settings.kwargs...)[1]

    # Apply a linear domain-correction factor [0,1]² ↦ [0,ρ]x[0,2π]
    return (2π * disk.radius) .* outerintegral
end

"""
    surfaceintegral(f, triangle::Meshes.Triangle, ::GaussKronrod)

Like [`surfaceintegral`](@ref) but integrates over the surface of a `triangle`
using a Gauss-Kronrod quadrature rule along each barycentric dimension of the
triangle.
"""
function surfaceintegral(
    f,
    triangle::Meshes.Ngon{3,Dim,T},
    settings::GaussKronrod
) where {Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Integrate the Barycentric triangle in (u,v)-space: (0,0), (0,1), (1,0)
    #   i.e. \int_{0}^{1} \int_{0}^{1-u} f(u,v) dv du
    innerintegral(u) = QuadGK.quadgk(v -> f(triangle(u,v)), 0, 1-u; settings.kwargs...)[1]
    outerintegral = QuadGK.quadgk(innerintegral, 0, 1; settings.kwargs...)[1]

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2.0 * area(triangle) .* outerintegral
end
