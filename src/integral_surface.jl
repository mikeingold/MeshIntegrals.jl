################################################################################
#                               Gauss-Legendre
################################################################################

function surfaceintegral(
    f,
    disk::Meshes.Ball{2,T},
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
    point(xi,xj) = disk(s(xi), t(xj))

    # Calculate weight-node product with curvilinear correction
    g(((wi,wj), (xρ,xϕ))) = wi * wj * f(point(xρ,xϕ)) * disk.radius * s(xρ)

    # Calculate 2D Gauss-Legendre integral of f over parametric coordinates [-1,1]²
    # Apply curvilinear domain-correction factor [-1,1]² ↦ [0,1]² ↦ [0,ρ]x[0,2π]
    return (0.25 * area(disk)) .* sum(g, zip(wws,xxs))
end

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
    disk::Meshes.Ball{2,T},
    settings::GaussKronrod
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    # Map parametric ρ₀ [0,1] ↦ ρ [0, disk.radius]
    ρ(ρ₀) = disk.radius * ρ₀

    # Integrate the box in parametric (ρ₀,φ₀)-space [0,1]
    integrand(ρ₀,φ₀) = f(disk(ρ₀,φ₀)) * ρ(ρ₀)
    innerintegral(φ₀) = QuadGK.quadgk(ρ₀ -> integrand(ρ₀,φ₀), 0, 1; settings.kwargs...)[1]
    outerintegral = QuadGK.quadgk(φ₀ -> innerintegral(φ₀), 0, 1; settings.kwargs...)[1]

    # Apply a linear domain-correction factor [0,1]² ↦ [0,ρ]x[0,2π]
    return (2π * disk.radius) .* outerintegral
end

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

    # Map parametric ρ₀ [0,1] ↦ ρ [0, disk.radius]
    ρ(ρ₀) = disk.radius * ρ₀

    # Integrate the box in parametric (ρ₀,φ₀)-space [0,1]
    integrand(ρ₀,φ₀) = f(disk(ρ₀,φ₀)) * ρ(ρ₀)
    innerintegral(φ₀) = QuadGK.quadgk(ρ₀ -> integrand(ρ₀,φ₀), 0, 1; settings.kwargs...)[1]
    outerintegral = QuadGK.quadgk(φ₀ -> innerintegral(φ₀), 0, 1; settings.kwargs...)[1]

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


################################################################################
#                               HCubature
################################################################################

function surfaceintegral(
    f,
    disk::Meshes.Ball{2,T},
    settings::HAdaptiveCubature
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    # Map parametric ρ₀ [0,1] ↦ ρ [0, disk.radius]
    ρ(ρ₀) = disk.radius * ρ₀

    # Integrate the box in parametric (ρ₀,φ₀)-space [0,1]
    integrand(ρ₀,φ₀) = f(disk(ρ₀,φ₀)) * ρ(ρ₀)
    integrand(ρφ) = integrand(ρφ[1], ρφ[2])
    intval = hcubature(integrand, [0,0], [1,1], settings.kwargs...)[1]

    # Apply a linear domain-correction factor [0,1]² ↦ [0,ρ]x[0,2π]
    return (2π * disk.radius) .* intval
end

function surfaceintegral(
    f,
    box::Meshes.Box{2,T},
    settings::HAdaptiveCubature
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    # Integrate the box in parametric (u,v)-space
    integrand(u,v) = f(box(u,v))
    integrand(uv) = integrand(uv[1], uv[2])
    intval = hcubature(integrand, [0,0], [1,1], settings.kwargs...)[1]

    # Apply a linear domain-correction factor 1 ↦ area(box)
    return area(box) .* intval
end

function surfaceintegral(
    f,
    disk::Meshes.Disk{T},
    settings::HAdaptiveCubature
) where {T}
    # Validate the provided integrand function
    # A Disk is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    # Map parametric ρ₀ [0,1] ↦ ρ [0, disk.radius]
    ρ(ρ₀) = disk.radius * ρ₀

    # Integrate the box in parametric (ρ₀,φ₀)-space [0,1]
    integrand(ρ₀,φ₀) = f(disk(ρ₀,φ₀)) * ρ(ρ₀)
    integrand(ρφ) = integrand(ρφ[1], ρφ[2])
    intval = hcubature(integrand, [0,0], [1,1], settings.kwargs...)[1]

    # Apply a linear domain-correction factor [0,1]² ↦ [0,ρ]x[0,2π]
    return (2π * disk.radius) .* intval
end

"""
    surfaceintegral(f, triangle::Meshes.Triangle, ::GaussKronrod)

Like [`surfaceintegral`](@ref) but integrates over the surface of a `triangle`
by transforming the triangle into a polar-barycentric coordinate system and
using an h-adaptive cubature rule.
"""
function surfaceintegral(
    f,
    triangle::Meshes.Ngon{3,Dim,T},
    settings::HAdaptiveCubature
) where {Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Integrate in Barycentric-polar-space with transformed R
    function integrand(Rφ)
        R,φ = Rφ
        a,b = sincos(φ)
        u = R * (1 - a/(a+b))
        v = R * (1 - b/(a+b))
        return f(triangle(u,v)) * R / (2*a*b)
    end

    return hcubature(integrand, [0,0], [1,π/2], settings.kwargs...)[1]
end
