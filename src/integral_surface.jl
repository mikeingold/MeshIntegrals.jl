################################################################################
#                               Gauss-Legendre
################################################################################

# Generalized method
function _integral_2d(f, geometry2d, settings::GaussLegendre)
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

function integral(
    f,
    disk::Meshes.Ball{2,T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    return _integral_2d(f, disk, settings)
end

function integral(
    f,
    box::Meshes.Box{2,T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    return _integral_2d(f, box, settings)
end

function integral(
    f,
    cyl::Meshes.CylinderSurface{T},
    settings::GaussLegendre
) where {T}
    error("Integrating a CylinderSurface{T} with GaussLegendre not supported.")
    # TODO Planned to support in the future
    # Waiting for resolution on whether CylinderSurface includes the terminating disks
    # on its surface by definition, and whether there will be parametric function to
    # generate those.
end

function integral(
    f,
    disk::Meshes.Disk{T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    # A Disk is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    return _integral_2d(f, disk, settings)
end

"""
    integral(f, triangle::Meshes.Triangle, ::GaussLegendre)

Like [`integral`](@ref) but integrates over the surface of a `triangle`
by transforming the triangle into a polar-barycentric coordinate system and
using a Gauss-Legendre quadrature rule along each barycentric dimension of the
triangle.
"""
function integral(
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

function integral(
    f,
    paraboloid::Meshes.ParaboloidSurface{T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    # A ParaboloidSurface is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    return _integral_2d(f, paraboloid, settings)
end

function integral(
    f,
    sphere::Meshes.Sphere{3,T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,3,T)

    return _integral_2d(f, sphere, settings)
end

function integral(
    f,
    torus::Meshes.Torus{T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,3,T)

    return _integral_2d(f, torus, settings)
end


################################################################################
#                               Gauss-Kronrod
################################################################################

# Generalized method
function _integral_2d(f, geometry2d, settings::GaussKronrod)
    function paramfactor(uv)
        J = jacobian(geometry2d, uv)
        return norm(J[1] × J[2])
    end

    integrand(u,v) = f(geometry2d(u,v)) * paramfactor([u,v])
    innerintegral(v) = QuadGK.quadgk(u -> integrand(u,v), 0, 1; settings.kwargs...)[1]
    return QuadGK.quadgk(v -> innerintegral(v), 0, 1; settings.kwargs...)[1]
end

function integral(
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

function integral(
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

function integral(
    f,
    cyl::Meshes.CylinderSurface{T},
    settings::GaussKronrod
) where {T}
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
    integral(f, triangle::Meshes.Triangle, ::GaussKronrod)

Like [`integral`](@ref) but integrates over the surface of a `triangle` using nested
Gauss-Kronrod quadrature rules along each barycentric dimension of the triangle.
"""
function integral(
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

function integral(
    f,
    paraboloid::Meshes.ParaboloidSurface{T},
    settings::GaussKronrod
) where {T}
    # Validate the provided integrand function
    # A Torus is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    return _integral_2d(f, paraboloid, settings)
end

function integral(
    f,
    sphere::Meshes.Sphere{3,T},
    settings::GaussKronrod
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,3,T)

    return _integral_2d(f, sphere, settings)
end

function integral(
    f,
    torus::Meshes.Torus{T},
    settings::GaussKronrod
) where {T}
    # Validate the provided integrand function
    # A Torus is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    return _integral_2d(f, torus, settings)
end


################################################################################
#                               HCubature
################################################################################

# Generalized method
function _integral_2d(f, geometry2d, settings::HAdaptiveCubature)
    function paramfactor(uv)
        J = jacobian(geometry2d, uv)
        return norm(J[1] × J[2])
    end

    integrand(uv) = paramfactor(uv) * f(geometry2d(uv...))
    return hcubature(integrand, [0,0], [1,1]; settings.kwargs...)[1]
end

function integral(
    f,
    disk::Meshes.Ball{2,T},
    settings::HAdaptiveCubature
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    return _integral_2d(f, disk, settings)
end

function integral(
    f,
    box::Meshes.Box{2,T},
    settings::HAdaptiveCubature
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,2,T)

    return _integral_2d(f, box, settings)
end

function integral(
    f,
    cyl::Meshes.CylinderSurface{T},
    settings::HAdaptiveCubature
) where {T}
    error("Integrating a CylinderSurface{T} with HAdaptiveCubature not supported.")
    # TODO Planned to support in the future
end

function integral(
    f,
    disk::Meshes.Disk{T},
    settings::HAdaptiveCubature
) where {T}
    # Validate the provided integrand function
    # A Disk is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    return _integral_2d(f, disk, settings)
end

"""
    integral(f, triangle::Meshes.Triangle, ::GaussKronrod)

Like [`integral`](@ref) but integrates over the surface of a `triangle` by
transforming the triangle into a polar-barycentric coordinate system and using
an h-adaptive cubature rule.
"""
function integral(
    f,
    triangle::Meshes.Ngon{3,Dim,T},
    settings::HAdaptiveCubature
) where {Dim, T}
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

function integral(
    f,
    paraboloid::Meshes.ParaboloidSurface{T},
    settings::HAdaptiveCubature
) where {T}
    # Validate the provided integrand function
    # A ParaboloidSurface is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    return _integral_2d(f, paraboloid, settings)
end

function integral(
    f,
    sphere::Meshes.Sphere{3,T},
    settings::HAdaptiveCubature
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,3,T)

    return _integral_2d(f, sphere, settings)
end

function integral(
    f,
    torus::Meshes.Torus{T},
    settings::HAdaptiveCubature
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,3,T)

    return _integral_2d(f, torus, settings)
end
