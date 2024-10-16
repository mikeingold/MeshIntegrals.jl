# This section tests for:
# - All supported combinations of integral(f, ::Geometry, ::IntegrationAlgorithm) produce accurate results
# - Invalid applications of integral aliases (e.g. lineintegral) produce a descriptive error

@testitem "Meshes.Ball 2D" setup=[Setup] begin
    origin = Point(0, 0)
    radius = 2.8
    ball = Ball(origin, radius)

    function f(p::P) where {P <: Meshes.Point}
        ur = hypot(p.coords.x, p.coords.y)
        r = ustrip(u"m", ur)
        exp(-r^2)
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = (π - π * exp(-radius^2)) * u"m^2"
    @test integral(f, ball, GaussLegendre(100)) ≈ sol
    @test integral(f, ball, GaussKronrod()) ≈ sol
    @test integral(f, ball, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, ball, GaussLegendre(100)) ≈ vsol
    @test integral(fv, ball, GaussKronrod()) ≈ vsol
    @test integral(fv, ball, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, ball)
    @test surfaceintegral(f, ball) ≈ sol
    @test_throws "not supported" volumeintegral(f, ball)
end

@testitem "Meshes.Ball 3D" setup=[Setup] begin
    # using SpecialFunctions: erf

    center = Point(1, 2, 3)
    radius = 2.8u"m"
    ball = Ball(center, radius)

    function f(p::P) where {P <: Meshes.Point}
        ur = hypot(p.coords.x, p.coords.y)
        r = ustrip(u"m", ur)
        # exp(-r^2)
        # 1 / r
        1.0
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    # r = ustrip(u"m", radius)
    # sol = (π^(3/2) * erf(r) - 2π * exp(-r^2) * r) * u"m^3"   # for f(p) = exp(-r^2)
    # sol = 2π * radius^2 * u"m"     # for f(p) = 1/r
    sol = (4 / 3) * π * radius^3   # for f(p) = 1
    @test integral(f, ball, GaussLegendre(100)) ≈ sol
    @test_throws "not supported" integral(f, ball, GaussKronrod())≈sol
    @test integral(f, ball, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, ball, GaussLegendre(100)) ≈ vsol
    @test_throws "not supported" integral(fv, ball, GaussKronrod())≈vsol
    @test integral(fv, ball, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, ball)
    @test_throws "not supported" surfaceintegral(f, ball)
    @test volumeintegral(f, ball) ≈ sol
end

@testitem "Meshes.BezierCurve" setup=[Setup] begin
    curve = BezierCurve(
        [Point(t * u"m", sin(t) * u"m", 0.0u"m") for t in range(-pi, pi, length = 361)]
    )

    function f(p::P) where {P <: Meshes.Point}
        ux = ustrip(p.coords.x)
        (1 / sqrt(1 + cos(ux)^2)) * u"Ω/m"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 2π * u"Ω"
    @test integral(f, curve, GaussLegendre(100))≈sol rtol=0.5e-2
    @test integral(f, curve, GaussKronrod())≈sol rtol=0.5e-2
    @test integral(f, curve, HAdaptiveCubature())≈sol rtol=0.5e-2

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, curve, GaussLegendre(100))≈vsol rtol=0.5e-2
    @test integral(fv, curve, GaussKronrod())≈vsol rtol=0.5e-2
    @test integral(fv, curve, HAdaptiveCubature())≈vsol rtol=0.5e-2

    # Integral aliases
    @test lineintegral(f, curve)≈sol rtol=0.5e-2
    @test_throws "not supported" surfaceintegral(f, curve)
    @test_throws "not supported" volumeintegral(f, curve)

    # Check Bezier-specific jacobian bounds
    @test_throws DomainError jacobian(curve, [1.1])
end

@testitem "Meshes.Box 1D" setup=[Setup] begin
    a = π
    box = Box(Point(0), Point(a))

    function f(p::P) where {P <: Meshes.Point}
        t = ustrip(p.coords.x)
        sqrt(a^2 - t^2) * u"Ω/m"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = π * a^2 / 4 * u"Ω"
    @test integral(f, box, GaussLegendre(100))≈sol rtol=1e-6
    @test integral(f, box, GaussKronrod()) ≈ sol
    @test integral(f, box, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, box, GaussLegendre(100))≈vsol rtol=1e-6
    @test integral(fv, box, GaussKronrod()) ≈ vsol
    @test integral(fv, box, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, box) ≈ sol
    @test_throws "not supported" surfaceintegral(f, box)
    @test_throws "not supported" volumeintegral(f, box)
end

@testitem "Meshes.Box 2D" setup=[Setup] begin
    a = π
    box = Box(Point(0, 0), Point(a, a))

    function f(p::P) where {P <: Meshes.Point}
        x, y = ustrip.((p.coords.x, p.coords.y))
        (sqrt(a^2 - x^2) + sqrt(a^2 - y^2)) * u"Ω/m^2"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 2a * (π * a^2 / 4) * u"Ω"
    @test integral(f, box, GaussLegendre(100))≈sol rtol=1e-6
    @test integral(f, box, GaussKronrod()) ≈ sol
    @test integral(f, box, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, box, GaussLegendre(100))≈vsol rtol=1e-6
    @test integral(fv, box, GaussKronrod()) ≈ vsol
    @test integral(fv, box, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, box)
    @test surfaceintegral(f, box) ≈ sol
    @test_throws "not supported" volumeintegral(f, box)

    # Test jacobian with wrong number of parametric coordinates
    @test_throws ArgumentError jacobian(box, zeros(3))
end

@testitem "Meshes.Box 3D" setup=[Setup] begin
    a = π
    box = Box(Point(0, 0, 0), Point(a, a, a))

    function f(p::P) where {P <: Meshes.Point}
        x, y, z = ustrip.((p.coords.x, p.coords.y, p.coords.z))
        (sqrt(a^2 - x^2) + sqrt(a^2 - y^2) + sqrt(a^2 - z^2)) * u"Ω/m^3"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 3a^2 * (π * a^2 / 4) * u"Ω"
    @test integral(f, box, GaussLegendre(100))≈sol rtol=1e-6
    @test_throws "not supported" integral(f, box, GaussKronrod())
    @test integral(f, box, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, box, GaussLegendre(100))≈vsol rtol=1e-6
    @test_throws "not supported" integral(fv, box, GaussKronrod())
    @test integral(fv, box, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, box)
    @test_throws "not supported" surfaceintegral(f, box)
    @test volumeintegral(f, box) ≈ sol
end

@testitem "Meshes.Box 4D" tags=[:extended] setup=[Setup] begin
    a = π
    box = Box(Point(0, 0, 0, 0), Point(a, a, a, a))

    function f(p::P) where {P <: Meshes.Point}
        x1, x2, x3, x4 = ustrip.(to(p).coords)
        σ(x) = sqrt(a^2 - x^2)
        (σ(x1) + σ(x2) + σ(x3) + σ(x4)) * u"Ω/m^4"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 4a^3 * (π * a^2 / 4) * u"Ω"
    @test integral(f, box, GaussLegendre(100))≈sol rtol=1e-6
    @test_throws "not supported" integral(f, box, GaussKronrod())
    @test integral(f, box, HAdaptiveCubature(rtol = 1e-6))≈sol rtol=1e-6

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, box, GaussLegendre(100))≈vsol rtol=1e-6
    @test_throws "not supported" integral(fv, box, GaussKronrod())
    @test integral(fv, box, HAdaptiveCubature(rtol = 1e-6))≈vsol rtol=1e-6

    # Integral aliases
    @test_throws "not supported" lineintegral(f, box)
    @test_throws "not supported" surfaceintegral(f, box)
    @test_throws "not supported" volumeintegral(f, box)
end

@testitem "Meshes.Circle" setup=[Setup] begin
    center = Point(0, 3, 0)
    n̂ = Vec(1/2, 1/2, sqrt(2)/2)
    plane = Plane(center, n̂)
    radius = 4.4
    circle = Circle(plane, radius)

    function f(p::P) where {P <: Meshes.Point}
        offset = p - center
        ur = hypot(offset.coords...)
        r = ustrip(u"m", ur)
        exp(-r^2)
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 2π * radius * exp(-radius^2) * u"m"
    @test integral(f, circle, GaussLegendre(100)) ≈ sol
    @test integral(f, circle, GaussKronrod()) ≈ sol
    @test integral(f, circle, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, circle, GaussLegendre(100)) ≈ vsol
    @test integral(fv, circle, GaussKronrod()) ≈ vsol
    @test integral(fv, circle, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, circle) ≈ sol
    @test_throws "not supported" surfaceintegral(f, circle)
    @test_throws "not supported" volumeintegral(f, circle)
end

@testitem "Meshes.Cone" setup=[Setup] begin
    r = 2.5u"m"
    h = 2.5u"m"
    origin = Point(0, 0, 0)
    xy_plane = Plane(origin, Vec(0, 0, 1))
    base = Disk(xy_plane, r)
    apex = Point(0.0u"m", 0.0u"m", h)
    cone = Cone(base, apex)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    _volume_cone_rightcircular(h, r) = π * r^2 * h / 3

    # Scalar integrand
    sol = _volume_cone_rightcircular(r, h)
    @test integral(f, cone, GaussLegendre(100)) ≈ sol
    @test_throws "not supported" integral(f, cone, GaussKronrod())
    @test integral(f, cone, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, cone, GaussLegendre(100)) ≈ vsol
    @test_throws "not supported" integral(fv, cone, GaussKronrod())
    @test integral(fv, cone, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, cone)
    @test_throws "not supported" surfaceintegral(f, cone)
    @test volumeintegral(f, cone) ≈ sol
end

@testitem "Meshes.ConeSurface" setup=[Setup] begin
    r = 2.5u"m"
    h = 2.5u"m"
    origin = Point(0, 0, 0)
    xy_plane = Plane(origin, Vec(0, 0, 1))
    base = Disk(xy_plane, r)
    apex = Point(0.0u"m", 0.0u"m", h)
    cone = ConeSurface(base, apex)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    _area_cone_rightcircular(h, r) = (π * r^2) + (π * r * hypot(h, r))

    # Scalar integrand
    sol = _area_cone_rightcircular(h, r)
    @test integral(f, cone, GaussLegendre(100))≈sol rtol=1e-6
    @test integral(f, cone, GaussKronrod())≈sol rtol=1e-6
    @test integral(f, cone, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, cone, GaussLegendre(100))≈vsol rtol=1e-6
    @test integral(fv, cone, GaussKronrod())≈vsol rtol=1e-6
    @test integral(fv, cone, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, cone)
    @test surfaceintegral(f, cone) ≈ sol
    @test_throws "not supported" volumeintegral(f, cone)
end

@testitem "Meshes.Cylinder" setup=[Setup] begin
    pt_w = Point(-1, 0, 0)
    pt_e = Point(1, 0, 0)
    cyl = Cylinder(pt_e, pt_w, 2.5)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(cyl)
    @test integral(f, cyl, GaussLegendre(100)) ≈ sol
    @test_throws "not supported" integral(f, cyl, GaussKronrod())
    @test integral(f, cyl, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, cyl, GaussLegendre(100)) ≈ vsol
    @test_throws "not supported" integral(fv, cyl, GaussKronrod())
    @test integral(fv, cyl, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, cyl)
    @test_throws "not supported" surfaceintegral(f, cyl)
    @test volumeintegral(f, cyl) ≈ sol
end

@testitem "Meshes.CylinderSurface" setup=[Setup] begin
    pt_w = Point(-1, 0, 0)
    pt_e = Point(1, 0, 0)
    cyl = CylinderSurface(pt_e, pt_w, 2.5)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(cyl)
    @test integral(f, cyl, GaussLegendre(100)) ≈ sol
    @test integral(f, cyl, GaussKronrod()) ≈ sol
    @test integral(f, cyl, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, cyl, GaussLegendre(100)) ≈ vsol
    @test integral(fv, cyl, GaussKronrod()) ≈ vsol
    @test integral(fv, cyl, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, cyl)
    @test surfaceintegral(f, cyl) ≈ sol
    @test_throws "not supported" volumeintegral(f, cyl)
end

@testitem "Meshes.Disk" setup=[Setup] begin
    center = Point(1, 2, 3)
    n̂ = Vec(1/2, 1/2, sqrt(2)/2)
    plane = Plane(center, n̂)
    radius = 2.5
    disk = Disk(plane, radius)

    function f(p::P) where {P <: Meshes.Point}
        offset = p - center
        ur = hypot(offset.coords...)
        r = ustrip(u"m", ur)
        exp(-r^2)
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = (π - π * exp(-radius^2)) * u"m^2"
    @test integral(f, disk, GaussLegendre(100)) ≈ sol
    @test integral(f, disk, GaussKronrod()) ≈ sol
    @test integral(f, disk, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, disk, GaussLegendre(100)) ≈ vsol
    @test integral(fv, disk, GaussKronrod()) ≈ vsol
    @test integral(fv, disk, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, disk)
    @test surfaceintegral(f, disk) ≈ sol
    @test_throws "not supported" volumeintegral(f, disk)
end

@testitem "Meshes.Ellipsoid" setup=[Setup] begin
    origin = Point(0, 0, 0)
    radii = (1.0, 2.0, 0.5)
    ellipsoid = Ellipsoid(radii, origin)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Tolerances are higher due to `measure` being only an approximation
    # Scalar integrand
    sol = Meshes.measure(ellipsoid)
    @test integral(f, ellipsoid, GaussLegendre(100))≈sol rtol=1e-2
    @test integral(f, ellipsoid, GaussKronrod())≈sol rtol=1e-2
    @test integral(f, ellipsoid, HAdaptiveCubature())≈sol rtol=1e-2

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, ellipsoid, GaussLegendre(100))≈vsol rtol=1e-2
    @test integral(fv, ellipsoid, GaussKronrod())≈vsol rtol=1e-2
    @test integral(fv, ellipsoid, HAdaptiveCubature())≈vsol rtol=1e-2

    # Integral aliases
    @test_throws "not supported" lineintegral(f, ellipsoid)
    @test surfaceintegral(f, ellipsoid)≈sol rtol=1e-2
    @test_throws "not supported" volumeintegral(f, ellipsoid)
end

@testitem "Meshes.FrustumSurface" setup=[Setup] begin
    # Create a frustum whose radius halves at the top,
    # i.e. the bottom half of a cone by height
    r_bot = 2.5u"m"
    r_top = 1.25u"m"
    cone_h = 2π * u"m"
    origin = Point(0, 0, 0)
    z = Vec(0, 0, 1)
    plane_bot = Plane(origin, z)
    disk_bot = Disk(plane_bot, r_bot)
    center_top = Point(0.0u"m", 0.0u"m", 0.5cone_h)
    plane_top = Plane(center_top, z)
    disk_top = Disk(plane_top, r_top)
    frustum = FrustumSurface(disk_bot, disk_top)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    _area_base(r) = π * r^2
    _area_cone_walls(h, r) = π * r * hypot(h, r)

    # Scalar integrand
    sol = let
        area_walls_projected = _area_cone_walls(cone_h, r_bot)
        area_walls_missing = _area_cone_walls(0.5cone_h, r_top)
        area_walls = area_walls_projected - area_walls_missing
        area_walls + _area_base(r_top) + _area_base(r_bot)
    end
    @test integral(f, frustum, GaussLegendre(100))≈sol rtol=1e-6
    @test integral(f, frustum, GaussKronrod())≈sol rtol=1e-6
    @test integral(f, frustum, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, frustum, GaussLegendre(100))≈vsol rtol=1e-6
    @test integral(fv, frustum, GaussKronrod())≈vsol rtol=1e-6
    @test integral(fv, frustum, HAdaptiveCubature()) ≈ vsol
end

@testitem "Meshes.Hexahedron" setup=[Setup] begin
    hexahedron = Hexahedron(Point(0, 0, 0), Point(2, 0, 0), Point(2, 2, 0),
        Point(0, 2, 0), Point(0, 0, 2), Point(1, 0, 2), Point(1, 1, 2), Point(0, 1, 2))

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(hexahedron)
    @test integral(f, hexahedron, GaussLegendre(100)) ≈ sol
    @test_throws "not supported" integral(f, hexahedron, GaussKronrod())≈sol
    @test integral(f, hexahedron, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, hexahedron, GaussLegendre(100)) ≈ vsol
    @test_throws "not supported" integral(fv, hexahedron, GaussKronrod())≈vsol
    @test integral(fv, hexahedron, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, hexahedron)
    @test_throws "not supported" surfaceintegral(f, hexahedron)
    @test volumeintegral(f, hexahedron) ≈ sol
end

@testitem "Meshes.Line" setup=[Setup] begin
    a = Point(0.0u"m", 0.0u"m", 0.0u"m")
    b = Point(1.0u"m", 1.0u"m", 1.0u"m")
    line = Line(a, b)

    function f(p::P) where {P <: Meshes.Point}
        ur = hypot(p.coords.x, p.coords.y, p.coords.z)
        r = ustrip(u"m", ur)
        exp(-r^2)
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = sqrt(π) * u"m"
    @test integral(f, line, GaussLegendre(100)) ≈ sol
    @test integral(f, line, GaussKronrod()) ≈ sol
    @test integral(f, line, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, line, GaussLegendre(100)) ≈ vsol
    @test integral(fv, line, GaussKronrod()) ≈ vsol
    @test integral(fv, line, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, line) ≈ sol
    @test_throws "not supported" surfaceintegral(f, line)
    @test_throws "not supported" volumeintegral(f, line)
end

@testitem "Meshes.ParaboloidSurface" setup=[Setup] begin
    origin = Point(0, 0, 0)
    parab = ParaboloidSurface(origin, 2.5, 4.15)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(parab)
    @test integral(f, parab, GaussLegendre(100)) ≈ sol
    @test integral(f, parab, GaussKronrod()) ≈ sol
    @test integral(f, parab, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, parab, GaussLegendre(100)) ≈ vsol
    @test integral(fv, parab, GaussKronrod()) ≈ vsol
    @test integral(fv, parab, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, parab)
    @test surfaceintegral(f, parab) ≈ sol
    @test_throws "not supported" volumeintegral(f, parab)
end

@testitem "ParametrizedCurve" setup=[Setup] begin
    # ParametrizedCurve has been added in Meshes v0.51.20
    # If the version is specified as minimal compat bound in the Project.toml, the downgrade test fails
    if pkgversion(Meshes) >= v"0.51.20"
        using CoordRefSystems: Polar
        using LinearAlgebra: norm

        # Parameterize a circle centered on origin with specified radius
        radius = 4.4
        curve_cart = ParametrizedCurve(
            t -> Point(radius * cos(t), radius * sin(t)), (0.0, 2π))
        curve_polar = ParametrizedCurve(t -> Point(Polar(radius, t)), (0.0, 2π))

        function f(p::P) where {P <: Meshes.Point}
            ur = norm(to(p))
            r = ustrip(u"m", ur)
            exp(-r^2)
        end
        fv(p) = fill(f(p), 3)

        # Scalar integrand
        sol = 2π * radius * exp(-radius^2) * u"m"
        @test integral(f, curve_cart, GaussLegendre(100)) ≈ sol
        @test integral(f, curve_cart, GaussKronrod()) ≈ sol
        @test integral(f, curve_cart, HAdaptiveCubature()) ≈ sol
        @test integral(f, curve_polar, GaussLegendre(100)) ≈ sol
        @test integral(f, curve_polar, GaussKronrod()) ≈ sol
        @test integral(f, curve_polar, HAdaptiveCubature()) ≈ sol

        # Vector integrand
        vsol = fill(sol, 3)
        @test integral(fv, curve_cart, GaussLegendre(100)) ≈ vsol
        @test integral(fv, curve_cart, GaussKronrod()) ≈ vsol
        @test integral(fv, curve_cart, HAdaptiveCubature()) ≈ vsol
        @test integral(fv, curve_polar, GaussLegendre(100)) ≈ vsol
        @test integral(fv, curve_polar, GaussKronrod()) ≈ vsol
        @test integral(fv, curve_polar, HAdaptiveCubature()) ≈ vsol

        # Integral aliases
        @test lineintegral(f, curve_cart) ≈ sol
        @test_throws "not supported" surfaceintegral(f, curve_cart)
        @test_throws "not supported" volumeintegral(f, curve_cart)
        @test lineintegral(f, curve_polar) ≈ sol
        @test_throws "not supported" surfaceintegral(f, curve_polar)
        @test_throws "not supported" volumeintegral(f, curve_polar)
    end
end

@testitem "Meshes.Plane" setup=[Setup] begin
    p = Point(0.0u"m", 0.0u"m", 0.0u"m")
    v = Vec(0.0u"m", 0.0u"m", 1.0u"m")
    plane = Plane(p, v)

    function f(p::P) where {P <: Meshes.Point}
        ur = hypot(p.coords.x, p.coords.y, p.coords.z)
        r = ustrip(u"m", ur)
        exp(-r^2)
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = π * u"m^2"
    @test integral(f, plane, GaussLegendre(100)) ≈ sol
    @test integral(f, plane, GaussKronrod()) ≈ sol
    @test integral(f, plane, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, plane, GaussLegendre(100)) ≈ vsol
    @test integral(fv, plane, GaussKronrod()) ≈ vsol
    @test integral(fv, plane, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, plane)
    @test surfaceintegral(f, plane) ≈ sol
    @test_throws "not supported" volumeintegral(f, plane)
end

@testitem "Meshes.Quadrangle" setup=[Setup] begin
    using SpecialFunctions: erf
    quadrangle = Quadrangle((-1.0, 0.0), (-1.0, 1.0), (1.0, 1.0), (1.0, 0.0))

    function f(p::P) where {P <: Meshes.Point}
        ur = hypot(p.coords.x, p.coords.y)
        r = ustrip(u"m", ur)
        exp(-r^2)
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 0.5 * π * erf(1)^2 * u"m^2"
    @test integral(f, quadrangle, GaussLegendre(100)) ≈ sol
    @test integral(f, quadrangle, GaussKronrod()) ≈ sol
    @test integral(f, quadrangle, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, quadrangle, GaussLegendre(100)) ≈ vsol
    @test integral(fv, quadrangle, GaussKronrod()) ≈ vsol
    @test integral(fv, quadrangle, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, quadrangle)
    @test surfaceintegral(f, quadrangle) ≈ sol
    @test_throws "not supported" volumeintegral(f, quadrangle)
end

@testitem "Meshes.Ray" setup=[Setup] begin
    a = Point(0.0u"m", 0.0u"m", 0.0u"m")
    v = Vec(1.0u"m", 1.0u"m", 1.0u"m")
    ray = Ray(a, v)

    function f(p::P) where {P <: Meshes.Point}
        ur = hypot(p.coords.x, p.coords.y, p.coords.z)
        r = ustrip(u"m", ur)
        exp(-r^2)
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = sqrt(π) / 2 * u"m"
    @test integral(f, ray, GaussLegendre(100)) ≈ sol
    @test integral(f, ray, GaussKronrod()) ≈ sol
    @test integral(f, ray, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, ray, GaussLegendre(100)) ≈ vsol
    @test integral(fv, ray, GaussKronrod()) ≈ vsol
    @test integral(fv, ray, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, ray) ≈ sol
    @test_throws "not supported" surfaceintegral(f, ray)
    @test_throws "not supported" volumeintegral(f, ray)
end

@testitem "Meshes.Ring" setup=[Setup] begin
    pt_a = Point(0.0u"m", 0.0u"m", 0.0u"m")
    pt_b = Point(1.0u"m", 0.0u"m", 0.0u"m")
    pt_c = Point(1.0u"m", 1.0u"m", 0.0u"m")
    pt_d = Point(1.0u"m", 1.0u"m", 1.0u"m")
    rope = Ring(pt_a, pt_b, pt_c, pt_d, pt_c, pt_b)

    function f(p::P) where {P <: Meshes.Point}
        x, y, z = (p.coords.x, p.coords.y, p.coords.z)
        (x + 2y + 3z) * u"A/m^2"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 14.0u"A"
    @test integral(f, rope, GaussLegendre(100)) ≈ sol
    @test integral(f, rope, GaussKronrod()) ≈ sol
    @test integral(f, rope, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, rope, GaussLegendre(100)) ≈ vsol
    @test integral(fv, rope, GaussKronrod()) ≈ vsol
    @test integral(fv, rope, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, rope) ≈ sol
    @test_throws "not supported" surfaceintegral(f, rope)
    @test_throws "not supported" volumeintegral(f, rope)
end

@testitem "Meshes.Rope" setup=[Setup] begin
    pt_a = Point(0.0u"m", 0.0u"m", 0.0u"m")
    pt_b = Point(1.0u"m", 0.0u"m", 0.0u"m")
    pt_c = Point(1.0u"m", 1.0u"m", 0.0u"m")
    pt_d = Point(1.0u"m", 1.0u"m", 1.0u"m")
    rope = Rope(pt_a, pt_b, pt_c, pt_d)

    function f(p::P) where {P <: Meshes.Point}
        x, y, z = (p.coords.x, p.coords.y, p.coords.z)
        (x + 2y + 3z) * u"A/m^2"
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 7.0u"A"
    @test integral(f, rope, GaussLegendre(100)) ≈ sol
    @test integral(f, rope, GaussKronrod()) ≈ sol
    @test integral(f, rope, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, rope, GaussLegendre(100)) ≈ vsol
    @test integral(fv, rope, GaussKronrod()) ≈ vsol
    @test integral(fv, rope, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, rope) ≈ sol
    @test_throws "not supported" surfaceintegral(f, rope)
    @test_throws "not supported" volumeintegral(f, rope)
end

@testitem "Meshes.Segment" setup=[Setup] begin
    # Connect a line segment from the origin to an arbitrary point on the unit sphere
    φ, θ = (7pi / 6, pi / 3)  # Arbitrary spherical angles
    pt_a = Point(0.0u"m", 0.0u"m", 0.0u"m")
    pt_b = Point(sin(θ) * cos(φ) * u"m", sin(θ) * sin(φ) * u"m", cos(θ) * u"m")
    segment = Segment(pt_a, pt_b)

    a, b = (7.1, 4.6)  # arbitrary constants > 0

    function f(p::P) where {P <: Meshes.Point}
        ur = hypot(p.coords.x, p.coords.y, p.coords.z)
        r = ustrip(u"m", ur)
        exp(r * log(a) + (1 - r) * log(b))
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = (a - b) / (log(a) - log(b)) * u"m"
    @test integral(f, segment, GaussLegendre(100)) ≈ sol
    @test integral(f, segment, GaussKronrod()) ≈ sol
    @test integral(f, segment, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, segment, GaussLegendre(100)) ≈ vsol
    @test integral(fv, segment, GaussKronrod()) ≈ vsol
    @test integral(fv, segment, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, segment) ≈ sol
    @test_throws "not supported" surfaceintegral(f, segment)
    @test_throws "not supported" volumeintegral(f, segment)
end

@testitem "Meshes.Sphere 2D" setup=[Setup] begin
    origin = Point(0, 0)
    radius = 4.4
    sphere = Sphere(origin, radius)

    function f(p::P) where {P <: Meshes.Point}
        ur = hypot(p.coords.x, p.coords.y)
        r = ustrip(u"m", ur)
        exp(-r^2)
    end
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = 2π * radius * exp(-radius^2) * u"m"
    @test integral(f, sphere, GaussLegendre(100)) ≈ sol
    @test integral(f, sphere, GaussKronrod()) ≈ sol
    @test integral(f, sphere, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, sphere, GaussLegendre(100)) ≈ vsol
    @test integral(fv, sphere, GaussKronrod()) ≈ vsol
    @test integral(fv, sphere, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test lineintegral(f, sphere) ≈ sol
    @test_throws "not supported" surfaceintegral(f, sphere)
    @test_throws "not supported" volumeintegral(f, sphere)
end

@testitem "Meshes.Sphere 3D" setup=[Setup] begin
    origin = Point(0, 0, 0)
    sphere = Sphere(origin, 4.4)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(sphere)
    @test integral(f, sphere, GaussLegendre(100)) ≈ sol
    @test integral(f, sphere, GaussKronrod()) ≈ sol
    @test integral(f, sphere, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, sphere, GaussLegendre(100)) ≈ vsol
    @test integral(fv, sphere, GaussKronrod()) ≈ vsol
    @test integral(fv, sphere, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, sphere)
    @test surfaceintegral(f, sphere) ≈ sol
    @test_throws "not supported" volumeintegral(f, sphere)
end

@testitem "Meshes.Tetrahedron" tags=[:extended] setup=[Setup] begin
    pt_n = Point(0, 1, 0)
    pt_w = Point(-1, 0, 0)
    pt_e = Point(1, 0, 0)
    ẑ = Vec(0, 0, 1)
    tetrahedron = Tetrahedron(pt_n, pt_w, pt_e, pt_n + ẑ)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(tetrahedron)
    @test_throws "not supported" integral(f, tetrahedron, GaussLegendre(100))
    @test integral(f, tetrahedron, GaussKronrod()) ≈ sol
    @test_throws "not supported" integral(f, tetrahedron, HAdaptiveCubature())

    # Vector integrand
    vsol = fill(sol, 3)
    @test_throws "not supported" integral(fv, tetrahedron, GaussLegendre(100))≈vsol
    @test integral(fv, tetrahedron, GaussKronrod()) ≈ vsol
    @test_throws "not supported" integral(fv, tetrahedron, HAdaptiveCubature())≈vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, tetrahedron)
    @test_throws "not supported" surfaceintegral(f, tetrahedron)
    @test volumeintegral(f, tetrahedron, GaussKronrod()) ≈ sol
end

@testitem "Meshes.Torus" setup=[Setup] begin
    origin = Point(0, 0, 0)
    ẑ = Vec(0, 0, 1)
    torus = Torus(origin, ẑ, 3.5, 1.25)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(torus)
    @test integral(f, torus, GaussLegendre(100)) ≈ sol
    @test integral(f, torus, GaussKronrod()) ≈ sol
    @test integral(f, torus, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, torus, GaussLegendre(100)) ≈ vsol
    @test integral(fv, torus, GaussKronrod()) ≈ vsol
    @test integral(fv, torus, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, torus)
    @test surfaceintegral(f, torus) ≈ sol
    @test_throws "not supported" volumeintegral(f, torus)
end

@testitem "Meshes.Triangle" setup=[Setup] begin
    pt_n = Point(0, 1, 0)
    pt_w = Point(-1, 0, 0)
    pt_e = Point(1, 0, 0)
    triangle = Triangle(pt_e, pt_n, pt_w)

    f(p) = 1.0
    fv(p) = fill(f(p), 3)

    # Scalar integrand
    sol = Meshes.measure(triangle)
    @test integral(f, triangle, GaussLegendre(100)) ≈ sol
    @test integral(f, triangle, GaussKronrod()) ≈ sol
    @test integral(f, triangle, HAdaptiveCubature()) ≈ sol

    # Vector integrand
    vsol = fill(sol, 3)
    @test integral(fv, triangle, GaussLegendre(100)) ≈ vsol
    @test integral(fv, triangle, GaussKronrod()) ≈ vsol
    @test integral(fv, triangle, HAdaptiveCubature()) ≈ vsol

    # Integral aliases
    @test_throws "not supported" lineintegral(f, triangle)
    @test surfaceintegral(f, triangle) ≈ sol
    @test_throws "not supported" volumeintegral(f, triangle)
end
