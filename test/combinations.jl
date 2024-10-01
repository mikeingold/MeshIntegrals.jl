# This section tests for:
# - All supported combinations of integral(f, ::Geometry, ::IntegrationAlgorithm) produce accurate results
# - Invalid applications of integral aliases (e.g. lineintegral) produce a descriptive error

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
    sol = π * a^3 / 2 * u"Ω"
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
end

@testitem "Meshes.Box 4D" setup=[Setup] begin
    box = Box(Point(zeros(4)...), Point(ones(4)...))

    f = p -> one(FP)

    # Test for currently-unsupported >3D differentials
    @test integral(f, box)≈1.0u"m^4" broken=true
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
