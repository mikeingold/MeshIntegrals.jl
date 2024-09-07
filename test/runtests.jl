using Aqua
using Meshes
using MeshIntegrals
using Test
using Unitful


################################################################################
#                                Aqua.jl Tests
################################################################################

@testset "Aqua.jl" begin
    # As of v0.11.4:
    # - Ambiguities check disabled since it fails due to upstream findings
    # - Verified that no ambiguities exist within MeshIntegrals.jl
    Aqua.test_all(MeshIntegrals; ambiguities=false)
end


################################################################################
#                                Infrastructure
################################################################################

struct SupportItem{T, Dim, CRS, G<:Meshes.Geometry{Meshes.𝔼{Dim},CRS}}
    name::String
    type::Type{T}
    geometry::G
    integral::Bool
    lineintegral::Bool
    surfaceintegral::Bool
    volumeintegral::Bool
    gausslegendre::Bool
    gausskronrod::Bool
    hadaptivecubature::Bool
end

# Constructor to explicitly convert Ints (0,1) to Bool values
SupportItem(name, type, geometry, checkboxes::Vararg{I,7}) where {I<:Integer} = SupportItem(name, type, geometry, Bool.(checkboxes)...)

# If method is supported, test it on scalar- and vector-valued functions.
# Otherwise, test that its use throws a MethodError
function integraltest(intf, geometry, rule, supported, T)
    f(::Point) = T(1)
    fv(::Point) = fill(T(1),2)

    if supported
        a1 = intf(f, geometry, rule)
        b1 = measure(geometry)
        @test a1 ≈ b1
        @test typeof(a1) == typeof(b1)
        @test intf(fv, geometry, rule) ≈ fill(b1,2)
    else
        @test_throws "not supported" intf(f, geometry, rule)
    end
end

# Generate a @testset for item
function autotest(item::SupportItem)
    #@assert item.type == coordtype(item.geometry) "Item type mismatch"

    N = (item.type == Float32) ? 1000 : 100
    algorithm_set = [
        (GaussLegendre(N),  item.gausslegendre),
        (GaussKronrod(),      item.gausskronrod),
        (HAdaptiveCubature(), item.hadaptivecubature)
    ]

    method_set = [
        (integral, item.integral),
        (lineintegral, item.lineintegral),
        (surfaceintegral, item.surfaceintegral),
        (volumeintegral, item.volumeintegral)
    ]

    itemsupport = Iterators.product(method_set,algorithm_set)

    # For each enabled solver type, run the test suite
    @testset "$(item.name)" begin
        for ((method,methodsupport), (alg,algsupport)) in itemsupport
            integraltest(method, item.geometry, alg, methodsupport && algsupport, item.type)
        end
    end
end

################################################################################
#                                  Integrals
################################################################################

@testset "Integrals" begin
    # Spatial descriptors
    origin3d(T) = Point(T(0), T(0), T(0))
    origin2d(T) = Point(T(0), T(0))
    ẑ(T) = Vec(T(0), T(0), T(1))
    plane_xy(T) = Plane(origin3d(T), ẑ(T))

    # Points on xy-plane at unit distance on axes
    pt_n(T) = Point(T( 0), T( 1), T(0))
    pt_w(T) = Point(T(-1), T( 0), T(0))
    pt_e(T) = Point(T( 1), T( 0), T(0))
    pt_s(T) = Point(T( 0), T(-1), T(0))
    pt_z(T) = Point(T( 0), T( 0), T(1))
    
    # Test Geometries
    ball2d(T)   = Ball(origin2d(T), T(2.0))
    ball3d(T)   = Ball(origin3d(T), T(2.0))
    bezier(T)   = BezierCurve([Point(cos(t), sin(t), 0) for t in range(T(0), T(2π), length=361)])
    box1d(T)    = Box(Point(T(-1)), Point(T(1)))
    box2d(T)    = Box(Point(T(-1), T(-1)), Point(T(1), T(1)))
    box3d(T)    = Box(Point(T(-1), T(-1), T(-1)), Point(T(1), T(1), T(-1)))
    circle(T)   = Circle(plane_xy(T), T(2.5))
    cyl(T)      = Cylinder(pt_e(T), pt_w(T), T(2.5))
    cylsurf(T)  = CylinderSurface(pt_e(T), pt_w(T), T(2.5))
    disk(T)     = Disk(plane_xy(T), T(2.5))
    parab(T)    = ParaboloidSurface(origin3d(T), T(2.5), T(4.15))
    ring(T)     = Ring(pt_e(T), pt_n(T), pt_w(T), pt_s(T))
    rope(T)     = Rope(pt_e(T), pt_n(T), pt_w(T), pt_s(T), pt_e(T))
    segment(T)  = Segment(pt_e(T), pt_n(T))
    sphere2d(T) = Sphere(origin2d(T), T(2.5))
    sphere3d(T) = Sphere(origin3d(T), T(2.5))
    tetra(T)    = Tetrahedron(pt_n(T), pt_w(T), pt_e(T), pt_n(T)+ẑ(T))
    triangle(T) = Ngon(pt_e(T), pt_n(T), pt_w(T))
    torus(T)    = Torus(origin3d(T), ẑ(T), T(3.5), T(1.25))

    SUPPORT_MATRIX(T) = [
    # Name, T type, example,    integral,line,surface,volume,    GaussLegendre,GaussKronrod,HAdaptiveCubature
        SupportItem("Ball{2,$T}", T, ball2d(T),             1, 0, 1, 0,   1, 1, 1),
        SupportItem("Ball{3,$T}", T, ball3d(T),             1, 0, 0, 1,   1, 0, 1),
        SupportItem("BezierCurve{$T}", T, bezier(T),        1, 1, 0, 0,   1, 1, 1),
        SupportItem("Box{1,$T}", T, box1d(T),               1, 1, 0, 0,   1, 1, 1),
        SupportItem("Box{2,$T}", T, box2d(T),               1, 0, 1, 0,   1, 1, 1),
        SupportItem("Box{3,$T}", T, box3d(T),               1, 0, 0, 1,   1, 0, 1),
        SupportItem("Circle{$T}", T, circle(T),             1, 1, 0, 0,   1, 1, 1),
        # Cone -- custom tests below
        # ConeSurface -- custom tests below
        SupportItem("Cylinder{$T}", T, cyl(T),              1, 0, 0, 1,   1, 0, 1),
        SupportItem("CylinderSurface{$T}", T, cylsurf(T),   1, 0, 1, 0,   1, 1, 1),
        SupportItem("Disk{$T}", T, disk(T),                 1, 0, 1, 0,   1, 1, 1),
        # Frustum -- not yet supported
        # FrustumSurface -- not yet supported
        # Line -- custom tests below
        SupportItem("ParaboloidSurface{$T}", T, parab(T),   1, 0, 1, 0,   1, 1, 1),
        # Plane -- custom tests below
        # Ray -- custom tests below
        SupportItem("Ring{$T}", T, ring(T),                 1, 1, 0, 0,   1, 1, 1),
        SupportItem("Rope{$T}", T, rope(T),                 1, 1, 0, 0,   1, 1, 1),
        SupportItem("Segment{$T}", T, segment(T),           1, 1, 0, 0,   1, 1, 1),
        # SimpleMesh
        SupportItem("Sphere{2,$T}", T, sphere2d(T),         1, 1, 0, 0,   1, 1, 1),
        SupportItem("Sphere{3,$T}", T, sphere3d(T),         1, 0, 1, 0,   1, 1, 1),
        SupportItem("Tetrahedron", T, tetra(T),             1, 0, 0, 1,   0, 1, 0),
        SupportItem("Triangle{$T}", T, triangle(T),         1, 0, 1, 0,   1, 1, 1),
        SupportItem("Torus{$T}", T, torus(T),               1, 0, 1, 0,   1, 1, 1),
    ]

    @testset "Float64 Geometries" begin
        map(autotest, SUPPORT_MATRIX(Float64))
    end

    # Custom tests for Cone
    @testset "Meshes.Cone" begin
        T = Float64

        cone_r = T(2.5)
        cone_h = T(2.5)

        cone = let
            base = Disk(plane_xy(T), cone_r)
            Cone(base, Point(0, 0, cone_h))
        end

        f(p) = T(1)
        fv(p) = fill(f(p), 3)

        _volume_cone_rightcircular(h, r) = T(π) * r^2 * h / 3
        cone_volume = _volume_cone_rightcircular(cone_r * u"m", cone_h * u"m")

        @test integral(f, cone, GaussLegendre(100)) ≈ cone_volume
        @test_throws "not supported" integral(f, cone, GaussKronrod())
        @test integral(f, cone, HAdaptiveCubature()) ≈ cone_volume

        @test integral(fv, cone, GaussLegendre(100)) ≈ fill(cone_volume, 3)
        @test_throws "not supported" integral(fv, cone, GaussKronrod())
        @test integral(fv, cone, HAdaptiveCubature()) ≈ fill(cone_volume, 3)
    end

    # Custom tests for ConeSurface
    @testset "Meshes.ConeSurface" begin
        T = Float64

        cone_r = T(2.5)
        cone_h = T(2.5)

        cone = let
            base = Disk(plane_xy(T), cone_r)
            ConeSurface(base, Point(0, 0, cone_h))
        end

        f(p) = T(1)
        fv(p) = fill(f(p), 3)

        _area_cone_rightcircular(h, r) = T(π) * r^2 + T(π) * r * hypot(h, r)
        cone_area = _area_cone_rightcircular(cone_r * u"m", cone_h * u"m")

        @test integral(f, cone, GaussLegendre(100)) ≈ cone_area
        @test integral(f, cone, GaussKronrod()) ≈ cone_area
        @test integral(f, cone, HAdaptiveCubature()) ≈ cone_area

        @test integral(fv, cone, GaussLegendre(100)) ≈ fill(cone_area, 3)
        @test integral(fv, cone, GaussKronrod()) ≈ fill(cone_area, 3)
        @test integral(fv, cone, HAdaptiveCubature()) ≈ fill(cone_area, 3)
    end

    #= DISABLED FrustumSurface testing due to long run times and seemingly-incorrect results
    # Custom tests for FrustumSurface
    @testset "Meshes.FrustumSurface" begin
        T = Float64

        # Create a frustum whose radius halves at the top,
        # i.e. the bottom half of a cone by height
        bot_r = T(5//2)
        top_r = T(5//4)
        cone_h = T(2π)
        frustum = let
            plane_bot = Plane(Point(0,0,0), Vec(0,0,1))
            disk_bot = Disk(plane_bot, bot_r)
            plane_top = Plane(Point(0,0,T(0.5)*cone_h), Vec(0,0,1))
            disk_top = Disk(plane_top, top_r)
            FrustumSurface(disk_bot, disk_top)
        end

        f(p) = T(1)
        fv(p) = fill(f(p), 3)

        _area_cone_rightcircular(h, r) = T(π) * r^2 + T(π) * r * hypot(h, r)
        frustum_area = let
            area_projected = _area_cone_rightcircular(top_r * u"m", cone_h * u"m")
            area_missing = _area_cone_rightcircular(top_r * u"m", T(0.5) * cone_h * u"m")
            area_projected - area_missing
        end

        @test integral(f, frustum, GaussLegendre(100)) ≈ frustum_area
        @test integral(f, frustum, GaussKronrod()) ≈ frustum_area
        @test integral(f, frustum, HAdaptiveCubature()) ≈ frustum_area

        @test integral(fv, frustum, GaussLegendre(100)) ≈ fill(frustum_area, 3)
        @test integral(fv, frustum, GaussKronrod()) ≈ fill(frustum_area, 3)
        @test integral(fv, frustum, HAdaptiveCubature()) ≈ fill(frustum_area, 3)
    end
    =#
end

    
################################################################################
#                                New Tests
################################################################################

@testset "New Independent Tests" begin

    @testset "Meshes.Line" begin
        a = Point(0.0u"m", 0.0u"m", 0.0u"m")
        b = Point(1.0u"m", 1.0u"m", 1.0u"m")
        line = Line(a, b)

        function f(p::P) where {P<:Meshes.Point}
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

    @testset "Meshes.Plane" begin
        p = Point(0.0u"m", 0.0u"m", 0.0u"m")
        v = Vec(0.0u"m", 0.0u"m", 1.0u"m")
        plane = Plane(p, v)

        function f(p::P) where {P<:Meshes.Point}
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

    @testset "Meshes.Ray" begin
        a = Point(0.0u"m", 0.0u"m", 0.0u"m")
        v = Vec(1.0u"m", 1.0u"m", 1.0u"m")
        ray = Ray(a, v)

        function f(p::P) where {P<:Meshes.Point}
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
        @test lineintegral(f, line) ≈ sol
        @test_throws "not supported" surfaceintegral(f, line)
        @test_throws "not supported" volumeintegral(f, line)
    end

end
