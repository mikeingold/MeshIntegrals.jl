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
#                          Automatic test generation
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
    

@testset verbose=true showtiming=true "Integrals" begin
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
        SupportItem("Cylinder{$T}", T, cyl(T),              1, 0, 0, 1,   1, 0, 1),
        SupportItem("CylinderSurface{$T}", T, cylsurf(T),   1, 0, 1, 0,   1, 1, 1),
        SupportItem("Disk{$T}", T, disk(T),                 1, 0, 1, 0,   1, 1, 1),
        # Frustum -- not yet supported
        SupportItem("ParaboloidSurface{$T}", T, parab(T),   1, 0, 1, 0,   1, 1, 1),
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
end

    
################################################################################
#                                New Tests
################################################################################

@testset verbose=true showtiming=true "Function-Geometry-Algorithm Combinations" begin
# This section tests for:
# - All supported combinations of integral(f, ::Geometry, ::IntegrationAlgorithm) produce accurate results
# - Invalid applications of integral aliases (e.g. lineintegral) produce a descriptive error

    @testset "Meshes.Cone" begin
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

    @testset "Meshes.ConeSurface" begin
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
        @test integral(f, cone, GaussLegendre(100)) ≈ sol
        @test integral(f, cone, GaussKronrod()) ≈ sol
        @test integral(f, cone, HAdaptiveCubature()) ≈ sol

        # Vector integrand
        vsol = fill(sol, 3)
        @test integral(fv, cone, GaussLegendre(100)) ≈ vsol
        @test integral(fv, cone, GaussKronrod()) ≈ vsol
        @test integral(fv, cone, HAdaptiveCubature()) ≈ vsol

        # Integral aliases
        @test_throws "not supported" lineintegral(f, cone)
        @test surfaceintegral(f, cone) ≈ sol
        @test_throws "not supported" volumeintegral(f, cone)
    end

    @testset "Meshes.FrustumSurface" begin
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
        @test integral(f, frustum, GaussLegendre(100)) ≈ sol
        @test integral(f, frustum, GaussKronrod()) ≈ sol
        @test integral(f, frustum, HAdaptiveCubature()) ≈ sol

        # Vector integrand
        vsol = fill(sol, 3)
        @test integral(fv, frustum, GaussLegendre(100)) ≈ vsol
        @test integral(fv, frustum, GaussKronrod()) ≈ vsol
        @test integral(fv, frustum, HAdaptiveCubature()) ≈ vsol
    end

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
        @test lineintegral(f, ray) ≈ sol
        @test_throws "not supported" surfaceintegral(f, ray)
        @test_throws "not supported" volumeintegral(f, ray)
    end

    @testset "Meshes.Ring" begin
        pt_a = Point(0.0u"m", 0.0u"m", 0.0u"m")
        pt_b = Point(1.0u"m", 0.0u"m", 0.0u"m")
        pt_c = Point(1.0u"m", 1.0u"m", 0.0u"m")
        pt_d = Point(1.0u"m", 1.0u"m", 1.0u"m")
        rope = Ring(pt_a, pt_b, pt_c, pt_d, pt_c, pt_b)

        function f(p::P) where {P<:Meshes.Point}
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

    @testset "Meshes.Rope" begin
        pt_a = Point(0.0u"m", 0.0u"m", 0.0u"m")
        pt_b = Point(1.0u"m", 0.0u"m", 0.0u"m")
        pt_c = Point(1.0u"m", 1.0u"m", 0.0u"m")
        pt_d = Point(1.0u"m", 1.0u"m", 1.0u"m")
        rope = Rope(pt_a, pt_b, pt_c, pt_d)

        function f(p::P) where {P<:Meshes.Point}
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

    @testset "Meshes.Segment" begin
        # Connect a line segment from the origin to an arbitrary point on the unit sphere
        phi, theta = (7pi/6, pi/3)  # Arbitrary spherical angles
        pt_a = Point(0.0u"m", 0.0u"m", 0.0u"m")
        pt_b = Point(sin(theta)*cos(phi)*u"m", sin(theta)*sin(phi)*u"m", cos(theta)*u"m")
        segment = Segment(pt_a, pt_b)

        a, b = (7.1, 4.6)  # arbitrary constants > 0

        function f(p::P) where {P<:Meshes.Point}
            ur = hypot(p.coords.x, p.coords.y, p.coords.z)
            r = ustrip(u"m", ur)
            exp(r*log(a) + (1-r)*log(b))
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

end

@testset verbose=true showtiming=true "Alternate Floating Point Types" begin
# For integral(f, geometry, settings, FP) where FP is not Float64, ensure results
# have expected level of accuracy and are produce results in appropriate type
    
    @testset "$FP" for FP in (Float32,)    # TODO BigFloat tests
        # Rectangular volume with unit integrand
        f = p -> one(FP)
        box1d = Box(Point(fill(zero(FP)*u"m", 1)...), Point(fill(one(FP)*u"m", 1)...))
        box2d = Box(Point(fill(zero(FP)*u"m", 2)...), Point(fill(one(FP)*u"m", 2)...))
        box3d = Box(Point(fill(zero(FP)*u"m", 3)...), Point(fill(one(FP)*u"m", 3)...))
        box4d = Box(Point(fill(zero(FP)*u"m", 4)...), Point(fill(one(FP)*u"m", 4)...))

        # Check HCubature integrals (same method invoked for all dimensions)
        int_HC = integral(f, box1d, HAdaptiveCubature(), FP)
        @test int_HC ≈ 1.0u"m"    atol=0.01u"m"
        @test typeof(int_HC.val) == FP    broken=true

        # Check Gauss-Legendre integral in 1D
        int_GL_1D = integral(f, box1d, GaussLegendre(100), FP)
        @test int_GL_1D ≈ 1.0u"m"     atol=0.01u"m"
        @test typeof(int_GL_1D.val) == FP    broken=true

        # Check Gauss-Legendre integral in 2D
        int_GL_2D = integral(f, box2d, GaussLegendre(100), FP)
        @test int_GL_2D ≈ 1.0u"m^2"   atol=0.02u"m^2"
        @test typeof(int_GL_2D.val) == FP    broken=true

        # Check Gauss-Legendre integral in 3D
        int_GL_3D = integral(f, box3d, GaussLegendre(100), FP)
        @test int_GL_3D ≈ 1.0u"m^3"   atol=0.03u"m^3"
        @test typeof(int_GL_3D.val) == FP    broken=true
    end

    @testset "Integral Aliases" begin
        f = p -> one(Float32)
        box4d = Box(Point(fill(0.0f0u"m", 4)...), Point(fill(1.0f0u"m", 4)...))
            
        # Check alias functions for accuracy
        @test lineintegral(f, box1d, GaussLegendre(100), FP)    ≈ 1.0f0u"m"     atol=0.01f0u"m"
        @test surfaceintegral(f, box2d, GaussLegendre(100), FP) ≈ 1.0f0u"m^2"   atol=0.02f0u"m^2"
        @test volumeintegral(f, box3d, GaussLegendre(100), FP)  ≈ 1.0f0u"m^3"   atol=0.03f0u"m^3"

        # Check for unsupported use of alias functions
        @test_throws "not supported" lineintegral(f, box4d, HAdaptiveCubature(), FP)
        @test_throws "not supported" surfaceintegral(f, box4d, HAdaptiveCubature(), FP)
        @test_throws "not supported" volumeintegral(f, box4d, HAdaptiveCubature(), FP)
    end

end
