using Aqua
using MeshIntegrals
using Meshes
using Test
using Unitful

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
    conesurf(T) = ConeSurface(disk(T), pt_z(T))
    cyl(T)      = Cylinder(pt_e(T), pt_w(T), T(2.5))
    cylsurf(T)  = CylinderSurface(pt_e(T), pt_w(T), T(2.5))
    disk(T)     = Disk(plane_xy(T), T(2.5))
    frusurf(T)  = FrustumSurface(Disk(plane_xy(T),T(1.2)), Disk(Plane(Point(T(0),T(0),T(π)),ẑ(T)),T(2.5)))
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
        # Ball{Dim,T}
        SupportItem("BezierCurve{$T}", T, bezier(T),        1, 1, 0, 0,   1, 1, 1),
        SupportItem("Box{1,$T}", T, box1d(T),               1, 1, 0, 0,   1, 1, 1),
        SupportItem("Box{2,$T}", T, box2d(T),               1, 0, 1, 0,   1, 1, 1),
        SupportItem("Box{3,$T}", T, box3d(T),               1, 0, 0, 1,   1, 0, 1),
        # Box{Dim,T}
        SupportItem("Circle{$T}", T, circle(T),             1, 1, 0, 0,   1, 1, 1),
        # Cone -- custom tests below
        SupportItem("ConeSurface{$T}", T, conesurf(T),      1, 0, 1, 0,   1, 1, 1),
        SupportItem("Cylinder{$T}", T, cyl(T),              1, 0, 0, 1,   1, 0, 1),
        SupportItem("CylinderSurface{$T}", T, cylsurf(T),   1, 0, 1, 0,   1, 1, 1),
        SupportItem("Disk{$T}", T, disk(T),                 1, 0, 1, 0,   1, 1, 1),
        # Frustum
        SupportItem("FrustumSurface{$T}", T, frusurf(T),    1, 0, 1, 0,   1, 1, 1),
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

    @testset "Float32 Geometries" begin
        # TODO temp disabled, see Issue #33
        #map(autotest, SUPPORT_MATRIX(Float64))
    end

    # Custom tests for Line (no measure available for reference)
    @testset "Meshes.Line" begin
        line = Line(pt_e(Float64), pt_w(Float64))

        function f(p::P) where {P<:Meshes.Point}
            x = ustrip(u"m", p.coords.x)
            y = ustrip(u"m", p.coords.y)
            z = ustrip(u"m", p.coords.z)
            exp(-x^2)
        end
        fv(p) = fill(f(p),3)

        @test integral(f, line, GaussLegendre(100)) ≈ sqrt(π)*u"m"
        @test integral(f, line, GaussKronrod()) ≈ sqrt(π)*u"m"
        @test integral(f, line, HAdaptiveCubature()) ≈ sqrt(π)*u"m"

        @test integral(fv, line, GaussLegendre(100)) ≈ fill(sqrt(π)*u"m",3)
        @test integral(fv, line, GaussKronrod()) ≈ fill(sqrt(π)*u"m",3)
        @test integral(fv, line, HAdaptiveCubature()) ≈ fill(sqrt(π)*u"m",3)
    end

    # Custom tests for Ray (no measure available for reference)
    @testset "Meshes.Ray" begin
        ray = Ray(origin3d(Float64), ẑ(Float64))

        function f(p::P) where {P<:Meshes.Point}
            x = ustrip(u"m", p.coords.x)
            y = ustrip(u"m", p.coords.y)
            z = ustrip(u"m", p.coords.z)
            2 * exp(-z^2)
        end
        fv(p) = fill(f(p),3)

        @test integral(f, ray, GaussLegendre(100)) ≈ sqrt(π)*u"m"
        @test integral(f, ray, GaussKronrod()) ≈ sqrt(π)*u"m"
        @test integral(f, ray, HAdaptiveCubature()) ≈ sqrt(π)*u"m"

        @test integral(fv, ray, GaussLegendre(100)) ≈ fill(sqrt(π)*u"m",3)
        @test integral(fv, ray, GaussKronrod()) ≈ fill(sqrt(π)*u"m",3)
        @test integral(fv, ray, HAdaptiveCubature()) ≈ fill(sqrt(π)*u"m",3)
    end

    # Custom tests for Plane (no measure available for reference)
    @testset "Meshes.Plane" begin
        plane = Plane(origin3d(Float64), ẑ(Float64))

        function f(p::P) where {P<:Meshes.Point}
            x = ustrip(u"m", p.coords.x)
            y = ustrip(u"m", p.coords.y)
            z = ustrip(u"m", p.coords.z)
            exp(-x^2 - y^2)
        end
        fv(p) = fill(f(p),3)

        @test integral(f, plane, GaussLegendre(100)) ≈ π*u"m^2"
        @test integral(f, plane, GaussKronrod()) ≈ π*u"m^2"
        @test integral(f, plane, HAdaptiveCubature()) ≈ π*u"m^2"

        @test integral(fv, plane, GaussLegendre(100)) ≈ fill(π*u"m^2",3)
        @test integral(fv, plane, GaussKronrod()) ≈ fill(π*u"m^2",3)
        @test integral(fv, plane, HAdaptiveCubature()) ≈ fill(π*u"m^2",3)
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
        @test integral(f, cone, GaussKronrod()) ≈ cone_volume
        @test integral(f, cone, HAdaptiveCubature()) ≈ cone_volume

        @test integral(fv, cone, GaussLegendre(100)) ≈ fill(cone_volume, 3)
        @test integral(fv, cone, GaussKronrod()) ≈ fill(cone_volume, 3)
        @test integral(fv, cone, HAdaptiveCubature()) ≈ fill(cone_volume, 3)
    end
end

################################################################################
#                                Aqua.jl Tests
################################################################################

@testset "Aqua.jl" begin
    # As of v0.11.4 -- Ambiguities check disabled since it fails due to upstream findings
    #   Verified that no ambiguities exist within MeshIntegrals.jl
    Aqua.test_all(MeshIntegrals; ambiguities=false)
end
