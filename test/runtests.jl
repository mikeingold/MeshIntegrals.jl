using MeshIntegrals
using Meshes
using Test

################################################################################
#                                Infrastructure
################################################################################

struct SupportItem{Dim, T, G<:Meshes.Geometry{Dim,T}}
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
    @assert item.type == coordtype(item.geometry) "Item type mismatch"

    algorithm_set = [
        (GaussLegendre(100),  item.gausslegendre),
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
    
    # Test Geometries
    ball2d(T)   = Ball{2,T}(origin2d(T), T(2.0))
    ball3d(T)   = Ball{3,T}(origin3d(T), T(2.0))
    bezier(T)   = BezierCurve([Point(cos(t), sin(t), 0) for t in range(T(0), T(2π), length=361)])
    box1d(T)    = Box{1,T}(Point(T(-1)), Point(T(1)))
    box2d(T)    = Box{2,T}(Point(T(-1), T(-1)), Point(T(1), T(1)))
    box3d(T)    = Box{3,T}(Point(T(-1), T(-1), T(-1)), Point(T(1), T(1), T(-1)))
    circle(T)   = Circle(plane_xy(T), T(2.5))
    cyl(T)      = Cylinder(pt_e(T), pt_w(T), T(2.5))
    cylsurf(T)  = CylinderSurface(pt_e(T), pt_w(T), T(2.5))
    disk(T)     = Disk(plane_xy(T), T(2.5))
    parab(T)    = ParaboloidSurface(origin3d(T), T(2.5), T(4.15))
    ring(T)     = Ring(pt_e(T), pt_n(T), pt_w(T), pt_s(T))
    rope(T)     = Rope(pt_e(T), pt_n(T), pt_w(T), pt_s(T), pt_e(T))
    segment(T)   = Segment(pt_e(T), pt_n(T))
    sphere2d(T) = Sphere(origin2d(T), T(2.5))
    sphere3d(T) = Sphere(origin3d(T), T(2.5))
    triangle(T) = Ngon(pt_e(T), pt_n(T), pt_w(T))
    torus(T)    = Torus(origin3d(T), ẑ(T), T(3.5), T(1.25))

    SUPPORT_MATRIX(T) = [
    # Name, T type, example,    integral,line,surface,volume,    GaussLegendre,GaussKronrod,HAdaptiveCubature
        SupportItem("Ball{2,T}", T, ball2d(T),          1, 0, 1, 0,   1, 1, 1),
        SupportItem("Ball{3,T}", T, ball3d(T),          1, 0, 0, 1,   1, 0, 1),
        # Ball{Dim,T}
        SupportItem("BezierCurve", T, bezier(T),        1, 1, 0, 0,   1, 1, 1),
        SupportItem("Box{1,T}", T, box1d(T),            1, 1, 0, 0,   1, 1, 1),
        SupportItem("Box{2,T}", T, box2d(T),            1, 0, 1, 0,   1, 1, 1),
        SupportItem("Box{3,T}", T, box3d(T),            1, 0, 0, 1,   1, 0, 1),
        # Box{Dim,T}
        SupportItem("Circle", T, circle(T),             1, 1, 0, 0,   1, 1, 1),
        # Cone
        # ConeSurface
        SupportItem("Cylinder", T, cyl(T),              1, 0, 0, 1,   1, 0, 1),
        SupportItem("CylinderSurface", T, cylsurf(T),   1, 0, 1, 0,   0, 1, 0),
        SupportItem("Disk", T, disk(T),                 1, 0, 1, 0,   1, 1, 1),
        # Frustum
        # FrustumSurface
        # Line -- custom tests below
        SupportItem("ParaboloidSurface", T, parab(T),   1, 0, 1, 0,   1, 1, 1),
        # Plane -- custom tests below
        # Ray
        SupportItem("Ring", T, ring(T),                 1, 1, 0, 0,   1, 1, 1),
        SupportItem("Rope", T, rope(T),                 1, 1, 0, 0,   1, 1, 1),
        SupportItem("Segment", T, segment(T),           1, 1, 0, 0,   1, 1, 1),
        # SimpleMesh
        SupportItem("Sphere{2,T}", T, sphere2d(T),      1, 1, 0, 0,   1, 1, 1),
        SupportItem("Sphere{3,T}", T, sphere3d(T),      1, 0, 1, 0,   1, 1, 1),
        SupportItem("Triangle", T, triangle(T),         1, 0, 1, 0,   1, 1, 1),
        SupportItem("Torus", T, torus(T),               1, 0, 1, 0,   1, 1, 1),
    ]

    # Run all integral tests
    #map(autotest, SUPPORT_MATRIX(Float32))  TODO re-enable
    map(autotest, SUPPORT_MATRIX(Float64))

    # Custom tests for Line (no measure available for reference)
    @testset "Meshes.Line" begin
        line = Line(pt_e, pt_w)

        function f(p::Point{3,T}) where {T}
            x, y, z = p.coords
            exp(-x^2)
        end
        fv(p) = fill(f(p),3)

        @test integral(f, line, GaussLegendre(100)) ≈ sqrt(π)
        @test integral(f, line, GaussKronrod()) ≈ sqrt(π)
        @test integral(f, line, HAdaptiveCubature()) ≈ sqrt(π)

        @test integral(fv, line, GaussLegendre(100)) ≈ fill(sqrt(π),3)
        @test integral(fv, line, GaussKronrod()) ≈ fill(sqrt(π),3)
        @test integral(fv, line, HAdaptiveCubature()) ≈ fill(sqrt(π),3)
    end

    # Custom tests for Plane (no measure available for reference)
    @testset "Meshes.Plane" begin
        plane = Plane(origin3d, ẑ)

        function f(p::Point{3,T}) where {T}
            x, y, z = p.coords
            exp(-x^2 - y^2)
        end
        fv(p) = fill(f(p),3)

        @test integral(f, plane, GaussLegendre(100)) ≈ π
        @test integral(f, plane, GaussKronrod()) ≈ π
        @test integral(f, plane, HAdaptiveCubature()) ≈ π

        @test integral(fv, plane, GaussLegendre(100)) ≈ fill(π,3)
        @test integral(fv, plane, GaussKronrod()) ≈ fill(π,3)
        @test integral(fv, plane, HAdaptiveCubature()) ≈ fill(π,3)
    end
end
