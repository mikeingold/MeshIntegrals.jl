using MeshIntegrals
using Meshes
using Test

#using DynamicQuantities
#using Unitful

################################################################################
#                                Infrastructure
################################################################################

struct SupportItem
    name::String
    geometry
    integral::Bool
    lineintegral::Bool
    surfaceintegral::Bool
    volumeintegral::Bool
    gausslegendre::Bool
    gausskronrod::Bool
    hadaptivecubature::Bool
end

# If method is supported, test it on scalar- and vector-valued functions.
# Otherwise, test that its use throws a MethodError
function integraltest(intf, geometry, rule, supported)
    f(::Point) = 1.0
    fv(::Point) = fill(1.0,3)

    if supported
        @test intf(f, geometry, rule) ≈ measure(geometry)
        @test intf(fv, geometry, rule) ≈ fill(measure(geometry),3)
    else
        @test_throws "not supported" intf(f, geometry, rule)
    end
end

# Generate a @testset for item
function autotest(item::SupportItem)
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
            integraltest(method, item.geometry, alg, methodsupport && algsupport)
        end
    end
end

################################################################################
#                                  Integrals
################################################################################

@testset "Integrals" begin
    # Spatial descriptors
    origin3d = Point(0.0, 0.0, 0.0)
    origin2d = Point(0.0, 0.0)
    ẑ = Vec(0.0, 0.0, 1.0)
    plane_xy = Plane(origin3d,ẑ)

    # Points on unit circle at axes
                   pt_n = Point( 0.0,  1.0, 0.0)
    pt_w = Point(-1.0,  0.0, 0.0);    pt_e = Point( 1.0,  0.0, 0.0)
                   pt_s = Point( 0.0, -1.0, 0.0)  
    
    # Test Geometries
    ball2d = Ball(origin2d, 2.0)
    ball3d = Ball(origin3d, 2.0)
    bezier = BezierCurve([Point(cos(t), sin(t), 0.0) for t in range(0, 2pi, length=361)])
    box1d = Box(Point(-1.0), Point(1.0))
    box2d = Box(Point(-1.0,-1.0), Point(1.0,1.0))
    box3d = Box(Point(-1.0,-1.0,-1.0), Point(1.0,1.0,1.0))
    circle = Circle(plane_xy, 2.5)
    disk = Disk(plane_xy, 2.5)
    line_ne = Line(pt_e, pt_n)
    ring_rect = Ring(pt_e, pt_n, pt_w, pt_s)
    rope_rect = Rope(pt_e, pt_n, pt_w, pt_s, pt_e)
    seg_ne = Segment(pt_e, pt_n)
    sphere2d = Sphere(origin2d, 2.5)
    sphere3d = Sphere(origin3d, 2.5)
    triangle = Ngon(pt_e, pt_n, pt_w)
    cylsurf = CylinderSurface(pt_e, pt_w, 2.5)   # TODO modify to a non-right-CylinderSurface when measure(c) supported in Meshes

    # TODO Custom tests: Line, Plane

    SUPPORT_MATRIX = [
    # Name, example,    integral,line,surface,volume,    GaussLegendre,GaussKronrod,HAdaptiveCubature
        SupportItem("BezierCurve", bezier,       1, 1, 0, 0,   1, 1, 1),
        SupportItem("Box{1,T}", box1d,           1, 1, 0, 0,   1, 1, 1),
        SupportItem("Circle", circle,            1, 1, 0, 0,   1, 1, 1),
        # Line -- custom test
        SupportItem("Ring", ring_rect,           1, 1, 0, 0,   1, 1, 1),
        SupportItem("Rope", rope_rect,           1, 1, 0, 0,   1, 1, 1),
        SupportItem("Segment", seg_ne,           1, 1, 0, 0,   1, 1, 1),
        SupportItem("Sphere{2,T}", sphere2d,     1, 1, 0, 0,   1, 1, 1),

        SupportItem("Ball{2,T}", ball2d,         1, 0, 1, 0,   1, 1, 1),
        SupportItem("Box{2,T}", box2d,           1, 0, 1, 0,   1, 1, 1),
        SupportItem("CylinderSurface", cylsurf,  1, 0, 1, 0,   0, 1, 0),
        SupportItem("Disk", disk,                1, 0, 1, 0,   1, 1, 1),
        # ParaboloidSurface -- not yet supported
        SupportItem("Sphere{3,T}", sphere3d,     1, 0, 1, 0,   1, 1, 1),
        SupportItem("Triangle", triangle,        1, 0, 1, 0,   1, 1, 1),
        # Torus -- not yet supported
        # SimpleMesh -- not yet supported

        SupportItem("Ball{3,T}", ball3d,         1, 0, 0, 1,   1, 0, 1),
        SupportItem("Box{3,T}", box3d,           1, 0, 0, 1,   1, 0, 1)
    ]

    # Run all integral tests
    map(autotest, SUPPORT_MATRIX)
end
