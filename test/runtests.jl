using MeshIntegrals
using Meshes
using Test

#using DynamicQuantities
#using Unitful

################################################################################
#                                Line Integrals
################################################################################

@testset "Integrals" begin
    # Spatial descriptors
    origin3d = Point(0.0, 0.0, 0.0)
    origin2d = Point(0.0, 0.0)
    ẑ = Vec(0.0, 0.0, 1.0)

    # Points on unit circle at axes
    pt_e = Point( 1.0,  0.0, 0.0)
    pt_n = Point( 0.0,  1.0, 0.0)
    pt_w = Point(-1.0,  0.0, 0.0)
    pt_s = Point( 0.0, -1.0, 0.0)

    # Line segments/paths oriented CCW between points
    seg_ne = Segment(pt_e, pt_n)
    line_ne = Line(pt_e, pt_n)
    ring_rect = Ring(pt_e, pt_n, pt_w, pt_s)
    rope_rect = Rope(pt_e, pt_n, pt_w, pt_s, pt_e)

    # Approximately circular trajectory CCW around the unit circle
    bezier = BezierCurve([Point(cos(t), sin(t), 0.0) for t in range(0, 2pi, length=361)])

    # Test Geometries
    ball2d = Ball(origin2d, 2.0)
    ball3d = Ball(origin3d, 2.0)
    box1d = Box(Point(-1.0), Point(1.0))
    box2d = Box(Point(-1.0,-1.0), Point(1.0,1.0))
    box3d = Box(Point(-1.0,-1.0,-1.0), Point(1.0,1.0,1.0))
    circle = Circle(Plane(origin3d,ẑ), 2.0)
    disk = Disk(Plane(origin3d,ẑ), 2.0)
    sphere2d = Sphere(origin2d, 2.0)
    sphere3d = Sphere(origin3d, 2.0)
    triangle = Ngon(pt_e, pt_n, pt_w)
    cylsurf = CylinderSurface(pt_e, pt_w, 2.5)
    # TODO add test for a non-right-cylinder surface when measure(c) is fixed in Meshes

    # TODO Custom tests: Line, Plane, CylinderSurface (non-right)?

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
            @test_throws MethodError intf(f, geometry, rule)
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

    SUPPORT_MATRIX = [
        SupportItem("BezierCurve", bezier,       false, true, false, false,   true, true, true),
        SupportItem("Box{1,T}", box1d,           false, true, false, false,   true, true, true),
        SupportItem("Circle", circle,            false, true, false, false,   true, true, true),
        # Line -- custom test
        SupportItem("Ring", ring_rect,           false, true, false, false,   true, true, true),
        SupportItem("Rope", rope_rect,           false, true, false, false,   true, true, true),
        SupportItem("Segment", seg_ne,           false, true, false, false,   true, true, true),
        SupportItem("Sphere{2,T}", sphere2d,     false, true, false, false,   true, true, true),

        SupportItem("Ball{2,T}", ball2d,         false, false, true, false,   true, true, true),
        SupportItem("Box{2,T}", box2d,           false, false, true, false,   true, true, true),
        SupportItem("CylinderSurface", cylsurf,  false, false, true, false,   false, true, false),
        SupportItem("Disk", disk,                false, false, true, false,   true, true, true),
        # ParaboloidSurface -- not yet supported
        SupportItem("Sphere{3,T}", sphere3d,     false, false, true, false,   true, true, true),
        SupportItem("Triangle", triangle,        false, false, true, false,   true, true, true),
        # Torus -- not yet supported
        # SimpleMesh -- not yet supported
        
        SupportItem("Ball{3,T}", ball3d,         false, false, false, true,   false, true, true),
        SupportItem("Box{3,T}", box3d,           false, false, false, true,   false, true, true)
    ]

    map(autotest, SUPPORT_MATRIX)
end
