# TODO Custom tests: Line, Plane, CylinderSurface (non-right)?

function integraltest_scalar(geometry, rule)
    f(::Point) = 1.0
    fv(::Point) = fill(1.0,3)
    @test integral(f, geometry, rule) ≈ measure(geometry)
    @test integral(fv, geometry, rule) ≈ fill(measure(geometry),3)
end

function integraltest_dimensionful(geometry, rule, Package)
    W = Package.W
    f(::Point) = 1.0W
    fv(::Point) = fill(1.0W,3)
    @test integral(f, geometry, rule) ≈ measure(geometry)*W
    @test integral(fv, geometry, rule) ≈ fill(measure(geometry)*W,3)
end

function integraltest_all(geometry, rule)
    integraltest_scalar(geometry, rule)
    # TODO enable these when ready
    #integraltest_dimensionful(geometry, rule, Unitful)
    #integraltest_dimensionful(geometry, rule, DynamicQuantities)
end

struct SupportItem
    geometry
    integralfunction::Function
    gausslegendre::Bool
    gausskronrod::Bool
    hadaptivecubature::Bool
end

function autotest(item::SupportItem)
    # For each enabled solver type, run the test suite
    @testset "$typeof(geometry)" begin
        item.gausslegendre     && integraltest_all(geometry, GaussLegendre(100))
        item.gausskronrod      && integraltest_all(geometry, GaussKronrod())
        item.hadaptivecubature && integraltest_scalar(geometry, HAdaptiveCubature())
    end
end

SUPPORT_MATRIX = [
    SupportItem(bezier, lineintegral, true, true, true),       # BezierCurve
    SupportItem(box1d, lineintegral, true, true, true),        # Box{1,T}
    SupportItem(circle, lineintegral, true, true, true),       # Circle
    # Line -- custom test
    SupportItem(ring_rect, lineintegral, true, true, true),    # Ring
    SupportItem(rope_rect, lineintegral, true, true, true),    # Rope
    SupportItem(seg_ne, lineintegral, true, true, true),       # Segment
    SupportItem(sphere2d, lineintegral, true, true, true),     # Sphere{2,T}

    SupportItem(ball2d, surfaceintegral, true, true, true),       # Ball{2,T}
    SupportItem(box2d, surfaceintegral, true, true, true),        # Box{2,T}
    SupportItem(cylsurf, surfaceintegral, false, true, false),    # CylinderSurface
    SupportItem(disk, surfaceintegral, true, true, true),         # Disk
    # ParaboloidSurface -- not yet supported
    SupportItem(sphere3d, surfaceintegral, true, true, true),    # Sphere{3,T}
    SupportItem(triangle, surfaceintegral, true, true, true),    # Triangle
    # Torus -- not yet supported
    # SimpleMesh -- not yet supported
    
    SupportItem(ball3d, volumeintegral, false, true, true),    # Ball{3,T}
    SupportItem(box3d, volumeintegral, false, true, true)      # Box{3,T}
]

foreach(autotest, SUPPORT_MATRIX)
