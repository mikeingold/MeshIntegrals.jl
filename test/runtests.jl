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
    cylsurf = CylinderSurface(pt_e, pt_w, 2.0)

    @testset "Errors Expected on Invalid Methods" begin
        f(::Point) = 1.0

        # lineintegrals should fail on these geometries with paramdims ≠ 1
        @test_throws MethodError lineintegral(f, ball2d)      # Ball{2,T}
        @test_throws MethodError lineintegral(f, ball3d)      # Ball{3,T}
        @test_throws MethodError lineintegral(f, box2d)       # Box{2,T}
        @test_throws MethodError lineintegral(f, box3d)       # Box{3,T}
        @test_throws MethodError lineintegral(f, cylsurf)     # CylinderSurface
        @test_throws MethodError lineintegral(f, disk)        # Disk
        @test_throws MethodError lineintegral(f, triangle)    # Ngon{3,Dim,T}
        # ParaboloidSurface
        @test_throws MethodError lineintegral(f, sphere3d)    # Sphere{3,T}
        # Torus

        # surfaceintegrals should fail on these geometries with paramdims ≠ 1
        @test_throws MethodError surfaceintegral(f, bezier)       # BezierCurve
        @test_throws MethodError surfaceintegral(f, ball3d)       # Ball{3,T}
        @test_throws MethodError surfaceintegral(f, box1d)        # Box{1,T}
        @test_throws MethodError surfaceintegral(f, box3d)        # Box{3,T}
        @test_throws MethodError surfaceintegral(f, circle)       # Circle
        @test_throws MethodError surfaceintegral(f, line_ne)      # Line
        @test_throws MethodError surfaceintegral(f, ring_rect)    # Ring
        @test_throws MethodError surfaceintegral(f, rope_rect)    # Rope
        @test_throws MethodError surfaceintegral(f, seg_ne)       # Segment
        @test_throws MethodError surfaceintegral(f, sphere2d)     # Sphere{2,T}
    end

    test_solvers = [
        ("Gauss-Legendre", GaussLegendre(100)),
        ("Gauss-Kronrod", GaussKronrod()),
        ("H-Adaptive Cubature", HAdaptiveCubature())
    ]

    for (name,rule) in test_solvers
        @testset "$name" begin
            @testset "Scalar-Valued Functions" begin
                f(::Point) = 1.0
                fLine(p::Point) = exp(-(p.coords[1])^2)
                fLineVal = sqrt(2pi)

                # Line Integrals
                @test lineintegral(f, bezier, rule) ≈ length(bezier)        # BezierCurve
                @test lineintegral(f, box1d, rule) ≈ length(box1d)          # Box{1,T}
                @test lineintegral(f, circle, rule) ≈ length(circle)        # Circle
                @test lineintegral(fLine, line_ne, rule) ≈ fLineVal         # Line
                @test lineintegral(f, ring_rect, rule) ≈ length(ring_rect)  # Ring
                @test lineintegral(f, rope_rect, rule) ≈ length(rope_rect)  # Rope
                @test lineintegral(f, seg_ne, rule) ≈ length(seg_ne)        # Segment
                @test lineintegral(f, sphere2d, rule) ≈ length(sphere2d)    # Sphere{2,T}

                # Surface Integrals
                @test isapprox(surfaceintegral(f, ball2d, rule), area(ball2d); rtol=1e-6)       # Ball{2,T}
                @test isapprox(surfaceintegral(f, box2d, rule), area(box2d); rtol=1e-6)         # Box{2,T}
                @test isapprox(surfaceintegral(f, disk, rule), area(disk); rtol=1e-6)           # Disk
                @test isapprox(surfaceintegral(f, triangle, rule), area(triangle); rtol=1e-6)   # Triangle
                @test isapprox(surfaceintegral(f, cylsurf, rule), area(cylsurf); rtol=1e-6)     # CylinderSurface
                    # TODO add test for a non-right-cylinder surface when measure(c) is fixed in Meshes

                # Volume Integrals (skip for GaussKronrod rules)
                if rule != GaussKronrod()
                    @test volumeintegral(f, box3d, rule) ≈ volume(box3d)     # Box{3,T}
                end
            end
            @testset "Vector-Valued Functions" begin
                f(::Point) = [1.0, 1.0, 1.0]
                fLine(p::Point) = fill(exp(-(p.coords[1])^2), 3)
                fLineVal = sqrt(2pi)

                # Line Integrals
                @test lineintegral(f, bezier, rule) ≈ fill(length(bezier),3)         # BezierCurve
                @test lineintegral(f, box1d, rule) ≈ fill(length(box1d),3)           # Box{1,T}
                @test lineintegral(f, circle, rule) ≈ fill(length(circle),3)         # Circle
                @test lineintegral(fLine, line_ne, rule) ≈ fill(fLineVal,3)          # Line
                @test lineintegral(f, ring_rect, rule) ≈ fill(length(ring_rect),3)   # Ring
                @test lineintegral(f, rope_rect, rule) ≈ fill(length(rope_rect),3)   # Rope
                @test lineintegral(f, seg_ne, rule) ≈ fill(length(seg_ne),3)         # Segment
                @test lineintegral(f, sphere2d, rule) ≈ fill(length(sphere2d),3)     # Sphere{2,T}

                # Surface Integrals
                @test isapprox(surfaceintegral(f, ball2d, rule), fill(area(ball2d),3); rtol=1e-6)      # Ball{2,T}
                @test isapprox(surfaceintegral(f, box2d, rule), fill(area(box2d),3); rtol=1e-6)        # Box{2,T}
                @test isapprox(surfaceintegral(f, disk, rule), fill(area(disk),3); rtol=1e-6)          # Disk
                @test isapprox(surfaceintegral(f, triangle, rule), fill(area(triangle),3); rtol=1e-6)  # Triangle
                @test isapprox(surfaceintegral(f, cylsurf, rule), fill(area(cylsurf),3); rtol=1e-6)    # CylinderSurface

                # Volume Integrals (skip for GaussKronrod rules)
                if rule != GaussKronrod()
                    @test volumeintegral(f, box3d, rule) ≈ fill(volume(box3d),3)     # Box{3,T}
                end
            end
        end
    end

    @testset "Contour Integrals on a Point{1,Complex}-Domain" begin
        fc(z::Complex) = 1/z
        fc(p::Point{1,ComplexF64}) = fc(p.coords[1])

        # Construct a unit circle on the complex domain
        unit_circle_complex = Meshes.BezierCurve([Point{1,ComplexF64}(cos(t) + sin(t)*im) for t in range(0,2pi,length=361)])

        # 2πi Res_{z=0}(1/z) = \int_C (1/z) dz
        # Res_{z=0}(1/z) = 1
        # ∴ \int_C (1/z) dz = 2πi
        @test lineintegral(fc, unit_circle_complex, GaussKronrod()) ≈ 2π*im
    end
end
