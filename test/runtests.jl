#using DynamicQuantities
using MeshIntegrals
using Meshes
#using Unitful
using Test

################################################################################
#                                Line Integrals
################################################################################

@testset "Integrals" begin
    # Points on unit circle at axes
    pt_e = Point( 1.0,  0.0, 0.0)
    pt_n = Point( 0.0,  1.0, 0.0)
    pt_w = Point(-1.0,  0.0, 0.0)
    pt_s = Point( 0.0, -1.0, 0.0)
    origin = Point(0.0, 0.0, 0.0)
    ẑ = Vec(0.0, 0.0, 1.0)

    # Line segments oriented CCW between points
    seg_ne = Segment(pt_e, pt_n)
    seg_nw = Segment(pt_n, pt_w)
    seg_sw = Segment(pt_w, pt_s)
    seg_se = Segment(pt_s, pt_e)

    # Rectangular trajectory CCW around the four points
    rect_traj_ring = Ring(pt_e, pt_n, pt_w, pt_s)
    rect_traj_rope = Rope(pt_e, pt_n, pt_w, pt_s, pt_e)

    # Approximately circular trajectory CCW around the unit circle
    unit_bezier = BezierCurve(
        [Point(cos(t), sin(t), 0.0) for t in range(0, 2pi, length=361)]
    )

    # Triangle on upper-half-plane
    triangle = Ngon(pt_e, pt_n, pt_w)

    # Unit circle/disk
    circle = Circle(Plane(origin,ẑ), 2.0)
    disk = Disk(Plane(origin,ẑ), 2.0)
    sphere2d = Sphere(Point(0.0,0.0), 2.0)

    # 2D Box on [-1,1]^2
    box2d = Box(Point(-1.0,-1.0), Point(1.0,1.0))

    for (name,rule) in [("Gauss-Legendre",GaussLegendre(10_000)), ("Gauss-Kronrod",GaussKronrod())]
        @testset "$name" begin
            @testset "Scalar-Valued Functions" begin
                f(::Point) = 1.0
                # Line Integrals
                @test lineintegral(f, seg_ne, rule) ≈ sqrt(2)                         # Meshes.Segment
                @test lineintegral(f, rect_traj_ring, rule) ≈ 4sqrt(2)                # Meshes.Ring
                @test lineintegral(f, rect_traj_rope, rule) ≈ 4sqrt(2)                # Meshes.Rope
                @test lineintegral(f, unit_bezier, rule) ≈ length(unit_bezier)        # Meshes.BezierCurve
                @test lineintegral(f, circle, rule) ≈ length(circle)                  # Meshes.Circle
                @test lineintegral(f, disk, rule) ≈ length(circle)                    # Meshes.Disk
                @test lineintegral(f, sphere2d, rule) ≈ length(sphere2d)              # Meshes.Sphere{2,T}
                @test lineintegral(f, triangle, rule) ≈ 2 + 2sqrt(2)                  # Meshes.Triangle
                @test lineintegral(f, box2d, rule) ≈ 8.0                              # Meshes.Box{2,T}

                # Surface Integrals
                @test isapprox(surfaceintegral(f, triangle, rule), 1.0; rtol=1e-3)    # Meshes.Triangle
                @test isapprox(surfaceintegral(f, box2d, rule), 4.0; rtol=1e-3)       # Meshes.Box{2,T}
                @test isapprox(surfaceintegral(f, circle, rule), area(disk); rtol=1e-3)     # Meshes.Circle
                @test isapprox(surfaceintegral(f, disk, rule), area(disk); rtol=1e-3)       # Meshes.Disk
                @test isapprox(surfaceintegral(f, sphere2d, rule), area(disk); rtol=1e-3)   # Meshes.Sphere{2,T}
            end
            @testset "Vector-Valued Functions" begin
                f(::Point) = [1.0, 1.0, 1.0]

                # Line Integrals
                @test lineintegral(f, seg_ne, rule) ≈ [sqrt(2), sqrt(2), sqrt(2)]              # Meshes.Segment
                @test lineintegral(f, rect_traj_ring, rule) ≈ [4sqrt(2), 4sqrt(2), 4sqrt(2)]   # Meshes.Ring
                @test lineintegral(f, rect_traj_rope, rule) ≈ [4sqrt(2), 4sqrt(2), 4sqrt(2)]   # Meshes.Rope
                @test lineintegral(f, unit_bezier, rule) ≈ length(unit_bezier) .* [1.0, 1.0, 1.0]    # Meshes.BezierCurve
                @test lineintegral(f, circle, rule) ≈ length(circle) .* [1.0, 1.0, 1.0]        # Meshes.Circle
                @test lineintegral(f, disk, rule) ≈ length(circle) .* [1.0, 1.0, 1.0]          # Meshes.Disk
                @test lineintegral(f, sphere2d, rule) ≈ length(sphere2d) .* [1.0, 1.0, 1.0]    # Meshes.Sphere{2,T}
                @test lineintegral(f, triangle, rule) ≈ (2 + 2sqrt(2)) .* [1.0, 1.0, 1.0]      # Meshes.Triangle
                @test lineintegral(f, box2d, rule) ≈ [8.0, 8.0, 8.0]                           # Meshes.Box{2,T}

                # Surface Integrals
                @test isapprox(surfaceintegral(f, triangle, rule), [1.0, 1.0, 1.0]; rtol=1e-3)   # Meshes.Triangle
                @test isapprox(surfaceintegral(f, box2d, rule), [4.0, 4.0, 4.0]; rtol=1e-3)      # Meshes.Box{2,T}
                @test isapprox(surfaceintegral(f, circle, rule), area(disk) .* [1.0, 1.0, 1.0]; rtol=1e-3)   # Meshes.Circle
                @test isapprox(surfaceintegral(f, disk, rule), area(disk) .* [1.0, 1.0, 1.0]; rtol=1e-3)     # Meshes.Disk
                @test isapprox(surfaceintegral(f, sphere2d, rule), area(disk) .* [1.0, 1.0, 1.0]; rtol=1e-3) # Meshes.Sphere{2,T}
            end
        end
    end

    @testset "Contour Integrals on a Point{1,Complex}-Domain" begin
        fc(z::Complex) = 1/z
        fc(p::Point{1,ComplexF64}) = fc(p.coords[1])

        # Construct a unit circle on the complex domain
        unit_circle_complex = Meshes.BezierCurve(
            [Point{1,ComplexF64}(cos(t) + sin(t)*im) for t in range(0,2pi,length=361)]
        )

        # 2πi Res_{z=0}(1/z) = \int_C (1/z) dz
        # Res_{z=0}(1/z) = 1
        # ∴ \int_C (1/z) dz = 2πi
        @test lineintegral(fc, unit_circle_complex, GaussKronrod()) ≈ 2π*im
    end
end

#= temp disabled
################################################################################
#                             Tests -- Unitful.jl
################################################################################

@testset "Integrate with Unitful.jl" begin
    m = Unitful.m
    Ω = Unitful.Ω

    # Points on unit circle at axes
    pt_e = Point( 1.0m,  0.0m, 0.0m)
    pt_n = Point( 0.0m,  1.0m, 0.0m)
    pt_w = Point(-1.0m,  0.0m, 0.0m)
    pt_s = Point( 0.0m, -1.0m, 0.0m)

    # Line segments oriented CCW between points
    seg_ne = Segment(pt_e, pt_n)
    seg_nw = Segment(pt_n, pt_w)
    seg_sw = Segment(pt_w, pt_s)
    seg_se = Segment(pt_s, pt_e)

    # Rectangular trajectory CCW around the four points
    rect_traj_segs = [seg_ne, seg_nw, seg_sw, seg_se]
    rect_traj_ring = Ring(pt_e, pt_n, pt_w, pt_s)
    rect_traj_rope = Rope(pt_e, pt_n, pt_w, pt_s, pt_e)

    # Approximately circular trajectory CCW around the unit-meter circle
    unit_circle = BezierCurve(
        [Point(cos(t)*m, sin(t)*m, 0.0m) for t in range(0, 2pi, length=361)]
    )

    @testset "Scalar-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = 1.0Ω/m
        @test lineintegral(f, seg_ne) ≈ sqrt(2)*Ω                         # Meshes.Segment
        @test lineintegral(f, rect_traj_segs) ≈ 4sqrt(2)*Ω                # Vector{::Meshes.Segment}
        @test lineintegral(f, rect_traj_ring) ≈ 4sqrt(2)*Ω                # Meshes.Ring
        @test lineintegral(f, rect_traj_rope) ≈ 4sqrt(2)*Ω                # Meshes.Rope
        @test isapprox(lineintegral(f, unit_circle), 2π*Ω; atol=0.15Ω)    # Meshes.BezierCurve
    end

    @testset "Vector-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = [1.0Ω/m, 1.0Ω/m, 1.0Ω/m]
        @test lineintegral(f, seg_ne) ≈ [sqrt(2), sqrt(2), sqrt(2)] .* Ω                  # Meshes.Segment
        @test lineintegral(f, rect_traj_segs)  ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)] .* Ω    # Vector{::Meshes.Segment}
        @test lineintegral(f, rect_traj_ring) ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)] .* Ω     # Meshes.Ring
        @test lineintegral(f, rect_traj_rope) ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)] .* Ω     # Meshes.Rope
        @test isapprox(lineintegral(f, unit_circle), [2π, 2π, 2π] .* Ω; atol=0.15Ω)    # Meshes.BezierCurve
    end
end

################################################################################
#                             Tests -- DynamicQuantities.jl
################################################################################

# TODO once these tests work identically to Unitful tests, consolidate them into
#   a single abstracted test loop
# for Package in (Unitful, DynamicQuantities)
#   @testset "Integrate with $Package"

@testset "Integrate with DynamicQuantities.jl" begin
    m = DynamicQuantities.m
    Ω = DynamicQuantities.Ω

    # Points on unit circle at axes
    pt_e = Point( 1.0m,  0.0m, 0.0m)
    pt_n = Point( 0.0m,  1.0m, 0.0m)
    pt_w = Point(-1.0m,  0.0m, 0.0m)
    pt_s = Point( 0.0m, -1.0m, 0.0m)

    # Line segments oriented CCW between points
    seg_ne = Segment(pt_e, pt_n)
    seg_nw = Segment(pt_n, pt_w)
    seg_sw = Segment(pt_w, pt_s)
    seg_se = Segment(pt_s, pt_e)

    # Rectangular trajectory CCW around the four points
    rect_traj_segs = [seg_ne, seg_nw, seg_sw, seg_se]
    rect_traj_ring = Ring(pt_e, pt_n, pt_w, pt_s)
    rect_traj_rope = Rope(pt_e, pt_n, pt_w, pt_s, pt_e)

    # Approximately circular trajectory CCW around the unit-meter circle
    unit_circle = BezierCurve(
        [Point(cos(t)*m, sin(t)*m, 0.0m) for t in range(0, 2pi, length=361)]
    )

    @testset "Scalar-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = 1.0Ω/m
        @test lineintegral(f, seg_ne) ≈ sqrt(2)*Ω                         # Meshes.Segment
        @test lineintegral(f, rect_traj_segs) ≈ 4sqrt(2)*Ω                # Vector{::Meshes.Segment}
        @test lineintegral(f, rect_traj_ring) ≈ 4sqrt(2)*Ω                # Meshes.Ring
        @test lineintegral(f, rect_traj_rope) ≈ 4sqrt(2)*Ω                # Meshes.Rope
        @test isapprox(lineintegral(f, unit_circle), 2π*Ω; atol=0.15)    # Meshes.BezierCurve
              # TODO change 0.15 => 0.15Ω once DynamicQuantities PR approved
    end

    @testset "Vector-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = [1.0Ω/m, 1.0Ω/m, 1.0Ω/m]
        @test all(isapprox.(lineintegral(f, seg_ne), sqrt(2)*Ω))                # Meshes.Segment
        @test all(isapprox.(lineintegral(f, rect_traj_segs), 4sqrt(2)*Ω))       # Vector{::Meshes.Segment}
        @test all(isapprox.(lineintegral(f, rect_traj_ring), 4sqrt(2)*Ω))       # Meshes.Ring
        @test all(isapprox.(lineintegral(f, rect_traj_rope), 4sqrt(2)*Ω))       # Meshes.Rope
        @test all(isapprox.(lineintegral(f, unit_circle), 2π*Ω; atol=0.15))    # Meshes.BezierCurve
             # TODO change 0.15 => 0.15Ω once DynamicQuantities PR approved
    end
end
=#