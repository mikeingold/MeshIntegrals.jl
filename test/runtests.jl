using DynamicQuantities
using LineIntegrals
using Meshes
using QuadGK
using Unitful
using Test

################################################################################
#                             Tests -- Integrals
################################################################################

@testset "Integrate" begin
    # Points on unit circle at axes
    pt_e = Point( 1.0,  0.0, 0.0)
    pt_n = Point( 0.0,  1.0, 0.0)
    pt_w = Point(-1.0,  0.0, 0.0)
    pt_s = Point( 0.0, -1.0, 0.0)

    # Line segments oriented CCW between points
    seg_ne = Segment(pt_e, pt_n)
    seg_nw = Segment(pt_n, pt_w)
    seg_sw = Segment(pt_w, pt_s)
    seg_se = Segment(pt_s, pt_e)

    # Rectangular trajectory CCW around the four points
    rect_traj_segs = [seg_ne, seg_nw, seg_sw, seg_se]
    rect_traj_ring = Ring(pt_e, pt_n, pt_w, pt_s)
    rect_traj_rope = Rope(pt_e, pt_n, pt_w, pt_s, pt_e)

    # Approximately circular trajectory CCW around the unit circle
    unit_circle = BezierCurve(
        [Point(cos(t), sin(t), 0.0) for t in range(0, 2pi, length=361)]
    )

    @testset "QuadGK Methods" begin
        f(::Point{Dim,T}) where {Dim,T} = 1.0
        @test quadgk(f, seg_ne)[1] ≈ sqrt(2)                        # Meshes.Segment
        @test quadgk(f, rect_traj_ring)[1] ≈ 4sqrt(2)               # Meshes.Ring
        @test quadgk(f, rect_traj_rope)[1] ≈ 4sqrt(2)               # Meshes.Rope
        @test isapprox(quadgk(f, unit_circle)[1], 2pi; atol=0.15)   # Meshes.BezierCurve
        @test quadgk(f, pt_e, pt_n, pt_w, pt_s, pt_e)[1] ≈ 4sqrt(2)    # Varargs of Meshes.Point
        @test quadgk(t -> exp(-t), 0, 10, 100)[1] ≈ 1    # Verify compatibility with generic Vararg{T} method
    end

    @testset "Caught Errors" begin
        # Catch wrong method signature: f(x,y,z) vs f(::Point)
        fvec(x,y,z) = x*y*z
        @test_throws ErrorException integral(fvec, seg_ne)          # Meshes.Segment
        @test_throws ErrorException integral(fvec, rect_traj_segs)  # Vector{::Meshes.Segment}
        @test_throws ErrorException integral(fvec, rect_traj_ring)  # Meshes.Ring
        @test_throws ErrorException integral(fvec, rect_traj_rope)  # Meshes.Rope
        @test_throws ErrorException integral(fvec, unit_circle)     # Meshes.BezierCurve
    end

    @testset "Scalar-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = 1.0
        @test integral(f, seg_ne) ≈ sqrt(2)                         # Meshes.Segment
        @test integral(f, rect_traj_segs) ≈ 4sqrt(2)                # Vector{::Meshes.Segment}
        @test integral(f, rect_traj_ring) ≈ 4sqrt(2)                # Meshes.Ring
        @test integral(f, rect_traj_rope) ≈ 4sqrt(2)                # Meshes.Rope
        @test isapprox(integral(f, unit_circle), 2pi; atol=0.15)    # Meshes.BezierCurve
    end

    @testset "Vector-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = [1.0, 1.0, 1.0]
        @test integral(f, seg_ne) ≈ [sqrt(2), sqrt(2), sqrt(2)]                # Meshes.Segment
        @test integral(f, rect_traj_segs) ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)]   # Vector{::Meshes.Segment}
        @test integral(f, rect_traj_ring) ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)]   # Meshes.Ring
        @test integral(f, rect_traj_rope) ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)]   # Meshes.Rope
        @test isapprox(integral(f, unit_circle), [2π, 2π, 2π]; atol=0.15)      # Meshes.BezierCurve
    end

    @testset "Results Consistent with QuadGK" begin
        # Test handling of real-valued functions
        fr(x) = exp(-x)
        fr(p::Point) = fr(p.coords[1])
        @test quadgk(fr, Point(0,0), Point(100,0))[1] ≈ quadgk(fr, 0, 100)[1]
    end

    @testset "Contour Integrals on the Complex-Domain" begin
        fc(z::Complex) = 1/z
        fc(p::Point{1,ComplexF64}) = fc(p.coords[1])

        # Construct a unit circle on the complex domain
        unit_circle_complex = Meshes.BezierCurve(
            [Point{1,ComplexF64}(cos(t) + sin(t)*im) for t in range(0,2pi,length=361)]
        )

        @test integral(fc, unit_circle_complex, n=1000) ≈ 2π*im
    end
end

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
        @test integral(f, seg_ne) ≈ sqrt(2)*Ω                         # Meshes.Segment
        @test integral(f, rect_traj_segs) ≈ 4sqrt(2)*Ω                # Vector{::Meshes.Segment}
        @test integral(f, rect_traj_ring) ≈ 4sqrt(2)*Ω                # Meshes.Ring
        @test integral(f, rect_traj_rope) ≈ 4sqrt(2)*Ω                # Meshes.Rope
        @test isapprox(integral(f, unit_circle), 2π*Ω; atol=0.15Ω)    # Meshes.BezierCurve
    end

    @testset "Vector-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = [1.0Ω/m, 1.0Ω/m, 1.0Ω/m]
        @test integral(f, seg_ne) ≈ [sqrt(2), sqrt(2), sqrt(2)] .* Ω                  # Meshes.Segment
        @test integral(f, rect_traj_segs)  ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)] .* Ω    # Vector{::Meshes.Segment}
        @test integral(f, rect_traj_ring) ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)] .* Ω     # Meshes.Ring
        @test integral(f, rect_traj_rope) ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)] .* Ω     # Meshes.Rope
        @test isapprox(integral(f, unit_circle), [2π, 2π, 2π] .* Ω; atol=0.15Ω)    # Meshes.BezierCurve
    end
end

################################################################################
#                             Tests -- DynamicQuantities.jl
################################################################################

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
        @test integral(f, seg_ne) ≈ sqrt(2)*Ω                         # Meshes.Segment
        @test integral(f, rect_traj_segs) ≈ 4sqrt(2)*Ω                # Vector{::Meshes.Segment}
        @test integral(f, rect_traj_ring) ≈ 4sqrt(2)*Ω                # Meshes.Ring
        @test integral(f, rect_traj_rope) ≈ 4sqrt(2)*Ω                # Meshes.Rope
        @test isapprox(integral(f, unit_circle), 2π*Ω; atol=0.15)    # Meshes.BezierCurve
              # TODO change 0.15 => 0.15Ω once DynamicQuantities PR approved
    end

    @testset "Vector-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = [1.0Ω/m, 1.0Ω/m, 1.0Ω/m]
        @test all(isapprox.(integral(f, seg_ne), sqrt(2)*Ω))                # Meshes.Segment
        @test all(isapprox.(integral(f, rect_traj_segs), 4sqrt(2)*Ω))       # Vector{::Meshes.Segment}
        @test all(isapprox.(integral(f, rect_traj_ring), 4sqrt(2)*Ω))       # Meshes.Ring
        @test all(isapprox.(integral(f, rect_traj_rope), 4sqrt(2)*Ω))       # Meshes.Rope
        @test all(isapprox.(integral(f, unit_circle), 2π*Ω; atol=0.15))    # Meshes.BezierCurve
             # TODO change 0.15 => 0.15Ω once DynamicQuantities PR approved
    end
end
