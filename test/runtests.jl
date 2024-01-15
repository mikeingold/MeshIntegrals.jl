using LineIntegrals
using Meshes
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

    @testset "Scalar-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = 1.0

        # integrate(f, ::Meshes.Segment)
        @test integrate(f, seg_ne) ≈ sqrt(2)

        # integrate(f, ::Vector{Meshes.Segment})
        @test integrate(f, rect_traj_segs) ≈ 4sqrt(2)

        # integate(f, ::Meshes.Ring)
        @test integrate(f, rect_traj_ring) ≈ 4sqrt(2)

        # integate(f, ::Meshes.Rope)
        @test integrate(f, rect_traj_rope) ≈ 4sqrt(2)

        # integrate(f, ::Meshes.BezierCurve)
        @test isapprox(integrate(f, unit_circle), 2pi; atol=0.15)
    end

    @testset "Vector-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = [1.0, 1.0, 1.0]

        # integrate(f, ::Meshes.Segment)
        @test integrate(f, seg_ne) ≈ [sqrt(2), sqrt(2), sqrt(2)]

        # integrate(f, ::Vector{Meshes.Segment})
        @test integrate(f, rect_traj_segs)  ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)]

        # integate(f, ::Meshes.Ring)
        @test integrate(f, rect_traj_ring) ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)]

        # integate(f, ::Meshes.Rope)
        @test integrate(f, rect_traj_rope) ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)]

        # integrate(f, ::Meshes.BezierCurve)
        @test isapprox(integrate(f, unit_circle), [2π, 2π, 2π]; atol=0.15)
    end
end

@testset "Integrate Unitful" begin
    # Points on unit circle at axes
    pt_e = Point( 1.0u"m",  0.0u"m", 0.0u"m")
    pt_n = Point( 0.0u"m",  1.0u"m", 0.0u"m")
    pt_w = Point(-1.0u"m",  0.0u"m", 0.0u"m")
    pt_s = Point( 0.0u"m", -1.0u"m", 0.0u"m")

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
        [Point(cos(t)*u"m", sin(t)*u"m", 0.0u"m") for t in range(0, 2pi, length=361)]
    )

    @testset "Scalar-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = 1.0u"Ω/m"

        # integrate(f, ::Meshes.Segment)
        @test integrate(f, seg_ne) ≈ sqrt(2)*u"Ω"

        # integrate(f, ::Vector{Meshes.Segment})
        @test integrate(f, rect_traj_segs) ≈ 4sqrt(2)*u"Ω"

        # integate(f, ::Meshes.Ring)
        @test integrate(f, rect_traj_ring) ≈ 4sqrt(2)*u"Ω"

        # integate(f, ::Meshes.Rope)
        @test integrate(f, rect_traj_rope) ≈ 4sqrt(2)*u"Ω"

        # integrate(f, ::Meshes.BezierCurve)
        @test isapprox(integrate(f, unit_circle), 2π*u"Ω"; atol=0.15u"Ω")
    end

    @testset "Vector-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = [1.0u"Ω/m", 1.0u"Ω/m", 1.0u"Ω/m"]

        # integrate(f, ::Meshes.Segment)
        @test integrate(f, seg_ne) ≈ [sqrt(2), sqrt(2), sqrt(2)] .* u"Ω"

        # integrate(f, ::Vector{Meshes.Segment})
        @test integrate(f, rect_traj_segs)  ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)] .* u"Ω"

        # integate(f, ::Meshes.Ring)
        @test integrate(f, rect_traj_ring) ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)] .* u"Ω"

        # integate(f, ::Meshes.Rope)
        @test integrate(f, rect_traj_rope) ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)] .* u"Ω"

        # integrate(f, ::Meshes.BezierCurve)
        @test isapprox(integrate(f, unit_circle), [2π, 2π, 2π] .* u"Ω"; atol=0.15u"Ω")
    end
end
