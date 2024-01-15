using LineIntegrals
using Meshes
using Test

################################################################################
#                             Tests -- Integrals
################################################################################

@testset "Integration" begin
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
    rect_traj = [seg_ne, seg_nw, seg_sw, seg_se]

    # Circular trajectory CCW around the unit circle
    unit_circle = BezierCurve(pt_e, pt_n, pt_w, pt_s, pt_e)

    @testset "Scalar-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = 1.0

        # integrate(f, ::Vector{Meshes.Segment})
        @test integrate(f, rect_traj) ≈ 4 * sqrt(2)

        # integrate(f, ::Meshes.BezierCurve)
        @test integrate(f, unit_circle) ≈ 2π
    end

    @testset "Vector-Valued Functions" begin
        f(::Point{Dim,T}) where {Dim,T} = [1.0, 1.0, 1.0]

        # integrate(f, ::Vector{Meshes.Segment})
        @test integrate(f, rect_traj)  ≈ 4 .* [sqrt(2), sqrt(2), sqrt(2)]

        # integrate(f, ::Meshes.BezierCurve)
        @test integrate(f, unit_circle) ≈ [2π, 2π, 2π]
    end
end
