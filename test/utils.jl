@testsnippet Utils begin
    using LinearAlgebra: norm
    using Meshes
    using MeshIntegrals
    using MeshIntegrals: _default_diff_method, _parametric, _units, _zeros, _ones
    using Unitful
    import Enzyme
end

@testitem "Utilities" setup=[Utils] begin
    # _KVector
    v = Meshes.Vec(3, 4)
    @test norm(MeshIntegrals._KVector(v)) ≈ 5.0u"m"

    # _units
    p = Point(1.0u"cm", 2.0u"mm", 3.0u"m")
    @test _units(p) == u"m"

    # _zeros
    @test _zeros(2) == (0.0, 0.0)
    @test _zeros(Float32, 2) == (0.0f0, 0.0f0)

    # _ones
    @test _ones(2) == (1.0, 1.0)
    @test _ones(Float32, 2) == (1.0f0, 1.0f0)
end

@testitem "Differentiation" setup=[Utils] begin
    # _default_diff_method
    sphere = Sphere(Point(0, 0, 0), 1.0)
    @test _default_diff_method(Meshes.Sphere, Float64) isa AutoEnzyme
    @test _default_diff_method(sphere, Float64) isa AutoEnzyme
    @test _default_diff_method(sphere, BigFloat) isa FiniteDifference

    # FiniteDifference
    @test FiniteDifference().ε ≈ 1e-6

    # Two-argument jacobian
    segment = Segment(Point(0), Point(1))
    @test MeshIntegrals.jacobian(segment, (0.5,)) == (Vec(1),)

    # Test jacobian with wrong number of parametric coordinates
    box = Box(Point(0, 0), Point(1, 1))
    @test_throws ArgumentError MeshIntegrals.jacobian(box, zeros(3), FiniteDifference())
    @test_throws ArgumentError MeshIntegrals.jacobian(box, zeros(3), AutoEnzyme())
end

@testitem "_ParametricGeometry" setup=[Utils] begin
    pt_n = Point(0, 3, 0)
    pt_w = Point(-7, 0, 0)
    pt_e = Point(8, 0, 0)
    ẑ = Vec(0, 0, 1)
    triangle = _parametric(Triangle(pt_n, pt_w, pt_e))
    tetrahedron = _parametric(Tetrahedron(pt_n, pt_w, pt_e, pt_n + ẑ))

    # Ensure error is thrown for out-of-bounds coordinate
    for ts in [(-1, 0), (0, -1), (2, 0), (0, 2)]
        @test_throws "not defined" triangle(ts...)
    end
    for ts in [(-1, 0, 0), (0, -1, 0), (0, 0, -1), (2, 0, 0), (0, 2, 0), (0, 0, 2)]
        @test_throws "not defined" tetrahedron(ts...)
    end
end
