@testitem "Utilities" setup=[Setup] begin
    using LinearAlgebra: norm
    using MeshIntegrals: _units, _zeros, _ones

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

@testitem "Differentiation" setup=[Setup] begin
    using MeshIntegrals: _default_diff_method

    # _default_diff_method
    sphere = Sphere(Point(0, 0, 0), 1.0)
    @test _default_diff_method(Meshes.Sphere) isa FiniteDifference
    @test _default_diff_method(sphere) isa FiniteDifference

    # FiniteDifference
    @test FiniteDifference().ε ≈ 1e-6

    # Test jacobian with wrong number of parametric coordinates
    box = Box(Point(0, 0), Point(1, 1))
    @test_throws ArgumentError jacobian(box, zeros(3))
end

@testitem "_ParametricGeometry" setup=[Setup] begin
    using MeshIntegrals: _parametric

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
