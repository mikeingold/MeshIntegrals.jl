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

@testitem "DifferentiationMethod" setup=[Setup] begin
    using MeshIntegrals: _has_analytical, _default_method, _guarantee_analytical

    #=
    # Test geometries
    sphere = Sphere(Point(0, 0, 0), 1.0)
    triangle = Triangle(Point(0, 0, 0), Point(0, 1, 0), Point(1, 0, 0))

    # _has_analytical of instances
    @test _has_analytical(sphere) == false
    @test _has_analytical(triangle) == true

    # _has_analytical of types
    @test _has_analytical(Meshes.Sphere) == false
    @test _has_analytical(Meshes.Triangle) == true
    =#

    # _default_method
    @test _default_method(Meshes.Sphere) isa FiniteDifference
    @test _default_method(sphere) isa FiniteDifference
    @test _default_method(Meshes.Triangle) isa Analytical
    @test _default_method(triangle) isa Analytical

    #=
    # _guarantee_analytical
    @test _guarantee_analytical(Meshes.Triangle, Analytical()) === nothing
    @test_throws "Analytical" _guarantee_analytical(Meshes.Triangle, FiniteDifference())
    =#

    # FiniteDifference
    @test FiniteDifference().ε ≈ 1e-6
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
