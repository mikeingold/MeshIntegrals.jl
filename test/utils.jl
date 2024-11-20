@testitem "Utilities" setup=[Setup] begin
    using LinearAlgebra: norm

    # _kvector
    v = Meshes.Vec(3, 4)
    @test norm(MeshIntegrals._kvector(v)) ≈ 5.0

    # _units
    p = Point(1.0u"cm", 2.0u"mm", 3.0u"m")
    @test MeshIntegrals._units(p) == u"m"
end

@testitem "DifferentiationMethod" setup=[Setup] begin
    using MeshIntegrals: _has_analytical, _default_method, _guarantee_analytical

    # _has_analytical of instances
    bezier = BezierCurve([Point(t, sin(t), 0.0) for t in range(-π, π, length = 361)])
    @test _has_analytical(bezier) == true
    sphere = Sphere(Point(0, 0, 0), 1.0)
    @test _has_analytical(sphere) == false

    # _has_analytical of types
    @test _has_analytical(Meshes.BezierCurve) == true
    @test _has_analytical(Meshes.Line) == true
    @test _has_analytical(Meshes.Plane) == true
    @test _has_analytical(Meshes.Ray) == true
    @test _has_analytical(Meshes.Sphere) == false
    @test _has_analytical(Meshes.Tetrahedron) == true
    @test _has_analytical(Meshes.Triangle) == true

    # _guarantee_analytical
    @test _guarantee_analytical(Meshes.Line, Analytical()) === nothing
    @test_throws "Analytical" _guarantee_analytical(Meshes.Line, FiniteDifference())

    # _default_method
    @test _default_method(Meshes.BezierCurve) isa Analytical
    @test _default_method(bezier) isa Analytical
    @test _default_method(Meshes.Sphere) isa FiniteDifference
    @test _default_method(sphere) isa FiniteDifference

    # FiniteDifference
    @test FiniteDifference().ε ≈ 1e-6

    # AutoEnzyme
    @test AutoEnzyme() isa AutoEnzyme
end
