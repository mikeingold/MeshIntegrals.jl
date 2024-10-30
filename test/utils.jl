@testitem "Utilities" setup=[Setup] begin
    using LinearAlgebra: norm

    # _kvector
    v = Meshes.Vec(3, 4)
    @test norm(MeshIntegrals._kvector(v)) ≈ 5.0

    # _units
    p = Point(1.0u"cm", 2.0u"mm", 3.0u"m")
    @test MeshIntegrals._units(p) == u"m"

    @testset "DifferentiationMethod" begin
        bezier = BezierCurve([Point(t, sin(t), 0.0) for t in range(-π, π, length = 361)])
        sphere = Sphere(Point(0, 0, 0), 1.0)
        
        @test MeshIntegrals.has_analytical(Meshes.BezierCurve) == true
        @test MeshIntegrals.has_analytical(bezier) == true
        @test MeshIntegrals.has_analytical(Meshes.Sphere) == false
        @test MeshIntegrals.has_analytical(sphere) == false

        @test MeshIntegrals.default_method(Meshes.BezierCurve) isa Analytic()
        @test MeshIntegrals.default_method(bezier) isa Analytic()
        @test MeshIntegrals.default_method(Meshes.Sphere) isa FiniteDifference()
        @test MeshIntegrals.default_method(sphere) isa FiniteDifference()

        @test FiniteDifference().ε ≈ 1e-6
    end
end
