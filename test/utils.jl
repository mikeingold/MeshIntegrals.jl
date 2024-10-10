@testitem "Utilities" setup=[Setup] begin
    using LinearAlgebra: norm

    # _kvector
    v = Meshes.Vec(3, 4)
    @test norm(MeshIntegrals._kvector(v)) ≈ 5.0

    # _units
    p = Point(1.0u"cm", 2.0u"mm", 3.0u"m")
    @test MeshIntegrals._units(p) == u"m"
end
