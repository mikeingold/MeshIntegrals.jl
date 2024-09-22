@testitem "Aqua.jl" begin
    using Aqua
    # As of v0.11.4:
    # - Ambiguities check disabled since it fails due to upstream findings
    # - Verified that no ambiguities exist within MeshIntegrals.jl
    Aqua.test_all(MeshIntegrals; ambiguities=false)
end
