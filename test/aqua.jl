@testitem "Aqua.jl" begin
    using Aqua
    # As of v0.11.4:
    # - Ambiguities check disabled since it fails due to upstream findings
    # - Verified that no ambiguities exist within MeshIntegrals.jl
    Aqua.test_all(MeshIntegrals; ambiguities = false)
end

@testitem "ExplicitImports.jl" begin
    using ExplicitImports: check_no_implicit_imports, check_no_stale_explicit_imports
    @test isnothing(check_no_implicit_imports(MeshIntegrals))
    @test isnothing(check_no_stale_explicit_imports(MeshIntegrals))
end
