@testitem "Aqua.jl" begin
    using Aqua
    # As of v0.11.4:
    # - Ambiguities check disabled since it fails due to upstream findings
    # - Verified that no ambiguities exist within MeshIntegrals.jl
    Aqua.test_all(MeshIntegrals; ambiguities = false)
end

@testitem "ExplicitImports.jl" begin
    using ExplicitImports: check_no_implicit_imports, check_no_stale_explicit_imports,
                           check_all_explicit_imports_via_owners,
                           check_all_explicit_imports_are_public,
                           check_all_qualified_accesses_via_owners,
                           check_no_self_qualified_accesses
    @test isnothing(check_no_implicit_imports(MeshIntegrals))
    @test isnothing(check_no_stale_explicit_imports(MeshIntegrals))
    @test isnothing(check_all_explicit_imports_via_owners(MeshIntegrals))
    @test isnothing(check_all_explicit_imports_are_public(MeshIntegrals))
    @test isnothing(check_all_qualified_accesses_via_owners(MeshIntegrals))
    @test isnothing(check_no_self_qualified_accesses(MeshIntegrals))
end
