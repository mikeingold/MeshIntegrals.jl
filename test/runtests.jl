using TestItemRunner
using TestItems

# For CI, run all tests not marked with the :extended tag
@run_package_tests filter=ti -> !(:extended in ti.tags) verbose=true

@testsnippet Setup begin
    using LinearAlgebra: norm
    using Meshes
    using Unitful

    # environment settings
    isCI = "CI" ∈ keys(ENV)

    # float settings
    AD = if isCI
        if ENV["AD"] == "FiniteDifference"
            FiniteDifference()
        elseif ENV["AD"] == "Enzyme"
            using Enzyme: Enzyme
            AutoEnzyme()
        end
    else
        using Enzyme: Enzyme
        AutoEnzyme()
    end

    @info "Running tests with autodiff backend $AD"

    for g in (:integral, :lineintegral, :surfaceintegral, :volumeintegral)
        newname = Symbol(g, "_test")
        @eval begin
            function $newname(f, geometry, rule; kwargs...)
                MeshIntegrals.$g(f, geometry, rule; diff_method = AD, kwargs...)
            end
        end
    end

    for g in (:lineintegral, :surfaceintegral, :volumeintegral)
        newname = Symbol(g, "_test")
        @eval begin
            function $newname(f, geometry; kwargs...)
                MeshIntegrals.$g(f, geometry; diff_method = AD, kwargs...)
            end
        end
    end
end
