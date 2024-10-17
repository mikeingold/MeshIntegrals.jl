using TestItemRunner
using TestItems

# For CI, run all tests not marked with the :extended tag
@run_package_tests filter=ti -> !(:extended in ti.tags) verbose=true

@testsnippet Setup begin
    using LinearAlgebra: norm
    using Meshes
    using Unitful
end
