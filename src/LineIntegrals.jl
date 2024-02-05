module LineIntegrals
    using FastGaussQuadrature
    using LinearAlgebra
    using Meshes
    using QuadGK

    include("integrate.jl")
    export integral
    # non-exported API: quadgk

    include("utils.jl")
    export derivative, unitdirection
end
