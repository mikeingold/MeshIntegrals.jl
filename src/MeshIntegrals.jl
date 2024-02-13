module MeshIntegrals
    using FastGaussQuadrature
    using HCubature
    using LinearAlgebra
    using Meshes
    using QuadGK

    include("integral.jl")
    include("lineintegral.jl")
    export GaussKronrod, GaussLegendre, HAdaptiveCubature
    export lineintegral, surfaceintegral, volumeintegral

    include("utils.jl")
    export derivative, unitdirection
end
