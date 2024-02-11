module MeshIntegrals
    using FastGaussQuadrature
    using HCubature
    using LinearAlgebra
    using Meshes
    using QuadGK

    include("integrate.jl")
    export lineintegral, surfaceintegral, volumeintegral

    include("utils.jl")
    export derivative, unitdirection
end
