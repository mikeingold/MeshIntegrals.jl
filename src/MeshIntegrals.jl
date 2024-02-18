module MeshIntegrals
    using FastGaussQuadrature
    using HCubature
    using LinearAlgebra
    using Meshes
    using QuadGK

    include("integral.jl")
    include("integral_line.jl")
    include("integral_surface.jl")
    include("integral_volume.jl")
    export GaussKronrod, GaussLegendre, HAdaptiveCubature
    export lineintegral, surfaceintegral, volumeintegral

    include("utils.jl")
    export derivative, unitdirection
end
