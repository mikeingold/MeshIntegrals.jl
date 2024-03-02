module MeshIntegrals
    using FastGaussQuadrature
    using HCubature
    using LinearAlgebra
    using Meshes
    using QuadGK

    include("utils.jl")
    export jacobian, derivative, unitdirection

    include("integral.jl")
    include("integral_line.jl")
    include("integral_surface.jl")
    include("integral_volume.jl")
    export GaussKronrod, GaussLegendre, HAdaptiveCubature
    export integral, lineintegral, surfaceintegral, volumeintegral
end
