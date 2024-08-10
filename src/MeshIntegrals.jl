module MeshIntegrals
    using LinearAlgebra
    using Meshes
    using Unitful

    import FastGaussQuadrature
    import HCubature
    import QuadGK

    include("utils.jl")
    export jacobian, derivative, unitdirection

    include("integral.jl")
    include("integral_aliases.jl")
    include("integral_line.jl")
    include("integral_surface.jl")
    include("integral_volume.jl")
    export GaussKronrod, GaussLegendre, HAdaptiveCubature
    export integral, lineintegral, surfaceintegral, volumeintegral
end
