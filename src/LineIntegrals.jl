module LineIntegrals
    using FastGaussQuadrature
    using LinearAlgebra
    using Meshes
    using QuadGK

    include("integrate.jl")
    export lineintegral, surfaceintegral, volumeintegral
    export quadgk_line, quadgk_surface

    include("utils.jl")
    export derivative, unitdirection
end
