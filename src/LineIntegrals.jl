module LineIntegrals
    using LinearAlgebra
    using Meshes
    using QuadGK

    #include("trajectory.jl")
    #export SurfacePathElement, SurfacePathSegment, SurfaceTrajectory

    include("integrate.jl")
    export integral

    include("utils.jl")
    export derivative, unitdirection
end
