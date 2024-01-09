module LineIntegrals
    using Meshes
    using QuadGK

    include("structs.jl")
    export SurfaceSegment, SurfaceTrajectory

    include("integrate.jl")
    export integrate

    include("utils.jl")
    export derivative
end
