struct SurfaceSegment{Dim,T}
    segment::Meshes.Segment{Dim,T}
    normal::Meshes.Vec{Dim,T}
end

struct SurfaceTrajectory{Dim,T}
    path::Vector{SurfaceSegment{Dim,T}}
end

function SurfaceTrajectory(segs, norms)
    # Convert all corresponding (Segments, Normals) into SurfaceSegments
    path = map(SurfaceSegment, Iterators.zip(segs,norms))
    return SurfaceTrajectory(path)
end
