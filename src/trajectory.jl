abstract type SurfacePathElement{Dim,T} end

"""
    SufacePathSegment

This represents a line segment that traverses a surface in space, defined by a
segment with start and end points and a vector pointing in the direction of
the local surface normal.
"""
struct SurfacePathSegment{Dim,T} <: SurfacePathElement{Dim,T}
    segment::Meshes.Segment{Dim,T}
    normal::Meshes.Vec{Dim,T}
end

# TODO struct SurfacePathBezierCurve{Dim,T} <: SurfacePathElement{Dim,T} end

"""
    SurfaceTrajectory

This represents a path that traverses a surface in space, where the path is
defined by an ordered collection of SurfacePathElements.
"""
struct SurfaceTrajectory{Dim,T}
    path::Vector{<:SurfacePathElement{Dim,T}}
end
# TODO fix bug where constructor on Vector{SurfacePathSegment} fails

# Additional constructor for separate vectors
function SurfaceTrajectory(
    segs::Vector{S}, norms::Vector{N}
) where {S<:SurfacePathElement, N<:Meshes.Vec}
    # Convert all corresponding (Segments, Normals) into SurfaceSegments
    path = map(SurfacePathSegment, Iterators.zip(segs,norms))
    return SurfaceTrajectory(path)
end
