"""
    derivative(b::BezierCurve, t)

Determine the vector derivative of a Bezier curve `b` for the point on the
curve parameterized by value `t`.
"""
function derivative(bz::Meshes.BezierCurve{Dim,T,V}, t) where {Dim,T,V}
    # Parameter t restricted to domain [0,1] by definition
    if t < 0 || t > 1
        throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
    end

    # Aliases
    P = bz.controls
    N = degree(bz)

    # For bz := ∑_{i=0}^{N-1} a_i B_i^N
    # bz' := N ∑_{i=0}^{N-1} b_i B_i^{N-1},
    # where b_i := a_{i+1}-a_i.
    BezierCurve(N*(P[2:end].-P[1:end-1]))(t)
end

"""
    unitdirection(b::BezierCurve, t)

Determine a unit vector pointing in the forward (t+) direction of a Bezier
curve `b` for a point on the curve parameterized by value `t`.
"""
function unitdirection(bz::Meshes.BezierCurve{Dim,T,V}, t) where {Dim,T,V}
    # Parameter t restricted to domain [0,1] by definition
    if t < 0 || t > 1
        throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
    end

    # Normalize the derivative of the curve
    u = derivative(bz,t)
    LinearAlgebra.normalize!(u)
    return u    # ::Vec{Dim,T}
end
