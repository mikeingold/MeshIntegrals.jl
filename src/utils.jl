################################################################################
#                     Derivatives, Tangents, Jacobians
################################################################################

"""
    jacobian(geometry, ts, ε=1e-6)

Calculate the Jacobian of a geometry at some parametric point `ts` using a simple
central-finite-difference approximation with step size `ε`.
"""
function jacobian(
    geometry,
    ts::AbstractVector{T};
    ε=1e-6
) where {T<:AbstractFloat}

    function ∂r_∂tn!(εv,n)
        # Construct a zero vector with ε in the n'th element 
        εv .= 0.0
        εv[n] = ε

        # Parametric points over which to perform finite difference approximation
        a = ts - εv
        b = ts + εv

        return (geometry(b...) - geometry(a...)) / 2ε
    end

    # Allocate a re-usable ε vector
    εv = similar(ts)
    ∂r_∂tn(n) = ∂r_∂tn!(εv,n)

    return map(∂r_∂tn, 1:length(ts))
end

"""
    derivative(b::BezierCurve, t)

Determine the vector derivative of a Bezier curve `b` for the point on the
curve parameterized by value `t`.
"""
function derivative(
    bz::Meshes.BezierCurve{Dim,T,V},
    t
) where {Dim,T,V}
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
function unitdirection(
    bz::Meshes.BezierCurve{Dim,T,V}, 
    t
) where {Dim,T,V}
    # Parameter t restricted to domain [0,1] by definition
    if t < 0 || t > 1
        throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
    end

    # Normalize the derivative of the curve
    u = derivative(bz,t)
    LinearAlgebra.normalize!(u)
    return u    # ::Vec{Dim,T}
end

################################################################################
#                               Internal Tools
################################################################################

# Validate that f has a method defined for f(::Point{Dim,T})
@inline function _validate_integrand(f,Dim,T)
    if hasmethod(f, (Point{Dim,T},))
        return nothing
    else
        error("The provided Function f must have a method f(::Meshes.Point{$Dim,$T})")
    end
end
