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
	
	# Generator for Bernstein polynomial functions
	B(i,n) = t -> binomial(n,i) * t^i * (1-t)^(n-i)
	
	# Derivative = N Σ_{i=0}^{N-1} sigma(i)
    #   P indices adjusted for Julia 1-based array indexing
    N = degree(bz)
    P = bz.controls
    sigma(i) = B(i,N-1)(t) .* (P[(i+1)+1] - P[(i)+1])
	return N .* sum(sigma, 0:N-1)    # ::Vec{Dim,T}
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
