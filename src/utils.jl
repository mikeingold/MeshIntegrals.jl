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

    # Ensure that this implementation is tractible: limited by ability to calculate
    #   binomial(N, N/2) without overflow. It's possible to extend this range by
    #   converting N to a BigInt, but this results in always returning BigFloat's.
    N <= 1028 || error("This algorithm overflows for curves with ⪆1000 control points.")

    # Generator for Bernstein polynomial functions
    B(i,n) = t -> binomial(Int128(n),i) * t^i * (1-t)^(n-i)

    # Derivative = N Σ_{i=0}^{N-1} sigma(i)
    #   P indices adjusted for Julia 1-based array indexing
    sigma(i) = B(i,N-1)(t) .* (P[(i+1)+1] - P[(i)+1])
    return N .* sum(sigma, 0:N-1)    # ::Vec{Dim,T}
end

#=
# Calculate the derivative of a `geometry` with 1 parametric dimension at `t`
function derivative(geometry, t)
    ε = 1e-6
    return (geometry(t+ε) - geometry(t-ε)) / 2ε
end

# Calculate the Jacobian of a `geometry` with 2 parametric dimensions at `(u,v)`
function jacobian2d(geometry, u, v)
    ε = 1e-6
    ∂r_∂u = (geometry(u+ε,v) - geometry(u-ε,v)) / 2ε
    ∂r_∂v = (geometry(u,v+ε) - geometry(u,v-ε)) / 2ε
    return [∂r_∂u, ∂r_∂v]
end

# Calculate the Jacobian of a `geometry` with 2 parametric dimensions at `(s,t,u)`
function jacobian2d(geometry, s, t, u)
    ε = 1e-6
    ∂r_∂s = (geometry(s+ε,t,u) - geometry(s-ε,t,u)) / 2ε
    ∂r_∂t = (geometry(s,t+ε,u) - geometry(s,t-ε,u)) / 2ε
    ∂r_∂u = (geometry(s,t,u+ε) - geometry(s,t,u-ε)) / 2ε
    return [∂r_∂s, ∂r_∂t, ∂r_∂u]
end
=#

# Calculate the Jacobian of a geometry at some parametric point using a simple
#   central-finite-difference approximation
function jacobian(geometry, ts::AbstractVector{T}; ε=1e-6) where {T<:AbstractFloat}

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

# Get the corners of a 2D Box
_corners(box::Meshes.Box{2,T}) where {T} = [box(0,0), box(1,0), box(1,1), box(0,1)]
