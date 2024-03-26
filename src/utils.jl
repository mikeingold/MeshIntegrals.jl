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
    # Get the partial derivative along the n'th axis via finite difference approximation
    #   where ts is the current parametric position (εv is a reusable buffer)
    function ∂r_∂tn!(εv, ts, n)
        if ts[n] < T(0.01)
            return ∂r_∂tn_right!(εv, ts, n)
        elseif T(0.99) < ts[n]
            return ∂r_∂tn_left!(εv, ts, n)
        else
            return ∂r_∂tn_central!(εv, ts, n)
        end
    end

    # Central finite difference
    function ∂r_∂tn_central!(εv, ts, n)
        εv .= T(0)
        εv[n] = ε
        a = ts - εv
        b = ts + εv
        return (geometry(b...) - geometry(a...)) / 2ε
    end

    # Left finite difference
    function ∂r_∂tn_left!(εv, ts, n)
        εv .= T(0)
        εv[n] = ε
        a = ts - εv
        b = ts
        return (geometry(b...) - geometry(a...)) / ε
    end

    # Right finite difference
    function ∂r_∂tn_right!(εv, ts, n)
        εv .= T(0)
        εv[n] = ε
        a = ts
        b = ts + εv
        return (geometry(b...) - geometry(a...)) / ε
    end

    # Allocate a re-usable ε vector
    εv = similar(ts)

    ∂r_∂tn(n) = ∂r_∂tn!(εv,ts,n)
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

# Calculate Gauss-Legendre nodes/weights and convert to type T
function _gausslegendre(T, n)
    xs, ws = FastGaussQuadrature.gausslegendre(n)
    return T.(xs), T.(ws)
end
