################################################################################
#                               API/DOCSTRINGS
################################################################################

"""
    integral(f, geometry; n=100)

Numerically integrate a given function `f(::Point)` along a 1D `geometry` using
a Gauss-Legendre quadrature of order `n`.

So long as `f` can be well-approximated by a polynomial of order `2n-1`, this
method should yield results with 16-digit accuracy in O(n) time. If `f` is know
to have some periodic content then `n` should (at a minimum) be greater than
the expected number of periods, e.g. `length(geometry)/lambda`.
"""
function integral end


################################################################################
#                              integral METHODS
################################################################################

# Validate that f has a method defined for f(::Point)
@inline function _validate_integrand(f,Dim,T)
    if !hasmethod(f, (Point{Dim,T},))
        error("The provided Function f must have a method f(::Point{$Dim,$T})")
    end

    return nothing
end

# Integrate f(::Point{Dim,T}) over a Segment
#   Allocations:
#     gausslegendre: 2n * sizeof(Float64)
function integral(
    f::F,
    segment::Meshes.Segment{Dim,T};
    n::Int64=100
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5

    # Wrapper function such that fx(x) ↦ f(segment(t))
    fx(x) = f(segment(t(x)))

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return 0.5 * length(segment) * sum(w .* fx(x) for (w,x) in zip(ws, xs))
end

# Integrate f(::Point{Dim,T}) over a BezierCurve
#   Allocations:
#     gausslegendre: 2n * sizeof(Float64)
function integral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V};
    n::Int64=100
) where {F<:Function, Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5

    # Wrapper function such that fx(x) ↦ f(segment(t))
    fx(x) = f(curve(t(x)))

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return 0.5 * length(curve) * sum(w .* fx(x) for (w,x) in zip(ws, xs))
end

# Integrate f(::Point{Dim,T}) over a Rope (an open Chain)
#   Allocations:
#     segments(rope)
#     integral(f, segment)
# TODO implement internal use _in-place methods with cached memory for performance
function integral(
    f::F,
    rope::Meshes.Rope{Dim,T};
    n::Int64=100
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    return sum(segment -> integral(f, segment; n=n), segments(rope))
end

# Integrate f(::Point{Dim,T}) over a Ring (a closed Chain)
#   Allocations:
#     segments(rope)
#     integral(f, segment)
function integral(
    f::F,
    ring::Meshes.Ring{Dim,T};
    n::Int64=100
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    return sum(segment -> integral(f, segment; n=n), segments(ring))
end

# Integrate f(::Point{Dim,T}) over an arbitrary geometry construct
function integral(
    f::F,
    path::Vector{<:Meshes.Geometry{Dim,T}};
    n::Int64=100
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    return sum(section -> integral(f, section; n=n), path)
end


################################################################################
#                       QuadGK.quadgk Methods
################################################################################

# Implement QuadGK.quadgk over a Meshes.Segment
function QuadGK.quadgk(
    f::F,
    segment::Meshes.Segment{Dim,T};
    kwargs...
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(segment)
    return quadgk(t -> len * f(segment(t)), 0, 1; kwargs...)
end

# Implement QuadGK.quadgk over a Meshes.Ring
function QuadGK.quadgk(
    f::F,
    ring::Meshes.Ring{Dim,T};
    kwargs...
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    chunks = map(segment -> quadgk(f, segment; kwargs...), segments(ring))
    return reduce(.+, chunks)
end

# Implement QuadGK.quadgk over a Meshes.Rope
function QuadGK.quadgk(
    f::F,
    rope::Meshes.Rope{Dim,T};
    kwargs...
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    chunks = map(segment -> quadgk(f, segment; kwargs...), segments(rope))
    return reduce(.+, chunks)
end

# Implement QuadGK.quadgk over a Meshes.BezierCurve
function QuadGK.quadgk(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V};
    kwargs...
) where {F<:Function,Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(curve)
    return quadgk(t -> len * f(curve(t)), 0, 1; kwargs...)
end

# Implement QuadGK/quadgk(f, pts...)
function QuadGK.quadgk(
    f,
    pts::Meshes.Point{Dim,T}...;
    kwargs...
) where {Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    rope = Meshes.Rope(pts...)
    return quadgk(f, rope; kwargs...)
end
