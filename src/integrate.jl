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

# Integrate f(::Point{Dim,T}) over a Segment
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
    point(x) = segment(t(x))

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return 0.5 * length(segment) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end

"""
    integral(f, curve; n=100, alg=Meshes.Horner())

Numerically integrate a given function `f(::Point)` along a BezierCurve using
a Gauss-Legendre quadrature of order `n`.

So long as `f` can be well-approximated by a polynomial of order `2n-1`, this
method should yield results with 16-digit accuracy in O(n) time. If `f` is know
to have some periodic content then `n` should (at a minimum) be greater than
the expected number of periods, e.g. `length(geometry)/lambda`.

By default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function integral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V};
    n::Int64=100,
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x) = curve(t(x), alg)

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return 0.5 * length(curve) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end

# Integrate f(::Point{Dim,T}) over a Rope (an open Chain)
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

"""
    quadgk(f, segment::Segment; kws...)

Like [`quadgk`](@ref), but integrates along the domain of a Segment. All standard
[`quadgk`](@ref) keyword arguments remain available.
"""
function QuadGK.quadgk(
    f::F,
    segment::Meshes.Segment{Dim,T};
    kwargs...
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(segment)
    point(t) = segment(t)
    return quadgk(t -> len * f(point(t)), 0, 1; kwargs...)
end

"""
    quadgk(f, curve::BezierCurve; alg=Horner(), kws...)

Like [`quadgk`](@ref), but integrates along the domain of a Bezier curve. All
standard [`quadgk`](@ref) keyword arguments remain available.

By default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function QuadGK.quadgk(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V};
    alg::Meshes.BezierEvalMethod=Meshes.Horner(),
    kwargs...
) where {F<:Function,Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(curve)
    point(t) = curve(t, alg)
    return quadgk(t -> len * f(point(t)), 0, 1; kwargs...)
end

"""
    quadgk(f, points::Point...; kws...)

Like [`quadgk`](@ref), but integrates along a domain befined by the linear
segments formed between a series of Points. All standard [`quadgk`](@ref)
keyword arguments remain available.
"""
function QuadGK.quadgk(
    f,
    pts::Meshes.Point{Dim,T}...;
    kwargs...
) where {Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Collect Points into a Rope, integrate that
    rope = Meshes.Rope(pts...)
    return quadgk(f, rope; kwargs...)
end

"""
    quadgk(f, ring::Ring; kws...)

Like [`quadgk`](@ref), but integrates along the domain of a Ring. All standard
[`quadgk`](@ref) keyword arguments remain available.
"""
function QuadGK.quadgk(
    f::F,
    ring::Meshes.Ring{Dim,T};
    kwargs...
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Partition the Ring into Segments, integrate each, sum results
    chunks = map(segment -> quadgk(f, segment; kwargs...), segments(ring))
    return reduce(.+, chunks)
end

"""
    quadgk(f, rope::Rope; kws...)

Like [`quadgk`](@ref), but integrates along the domain of a Rope. All standard
[`quadgk`](@ref) keyword arguments remain available.
"""
function QuadGK.quadgk(
    f::F,
    rope::Meshes.Rope{Dim,T};
    kwargs...
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Partition the Rope into Segments, integrate each, sum results
    chunks = map(segment -> quadgk(f, segment; kwargs...), segments(rope))
    return reduce(.+, chunks)
end
