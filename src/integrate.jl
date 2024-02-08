# Validate that f has a method defined for f(::Point{Dim,T})
@inline function _validate_integrand(f,Dim,T)
    if hasmethod(f, (Point{Dim,T},))
        return nothing
    else
        error("The provided Function f must have a method f(::Meshes.Point{$Dim,$T})")
    end
end

################################################################################
#                              lineintegral
################################################################################

"""
    lineintegral(f, geometry; n=100)

Numerically integrate a given function `f(::Point)` along a 1D `geometry` using
a Gauss-Legendre quadrature of order `n`.

So long as `f` can be well-approximated by a polynomial of order `2n-1`, this
method should yield results with 16-digit accuracy in O(n) time. If `f` is know
to have some periodic content then `n` should (at a minimum) be greater than
the expected number of periods, e.g. `length(geometry)/lambda`.
"""
function lineintegral end

# Integrate f(::Point{Dim,T}) along a Segment{Dim,T}
function lineintegral(
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
    lineintegral(f, curve::Meshes.BezierCurve; n=100, alg=Meshes.Horner())

Like [`lineintegral`](@ref) but integrates along the domain defined a `curve`.
By default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function lineintegral(
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

# Integrate f(::Point{Dim,T}) over a Rope{Dim,T} (an open Chain)
function lineintegral(
    f::F,
    rope::Meshes.Rope{Dim,T};
    n::Int64=100
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    return sum(segment -> lineintegral(f, segment; n=n), segments(rope))
end

# Integrate f(::Point{Dim,T}) over a Ring{Dim,T} (a closed Chain)
function lineintegral(
    f::F,
    ring::Meshes.Ring{Dim,T};
    n::Int64=100
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    return sum(segment -> lineintegral(f, segment; n=n), segments(ring))
end

# Integrate f(::Point{Dim,T}) over an arbitrary geometry in {Dim,T}
function lineintegral(
    f::F,
    path::Vector{<:Meshes.Geometry{Dim,T}};
    n::Int64=100
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    return sum(section -> lineintegral(f, section; n=n), path)
end


################################################################################
#                            surfaceintegral
################################################################################

"""
    surfaceintegral(f, geometry; n=100)

Numerically integrate a given function `f(::Point)` over a 2D surface `geometry`
using a Gauss-Legendre quadrature of order `n`.

So long as `f` can be well-approximated by a polynomial of order `2n-1`, this
method should yield results with 16-digit accuracy in O(n) time. If `f` is know
to have some periodic content then `n` should (at a minimum) be greater than
the expected number of periods, e.g. `length(geometry)/lambda`.
"""
function surfaceintegral end

"""
    surfaceintegral(f, triangle::Meshes.Triangle; n=100)

Like [`surfaceintegral`](@ref) but integrates over the surface of a `triangle`
using a Gauss-Legendre quadrature rule of order `n` along each Barycentric
dimension of the triangle.
"""
function surfaceintegral(
    f,
    triangle::Meshes.Ngon{3,Dim,T};
    n::Int64 = 100
) where {Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = gausslegendre(n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x1,x2) = triangle(t(x1), t(x2))

    # Determine output type of f at a Point inside the triangle
    # Define an applicable zero value
    fzero = zero(f(point(-0.5,-0.5)))

    # Calculate weight-node product
    function weightednode((w1,w2), (x1,x2))
        if 0.0 <= (w1 + w2) <= 1.0
            # Valid coordinate (inside triangle)
            return w1 * w2 * f(point(x1,x2))
        else
            # Invalid coordinate (outside triangle)
            return fzero
        end
    end

    # Calculate 2D Gauss-Legendre integral of f over Barycentric coordinates [-1,1]^2
    # Apply a linear domain-correction factor [-1,1]^2 ↦ area(triangle)
    return 0.25 * abs(signarea(triangle)) .* sum(weightednode, zip(wws,xxs))
end


################################################################################
#                             volumeintegral
################################################################################

"""
    volumeintegral(f, geometry; n=100)

Numerically integrate a given function `f(::Point)` throughout a volumetric
`geometry` using a Gauss-Legendre quadrature of order `n`.

So long as `f` can be well-approximated by a polynomial of order `2n-1`, this
method should yield results with 16-digit accuracy in O(n) time. If `f` is know
to have some periodic content then `n` should (at a minimum) be greater than
the expected number of periods, e.g. `length(geometry)/lambda`.
"""
function volumeintegral end


################################################################################
#                                quadgk_line
################################################################################

"""
    quadgk_line(f, geometry; kws...)

Numerically integrate a given function `f(::Point)` along a 1D `geometry` using
the h-adaptive Gauss-Kronrod quadrature rule from QuadGK.jl. All standard
[`QuadGK.quadgk`](@ref) keyword arguments are supported.
"""
function quadgk_line end

# Integrate f(::Point{Dim,T}) along a Segment{Dim,T}
function quadgk_line(
    f::F,
    segment::Meshes.Segment{Dim,T};
    kwargs...
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(segment)
    point(t) = segment(t)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; kwargs...)
end

"""
    quadgk_line(f, curve::BezierCurve; alg=Horner(), kws...)

Like [`quadgk_line`](@ref) but integrates along the domain defined a `curve`.
By default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function quadgk_line(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V};
    alg::Meshes.BezierEvalMethod=Meshes.Horner(),
    kwargs...
) where {F<:Function,Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(curve)
    point(t) = curve(t, alg)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; kwargs...)
end

# Integrate f(::Point{Dim,T}) over a Ring{Dim,T} (a closed Chain)
function quadgk_line(
    f::F,
    ring::Meshes.Ring{Dim,T};
    kwargs...
) where {F<:Function, Dim, T}
    # Partition the Ring into Segments, integrate each, sum results
    chunks = map(segment -> quadgk_line(f, segment; kwargs...), segments(ring))
    return reduce(.+, chunks)
end

# Integrate f(::Point{Dim,T}) over a Rope{Dim,T} (an open Chain)
function quadgk_line(
    f::F,
    rope::Meshes.Rope{Dim,T};
    kwargs...
) where {F<:Function, Dim, T}
    # Partition the Rope into Segments, integrate each, sum results
    chunks = map(segment -> quadgk_line(f, segment; kwargs...), segments(rope))
    return reduce(.+, chunks)
end

"""
    quadgk_line(f, points::Point...; kws...)

Like [`quadgk_line`](@ref), but integrates along a domain befined by the linear
segments formed between a series of Points.
"""
function quadgk_line(
    f,
    pts::Meshes.Point{Dim,T}...;
    kwargs...
) where {Dim, T}
    # Collect Points into a Rope, integrate that
    rope = Meshes.Rope(pts...)
    return quadgk_line(f, rope; kwargs...)
end


################################################################################
#                                quadgk_surface
################################################################################

"""
    quadgk_surface(f, geometry; kws...)

Numerically integrate a given function `f(::Point)` over the surface of a
`geometry` using the h-adaptive Gauss-Kronrod quadrature rule from QuadGK.jl.
All standard [`QuadGK.quadgk`](@ref) keyword arguments are supported.
"""
function quadgk_surface end

"""
    quadgk_surface(f, triangle::Meshes.Triangle; kws...)

Like [`quadgk_surface`](@ref), but integrates `f` over the surface of a `triangle`
using a nested integration in the Barycentric coordinate domain.
"""
function quadgk_surface(
    f,
    triangle::Meshes.Ngon{3,Dim,T};
    kwargs...
) where {Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Change of variables: u,v [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(u,v) = triangle(t(u), t(v))

    # Integrate the Barycentric triangle in (u,v)-space: (0,0), (0,1), (1,0)
    #   i.e. \int_{0}^{1} \int_{0}^{1-u} f(u,v) dv du
    innerintegral(u) = QuadGK.quadgk(v -> f(point(u,v)), 0, 1-u; kwargs...)
    outerintegral = QuadGK.quadgk(innerintegral, 0, 1; kwargs...)

    # Apply a linear domain-correction factor 0.5 ↦ area(triangle)
    return 2.0 * abs(signarea(triangle)) .* outerintegral
end
