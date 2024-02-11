################################################################################
#                            Gauss-Legendre
################################################################################

function lineintegral(
    f::F,
    segment::Meshes.Segment{Dim,T},
    settings::GaussLegendre
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x) = segment(t(x))

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return 0.5 * length(segment) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end

"""
    lineintegral(f, curve::Meshes.BezierCurve, ::GaussLegendre; alg=Meshes.Horner())

Like [`lineintegral`](@ref) but integrates along the domain defined a `curve`.
By default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function lineintegral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::GaussLegendre;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x) = curve(t(x), alg)

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return 0.5 * length(curve) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end

function lineintegral(
    f::F,
    rope::Meshes.Rope{Dim,T},
    settings::GaussLegendre
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    return sum(segment -> lineintegral(f, segment, settings), segments(rope))
end

function lineintegral(
    f::F,
    ring::Meshes.Ring{Dim,T},
    settings::GaussLegendre
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    return sum(segment -> lineintegral(f, segment, settings), segments(ring))
end

function lineintegral(
    f,
    pts::Meshes.Point{Dim,T}...,
    settings::GaussLegendre
) where {Dim, T}
    # Collect Points into a Rope, integrate that
    rope = Meshes.Rope(pts...)
    return lineintegral(f, rope, settings)
end


################################################################################
#                                Gauss-Kronrod
################################################################################

function lineintegral(
    f::F,
    segment::Meshes.Segment{Dim,T},
    settings::GaussKronrod
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(segment)
    point(t) = segment(t)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end

"""
    lineintegral(f, curve::BezierCurve, ::GaussKronrod; alg=Horner(), kws...)

Like [`lineintegral`](@ref) but integrates along the domain defined a `curve`.
By default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function lineintegral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::GaussKronrod;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function,Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(curve)
    point(t) = curve(t, alg)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end

function lineintegral(
    f::F,
    ring::Meshes.Ring{Dim,T};
    kwargs...
) where {F<:Function, Dim, T}
    # Partition the Ring into Segments, integrate each, sum results
    chunks = map(segment -> quadgk_line(f, segment; kwargs...), segments(ring))
    return reduce(.+, chunks)
end

function lineintegral(
    f,
    pts::Meshes.Point{Dim,T}...,
    settings::GaussKronrod
) where {Dim, T}
    # Collect Points into a Rope, integrate that
    rope = Meshes.Rope(pts...)
    return lineintegral(f, rope, settings)
end

function lineintegral(
    f::F,
    rope::Meshes.Rope{Dim,T},
    settings::GaussKronrod
) where {F<:Function, Dim, T}
    # Partition the Rope into Segments, integrate each, sum results
    chunks = map(segment -> lineintegral(f, segment, settings), segments(rope))
    return reduce(.+, chunks)
end
