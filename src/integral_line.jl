################################################################################
#                          Generalized 1D Methods
################################################################################

function _integral_1d(
    f,
    geometry,
    settings::GaussLegendre
)
    T = coordtype(geometry)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(T, settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = T(1/2) * x + T(1/2)

    # Integrate f along the geometry and apply a domain-correction factor for [-1,1] ↦ [0, 1]
    integrand((w,x)) = w * f(geometry(t(x))) * differential(geometry, [t(x)])
    return T(1/2) * sum(integrand, zip(ws, xs))
end

function _integral_1d(
    f,
    geometry,
    settings::GaussKronrod
)
    T = coordtype(geometry)
    integrand(t) = f(geometry(t)) * differential(geometry, [t])
    return QuadGK.quadgk(integrand, T(0), T(1); settings.kwargs...)[1]
end

function _integral_1d(
    f,
    geometry,
    settings::HAdaptiveCubature
)
    T = coordtype(geometry)
    integrand(t) = f(geometry(t[1])) * differential(geometry, t)
    return HCubature.hcubature(integrand, T[0], T[1]; settings.kwargs...)[1]
end


################################################################################
#                   Specialized Methods for BezierCurve
################################################################################

function lineintegral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::I;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, Dim, T, V, I<:IntegrationAlgorithm}
    return integral(f, curve, settings; alg=alg)
end

"""
    integral(f, curve::Meshes.BezierCurve, ::GaussLegendre; alg=Meshes.Horner())

Like [`integral`](@ref) but integrates along the domain defined a `curve`. By
default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function integral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::GaussLegendre;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(T, settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = T(1/2) * x + T(1/2)
    point(x) = curve(t(x), alg)

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return T(1/2) * length(curve) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end

"""
    integral(f, curve::BezierCurve, ::GaussKronrod; alg=Horner(), kws...)

Like [`integral`](@ref) but integrates along the domain defined a `curve`. By
default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function integral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::GaussKronrod;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function,Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(curve)
    point(t) = curve(t, alg)
    return QuadGK.quadgk(t -> len * f(point(t)), T(0), T(1); settings.kwargs...)[1]
end

"""
    integral(f, curve::BezierCurve, ::HAdaptiveCubature; alg=Horner(), kws...)

Like [`integral`](@ref) but integrates along the domain defined a `curve`. By
default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function integral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::HAdaptiveCubature;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function,Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(curve)
    point(t) = curve(t, alg)
    return hcubature(t -> len * f(point(t[1])), T[0], T[1]; settings.kwargs...)[1]
end

################################################################################
#                   Specialized Methods for Line
################################################################################

function integral(
    f::F,
    line::Meshes.Line{Dim,T},
    settings::GaussLegendre
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(T, settings.n)

    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Line(line.a, line.a + normalize(line.b - line.a))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f along the Line
    integrand(x) = f(line(t(x))) * t′(x)
    return sum(w .* integrand(x) for (w,x) in zip(ws, xs))
end

function integral(
    f::F,
    line::Meshes.Line{Dim,T},
    settings::GaussKronrod
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Line(line.a, line.a + normalize(line.b - line.a))

    # Integrate f along the Line
    return QuadGK.quadgk(t -> f(line(t)), T(-Inf), T(Inf); settings.kwargs...)[1]
end

function integral(
    f::F,
    line::Meshes.Line{Dim,T},
    settings::HAdaptiveCubature
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Line(line.a, line.a + normalize(line.b - line.a))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f along the Line
    integrand(x::AbstractVector) = f(line(t(x[1]))) * t′(x[1])
    return HCubature.hcubature(integrand, T[-1], T[1]; settings.kwargs...)[1]
end

################################################################################
#                   Specialized Methods for Ray
################################################################################

function integral(
    f::F,
    ray::Meshes.Ray{Dim,T},
    settings::GaussLegendre
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(T, settings.n)

    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, normalize(ray.v))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ [0,∞)
    t₁(x) = T(1/2) * x + T(1/2)
    t₂(x) = x / (1 - x^2)
    t₂′(x) = (1 + x^2) / (1 - x^2)^2
    t = t₂ ∘ t₁
    t′(x) = T(1/2) * t₂′(x)

    # Integrate f along the Ray
    integrand(x) = f(ray(t(x))) * t′(x)
    return sum(w .* integrand(x) for (w,x) in zip(ws, xs))
end

function integral(
    f::F,
    ray::Meshes.Ray{Dim,T},
    settings::GaussKronrod
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, normalize(ray.v))

    # Integrate f along the Ray
    return QuadGK.quadgk(t -> f(ray(t)), T(0), T(Inf); settings.kwargs...)[1]
end

function integral(
    f::F,
    ray::Meshes.Ray{Dim,T},
    settings::HAdaptiveCubature
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, normalize(ray.v))

    # Domain transformation: x ∈ [0,1] ↦ t ∈ [0,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2
    
    # Integrate f along the Ray
    integrand(x::AbstractVector) = f(ray(t(x[1]))) * t′(x[1])
    return HCubature.hcubature(integrand, T[0], T[1]; settings.kwargs...)[1]
end

################################################################################
#                    Specialized Methods for Ring, Rope
################################################################################

function integral(
    f::F,
    ring::Meshes.Ring,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(f, segment, settings), segments(ring))
end

function integral(
    f::F,
    rope::Meshes.Rope,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(f, segment, settings), segments(rope))
end