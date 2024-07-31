################################################################################
#                          Generalized 1D Methods
################################################################################

function _integral_1d(
    FP::Type{T} = Float64,
    f,
    geometry,
    settings::GaussLegendre
) where {T<:AbstractFloat}
    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = FP(1/2) * x + FP(1/2)

    # Integrate f along the geometry and apply a domain-correction factor for [-1,1] ↦ [0, 1]
    integrand((w,x)) = w * f(geometry(t(x))) * differential(geometry, [t(x)])
    return FP(1/2) * sum(integrand, zip(ws, xs))
end

function _integral_1d(
    FP::Type{T} = Float64,
    f,
    geometry,
    settings::GaussKronrod
)
    integrand(t) = f(geometry(t)) * differential(geometry, [t])
    return QuadGK.quadgk(integrand, FP(0), FP(1); settings.kwargs...)[1]
end

function _integral_1d(
    FP::Type{T} = Float64,
    f,
    geometry,
    settings::HAdaptiveCubature
)
    integrand(t) = f(geometry(t[1])) * differential(geometry, t)
    return HCubature.hcubature(integrand, FP[0], FP[1]; settings.kwargs...)[1]
end


################################################################################
#                   Specialized Methods for BezierCurve
################################################################################

function lineintegral(
    FP::Type{T} = Float64,
    f::F,
    curve::Meshes.BezierCurve,
    settings::I;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {T<:AbstractFloat, F<:Function, I<:IntegrationAlgorithm}
    return integral(FP, f, curve, settings; alg=alg)
end

function lineintegral(
    f::F,
    curve::Meshes.BezierCurve,
    settings::I;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, I<:IntegrationAlgorithm}
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
    FP::Type{T} = Float64,
    f::F,
    curve::Meshes.BezierCurve,
    settings::GaussLegendre;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,FP)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = FP(1/2) * x + FP(1/2)
    point(x) = curve(t(x), alg)

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return FP(1/2) * length(curve) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
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
    FP::Type{T} = Float64,
    f::F,
    curve::Meshes.BezierCurve,
    settings::GaussKronrod;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,FP)

    len = length(curve)
    point(t) = curve(t, alg)
    return QuadGK.quadgk(t -> len * f(point(t)), FP(0), FP(1); settings.kwargs...)[1]
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
    FP::Type{T} = Float64,
    f::F,
    curve::Meshes.BezierCurve,
    settings::HAdaptiveCubature;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function}
    # Validate the provided integrand function
    Dim = embeddim(curve)
    _validate_integrand(f,Dim,FP)

    len = length(curve)
    point(t) = curve(t, alg)
    return HCubature.hcubature(t -> len * f(point(t[1])), T[0], T[1]; settings.kwargs...)[1]
end

################################################################################
#                   Specialized Methods for Line
################################################################################

function integral(
    FP::Type{T} = Float64,
    f::F,
    line::Meshes.Line,
    settings::GaussLegendre
) where {F<:Function}
    # Validate the provided integrand function
    Dim = embeddim(line)
    _validate_integrand(f,Dim,FP)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, settings.n)

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
    FP::Type{T} = Float64,
    f::F,
    line::Meshes.Line,
    settings::GaussKronrod
) where {F<:Function}
    # Validate the provided integrand function
    Dim = embeddim(line)
    _validate_integrand(f,Dim,FP)

    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Line(line.a, line.a + normalize(line.b - line.a))

    # Integrate f along the Line
    return QuadGK.quadgk(t -> f(line(t)), FP(-Inf), FP(Inf); settings.kwargs...)[1]
end

function integral(
    FP::Type{T} = Float64,
    f::F,
    line::Meshes.Line,
    settings::HAdaptiveCubature
) where {F<:Function}
    # Validate the provided integrand function
    Dim = embeddim(line)
    _validate_integrand(f,Dim,FP)

    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Line(line.a, line.a + normalize(line.b - line.a))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f along the Line
    integrand(x::AbstractVector) = f(line(t(x[1]))) * t′(x[1])
    return HCubature.hcubature(integrand, FP[-1], FP[1]; settings.kwargs...)[1]
end

################################################################################
#                   Specialized Methods for Ray
################################################################################

function integral(
    FP::Type{T} = Float64,
    f::F,
    ray::Meshes.Ray,
    settings::GaussLegendre
) where {F<:Function}
    # Validate the provided integrand function
    Dim = embeddim(line)
    _validate_integrand(f,Dim,FP)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, settings.n)

    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, normalize(ray.v))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ [0,∞)
    t₁(x) = FP(1/2) * x + FP(1/2)
    t₁′(x) = FP(1/2)
    t₂(x) = x / (1 - x^2)
    t₂′(x) = (1 + x^2) / (1 - x^2)^2
    t = t₂ ∘ t₁
    t′(x) = t₂′(t₁(x)) * t₁′(x)

    # Integrate f along the Ray
    integrand(x) = f(ray(t(x))) * t′(x)
    return sum(w .* integrand(x) for (w,x) in zip(ws, xs))
end

function integral(
    FP::Type{T} = Float64,
    f::F,
    ray::Meshes.Ray,
    settings::GaussKronrod
) where {F<:Function}
    # Validate the provided integrand function
    Dim = embeddim(line)
    _validate_integrand(f,Dim,FP)

    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, normalize(ray.v))

    # Integrate f along the Ray
    return QuadGK.quadgk(t -> f(ray(t)), FP(0), FP(Inf); settings.kwargs...)[1]
end

function integral(
    FP::Type{T} = Float64,
    f::F,
    ray::Meshes.Ray,
    settings::HAdaptiveCubature
) where {F<:Function}
    # Validate the provided integrand function
    Dim = embeddim(line)
    _validate_integrand(f,Dim,FP)

    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, normalize(ray.v))

    # Domain transformation: x ∈ [0,1] ↦ t ∈ [0,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2
    
    # Integrate f along the Ray
    integrand(x::AbstractVector) = f(ray(t(x[1]))) * t′(x[1])
    return HCubature.hcubature(integrand, FP[0], FP[1]; settings.kwargs...)[1]
end

################################################################################
#                    Specialized Methods for Ring, Rope
################################################################################

function integral(
    FP::Type{T} = Float64,
    f::F,
    ring::Meshes.Ring,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(FP, f, segment, settings), segments(ring))
end

function integral(
    f::F,
    ring::Meshes.Ring,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(f, segment, settings), segments(ring))
end

function integral(
    FP::Type{T} = Float64,
    f::F,
    rope::Meshes.Rope,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(FP, f, segment, settings), segments(rope))
end

function integral(
    f::F,
    rope::Meshes.Rope,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(f, segment, settings), segments(rope))
end
