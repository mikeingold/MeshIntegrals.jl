################################################################################
#                          Generalized 1D Methods
################################################################################

function _integral_1d(
    f,
    geometry,
    settings::GaussLegendre,
    FP::Type{T} = Float64
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
    f,
    geometry,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    integrand(t) = f(geometry(t)) * differential(geometry, [t])
    return QuadGK.quadgk(integrand, FP(0), FP(1); settings.kwargs...)[1]
end

function _integral_1d(
    f,
    geometry,
    settings::HAdaptiveCubature,
)
    return _integral(f, geometry, settings)
end

function _integral_1d(
    f,
    geometry,
    settings::HAdaptiveCubature,
    FP::Type{T}
) where {T<:AbstractFloat}
    return _integral(f, geometry, settings, FP)
end


################################################################################
#                   Specialized Methods for BezierCurve
################################################################################

function lineintegral(
    f::F,
    curve::Meshes.BezierCurve,
    settings::I,
    FP::Type{T};
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, I<:IntegrationAlgorithm, T<:AbstractFloat}
    return integral(f, curve, settings, FP; alg=alg)
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
    f::F,
    curve::Meshes.BezierCurve,
    settings::GaussLegendre,
    FP::Type{T} = Float64;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, T<:AbstractFloat}
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
    f::F,
    curve::Meshes.BezierCurve,
    settings::GaussKronrod,
    FP::Type{T} = Float64;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, T<:AbstractFloat}
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
    f::F,
    curve::Meshes.BezierCurve,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, T<:AbstractFloat}
    len = length(curve)
    point(t) = curve(t, alg)
    integrand(t) = len * f(point(t[1]))

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = fill(FP(0.5),3)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, FP[0], FP[1]; settings.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end

################################################################################
#                   Specialized Methods for Line
################################################################################

function integral(
    f::F,
    line::Meshes.Line,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, settings.n)

    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Line(line.a, line.a + Meshes.unormalize(line.b - line.a))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f along the Line
    integrand(x) = f(line(t(x))) * t′(x)
    return sum(w .* integrand(x) for (w,x) in zip(ws, xs))
end

function integral(
    f::F,
    line::Meshes.Line,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Line(line.a, line.a + Meshes.unormalize(line.b - line.a))

    # Integrate f along the Line
    return QuadGK.quadgk(t -> f(line(t)), FP(-Inf), FP(Inf); settings.kwargs...)[1]
end

function integral(
    f::F,
    line::Meshes.Line,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Line(line.a, line.a + Meshes.unormalize(line.b - line.a))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f along the Line
    integrand(x::AbstractVector) = f(line(t(x[1]))) * t′(x[1])

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = FP[0.5]
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, FP[-1], FP[1]; settings.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end

################################################################################
#                   Specialized Methods for Ray
################################################################################

function integral(
    f::F,
    ray::Meshes.Ray,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, settings.n)

    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, Meshes.unormalize(ray.v))

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
    f::F,
    ray::Meshes.Ray,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, Meshes.unormalize(ray.v))

    # Integrate f along the Ray
    return QuadGK.quadgk(t -> f(ray(t)), FP(0), FP(Inf); settings.kwargs...)[1]
end

function integral(
    f::F,
    ray::Meshes.Ray,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, Meshes.unormalize(ray.v))

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
    f::F,
    ring::Meshes.Ring,
    settings::I,
    FP::Type{T}
) where {F<:Function, I<:IntegrationAlgorithm, T<:AbstractFloat}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(f, segment, settings, FP), segments(ring))
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
    f::F,
    rope::Meshes.Rope,
    settings::I,
    FP::Type{T}
) where {F<:Function, I<:IntegrationAlgorithm, T<:AbstractFloat}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(f, segment, settings, FP), segments(rope))
end

function integral(
    f::F,
    rope::Meshes.Rope,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(f, segment, settings), segments(rope))
end
