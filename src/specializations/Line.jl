################################################################################
#                   Specialized Methods for Line
################################################################################

function integral(
        f::F,
        line::Meshes.Line,
        rule::GaussLegendre,
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, rule.n)

    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Line(line.a, line.a + Meshes.unormalize(line.b - line.a))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f along the Line
    domainunits = _units(line(0))
    integrand(x) = f(line(t(x))) * t′(x) * domainunits
    return sum(w .* integrand(x) for (w, x) in zip(ws, xs))
end

function integral(
        f::F,
        line::Meshes.Line,
        rule::GaussKronrod,
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Line(line.a, line.a + Meshes.unormalize(line.b - line.a))

    # Integrate f along the Line
    domainunits = _units(line(0))
    integrand(t) = f(line(t)) * domainunits
    return QuadGK.quadgk(integrand, FP(-Inf), FP(Inf); rule.kwargs...)[1]
end

function integral(
        f::F,
        line::Meshes.Line,
        rule::HAdaptiveCubature,
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Line(line.a, line.a + Meshes.unormalize(line.b - line.a))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f along the Line
    domainunits = _units(line(0))
    integrand(x::AbstractVector) = f(line(t(x[1]))) * t′(x[1]) * domainunits

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = FP[0.5]
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, FP[-1], FP[1]; rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end
