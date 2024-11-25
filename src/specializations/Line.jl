################################################################################
#                   Specialized Methods for Line
#
# Why Specialized?
#   The Line geometry is a special case, representing a line of infinite length
#   that passes through two points. This requires another domain transformation
#   mapping from the typical parametric region [0,1] to an infinite one (-∞,∞).
################################################################################

function integral(
        f::F,
        line::Meshes.Line,
        rule::GaussLegendre;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    _guarantee_analytical(Meshes.Line, diff_method)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, rule.n)

    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Meshes.Line(line.a, line.a + Meshes.unormalize(line.b - line.a))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f along the Line
    differential(line, x) = t′(x) * _units(line(0))
    integrand(x) = f(line(t(x))) * differential(line, x)
    return sum(w .* integrand(x) for (w, x) in zip(ws, xs))
end

function integral(
        f::F,
        line::Meshes.Line,
        rule::GaussKronrod;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    _guarantee_analytical(Meshes.Line, diff_method)

    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Meshes.Line(line.a, line.a + Meshes.unormalize(line.b - line.a))

    # Integrate f along the Line
    domainunits = _units(line(0))
    integrand(t) = f(line(t)) * domainunits
    return QuadGK.quadgk(integrand, FP(-Inf), FP(Inf); rule.kwargs...)[1]
end

function integral(
        f::F,
        line::Meshes.Line,
        rule::HAdaptiveCubature;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    _guarantee_analytical(Meshes.Line, diff_method)

    # Normalize the Line s.t. line(t) is distance t from origin point
    line = Meshes.Line(line.a, line.a + Meshes.unormalize(line.b - line.a))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f along the Line
    differential(line, x) = t′(x) * _units(line(0))
    integrand(x::AbstractVector) = f(line(t(x[1]))) * differential(line, x[1])

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = (FP(0.5),)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, (-_one(FP),), (_one(FP),); rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end

################################################################################
#                               jacobian
################################################################################

_has_analytical(::Type{T}) where {T <: Meshes.Line} = true
