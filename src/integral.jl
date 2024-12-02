################################################################################
#                         Master Integral Function
################################################################################

"""
    integral(f, geometry[, rule]; diff_method=_default_method(geometry), FP=Float64)

Numerically integrate a given function `f(::Point)` over the domain defined by
a `geometry` using a particular numerical integration `rule` with floating point
precision of type `FP`.

# Arguments
- `f`: an integrand function, i.e. any callable with a method `f(::Meshes.Point)`
- `geometry`: some `Meshes.Geometry` that defines the integration domain
- `rule`: optionally, the `IntegrationRule` used for integration (by default
`GaussKronrod()` in 1D and `HAdaptiveCubature()` else)

# Keyword Arguments
- `diff_method::DifferentiationMethod = _default_method(geometry)`: the method to
use for calculating Jacobians that are used to calculate differential elements
- `FP = Float64`: the floating point precision desired.
"""
function integral end

# If only f and geometry are specified, select default rule
function integral(
        f,
        geometry::Geometry,
        rule::I = Meshes.paramdim(geometry) == 1 ? GaussKronrod() : HAdaptiveCubature();
        kwargs...
) where {I <: IntegrationRule}
    _integral(f, geometry, rule; kwargs...)
end

################################################################################
#                            Integral Workers
################################################################################

# GaussKronrod
function _integral(
        f,
        geometry,
        rule::GaussKronrod;
        FP::Type{T} = Float64,
        diff_method::DM = _default_diff_method(geometry)
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    # Implementation depends on number of parametric dimensions over which to integrate
    N = Meshes.paramdim(geometry)
    if N == 1
        integrand(t) = f(geometry(t)) * differential(geometry, (t,), diff_method)
        return QuadGK.quadgk(integrand, zero(FP), one(FP); rule.kwargs...)[1]
    elseif N == 2
        # Issue deprecation warning
        Base.depwarn("Use `HAdaptiveCubature` instead of \
                     `GaussKronrod` for surfaces.", :integral)

        # Nested integration
        integrand2D(u, v) = f(geometry(u, v)) * differential(geometry, (u, v), diff_method)
        ∫₁(v) = QuadGK.quadgk(u -> integrand2D(u, v), zero(FP), one(FP); rule.kwargs...)[1]
        return QuadGK.quadgk(v -> ∫₁(v), zero(FP), one(FP); rule.kwargs...)[1]
    else
        _error_unsupported_combination("geometry with more than two parametric dimensions",
            "GaussKronrod")
    end
end

# GaussLegendre
function _integral(
        f,
        geometry,
        rule::GaussLegendre;
        FP::Type{T} = Float64,
        diff_method::DM = _default_diff_method(geometry)
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    N = Meshes.paramdim(geometry)

    # Get Gauss-Legendre nodes and weights of type FP for a region [-1,1]ᴺ
    xs, ws = FastGaussQuadrature.gausslegendre(N)
    xsT = Iterators.map(FP, xs)
    wsT = Iterators.map(FP, ws)
    weights = Iterators.product(ntuple(Returns(wsT), N)...)
    nodes = Iterators.product(ntuple(Returns(xsT), N)...)

    # Domain transformation: x [-1,1] ↦ t [0,1]
    t(x) = (1 // 2) * x + (1 // 2)

    function integrand((weights, nodes))
        ts = t.(nodes)
        prod(weights) * f(geometry(ts...)) * differential(geometry, ts, diff_method)
    end

    return (1 // (2^N)) .* sum(integrand, zip(weights, nodes))
end

# HAdaptiveCubature
function _integral(
        f,
        geometry,
        rule::HAdaptiveCubature;
        FP::Type{T} = Float64,
        diff_method::DM = _default_diff_method(geometry)
) where {DM <: DifferentiationMethod, T <: AbstractFloat}
    N = Meshes.paramdim(geometry)

    integrand(ts) = f(geometry(ts...)) * differential(geometry, ts, diff_method)

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = _zeros(FP, N)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(ts) = Unitful.ustrip.(integrandunits, integrand(ts))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, _zeros(FP, N), _ones(FP, N); rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end
