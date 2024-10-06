################################################################################
#                         Master Integral Function
################################################################################

"""
    integral(f, geometry[, rule]; FP=Float64)

Numerically integrate a given function `f(::Point)` over the domain defined by
a `geometry` using a particular numerical integration `rule` with floating point
precision of type `FP`.

# Arguments
- `f`: an integrand function with a method `f(::Meshes.Point)`
- `geometry`: some `Meshes.Geometry` that defines the integration domain
- `rule`: optionally, the `IntegrationRule` used for integration (by default `GaussKronrod()` in 1D and `HAdaptiveCubature()` else)
- `FP`: optionally, the floating point precision desired (`Float64` by default)

Note that reducing `FP` below `Float64` will incur some loss of precision. By
contrast, increasing `FP` to e.g. `BigFloat` will typically increase precision
(at the expense of longer runtimes).
"""
function integral end

# If only f and geometry are specified, select default rule
function integral(
        f::F,
        geometry::G,
        rule::I = Meshes.paramdim(geometry) == 1 ? GaussKronrod() : HAdaptiveCubature();
        kwargs...
) where {F <: Function, G <: Meshes.Geometry, I <: IntegrationRule}
    _integral(f, geometry, rule; kwargs...)
end

################################################################################
#                    Generalized (n-Dimensional) Worker Methods
################################################################################

# GaussKronrod
function _integral(
        f,
        geometry,
        rule::GaussKronrod;
        kwargs...
)
    # Pass through to dim-specific workers in next section of this file
    N = Meshes.paramdim(geometry)
    if N == 1
        return _integral_gk_1d(f, geometry, rule; kwargs...)
    elseif N == 2
        return _integral_gk_2d(f, geometry, rule; kwargs...)
    elseif N == 3
        return _integral_gk_3d(f, geometry, rule; kwargs...)
    else
        error("Integrating this geometry type with GaussKronrod not supported.")
    end
end

# GaussLegendre
function _integral(
        f,
        geometry,
        rule::GaussLegendre;
        FP::Type{T} = Float64
) where {T <: AbstractFloat}
    N = Meshes.paramdim(geometry)

    # Get Gauss-Legendre nodes and weights for a region [-1,1]^N
    xs, ws = _gausslegendre(FP, rule.n)
    weights = Iterators.product(ntuple(Returns(ws), N)...)
    nodes = Iterators.product(ntuple(Returns(xs), N)...)

    # Domain transformation: x [-1,1] ↦ t [0,1]
    t(x) = FP(1 // 2) * x + FP(1 // 2)

    function integrand((weights, nodes))
        ts = t.(nodes)
        prod(weights) * f(geometry(ts...)) * differential(geometry, ts)
    end

    return FP(1 // (2^N)) .* sum(integrand, zip(weights, nodes))
end

# HAdaptiveCubature
function _integral(
        f,
        geometry,
        rule::HAdaptiveCubature;
        FP::Type{T} = Float64
) where {T <: AbstractFloat}
    N = Meshes.paramdim(geometry)

    integrand(ts) = f(geometry(ts...)) * differential(geometry, ts)

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = zeros(FP, N)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(ts) = Unitful.ustrip.(integrandunits, integrand(ts))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, zeros(FP, N), ones(FP, N); rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end

################################################################################
#                    Specialized GaussKronrod Methods
################################################################################

function _integral_gk_1d(
        f,
        geometry,
        rule::GaussKronrod;
        FP::Type{T} = Float64
) where {T <: AbstractFloat}
    integrand(t) = f(geometry(t)) * differential(geometry, (t,))
    return QuadGK.quadgk(integrand, zero(FP), one(FP); rule.kwargs...)[1]
end

function _integral_gk_2d(
        f,
        geometry2d,
        rule::GaussKronrod;
        FP::Type{T} = Float64
) where {T <: AbstractFloat}
    integrand(u, v) = f(geometry2d(u, v)) * differential(geometry2d, (u, v))
    ∫₁(v) = QuadGK.quadgk(u -> integrand(u, v), zero(FP), one(FP); rule.kwargs...)[1]
    return QuadGK.quadgk(v -> ∫₁(v), zero(FP), one(FP); rule.kwargs...)[1]
end

# Integrating volumes with GaussKronrod not supported by default
function _integral_gk_3d(
        f,
        geometry,
        rule::GaussKronrod;
        FP::Type{T} = Float64
) where {T <: AbstractFloat}
    error("Integrating this volume type with GaussKronrod not supported.")
end
