################################################################################
#                   Specialized Methods for BezierCurve
#
# Why Specialized?
#   The parametric function in Meshes.jl for BezierCurve accepts an argument
#   of type Meshes.BezierEvalMethod, in which the method for generating
#   parametric points along the curve is specified. These specialized methods
#   are essentially identical to the generic ones, except that they provide a
#   keyword argument and pass through the specified algorithm choice.
################################################################################

################################################################################
#                              integral
################################################################################
"""
    integral(f, curve::BezierCurve, rule = GaussKronrod();
             diff_method=Analytical(), FP=Float64, alg=Meshes.Horner())

Like [`integral`](@ref) but integrates along the domain defined a `curve`.

# Arguments
- `f`: an integrand function with a method `f(::Meshes.Point)`
- `curve`: a `BezierCurve` that defines the integration domain
- `rule`: optionally, the `IntegrationRule` used for integration (by default
`GaussKronrod()` in 1D and `HAdaptiveCubature()` else)

# Keyword Arguments
- `diff_method::DifferentiationMethod = Analytical()`: the method to use for
calculating Jacobians that are used to calculate differential elements
- `FP = Float64`: the floating point precision desired
- `alg = Meshes.Horner()`:  the method to use for parameterizing `curve`. Alternatively,
`alg=Meshes.DeCasteljau()` can be specified for increased accuracy, but comes with a
steep performance cost, especially for curves with a large number of control points.
"""
function integral(
        f::F,
        curve::Meshes.BezierCurve,
        rule::GaussLegendre;
        diff_method::DM = default_method(curve),
        FP::Type{T} = Float64,
        alg::Meshes.BezierEvalMethod = Meshes.Horner()
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, rule.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = FP(1 // 2) * x + FP(1 // 2)
    point(x) = curve(t(x), alg)
    integrand(x) = f(point(x)) * differential(curve, (t(x),), diff_method)

    # Integrate f along curve and apply domain-correction for [-1,1] ↦ [0, length]
    return FP(1 // 2) * sum(w .* integrand(x) for (w, x) in zip(ws, xs))
end

function integral(
        f::F,
        curve::Meshes.BezierCurve,
        rule::GaussKronrod;
        diff_method::DM = default_method(curve),
        FP::Type{T} = Float64,
        alg::Meshes.BezierEvalMethod = Meshes.Horner()
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    point(t) = curve(t, alg)
    integrand(t) = f(point(t)) * differential(curve, (t,), diff_method)
    return QuadGK.quadgk(integrand, zero(FP), one(FP); rule.kwargs...)[1]
end

function integral(
        f::F,
        curve::Meshes.BezierCurve,
        rule::HAdaptiveCubature;
        diff_method::DM = default_method(curve),
        FP::Type{T} = Float64,
        alg::Meshes.BezierEvalMethod = Meshes.Horner()
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    point(t) = curve(t, alg)
    integrand(ts) = f(point(only(ts))) * differential(curve, ts, diff_method)

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = zeros(FP, 1)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(ts) = Unitful.ustrip.(integrandunits, integrand(ts))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, zeros(FP, 1), ones(FP, 1); rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end

################################################################################
#                               jacobian
################################################################################

function jacobian(
        bz::Meshes.BezierCurve,
        ts::V,
        diff_method::Analytical
) where {V <: Union{AbstractVector, Tuple}}
    t = only(ts)
    # Parameter t restricted to domain [0,1] by definition
    if t < 0 || t > 1
        throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
    end

    # Aliases
    P = bz.controls
    N = Meshes.degree(bz)

    # Ensure that this implementation is tractible: limited by ability to calculate
    #   binomial(N, N/2) without overflow. It's possible to extend this range by
    #   converting N to a BigInt, but this results in always returning BigFloat's.
    N <= 1028 || error("This algorithm overflows for curves with ⪆1000 control points.")

    # Generator for Bernstein polynomial functions
    B(i, n) = t -> binomial(Int128(n), i) * t^i * (1 - t)^(n - i)

    # Derivative = N Σ_{i=0}^{N-1} sigma(i)
    #   P indices adjusted for Julia 1-based array indexing
    sigma(i) = B(i, N - 1)(t) .* (P[(i + 1) + 1] - P[(i) + 1])
    derivative = N .* sum(sigma, 0:(N - 1))

    return (derivative,)
end

has_analytical(::Type{Meshes.BezierCurve}) = true
