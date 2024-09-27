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

"""
    integral(f, curve::BezierCurve, ::GaussLegendre;
             FP=Float64, alg=Meshes.Horner())

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
        rule::GaussLegendre;
        FP::Type{T} = Float64,
        alg::Meshes.BezierEvalMethod = Meshes.Horner()
) where {F <: Function, T <: AbstractFloat}
    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, rule.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = FP(1 // 2) * x + FP(1 // 2)
    point(x) = curve(t(x), alg)
    integrand(x) = f(point(x)) * norm(derivative(curve, t(x)))

    # Integrate f along curve and apply domain-correction for [-1,1] ↦ [0, length]
    return FP(1 // 2) * sum(w .* integrand(x) for (w, x) in zip(ws, xs))
end

"""
    integral(f, curve::BezierCurve, ::GaussKronrod;
             FP=Float64, alg=Meshes.Horner())

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
        rule::GaussKronrod;
        FP::Type{T} = Float64,
        alg::Meshes.BezierEvalMethod = Meshes.Horner()
) where {F <: Function, T <: AbstractFloat}
    point(t) = curve(t, alg)
    integrand(t) = f(point(t)) * norm(derivative(curve, t))
    return QuadGK.quadgk(integrand, zero(FP), one(FP); rule.kwargs...)[1]
end

"""
    integral(f, curve::BezierCurve, ::HAdaptiveCubature;
             FP=Float64, alg=Meshes.Horner())

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
        rule::HAdaptiveCubature;
        FP::Type{T} = Float64,
        alg::Meshes.BezierEvalMethod = Meshes.Horner()
) where {F <: Function, T <: AbstractFloat}
    point(t) = curve(t, alg)
    integrand(t) = f(point(only(t))) *  norm(derivative(curve, only(t)))

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = zeros(FP, 3)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, zeros(FP, 1), ones(FP, 1); rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end
