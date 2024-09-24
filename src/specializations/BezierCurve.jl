################################################################################
#                   Specialized Methods for BezierCurve
################################################################################

function lineintegral(
        f::F,
        curve::Meshes.BezierCurve,
        settings::I,
        FP::Type{T};
        alg::Meshes.BezierEvalMethod = Meshes.Horner()
) where {F <: Function, I <: IntegrationAlgorithm, T <: AbstractFloat}
    return integral(f, curve, settings, FP; alg = alg)
end

function lineintegral(
        f::F,
        curve::Meshes.BezierCurve,
        settings::I;
        alg::Meshes.BezierEvalMethod = Meshes.Horner()
) where {F <: Function, I <: IntegrationAlgorithm}
    return integral(f, curve, settings; alg = alg)
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
        alg::Meshes.BezierEvalMethod = Meshes.Horner()
) where {F <: Function, T <: AbstractFloat}
    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = FP(1 / 2) * x + FP(1 / 2)
    point(x) = curve(t(x), alg)

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return FP(1 / 2) * length(curve) * sum(w .* f(point(x)) for (w, x) in zip(ws, xs))
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
        alg::Meshes.BezierEvalMethod = Meshes.Horner()
) where {F <: Function, T <: AbstractFloat}
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
        alg::Meshes.BezierEvalMethod = Meshes.Horner()
) where {F <: Function, T <: AbstractFloat}
    len = length(curve)
    point(t) = curve(t, alg)
    integrand(t) = len * f(point(t[1]))

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = fill(FP(0.5), 3)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, FP[0], FP[1]; settings.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end
