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

Like [`integral`](@ref) but integrates along the domain defined by `curve`.

# Arguments
- `f`: an integrand function, i.e. any callable with a method `f(::Meshes.Point)`
- `curve`: a `Meshes.BezierCurve` that defines the integration domain
- `rule = GaussKronrod()`: optionally, the `IntegrationRule` used for integration

# Keyword Arguments
- `diff_method::DifferentiationMethod = Analytical()`: the method to use for
calculating Jacobians that are used to calculate differential elements
- `FP = Float64`: the floating point precision desired
- `alg = Meshes.Horner()`:  the method to use for parametrizing `curve`. Alternatively,
`alg=Meshes.DeCasteljau()` can be specified for increased accuracy, but comes with a
steep performance cost, especially for curves with a large number of control points.
"""
function integral(
        f,
        curve::Meshes.BezierCurve,
        rule::IntegrationRule;
        alg::Meshes.BezierEvalMethod = Meshes.Horner(),
        kwargs...
)
    # Generate a _ParametricGeometry whose parametric function auto-applies the alg kwarg
    param_curve = _ParametricGeometry(_parametric(curve, alg), Meshes.paramdim(curve))

    # Integrate the _ParametricGeometry using the standard methods
    return _integral(f, param_curve, rule; kwargs...)
end

################################################################################
#                              Parametric
################################################################################

# Wrap (::BezierCurve)(t, ::BezierEvalMethod) into f(t) by embedding second argument
function _parametric(curve::Meshes.BezierCurve, alg::Meshes.BezierEvalMethod)
    return t -> curve(t, alg)
end
