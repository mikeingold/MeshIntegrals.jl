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
- `f`: an integrand function with a method `f(::Meshes.Point)`
- `curve`: a `Meshes.BezierCurve` that defines the integration domain
- `rule = GaussKronrod()`: optionally, the `IntegrationRule` used for integration

# Keyword Arguments
- `diff_method::DifferentiationMethod = Analytical()`: the method to use for
calculating Jacobians that are used to calculate differential elements
- `FP = Float64`: the floating point precision desired
- `alg = Meshes.Horner()`:  the method to use for parameterizing `curve`. Alternatively,
`alg=Meshes.DeCasteljau()` can be specified for increased accuracy, but comes with a
steep performance cost, especially for curves with a large number of control points.
"""
function integral(
        f::Function,
        curve::Meshes.BezierCurve,
        rule::IntegrationRule;
        alg::Meshes.BezierEvalMethod = Meshes.Horner(),
        kwargs...
)
    paramfunction(t) = _parametric(curve, t; alg = alg)
    param_curve = _ParametricGeometry(paramfunction, Meshes.BezierCurve, 1)
    return _integral(f, param_curve, rule; kwargs...)
end

################################################################################
#                              Parametric
################################################################################

function _parametric(curve::Meshes.BezierCurve, t; alg::Meshes.BezierEvalMethod)
    return curve(t, alg)
end

################################################################################
#                               jacobian
################################################################################

function jacobian(
        curve::_ParametricGeometry{M, C, Meshes.BezierCurve, F, Dim},
        ts::Union{AbstractVector{T}, Tuple{T, Vararg{T}}},
        ::Analytical
) where {M, C, F, Dim, T <: AbstractFloat}
    t = only(ts)
    # Parameter t restricted to domain [0,1] by definition
    if t < 0 || t > 1
        throw(DomainError(t, "b(t) is not defined for t outside [0, 1]."))
    end

    # Aliases
    P = curve.controls
    N = Meshes.degree(curve)

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

function _has_analytical(
    ::Type{_ParametricGeometry{M, C, Meshes.BezierCurve, F, Dim}}
) where {M, C, F, Dim}
    return true
end
