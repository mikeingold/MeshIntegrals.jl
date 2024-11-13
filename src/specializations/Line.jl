################################################################################
#                   Specialized Methods for Line
#
# Why Specialized?
#   The Line geometry is a special case, representing a line of infinite length
#   that passes through two points. This requires another domain transformation
#   mapping from the typical parametric region [0,1] to an infinite one (-∞,∞).
################################################################################

function integral(
        f::Function,
        line::Meshes.Line,
        rule::IntegrationRule;
        kwargs...
)
    paramfunction(t) = _parametric(line, t)
    param_line = _ParametricGeometry(paramfunction, line, 1)
    return _integral(f, param_line, rule; kwargs...)
end

################################################################################
#                              Parametric
################################################################################

function _parametric(line::Meshes.Line, t)
    f1(t) = t / (1 - t^2)
    f2(t) = 2t - 1
    return line((f1 ∘ f2)(t))
end
