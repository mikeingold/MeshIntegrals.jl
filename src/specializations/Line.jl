################################################################################
#                   Specialized Methods for Line
#
# Why Specialized?
#   The Line geometry is a special case, representing a line of infinite length
#   that passes through two points. This requires another domain transformation
#   mapping from the typical parametric region [0,1] to an infinite one (-∞,∞).
################################################################################

function integral(
        f,
        line::Meshes.Line,
        rule::IntegrationRule;
        kwargs...
)
    paramfunction(t) = _parametric(line, t)
    param_line = _ParametricGeometry(paramfunction, 1)
    return _integral(f, param_line, rule; kwargs...)
end

################################################################################
#                              Parametric
################################################################################

# Map [0, 1] ↦ (-∞, ∞)
function _parametric(line::Meshes.Line, t)
    # [-1, 1] ↦ (-∞, ∞)
    f1(t) = t / (1 - t^2)
    # [0, 1] ↦ [-1, 1]
    f2(t) = 2t - 1
    # Compose the two transforms
    return line((f1 ∘ f2)(t))
end
