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
    # Generate a _ParametricGeometry whose parametric function spans the domain [0, 1]
    param_line = _ParametricGeometry(_parametric(line), Meshes.paramdim(line))

    # Integrate the _ParametricGeometry using the standard methods
    return _integral(f, param_line, rule; kwargs...)
end

################################################################################
#                              Parametric
################################################################################

# Map argument domain from [0, 1] to (-∞, ∞) for (::Line)(t)
function _parametric(line::Meshes.Line)
    # [-1, 1] ↦ (-∞, ∞)
    f1(t) = t / (1 - t^2)
    # [0, 1] ↦ [-1, 1]
    f2(t) = 2t - 1
    # Compose the two transforms
    return t -> line((f1 ∘ f2)(t))
end
