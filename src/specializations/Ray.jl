################################################################################
#                   Specialized Methods for Ray
#
# Why Specialized?
#   The Ray geometry is a special case, representing a ray, originating at a point
#   and extending an infinite length in a particular direction. This requires
#   a domain transformation mapping from the typical parametric region [0,1] to
#   an infinite one (-∞,∞).
################################################################################

function integral(
        f,
        ray::Meshes.Ray,
        rule::IntegrationRule;
        kwargs...
)
    # Generate a _ParametricGeometry whose parametric function spans the domain [0, 1]
    param_ray = _ParametricGeometry(_parametric(ray), Meshes.paramdim(ray))

    # Integrate the _ParametricGeometry using the standard methods
    return _integral(f, param_ray, rule; kwargs...)
end

################################################################################
#                                 Parametric
################################################################################

# Map argument domain from [0, 1] to [0, ∞) for (::Ray)(t)
# f(t) = t / 1 - t^2)
# f'(t) = (t^2 + 1) / (1 - t^2)^2
function _parametric(ray::Meshes.Ray)
    f(t) = t / (1 - t^2)
    return t -> ray(f(t))
end
