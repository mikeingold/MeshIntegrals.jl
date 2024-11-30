################################################################################
#                      Specialized Methods for Plane
#
# Why Specialized?
#   The Plane geometry is a special case, representing a planar surface with
#   infinite extent along two basis vectors. This requires a pair of domain
#   transformations mapping from the typical parametric region [0,1]² to an
#   infinite one (-∞,∞)².
################################################################################

function integral(
        f,
        plane::Meshes.Plane,
        rule::IntegrationRule;
        kwargs...
)
    # Generate a _ParametricGeometry whose parametric function spans the domain [0, 1]²
    paramfunction(t1, t2) = _parametric(plane, t1, t2)

    # Integrate the _ParametricGeometry using the standard methods
    return _integral(f, param_plane, rule; kwargs...)
end

############################################################################################
#                                       Parametric
############################################################################################

# Map argument domain from [0, 1]² to (-∞, ∞)² for (::Plane)(t1, t2)
function _parametric(plane::Meshes.Plane, t1, t2)
    # [-1, 1] ↦ (-∞, ∞)
    f1(t) = t / (1 - t^2)
    # [0, 1] ↦ [-1, 1]
    f2(t) = 2t - 1
    # Compose the two transforms
    f = f1 ∘ f2
    return plane(f(t1), f(t2))
end
