############################################################################################
#                               Specialized Methods for Plane
#
# Why Specialized?
#   The Plane geometry is a special case, representing a planar surface with infinite extent
#   along two basis vectors. This requires a pair of domain transformations mapping from the
#   typical parametric region [0,1]² to an infinite one (-∞,∞)².
############################################################################################

function integral(
        f::Function,
        plane::Meshes.Plane,
        rule::IntegrationRule;
        kwargs...
)
    paramfunction(t1, t2) = _parametric(plane, t1, t2)
    param_plane = _ParametricGeometry(paramfunction, 2)
    return _integral(f, param_plane, rule; kwargs...)
end

############################################################################################
#                                       Parametric
############################################################################################

function _parametric(plane::Meshes.Plane, t1, t2)
    f1(t) = t / (1 - t^2)
    f2(t) = 2t - 1
    f = f1 ∘ f2
    return plane(f(t1), f(t2))
end
