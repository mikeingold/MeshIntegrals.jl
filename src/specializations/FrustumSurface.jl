################################################################################
#                   Specialized Methods for FrustumSurface
#
# Why Specialized?
#   The parametric function that Meshes.jl currently implements for FrustumSurface
#   only parameterizes the rounded walls, but this Geometry surface is defined as
#   including the truncated top and bottom surfaces as well. These methods simply
#   integrate both the walls and the ends and return the sum of the these integrals.
################################################################################

function integral(
        f,
        frust::Meshes.FrustumSurface,
        rule::I;
        kwargs...
) where {I <: IntegrationRule}
    # The generic method only parameterizes the sides
    sides = _integral(f, frust, rule; kwargs...)

    # Integrate the Disks at the top and bottom
    top = _integral(f, frust.top, rule; kwargs...)
    bottom = _integral(f, frust.bot, rule; kwargs...)

    return sides + top + bottom
end
