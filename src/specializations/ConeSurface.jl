################################################################################
#                      Specialized Methods for ConeSurface
#
# Why Specialized?
#   The parametric function that Meshes.jl currently implements for ConeSurface
#   only parametrizes the rounded walls of the cone, but this Geometry surface is
#   defined as including the circular base surface as well. These methods simply
#   integrate both the base and walls and return the sum of the two integrals.
################################################################################

function integral(
        f,
        cone::Meshes.ConeSurface,
        rule::I;
        kwargs...
) where {I <: IntegrationRule}
    # The generic method only parametrizes the sides
    sides = _integral(f, cone, rule; kwargs...)

    # Integrate the Disk at the base
    base = _integral(f, cone.base, rule; kwargs...)

    return sides + base
end
