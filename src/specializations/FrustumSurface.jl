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
        f::F,
        frust::Meshes.FrustumSurface,
        rule::I,
        FP::Type{T} = Float64
) where {F <: Function, I <: IntegrationRule, T <: AbstractFloat}
    # The generic method only parameterizes the sides
    sides = _integral(f, frust, rule, FP)

    # Integrate the Disks at the top and bottom
    top = _integral(f, frust.top, rule, FP)
    bottom = _integral(f, frust.bot, rule, FP)

    return sides + top + bottom
end
