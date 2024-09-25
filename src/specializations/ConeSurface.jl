################################################################################
#                      Specialized Methods for ConeSurface
#
# Why Specialized?
#   The parametric function that Meshes.jl currently implements for ConeSurface
#   only parameterizes the rounded walls of the cone, but this Geometry surface is
#   defined as including the circular base surface as well. These methods simply
#   integrate both the base and walls and return the sum of the two integrals.
################################################################################

function integral(
        f::F,
        cone::Meshes.ConeSurface,
        rule::I,
        FP::Type{T} = Float64
) where {F <: Function, I <: IntegrationRule, T <: AbstractFloat}
    # The generic method only parameterizes the sides
    sides = _integral(f, cone, rule, FP)

    # Integrate the Disk at the base
    base = _integral(f, cone.base, rule, FP)

    return sides + base
end
