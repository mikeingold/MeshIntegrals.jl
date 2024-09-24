################################################################################
#                      Specialized Methods for ConeSurface
################################################################################

function integral(
    f::F,
    cone::Meshes.ConeSurface,
    rule::I,
    FP::Type{T} = Float64
) where {F<:Function, I<:IntegrationRule, T<:AbstractFloat}
    # The generic method only parameterizes the sides
    sides = _integral(f, cone, rule, FP)

    # Integrate the Disk at the base
    base = _integral(f, cone.base, rule, FP)

    return sides + base
end
