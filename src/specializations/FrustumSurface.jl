################################################################################
#                   Specialized Methods for FrustumSurface
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
