################################################################################
#                   Specialized Methods for FrustumSurface
################################################################################

function integral(
    f::F,
    frust::Meshes.FrustumSurface,
    settings::I,
    FP::Type{T} = Float64
) where  {F<:Function, I<:IntegrationAlgorithm, T<:AbstractFloat}
    # The generic method only parameterizes the sides
    sides = _integral(f, frust, settings, FP)

    # Integrate the Disks at the top and bottom
    top = _integral(f, frust.top, settings, FP)
    bottom = _integral(f, frust.bot, settings, FP)

    return sides + top + bottom
end
