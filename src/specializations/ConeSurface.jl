################################################################################
#                      Specialized Methods for ConeSurface
################################################################################

function integral(
        f::F,
        cone::Meshes.ConeSurface,
        settings::I,
        FP::Type{T} = Float64
) where {F <: Function, I <: IntegrationAlgorithm, T <: AbstractFloat}
    # The generic method only parameterizes the sides
    sides = _integral(f, cone, settings, FP)

    # Integrate the Disk at the base
    base = _integral(f, cone.base, settings, FP)

    return sides + base
end
