################################################################################
#                      Specialized Methods for ConeSurface
################################################################################

function integral(
    f::F,
    cone::Meshes.ConeSurface,
    settings::I,
    FP::Type{T} = Float64
) where {F<:Function, I<:IntegrationAlgorithm, T<:AbstractFloat}
    # The generic method only parameterizes the sides
    sides = _integral_2d(f, cone, settings, FP)

    # Integrate the Disk at the base
    base = _integral_2d(f, cone.base, settings, FP)

    return sides + base
end

function integral(
    f::F,
    cone::Meshes.ConeSurface,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # The generic method only parameterizes the sides
    sides = _integral(f, cone, settings, FP)

    # Integrate the Disk at the base
    base = _integral(f, cone.base, settings, FP)

    return sides + base
end
