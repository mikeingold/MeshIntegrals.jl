################################################################################
#                  Specialized Methods for CylinderSurface
################################################################################

function integral(
    f::F,
    cyl::Meshes.CylinderSurface,
    settings::I,
    FP::Type{T} = Float64
) where {F<:Function, I<:IntegrationAlgorithm, T<:AbstractFloat}
    # The generic method only parameterizes the sides
    sides = _integral_2d(f, cyl, settings, FP)

    # Integrate the Disk at the top
    disk_top = Meshes.Disk(cyl.top, cyl.radius)
    top = _integral_2d(f, disk_top, settings, FP)

    # Integrate the Disk at the bottom
    disk_bottom = Meshes.Disk(cyl.bot, cyl.radius)
    bottom = _integral_2d(f, disk_bottom, settings, FP)

    return sides + top + bottom
end

function integral(
    f::F,
    cyl::Meshes.CylinderSurface,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # The generic method only parameterizes the sides
    sides = _integral(f, cyl, settings, FP)

    # Integrate the Disk at the top
    disk_top = Meshes.Disk(cyl.top, cyl.radius)
    top = _integral(f, disk_top, settings, FP)

    # Integrate the Disk at the bottom
    disk_bottom = Meshes.Disk(cyl.bot, cyl.radius)
    bottom = _integral(f, disk_bottom, settings, FP)

    return sides + top + bottom
end