################################################################################
#                  Specialized Methods for CylinderSurface
#
# Why Specialized?
#   The parametric function that Meshes.jl currently implements for CylinderSurface
#   only parameterizes the rounded walls, but this Geometry surface is defined as
#   including the top and bottom circular base surfaces as well. These methods
#   simply integrate the base and walls and return the sum of the three integrals.
################################################################################

function integral(
        f::F,
        cyl::Meshes.CylinderSurface,
        rule::I,
        FP::Type{T} = Float64
) where {F <: Function, I <: IntegrationRule, T <: AbstractFloat}
    # The generic method only parameterizes the sides
    sides = _integral(f, cyl, rule, FP)

    # Integrate the Disk at the top
    disk_top = Meshes.Disk(cyl.top, cyl.radius)
    top = _integral(f, disk_top, rule, FP)

    # Integrate the Disk at the bottom
    disk_bottom = Meshes.Disk(cyl.bot, cyl.radius)
    bottom = _integral(f, disk_bottom, rule, FP)

    return sides + top + bottom
end
