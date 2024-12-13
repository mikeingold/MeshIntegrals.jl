################################################################################
#                  Specialized Methods for CylinderSurface
#
# Why Specialized?
#   The parametric function that Meshes.jl currently implements for CylinderSurface
#   only parametrizes the rounded walls, but this Geometry surface is defined as
#   including the top and bottom circular base surfaces as well. These methods
#   simply integrate the base and walls and return the sum of the three integrals.
################################################################################

function integral(
        f,
        cyl::Meshes.CylinderSurface,
        rule::I;
        diff_method::DM = _default_diff_method(cyl),
        kwargs...
) where {I <: IntegrationRule, DM <: DifferentiationMethod}
    _check_diff_method_support(cyl, diff_method)

    # The generic method only parametrizes the sides
    sides = _integral(f, cyl, rule; diff_method = diff_method, kwargs...)

    # Integrate the Disk at the top
    disk_top = Meshes.Disk(cyl.top, cyl.radius)
    top = _integral(f, disk_top, rule; diff_method = diff_method, kwargs...)

    # Integrate the Disk at the bottom
    disk_bottom = Meshes.Disk(cyl.bot, cyl.radius)
    bottom = _integral(f, disk_bottom, rule; diff_method = diff_method, kwargs...)

    return sides + top + bottom
end
