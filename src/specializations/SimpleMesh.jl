################################################################################
#                      Specialized Methods for SimpleMesh
#
# Why Specialized?
#   The SimpleMesh <: Meshes.Domain type defines a meshed surface.
################################################################################

"""
    integral(f, mesh::SimpleMesh, rule = HAdaptiveCubature(); kwargs...)

Like [`integral`](@ref) but integrates throughout the surface defined by `mesh`.
The mesh is first partitioned into surface components, often `Triangle` or `Quadrangle`
primitives, and then the specified integration `rule` is applied independently to each.

# Arguments
- `f`: an integrand function, i.e. any callable with a method `f(::Meshes.Point)`
- `mesh`: a `SimpleMesh` that defines the integration domain
- `rule = HAdaptiveCubature()`: optionally, the `IntegrationRule` used for integration

# Keyword Arguments
See [`integral`](@ref)
"""
function integral(
        f,
        mesh::Meshes.SimpleMesh,
        rule::I = HAdaptiveCubature();
        kwargs...
) where {I <: IntegrationRule}
    # Partition the mesh and sum the integrals
    return sum(geometry -> _integral(f, geometry, rule; kwargs...), collect(mesh))
end
