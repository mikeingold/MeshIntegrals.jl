module MeshIntegralsEnzymeExt

using MeshIntegrals: MeshIntegrals, AutoEnzyme
using Meshes: Meshes
using Enzyme: Enzyme

function MeshIntegrals.jacobian(
        geometry::Meshes.Geometry,
        ts::Union{AbstractVector{T}, Tuple{T, Vararg{T}}},
        ::AutoEnzyme
) where {T <: AbstractFloat}
    return Meshes.to.(Enzyme.jacobian(Enzyme.Forward, geometry, ts...))
end

end
