module MeshIntegralsEnzymeExt

using MeshIntegrals: MeshIntegrals, AutoEnzyme
using Meshes: Geometry, to
using Enzyme: Enzyme

function Meshes.jacobian(
        geometry::Geometry,
        ts::Union{AbstractVector{T}, Tuple{T, Vararg{T}}},
        ::AutoEnzyme
) where {T <: AbstractFloat}
    return to.(Enzyme.jacobian(Enzyme.Forward, geometry, ts...))
end

end
