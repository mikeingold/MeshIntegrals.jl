module MeshIntegralsEnzymeExt

using MeshIntegrals: MeshIntegrals, AutoEnzyme
using Meshes: Meshes
using Enzyme: Enzyme

function MeshIntegrals.jacobian(
        geometry::Meshes.Geometry,
        ts::Union{AbstractVector{T}, Tuple{T, Vararg{T}}},
        ::AutoEnzyme
) where {T <: AbstractFloat}
    Dim = Meshes.paramdim(geometry)
    if Dim != length(ts)
        throw(ArgumentError("ts must have same number of dimensions as geometry."))
    end
    return Meshes.to.(Enzyme.jacobian(Enzyme.Forward, geometry, ts...))
end

end
