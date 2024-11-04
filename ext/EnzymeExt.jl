module EnzymeExt

using MeshIntegrals
using Meshes: Geometry
using Enzyme

function MeshIntegrals.jacobian(
        geometry::Geometry, ts::V, ::AutoEnzyme) where {V <: Union{AbstractVector, Tuple}}
    to.(Enzyme.jacobian(Enzyme.Forward, geometry, ts...))
end

end
