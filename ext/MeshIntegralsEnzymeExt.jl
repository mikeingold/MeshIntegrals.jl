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

# Supports all geometries except for a few
# See GitHub Issue #154 for more information.
supports_autoenzyme(::Type{<:Meshes.Geometry}) = true
supports_autoenzyme(::Type{<:Meshes.BezierCurve}) = false
supports_autoenzyme(::Type{<:Meshes.CylinderSurface}) = false
supports_autoenzyme(::Type{<:Meshes.Cylinder}) = false
supports_autoenzyme(::Type{<:Meshes.ParametrizedCurve}) = false

end
