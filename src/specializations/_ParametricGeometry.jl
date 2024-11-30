"""
    _ParametricGeometry <: Meshes.Primitive <: Meshes.Geometry

_ParametricGeometry is used internally in MeshIntegrals.jl to behave like a generic wrapper
for geometries with custom parametric functions. This type is used for transforming other
geometries to enable integration over the standard rectangular `[0,1]^n` domain.

Meshes.jl adopted a `ParametrizedCurve` type that performs a similar role as of `v0.51.20`,
but only supports geometries with one parametric dimension. Support is additionally planned
for more types that span surfaces and volumes, at which time this custom type will probably
no longer be required.

# Fields
- `fun::Function` - a parametric function: (ts...) -> Meshes.Point

# Type Structure
- `M <: Meshes.Manifold` - same usage as in `Meshes.Geometry{M, C}`
- `C <: CoordRefSystems.CRS` - same usage as in `Meshes.Geometry{M, C}`
- `F` - type of the callable integrand function
- `Dim` - number of parametric dimensions
"""
struct _ParametricGeometry{M <: Meshes.Manifold, C <: CRS, F, Dim} <:
       Meshes.Primitive{M, C}
    fun::F

    function _ParametricGeometry{M, C}(
            fun::F,
            Dim::Int64
    ) where {M <: Meshes.Manifold, C <: CRS, F}
        return new{M, C, F, Dim}(fun)
    end
end

"""
    _ParametricGeometry(fun, dims)

Construct a `_ParametricGeometry` using a provided parametric function `fun` for a geometry
with `dims` parametric dimensions.

# Arguments
- `fun::Function` - parametric function mapping `(ts...) -> Meshes.Point`
- `dims::Int64` - number of parametric dimensions, i.e. `length(ts)`
"""
function _ParametricGeometry(
        fun,
        Dim::Int64
)
    p = fun(_zeros(Dim)...)
    return _ParametricGeometry{Meshes.manifold(p), Meshes.crs(p)}(fun, Dim)
end

# Allow a _ParametricGeometry to be called like a Geometry
(g::_ParametricGeometry)(ts...) = g.fun(ts...)

Meshes.paramdim(::_ParametricGeometry{M, C, F, Dim}) where {M, C, F, Dim} = Dim

"""
    _parametric(geometry::G, ts...) where {G <: Meshes.Geometry}

Used in MeshIntegrals.jl for defining parametric functions that transform non-standard
geometries into a form that can be integrated over the standard rectangular [0,1]^n limits.
"""
function _parametric end
