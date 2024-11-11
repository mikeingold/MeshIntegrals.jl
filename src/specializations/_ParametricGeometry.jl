"""
    _ParametricGeometry <: Meshes.Geometry

_ParametricGeometry is used internally in MeshIntegrals.jl to behave like a generic wrapper
for geometries with custom parametric functions. This type is used for transforming other
geometries to enable integration over the standard rectangular `[0,1]^n` domain.

Meshes.jl adopted a ParametrizedCurve type that performs a similar role as of `v0.51.20`,
but only supports geometries with one parametric dimension. Support is additionally planned
for more types that span surfaces and volumes, at which time this custom type will probably
no longer be required.

# Fields
- `fun::Function` - a parametric function: (ts...) -> Meshes.Point
- `dims::Int64` - the geometry's number of parametric dimensions
"""
struct _ParametricGeometry{M <: Meshes.Manifold, C <: CRS, F <: Function, Dim} <:
       Meshes.Primitive{M, C}
    fun::F

    function _ParametricGeometry{M, C}(
            fun::F,
            dims::Int64
    ) where {M <: Meshes.Manifold, C <: CRS, F <: Function}
        return new{M, C, F, dims}(fun)
    end
end

function _ParametricGeometry(
        fun::F,
        dims::Int64
) where {F <: Function}
    p = fun(_zeros(dims)...)
    return _ParametricGeometry{Meshes.manifold(p), Meshes.crs(p)}(fun, dims)
end

# geometry(ts...) == geometry.fun(ts...)
(g::_ParametricGeometry{M, C, F, 1})(t) where {M, C, F} = g.fun(t)
(g::_ParametricGeometry{M, C, F, 2})(t1, t2) where {M, C, F} = g.fun(t1, t2)
(g::_ParametricGeometry{M, C, F, 3})(t1, t2, t3) where {M, C, F} = g.fun(t1, t2, t3)

Meshes.paramdim(::_ParametricGeometry{M, C, F, Dim}) where {M, C, F, Dim} = Dim

"""
    _parametric(geometry::G, ts...) where {G <: Meshes.Geometry}

Used in MeshIntegrals.jl for defining parametric functions that implement domain
transformations. Usages are specialized on 
"""
function _parametric end
