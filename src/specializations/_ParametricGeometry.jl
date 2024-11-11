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
struct _ParametricGeometry{M <: Meshes.Manifold, C <: CRS, F <: Function} <: Meshes.Primitive{M, C}
    fun::F
    dims::Int64

    function _ParametricGeometry{M, C}(
            fun::F,
            dims::Int64
    ) where {M <: Meshes.Manifold, C <: CRS, F <: Function}
        new{M, C, F}(fun, dims)
    end
end

function _ParametricGeometry(
        fun::F,
        dims::Int64
) where {F <: Function}
    p = fun(zeros(dims)...)
    _ParametricGeometry{Meshes.manifold(p), Meshes.crs(p)}(fun, dims)
end

(j::_ParametricGeometry)(t) = j.fun(t)
(j::_ParametricGeometry)(t1, t2) = j.fun(t1, t2)
(j::_ParametricGeometry)(t1, t2, t3) = j.fun(t1, t2, t3)

Meshes.paramdim(j::_ParametricGeometry) = j.dims
Meshes.paramdim(::Type{<:_ParametricGeometry}) = j.dims

"""
    _parametric(geometry, ts...)

Used for defining parametric functions for domain transformations.
"""
function _parametric end
