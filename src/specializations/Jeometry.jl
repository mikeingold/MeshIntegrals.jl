# _Jeometry is a fake Meshes.Geometry whose only purpose is to provide a wrapper
# for manually-defined parametric functions.
struct _Jeometry{M <: Meshes.Manifold, C <: CRS, F <: Function} <: Meshes.Primitive{M, C}
    fun::F
    dims::Int64

    function _Jeometry{M, C}(
        fun::F,
        dims::Int64
    ) where {M <: Meshes.Manifold, C <: CRS, F <: Function}
        new{M, C, F}(fun, dims)
    end
end

function _Jeometry(
    fun::F,
    dims::Int64
) where {F <: Function}
    p = fun(zeros(dims))
    _Jeometry{Meshes.manifold(p), Meshes.crs(p)}(fun, dims)
end

(j::_Jeometry)(t::T) where {T <: AbstractFloat} = j.fun(t)
(j::_Jeometry)(t1::T, t2::T) where {T <: AbstractFloat} = j.fun(t1, t2)
(j::_Jeometry)(t1::T, t2::T, t3::T) where {T <: AbstractFloat} = j.fun(t1, t2, t3)

Meshes.paramdim(j::_Jeometry) = j.dims
Meshes.paramdim(::Type{<:_Jeometry}) = j.dims
