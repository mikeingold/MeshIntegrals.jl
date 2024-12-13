################################################################################
#                          DifferentiationMethods
################################################################################

"""
    DifferentiationMethod

A category of types used to specify the desired method for calculating derivatives.
Derivatives are used to form Jacobian matrices when calculating the differential
element size throughout the integration region.

See also [`FiniteDifference`](@ref), [`AutoEnzyme`](@ref).
"""
abstract type DifferentiationMethod end

"""
    FiniteDifference(ε=1e-6)

Use to specify use of a finite-difference approximation method with a step size
of `ε` for calculating derivatives.
"""
struct FiniteDifference{T <: AbstractFloat} <: DifferentiationMethod
    ε::T
end

# If ε not specified, default to 1e-6
FiniteDifference() = FiniteDifference(1e-6)

"""
    AutoEnzyme()

Use to specify use of the Enzyme.jl for calculating derivatives.
"""
struct AutoEnzyme <: DifferentiationMethod end

# Future Support:
#   struct AutoZygote <: DifferentiationMethod end

################################################################################
#                                  Jacobian
################################################################################

"""
    jacobian(geometry, ts[, diff_method])

Calculate the Jacobian of a `geometry`'s parametric function at some point `ts`.
Optionally, direct the use of a particular differentiation method `diff_method`; by default
use analytic solutions where possible and finite difference approximations
otherwise.

# Arguments
- `geometry`: some `Meshes.Geometry` of N parametric dimensions
- `ts`: a parametric point specified as a vector or tuple of length N
- `diff_method`: the desired [`DifferentiationMethod`](@ref) to use
"""
function jacobian(
        geometry::G,
        ts::Union{AbstractVector{T}, Tuple{T, Vararg{T}}}
) where {G <: Geometry, T <: AbstractFloat}
    return jacobian(geometry, ts, _default_diff_method(G))
end

function jacobian(
        geometry::Geometry,
        ts::Union{AbstractVector{T}, Tuple{T, Vararg{T}}},
        diff_method::FiniteDifference
) where {T <: AbstractFloat}
    Dim = Meshes.paramdim(geometry)
    if Dim != length(ts)
        throw(ArgumentError("ts must have same number of dimensions as geometry."))
    end

    # Get the partial derivative along the n'th axis via finite difference
    #   approximation, where ts is the current parametric position
    function ∂ₙr(ts, n, ε)
        # Build left/right parametric coordinates with non-allocating iterators
        left = Iterators.map(((i, t),) -> i == n ? t - ε : t, enumerate(ts))
        right = Iterators.map(((i, t),) -> i == n ? t + ε : t, enumerate(ts))
        # Select orientation of finite-diff
        if ts[n] < 0.01
            # Right
            return (geometry(right...) - geometry(ts...)) / ε
        elseif 0.99 < ts[n]
            # Left
            return (geometry(ts...) - geometry(left...)) / ε
        else
            # Central
            return (geometry(right...) - geometry(left...)) / 2ε
        end
    end

    return ntuple(n -> ∂ₙr(ts, n, T(diff_method.ε)), Dim)
end

################################################################################
#                          Differential Elements
################################################################################

"""
    differential(geometry, ts[, diff_method])

Calculate the differential element (length, area, volume, etc) of the parametric
function for `geometry` at arguments `ts`. Optionally, direct the use of a
particular differentiation method `diff_method`; by default use analytic solutions where
possible and finite difference approximations otherwise.

# Arguments
- `geometry`: some `Meshes.Geometry` of N parametric dimensions
- `ts`: a parametric point specified as a vector or tuple of length N
- `diff_method`: the desired [`DifferentiationMethod`](@ref) to use
"""
function differential(
        geometry::G,
        ts::Union{AbstractVector{T}, Tuple{T, Vararg{T}}},
        diff_method::DifferentiationMethod = _default_diff_method(G)
) where {G <: Geometry, T <: AbstractFloat}
    J = Iterators.map(_KVector, jacobian(geometry, ts, diff_method))
    return LinearAlgebra.norm(foldl(∧, J))
end
