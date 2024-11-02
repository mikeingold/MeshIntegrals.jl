################################################################################
#                          DifferentiationMethods
################################################################################

abstract type DifferentiationMethod end

"""
    FiniteDifference(ε=1e-6)

Use a finite-difference approximation method to calculate derivatives with a
step size of `ε`.
"""
struct FiniteDifference{T <: AbstractFloat} <: DifferentiationMethod
    ε::T
end

# If ε not specified, default to 1e-6
FiniteDifference() = FiniteDifference(1e-6)

struct Analytical <: DifferentiationMethod end

# Future Support: (maybe define in package extensions?)
#   struct AutoEnzyme <: DifferentiationMethod end
#   struct AutoZygote <: DifferentiationMethod end

################################################################################
#                                  Jacobian
################################################################################

"""
    jacobian(geometry, ts[, diff_method])

Calculate the Jacobian of a `geometry`'s parametric function at some point `ts`.
Optionally, direct the use of a particular `differentiation method`; by default
use analytic solutions where possible and finite difference approximations
otherwise.

# Arguments
- `geometry`: some `Meshes.Geometry` of N parametric dimensions
- `ts`: a parametric point specified as a vector or tuple of length N
- `diff_method`: the desired `DifferentiationMethod` to use
"""
function jacobian end

function jacobian(
        geometry::G,
        ts::V
) where {G <: Geometry, V <: Union{AbstractVector, Tuple}}
    return jacobian(geometry, ts, _default_method(G))
end

function jacobian(
        geometry::Geometry,
        ts::V,
        diff_method::FiniteDifference
) where {V <: Union{AbstractVector, Tuple}}
    Dim = Meshes.paramdim(geometry)
    if Dim != length(ts)
        throw(ArgumentError("ts must have same number of dimensions as geometry."))
    end

    T = eltype(ts)
    ε = T(diff_method.ε)

    # Get the partial derivative along the n'th axis via finite difference
    #   approximation, where ts is the current parametric position
    function ∂ₙr(ts, n)
        # Build left/right parametric coordinates with non-allocating iterators 
        left = Iterators.map(((i, t),) -> i == n ? t - ε : t, enumerate(ts))
        right = Iterators.map(((i, t),) -> i == n ? t + ε : t, enumerate(ts))
        # Select orientation of finite-diff
        if ts[n] < T(0.01)
            # Right
            return (geometry(right...) - geometry(ts...)) / ε
        elseif T(0.99) < ts[n]
            # Left
            return (geometry(ts...) - geometry(left...)) / ε
        else
            # Central
            return (geometry(right...) - geometry(left...)) / 2ε
        end
    end

    return ntuple(n -> ∂ₙr(ts, n), Dim)
end

################################################################################
#                          Differential Elements
################################################################################

"""
    differential(geometry, ts[, diff_method])

Calculate the differential element (length, area, volume, etc) of the parametric
function for `geometry` at arguments `ts`. Optionally, direct the use of a
particular `differentiation method`; by default use analytic solutions where
possible and finite difference approximations otherwise.

# Arguments
- `geometry`: some `Meshes.Geometry` of N parametric dimensions
- `ts`: a parametric point specified as a vector or tuple of length N
- `diff_method`: the desired `DifferentiationMethod` to use
"""
function differential(
        geometry::G,
        ts::V,
        diff_method::DifferentiationMethod = default_method(G)
) where {G <: Geometry, V <: Union{AbstractVector, Tuple}}
    # Calculate the Jacobian, convert Vec -> KVector
    J = jacobian(geometry, ts, diff_method)
    J_kvecs = Iterators.map(_kvector, J)

    # Extract units from Geometry type
    Dim = Meshes.paramdim(geometry)
    units = _units(geometry)^Dim

    # Return norm of the exterior products
    element = foldl(∧, J_kvecs)
    return LinearAlgebra.norm(element) * units
end
