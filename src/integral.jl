################################################################################
#                         Integration Algorithms
################################################################################

abstract type IntegrationAlgorithm end

"""
    GaussKronrod(kwargs...)

Numerically integrate using the h-adaptive Gauss-Kronrod quadrature rule implemented
by QuadGK.jl. All standard `QuadGK.quadgk` keyword arguments are supported.
"""
struct GaussKronrod <: IntegrationAlgorithm
    kwargs
    GaussKronrod(; kwargs...) = new(kwargs)
end

"""
    GaussLegendre(n)

Numerically integrate using an `n`'th-order Gauss-Legendre quadrature rule. nodes
and weights are efficiently calculated using FastGaussQuadrature.jl.

So long as the integrand function can be well-approximated by a polynomial of
order `2n-1`, this method should yield results with 16-digit accuracy in `O(n)`
time. If the function is know to have some periodic content, then `n` should
(at a minimum) be greater than the expected number of periods over the geometry,
e.g. `length(geometry)/lambda`.
"""
struct GaussLegendre <: IntegrationAlgorithm
    n::Int64
end

"""
    GaussKronrod(kwargs...)

Numerically integrate areas and surfaces using the h-adaptive cubature rule
implemented by HCubature.jl. All standard `HCubature.hcubature` keyword
arguments are supported.
"""
struct HAdaptiveCubature <: IntegrationAlgorithm
    kwargs
    HAdaptiveCubature(; kwargs...) = new(kwargs)
end

################################################################################
#                         Master Integral Function
################################################################################

"""
    integral(f, geometry)
    integral(f, geometry, algorithm)
    integral(f, geometry, algorithm, FP)

Numerically integrate a given function `f(::Point)` over the domain defined by
a `geometry` using a particular `integration algorithm` with floating point
precision of type `FP`.

# Arguments
- `f`: an integrand function with a method `f(::Meshes.Point)`
- `geometry`: some `Meshes.Geometry` that defines the integration domain
- `algorithm`: optionally, the `IntegrationAlgorithm` settings to integration with
- `FP`: optionally, the floating point precision desired (`Float64` by default)

Note that reducing `FP` below `Float64` will incur some loss of precision.
"""
function integral end

# If only f and geometry are specified, default to HAdaptiveCubature
function integral(
    f::F,
    geometry::G
) where {F<:Function, G<:Meshes.Geometry}
    _integral(f, geometry, HAdaptiveCubature())
end

# If algorithm is HAdaptiveCubature, use the generalized n-dim solution
function integral(
    f::F,
    geometry::G,
    settings::HAdaptiveCubature
) where {F<:Function, G<:Meshes.Geometry}
    _integral(f, geometry, settings)
end

# If algorithm is HAdaptiveCubature, and FP specified, use the generalized n-dim solution
function integral(
    f::F,
    geometry::G,
    settings::HAdaptiveCubature,
    FP::Type{T}
) where {F<:Function, G<:Meshes.Geometry, T<:AbstractFloat}
    _integral(f, geometry, settings, FP)
end

# If algorithm is not HAdaptiveCubature, specialize on number of dimensions
function integral(
    f::F,
    geometry::G,
    settings::I
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationAlgorithm}
    # Run the appropriate integral type
    Dim = Meshes.paramdim(geometry)
    if Dim == 1
        return _integral_1d(f, geometry, settings)
    elseif Dim == 2
        return _integral_2d(f, geometry, settings)
    elseif Dim == 3
        return _integral_3d(f, geometry, settings)
    end
end

# If algorithm is not HAdaptiveCubature, and FP specified, specialize on number of dimensions
function integral(
    f::F,
    geometry::G,
    settings::I,
    FP::Type{T}
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationAlgorithm, T<:AbstractFloat}
    # Run the appropriate integral type
    Dim = Meshes.paramdim(geometry)
    if Dim == 1
        return _integral_1d(f, geometry, settings, FP)
    elseif Dim == 2
        return _integral_2d(f, geometry, settings, FP)
    elseif Dim == 3
        return _integral_3d(f, geometry, settings, FP)
    end
end


################################################################################
#                    Generalized (n-Dimensional) Worker Methods
################################################################################

# General solution for HAdaptiveCubature methods
function _integral(
    f,
    geometry,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    Dim = Meshes.paramdim(geometry)

    integrand(t) = f(geometry(t...)) * differential(geometry, t)

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = fill(FP(0.5),Dim)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, zeros(FP,Dim), ones(FP,Dim); settings.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end
