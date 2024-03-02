# Validate that f has a method defined for f(::Point{Dim,T})
@inline function _validate_integrand(f,Dim,T)
    if hasmethod(f, (Point{Dim,T},))
        return nothing
    else
        error("The provided Function f must have a method f(::Meshes.Point{$Dim,$T})")
    end
end

################################################################################
#                         Integration Algorithms
################################################################################

abstract type IntegrationAlgorithm end

"""
    GaussKronrod(kwargs...)

Numerically integrate using the h-adaptive Gauss-Kronrod quadrature rule implemented
by QuadGK.jl. All standard [`QuadGK.quadgk`](@ref) keyword arguments are supported.
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
implemented by HCubature.jl. All standard [`HCubature.hcubature`](@ref) keyword
arguments are supported.
"""
struct HAdaptiveCubature <: IntegrationAlgorithm
    kwargs
    HAdaptiveCubature(; kwargs...) = new(kwargs)
end

################################################################################
#                        Integral Function API
################################################################################

"""
    integral(f, geometry, algorithm::IntegrationAlgorithm)

Numerically integrate a given function `f(::Point)` over the domain defined by
a `geometry` using a particular `integration algorithm`.
"""
function integral(
    f::F,
    geometry::Meshes.Geometry{Dim,T},
    settings::I=HAdaptiveCubature()
) where {F<:Function, Dim, T, I<:IntegrationAlgorithm}
    # Validate that the provided function has an appropriate f(::Point{Dim,T}) method
    _validate_integrand(f, Dim, T)

    # Run the appropriate integral type
    dim_param = paramdim(geometry)
    if dim_param == 1
        return _integral_1d(f, geometry, settings)
    elseif dim_param == 2
        return _integral_2d(f, geometry, settings)
    elseif dim_param == 3
        return _integral_3d(f, geometry, settings)
    end
end

"""
    lineintegral(f, geometry, algorithm::IntegrationAlgorithm=GaussKronrod)

Numerically integrate a given function `f(::Point)` along a line-like `geometry`
using a particular `integration algorithm`.

Algorithm types available:
- GaussKronrod
- GaussLegendre
- HAdaptiveCubature
"""
function lineintegral(
    f::F,
    geometry::G,
    settings::I=GaussKronrod()
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationAlgorithm}
    dim = paramdim(geometry)
    if dim == 1
        return _integral_1d(f, geometry, settings)
    else
        error("Performing a line integral on a geometry with $dim parametric dimensions not supported.")
    end
end

"""
    surfaceintegral(f, geometry, algorithm::IntegrationAlgorithm=HAdaptiveCubature)

Numerically integrate a given function `f(::Point)` over a surface `geometry`
using a particular `integration algorithm`.
"""
function surfaceintegral(
    f::F,
    geometry::G,
    settings::I=HAdaptiveCubature()
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationAlgorithm}
    dim = paramdim(geometry)
    if dim == 2
        return _integral_2d(f, geometry, settings)
    else
        error("Performing a surface integral on a geometry with $dim parametric dimensions not supported.")
    end
end

"""
    volumeintegral(f, geometry, algorithm::IntegrationAlgorithm=HAdaptiveCubature)

Numerically integrate a given function `f(::Point)` throughout a volumetric
`geometry` using a particular `integration algorithm`.

"""
function volumeintegral(
    f::F,
    geometry::G,
    settings::I=HAdaptiveCubature()
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationAlgorithm}
    dim = paramdim(geometry)
    if dim == 3
        return _integral_3d(f, geometry, settings)
    else
        error("Performing a volume integral on a geometry with $dim parametric dimensions not supported.")
    end
end
