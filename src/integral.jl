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
    lineintegral(f, geometry, algorithm::IntegrationAlgorithm=GaussKronrod)

Numerically integrate a given function `f(::Point)` along a line-like `geometry`
using a particular `integration algorithm`.

Algorithm types available:
- GaussKronrod
- GaussLegendre
"""
function lineintegral(
    f::F,
    geometry::G
) where {F<:Function, G<:Meshes.Geometry}
    return lineintegral(f, geometry, GaussKronrod())
end

"""
    surfaceintegral(f, geometry, algorithm::IntegrationAlgorithm=HAdaptiveCubature)

Numerically integrate a given function `f(::Point)` over a surface `geometry`
using a particular `integration algorithm`.
"""
function surfaceintegral(
    f::F,
    geometry::G
) where {F<:Function, G<:Meshes.Geometry}
    return surfaceintegral(f, geometry, HAdaptiveCubature())
end

"""
    volumeintegral(f, geometry, algorithm::IntegrationAlgorithm=HAdaptiveCubature)

Numerically integrate a given function `f(::Point)` throughout a volumetric
`geometry` using a particular `integration algorithm`.

"""
function volumeintegral(
    f::F,
    geometry::G
) where {F<:Function, G<:Meshes.Geometry}
    return volumeintegral(f, geometry, HAdaptiveCubature())
end
