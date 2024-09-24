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
  kwargs::Any
  GaussKronrod(; kwargs...) = new(kwargs)
end

"""
    GaussLegendre(n)

Numerically integrate using an `n`'th-order Gauss-Legendre quadrature rule. Nodes
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
    HAdaptiveCubature(kwargs...)

Numerically integrate areas and surfaces using the h-adaptive cubature rule
implemented by HCubature.jl. All standard `HCubature.hcubature` keyword
arguments are supported.
"""
struct HAdaptiveCubature <: IntegrationAlgorithm
  kwargs::Any
  HAdaptiveCubature(; kwargs...) = new(kwargs)
end
