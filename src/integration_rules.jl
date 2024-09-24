################################################################################
#                           Integration Rules
################################################################################

abstract type IntegrationRule end

"""
    GaussKronrod(kwargs...)

The h-adaptive Gauss-Kronrod quadrature rule implemented by QuadGK.jl. All standard
`QuadGK.quadgk` keyword arguments are supported. This rule works natively for one
dimensional geometries; some two- and three-dimensional geometries are additionally
supported using nested integral solvers with the specified `kwarg` settings.
"""
struct GaussKronrod <: IntegrationRule
    kwargs
    GaussKronrod(; kwargs...) = new(kwargs)
end

"""
    GaussLegendre(n)

An `n`'th-order Gauss-Legendre quadrature rule. Nodes and weights are
efficiently calculated using FastGaussQuadrature.jl.

So long as the integrand function can be well-approximated by a polynomial of
order `2n-1`, this method should yield results with 16-digit accuracy in `O(n)`
time. If the function is know to have some periodic content, then `n` should
(at a minimum) be greater than the expected number of periods over the geometry,
e.g. `length(geometry)/λ`.
"""
struct GaussLegendre <: IntegrationRule
    n::Int64
end

"""
    HAdaptiveCubature(kwargs...)

The h-adaptive cubature rule implemented by HCubature.jl. All standard
`HCubature.hcubature` keyword arguments are supported.
"""
struct HAdaptiveCubature <: IntegrationRule
    kwargs
    HAdaptiveCubature(; kwargs...) = new(kwargs)
end
