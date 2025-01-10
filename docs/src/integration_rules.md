# Integration Rules

An **integration rule** is a method or algorithm used to numerically calculate the value of an integral. When an integral is calculated using the two-argument form `integral(f, geometry)`, **MeshIntegrals.jl** will automatically select an integration rule. Default rules are generally well-behaved, but may not be optimal for every problem.

**MeshIntegrals.jl** defines an abstract type `IntegrationRule` with sub-types representing the various integration rules supported by this package. These rules can be specified when calculating an integral using the three-argument form `integral(f, geometry, rule)`.

## Gauss-Kronrod

The `GaussKronrod` type is used to specify an adaptive Gauss-Kronrod quadrature rule, as implemented by [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl).

```julia
rule = GaussKronrod()
```

All standard `QuadGK.quadgk` keyword-argument options are supported. These can be specified when constructing the rule, where the `kwargs` in `GaussKronrod(kwargs...)` is equivalent to `quadgk(f, a, b; kwargs...)`, e.g.:
```julia
rule = GaussKronrod(order = 5, rtol = 1e-4)
```

## Gauss-Legendre

The `GaussLegendre` type is used to specify a Gauss-Legendre quadrature rule, as implemented by [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl).

```julia
rule = GaussLegendre(N)
```

## H-Adaptive Cubature

The `HAdaptiveCubature` type is used to specify an h-adaptive cubature rule, as implemented by [HCubature.jl](https://github.com/JuliaMath/HCubature.jl).

```julia
rule = HAdaptiveCubature()
```

All standard `HCubature.hcubature` keyword-argument options are supported. These can be specified when constructing the rule, where the `kwargs` in `HAdaptiveCubature(kwargs...)` is equivalent to `hcubature(f, a, b; kwargs...)`, e.g.:
```julia
rule = HAdaptiveCubature(order = 5, rtol = 1e-4)
```
