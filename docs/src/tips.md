# Tips

## General Usage

### Make use of the Enzyme extension

!!! note
    See [How it Works](@ref howitworks) for a more thorough explanation of how **MeshIntegrals.jl** uses differential forms to solve integral problems.

**MeshIntegrals.jl** uses differential forms to solve integral problems. At every sampled point in the geometry that forms the integration domain a differential element magnitude must be calculated, which involves calculating the `jacobian` of the geometry's parametric function. **MeshIntegrals.jl** includes multiple method options for calculating this `jacobian` which are selectable via `integral`s keyword argument `diff_method`.

The default/fallback `jacobian` method uses a finite-difference approximation method, which is reasonably performant and compatible with all geometry types. However, this method can potentially introduce a small amount of approximation error into solutions. This method can be explicitly selected via the keyword argument:
```julia
integral(f, geometry; diff_method = FiniteDifference())
```

**MeshIntegrals.jl** includes an extension for [**Enzyme.jl**](https://github.com/EnzymeAD/Enzyme.jl), using it to implement an [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation)-based `jacobian` method. This method is typically faster and more accurate, but it is not currently compatible with all geometries. This extension will automatically be loaded when **Enzyme.jl** is present in the active Julia environment. When loaded, this method will automatically be used by `integral` whenever it is compatible with the given geometry. This method can also be explicitly selected via the keyword argument:
```julia
using Enzyme

integral(f, geometry; diff_method = AutoEnzyme())
```

## Performance

### Using explicit tolerance settings with adaptive integration rules

!!! note
    This tip is applicable to the use of `GaussKronrod` and `HAdaptiveCubature` integration rules.

The `GaussKronrod` and `HAdaptiveCuabture` integration rules both make use of adaptive solvers that use relative tolerance (`rtol`) settings to determine when the result is sufficiently precise to return a solution. Under some circumstances these solvers may struggle to converge on a precise enough solution, resulting in extended compute times. This commonly occurs for integrand functions that are numerically unstable or whose solution is zero or near-zero. In such a case, performance can often be improved by explicitly providing an absolute tolerance setting (`atol`).

This can be observed by benchmarking the integration of an integral problem whose true solution is exactly zero.
```math
\int_{-1\text{m}}^{1\text{m}} x \left[\tfrac{\text{N}}{\text{m}}\right] ~ \text{d}x
= 0 \,\text{Nm}
```

```@repl tip_tolerances
using BenchmarkTools, Meshes, MeshIntegrals, Unitful
segment = Segment(Point(-1u"m"), Point(1u"m"))
f(p::Meshes.Point) = p.coords.x * u"N/m"
```

Calculating this problem with a default `GaussKronrod()` rule produces a solution that is very close to zero, on the order of $10^{-11}$. However, the result took a non-trivial amount of time due to a large number of memory allocatons for such a trivial problem.
```@repl tip_tolerances
@btime integral(f, segment, GaussKronrod())
```

Providing an explicit `atol` setting produces a solution that happens to be just as accurate but returns much faster.
```@repl tip_tolerances
@btime integral(f, segment, GaussKronrod(atol = 1e-8u"N*m"))
```
