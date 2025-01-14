# Tips

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

```julia
julia> using BenchmarkTools, Meshes, MeshIntegrals, Unitful

julia> segment = Segment(Point(-1u"m"), Point(1u"m"))
Segment
├─ Point(x: -1.0 m)
└─ Point(x: 1.0 m)

julia> f(p::Meshes.Point) = p.coords.x * u"N/m"
f (generic function with 1 method)
```

Calculating this problem with a default `GaussKronrod()` rule produces a solution that is very close to zero, on the order of $10^{-11}$. However, the result took a non-trivial amount of time due to a large number of memory allocatons for such a trivial problem.
```julia
julia> @btime integral(f, segment, GaussKronrod())
  154.285 ms (28 allocations: 33.70 MiB)
2.3298940582252377e-11 m N
```

Providing an explicit `atol` setting produces a solution that happens to be just as accurate but returns much faster.
```julia
julia> @btime integral(f, segment, GaussKronrod(atol = 1e-8u"N*m"))
  846.377 ns (9 allocations: 208 bytes)
1.1316035800104648e-11 m N
```
