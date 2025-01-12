# Tips

## Performance

### Use an `atol` setting for integrals with near-zero solutions

By default, the solvers used for `GaussKronrod` and `HAdaptiveCuabture` integration use a relative tolerance (`rtol`) setting to determine when a solution is sufficiently accurate. For integrals whose true solution equal exactly or very nearly zero, the solver can struggle to meet this relative error tolerance and significantly lengthen compute times. In such a case, performance can often be improved by explicitly providing an absolute tolerance setting (`atol`).

```julia
julia> using BenchmarkTools, Meshes, MeshIntegrals, Unitful

julia> segment = Segment(Point(-1u"m"), Point(1u"m"))
Segment
├─ Point(x: -1.0 m)
└─ Point(x: 1.0 m)

julia> f(p::Meshes.Point) = p.coords.x * u"N/m"
f (generic function with 1 method)

julia> @btime integral(f, segment, GaussKronrod())
  154.285 ms (28 allocations: 33.70 MiB)
2.3298940582252377e-11 m N

julia> @btime integral(f, segment, GaussKronrod(atol = 1e-8u"N*m"))
  846.377 ns (9 allocations: 208 bytes)
1.1316035800104648e-11 m N
```
