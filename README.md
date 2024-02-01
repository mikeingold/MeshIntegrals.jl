# LineIntegrals.jl

This package implements methods for computing line integrals along geometric 1-Dim polytopes
from [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl). Two sets of methods are
currently implemented:
- `LineIntegrals.integral(f, geometry)` uses Gauss-Legendre quadratures from [**FastGaussQuadrature.jl**](https://github.com/JuliaApproximation/FastGaussQuadrature.jl)
- `QuadGK.quadgk(f, geometry)` is exported using the adaptive Gauss-Kronrod quadrature from [**QuadGK.jl**](https://github.com/JuliaMath/QuadGK.jl)

All methods are verified to work with
- Meshes.jl geometries with **Unitful.jl** coordinate types, e.g. `Point(1.0u"m", 2.0u"m")`
- Meshes.jl geometries with **DynamicQuantities.jl** coordinate types, e.g. `Point(1.0u"m", 2.0u"m")`
- Any `f(::Meshes.Point)` that maps to a value type that QuadGK can integrate, including:
    - Real or complex-valued scalars
    - Real or complex-valued vectors
    - Dimensionful scalars or vectors from Unitful.jl
    - Dimensionful scalars or vectors from DynamicQuantities.jl

Implements `QuadGK.quadgk` methods for
- `quadgk(f, ::Meshes.Point...) `
- `quadgk(f, ::Meshes.Segment)`
- `quadgk(f, ::Meshes.Ring)`
- `quadgk(f, ::Meshes.Rope)`
- `quadgk(f, ::Meshes.BezierCurve)`

## Example Usage

```julia
using BenchmarkTools
using Meshes
using QuadGK
using LineIntegrals

# Construct a path that approximates a unit circle on the xy-plane
#   embedded in 3D space using a Bezier curve
unit_circle = BezierCurve(
    [Point(cos(t), sin(t), 0.0) for t in range(0,2pi,length=361)]
)

# Real function
fr(x,y,z) = abs(x + y)
fr(p) = fr(p.coords...)

@btime integral(fr, unit_circle)  # default n=100
# 9.970 ms (18831 allocations: 78.40 MiB)
# 5.55240987912083

@btime integral(fr, unit_circle, n=10_000)
# 16.932 ms (18835 allocations: 78.69 MiB)
# 5.551055240210768

@btime quadgk(fr, unit_circle)
# 44.874 ms (78229 allocations: 331.95 MiB)
# (5.551055333711397, 1.1102230246251565e-16)
```

# Work in Progress

- Need to troubleshoot functions of complex variables
    - Currently failing a test. Not sure if math error on my part or calculating incorrectly.