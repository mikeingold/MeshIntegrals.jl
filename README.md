# MeshIntegrals.jl

This package implements methods for computing integrals over geometric polytopes
from [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl).

Using Gauss-Legendre quadrature rules from [**FastGaussQuadrature.jl**](https://github.com/JuliaApproximation/FastGaussQuadrature.jl):
- Line integrals
    - `lineintegral(f, ::Meshes.Point...)`
    - `lineintegral(f, ::Meshes.Segment)`
    - `lineintegral(f, ::Meshes.Ring)`
    - `lineintegral(f, ::Meshes.Rope)`
    - `lineintegral(f, ::Meshes.BezierCurve)`
- Surface integrals
    - `surfaceintegral(g, ::Meshes.Triangle)`

Using the h-adaptive Gauss-Kronrod quadrature rules from [**QuadGK.jl**](https://github.com/JuliaMath/QuadGK.jl):
- Line integrals
    - `quadgk_line(f, ::Meshes.Point...)`
    - `quadgk_line(f, ::Meshes.Segment)`
    - `quadgk_line(f, ::Meshes.Ring)`
    - `quadgk_line(f, ::Meshes.Rope)`
    - `quadgk_line(f, ::Meshes.BezierCurve)`
- Surface integrals
    - `quadgk_surface(g, ::Meshes.Triangle)`

Methods are tested to ensure compatibility with
- Meshes.jl geometries with **Unitful.jl** coordinate types, e.g. `Point(1.0u"m", 2.0u"m")`
- Meshes.jl geometries with **DynamicQuantities.jl** coordinate types, e.g. `Point(1.0u"m", 2.0u"m")`
- Any `f(::Meshes.Point{Dim,<:Real})` that maps to a value type that **QuadGK.jl** can integrate, including:
    - Real or complex-valued scalars
    - Real or complex-valued vectors
    - Dimensionful scalars or vectors from Unitful.jl
    - Dimensionful scalars or vectors from DynamicQuantities.jl

## Example Usage

```julia
using BenchmarkTools
using Meshes
using LineIntegrals

# Construct a path that approximates a unit circle on the xy-plane
#   embedded in 3D space using a Bezier curve
unit_circle = BezierCurve(
    [Point(cos(t), sin(t), 0.0) for t in range(0,2pi,length=361)]
)

# Real function
fr(x,y,z) = abs(x + y)
fr(p) = fr(p.coords...)

@btime lineintegral(fr, unit_circle)  # default n=100
    # 9.970 ms (18831 allocations: 78.40 MiB)
    # 5.55240987912083

@btime lineintegral(fr, unit_circle, n=10_000)
    # 16.932 ms (18835 allocations: 78.69 MiB)
    # 5.551055240210768

@btime quadgk_line(fr, unit_circle)
    # 9.871 ms (18829 allocations: 78.40 MiB)
    # (5.551055333711397, 1.609823385706477e-15)
```

# Plans and Work in Progress

- Register in General
- Implement Aqua.jl tests
- Implement Documenter docs
- Implement methods
    - `Meshes.Circle`: `surfaceintegral`, `quadgk_surface`
    - `Meshes.Box{Dim,T}`: `surfaceintegral where {2,T}`, `volumeintegral where {>=3,T}`
    - `Meshes.Ball`: `surfaceintegral`, `volumeintegral`
    - `Meshes.Sphere`: `surfaceintegral`, `volumeintegral`