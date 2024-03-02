# MeshIntegrals.jl

This package implements methods for numerically-computing integrals over geometric polytopes
from [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl) using:
- Gauss-Legendre quadrature rules from [**FastGaussQuadrature.jl**](https://github.com/JuliaApproximation/FastGaussQuadrature.jl): `GaussLegendre(n)`
- H-adaptive Gauss-Kronrod quadrature rules from [**QuadGK.jl**](https://github.com/JuliaMath/QuadGK.jl): `GaussKronrod(kwargs...)`
- H-adaptive cubature rules from [**HCubature.jl**](https://github.com/JuliaMath/HCubature.jl): `HAdaptiveCubature(kwargs...)`

Functions available:
- `integral(f, geometry, ::IntegrationAlgorithm)`: integrates a function `f` over a domain defined by `geometry` using a particular
- `lineintegral`, `surfaceintegral`, and `volumeintegral` are available as aliases for `integral` that first verify that `geometry` has the appropriate number of parametric dimensions

Methods are tested to ensure compatibility with
- Meshes.jl geometries with **Unitful.jl** coordinate types, e.g. `Point(1.0u"m", 2.0u"m")`
- Meshes.jl geometries with **DynamicQuantities.jl** coordinate types, e.g. `Point(1.0u"m", 2.0u"m")`
- Any `f(::Meshes.Point{Dim,<:Real})` that maps to a value type that **QuadGK.jl** can integrate, including:
    - Real or complex-valued scalars
    - Real or complex-valued vectors
    - Dimensionful scalars or vectors from Unitful.jl
    - Dimensionful scalars or vectors from DynamicQuantities.jl

# Support Matrix

| Symbol | Meaning |
|--------|---------|
| :white_check_mark: | Implemented, passes tests |
| :yellow_square: | Implemented, not yet validated |
| :x: | Not yet implemented |

### Integral
| Geometry | Gauss-Legendre | Gauss-Kronrod | H-Adaptive Cubature |
|----------|----------------|---------------|---------------------|
| `Meshes.Box{Dim>3,T}` | :x: | :x: | :x: |

### Line Integral
| Geometry | Gauss-Legendre | Gauss-Kronrod | H-Adaptive Cubature |
|----------|----------------|---------------|---------------------|
| `Meshes.BezierCurve` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.Box{1,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.Circle` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.Line` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.Ring` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.Rope` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.Segment` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.Sphere{2,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |

### Surface Integral
| Geometry | Gauss-Legendre | Gauss-Kronrod | H-Adaptive Cubature |
|----------|----------------|---------------|-------------------|
| `Meshes.Ball{2,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.Box{2,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.CylinderSurface` | :x: | :white_check_mark: | :x: |
| `Meshes.Disk` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.ParaboloidSurface` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.Sphere{3,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.Triangle` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.Torus` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Meshes.SimpleMesh` | :x: | :x: | :x: |

### Volume Integral
| Geometry | Gauss-Legendre | H-Adaptive Cubature |
|----------|----------------|---------------|
| `Meshes.Ball{3,T}` | :white_check_mark: | :white_check_mark: |
| `Meshes.Box{3,T}` | :white_check_mark: | :white_check_mark: |

# Example Usage

```julia
using BenchmarkTools
using Meshes
using MeshIntegrals

# Construct a path that approximates a unit circle on the xy-plane
#   embedded in 3D space using a Bezier curve
unit_circle = BezierCurve(
    [Point(cos(t), sin(t), 0.0) for t in range(0,2pi,length=361)]
)

# Real function
fr(x,y,z) = abs(x + y)
fr(p) = fr(p.coords...)

@btime lineintegral(fr, unit_circle)
    # 9.970 ms (18831 allocations: 78.40 MiB)
    # 5.55240987912083

@btime lineintegral(fr, unit_circle, GaussLegendre(10_000))
    # 16.932 ms (18835 allocations: 78.69 MiB)
    # 5.551055240210768

@btime lineintegral(fr, unit_circle, GaussKronrod())
    # 9.871 ms (18829 allocations: 78.40 MiB)
    # (5.551055333711397, 1.609823385706477e-15)
```

# Plans and Work in Progress

- TODO
    - Update Example Usage and benchmarks
    - Re-implement all tests for Unitful compatibility
    - Implement Aqua.jl tests
    - Improve README documentation in lieu of Documenter
    - Implement Monte Carlo integration methods
    - Continue working to generalize and consolidate methods

- Longer term plans
    - Assess impact of upcoming Meshes.jl CRS refactor
    - Once functionality is established and stable, evaluate transition to JuliaGeometry organization or direct absorption into Meshes.jl
