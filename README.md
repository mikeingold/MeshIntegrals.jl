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

# Example Usage

```julia
using Meshes
using MeshIntegrals

# Define a unit circle on the xy-plane
origin = Point(0,0,0)
ẑ = Vec(0,0,1)
xy_plane = Plane(origin,ẑ)
unit_circle_xy = Circle(xy_plane, 1.0)

# Approximate unit_circle_xy with a high-order Bezier curve
unit_circle_bz = BezierCurve(
    [Point(cos(t), sin(t), 0.0) for t in range(0,2pi,length=361)]
)

# A Real-valued function
f(x,y,z) = abs(x + y)
f(p) = f(p.coords...)

integral(f, unit_circle_xy, GaussKronrod())
    # 56.500 μs (1819 allocations: 100.95 KiB)
    # ans == 5.656854249502878

integral(f, unit_circle_bz, GaussKronrod())
    # 9.638 ms (18830 allocations: 78.40 MiB)
    # ans = 5.551055333711397
```

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
| `Meshes.Cylinder{T}` | :white_check_mark: | :white_check_mark: |
