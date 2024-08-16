# MeshIntegrals.jl

[![Build Status](https://github.com/mikeingold/MeshIntegrals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mikeingold/MeshIntegrals.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

This package implements methods for numerically-computing integrals over geometric polytopes
from [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl) using the following `::IntegrationAlgorithms`:
- Gauss-Legendre quadrature rules from [**FastGaussQuadrature.jl**](https://github.com/JuliaApproximation/FastGaussQuadrature.jl): `GaussLegendre(n)`
- H-adaptive Gauss-Kronrod quadrature rules from [**QuadGK.jl**](https://github.com/JuliaMath/QuadGK.jl): `GaussKronrod(kwargs...)`
- H-adaptive cubature rules from [**HCubature.jl**](https://github.com/JuliaMath/HCubature.jl): `HAdaptiveCubature(kwargs...)`

Functions available:
- `integral(f, ::Geometry, ::IntegrationAlgorithm)`: integrates a function `f` over a domain defined by `geometry` using a particular `::IntegrationAlgorithm`
- `lineintegral`, `surfaceintegral`, and `volumeintegral` are available as aliases for `integral` that first verify that `geometry` has the appropriate number of parametric dimensions

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
f(x, y, z) = abs(x + y)
f(p) = f(to(p)...)

integral(f, unit_circle_xy, GaussKronrod())
    # 0.000170 seconds (5.00 k allocations: 213.531 KiB)
    # ans == 5.656854249525293 m^2

integral(f, unit_circle_bz, GaussKronrod())
    # 0.017122 seconds (18.93 k allocations: 78.402 MiB)
    # ans = 5.551055333711397 m^2
```

# Support Matrix

| Symbol | Meaning |
|--------|---------|
| :white_check_mark: | Implemented, passes tests |
| :x: | Planned but not yet implemented |
| :warning: | Unable to implement: parameterization not available (see [Issue #28](https://github.com/mikeingold/MeshIntegrals.jl/issues/28)) |

### Integral
| Geometry | Gauss-Legendre | Gauss-Kronrod | H-Adaptive Cubature |
|----------|----------------|---------------|---------------------|
| `Ball{2,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Ball{3,T}` | :white_check_mark: | :x: | :white_check_mark: |
| `Ball{Dim,T}` | :warning: | :warning: | :warning: |
| `BezierCurve{Dim,T,V}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Box{1,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Box{2,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Box{3,T}` | :white_check_mark: | :x: | :white_check_mark: |
| `Box{Dim,T}` | :x: | :x: | :x: |
| `Circle{Dim,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Cone{T}` | :warning: | :warning: | :warning: |
| `ConeSurface{T}` | :x: | :x: | :x: |
| `Cylinder{T}` | :white_check_mark: | :x: | :white_check_mark: |
| `CylinderSurface{T}` | :x: | :white_check_mark: | :x: |
| `Disk{T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Frustum{T}` | :warning: | :warning: | :warning: |
| `FrustumSurface{T}` | :warning: | :warning: | :warning: |
| `Line{Dim,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `ParaboloidSurface{T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Plane{T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Ray{Dim,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Ring{Dim,T,V}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Rope{Dim,T,V}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Segment{Dim,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `SimpleMesh{Dim,T,V}` | :x: | :x: | :x: |
| `Sphere{2,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Sphere{3,T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Tetrahedron{3,T}` | :x: | :white_check_mark: | :x: |
| `Triangle{T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Torus{T}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
