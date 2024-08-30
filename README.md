# MeshIntegrals.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![ColPrac](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet?style=flat-square)](https://github.com/SciML/ColPrac)

[![Build Status](https://github.com/mikeingold/MeshIntegrals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mikeingold/MeshIntegrals.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/mikeingold/MeshIntegrals.jl/graph/badge.svg)](https://codecov.io/gh/mikeingold/MeshIntegrals.jl)
[![Coveralls](https://coveralls.io/repos/github/mikeingold/MeshIntegrals.jl/badge.svg?branch=main)](https://coveralls.io/github/mikeingold/MeshIntegrals.jl?branch=main)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


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
| :stop_sign: | Not supported |

### Integral
| Geometry | Gauss-Legendre | Gauss-Kronrod | H-Adaptive Cubature |
|----------|----------------|---------------|---------------------|
| `Ball` in `𝔼{2}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Ball` in `𝔼{3}` | :white_check_mark: | :stop_sign: | :white_check_mark: |
| `BezierCurve` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Box` in `𝔼{1}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Box` in `𝔼{2}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Box` in `𝔼{3}` | :white_check_mark: | :stop_sign: | :white_check_mark: |
| `Circle` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Cone` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `ConeSurface` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Cylinder` | :white_check_mark: | :stop_sign: | :white_check_mark: |
| `CylinderSurface` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Disk` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Frustum` | :warning: | :warning: | :warning: |
| `FrustumSurface` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Line` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `ParaboloidSurface` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Plane` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Ray` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Ring` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Rope` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Segment` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `SimpleMesh` | :x: | :x: | :x: |
| `Sphere` in `𝔼{2}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Sphere` in `𝔼{3}` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Tetrahedron` in `𝔼{3}` | :x: | :white_check_mark: | :x: |
| `Triangle` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| `Torus` | :white_check_mark: | :white_check_mark: | :white_check_mark: |
