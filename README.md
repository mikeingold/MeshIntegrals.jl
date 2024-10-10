# MeshIntegrals.jl

[![Docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mikeingold.github.io/MeshIntegrals.jl/dev/)
[![Docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mikeingold.github.io/MeshIntegrals.jl/stable/)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![ColPrac](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet?style=flat-square)](https://github.com/SciML/ColPrac)

[![Build Status](https://github.com/mikeingold/MeshIntegrals.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mikeingold/MeshIntegrals.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/mikeingold/MeshIntegrals.jl/graph/badge.svg)](https://codecov.io/gh/mikeingold/MeshIntegrals.jl)
[![Coveralls](https://coveralls.io/repos/github/mikeingold/MeshIntegrals.jl/badge.svg?branch=main)](https://coveralls.io/github/mikeingold/MeshIntegrals.jl?branch=main)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


**MeshIntegrals.jl** uses differential forms to enable fast and easy numerical integration of arbitrary integrand functions over domains defined via [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl) geometries. This is achieved using:
- Gauss-Legendre quadrature rules from [**FastGaussQuadrature.jl**](https://github.com/JuliaApproximation/FastGaussQuadrature.jl): `GaussLegendre(n)`
- H-adaptive Gauss-Kronrod quadrature rules from [**QuadGK.jl**](https://github.com/JuliaMath/QuadGK.jl): `GaussKronrod(kwargs...)`
- H-adaptive cubature rules from [**HCubature.jl**](https://github.com/JuliaMath/HCubature.jl): `HAdaptiveCubature(kwargs...)`

These solvers have support for integrand functions that produce scalars, vectors, and [**Unitful.jl**](https://github.com/PainterQubits/Unitful.jl) `Quantity` types. While HCubature.jl does not natively support `Quantity` type integrands, this package provides a compatibility layer to enable this feature.

## Usage

```julia
integral(f, geometry)
```
Performs a numerical integration of some integrand function `f(p::Meshes.Point)` over the domain specified by `geometry`. A default integration method will be automatically selected according to the geometry: `GaussKronrod()` for 1D, and `HAdaptiveCubature()` for all others.

```julia
integral(f, geometry, rule, FP=Float64)
```
Performs a numerical integration of some integrand function `f(p::Meshes.Point)` over the domain specified by `geometry` using the specified integration rule, e.g. `GaussKronrod()`.

Optionally, a fourth argument can be provided to specify the floating point precision level desired. This setting can be manipulated if your integrand function produces outputs with alternate floating point precision (e.g. `Float16`, `BigFloat`, etc) AND you'd prefer to avoid implicit type promotions.

```julia
lineintegral(f, geometry)
surfaceintegral(f, geometry)
volumeintegral(f, geometry)
```
Alias functions are provided for convenience. These are simply wrappers for `integral` that first validate that the provided `geometry` has the expected number of parametric/manifold dimensions. Like with `integral` in the examples above, the `rule` can also be specified as a third-argument.
- `lineintegral` for curve-like geometries or polytopes (e.g. `Segment`, `Ray`, `BezierCurve`, `Rope`, etc)
- `surfaceintegral` for surfaces (e.g. `Disk`, `Sphere`, `CylinderSurface`, etc)
- `volumeintegral` for (3D) volumes (e.g. `Ball`, `Cone`, `Torus`, etc)

# Demo

```julia
using Meshes
using MeshIntegrals
using Unitful

# Define a path that approximates a sine-wave on the xy-plane
mypath = BezierCurve(
    [Point(t*u"m", sin(t)*u"m", 0.0u"m") for t in range(-pi, pi, length=361)]
)

# Map f(::Point) -> f(x, y, z) in unitless coordinates
f(p::Meshes.Point) = f(ustrip(to(p))...)

# Integrand function in units of Ohms/meter
f(x, y, z) = (1 / sqrt(1 + cos(x)^2)) * u"Ω/m"

integral(f, mypath)
# -> Approximately 2*Pi Ω
```
