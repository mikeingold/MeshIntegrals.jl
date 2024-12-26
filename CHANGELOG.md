# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]

No changes yet since v0.16.0.


## [0.16.0] - 2024-12-14

### Added

- Added a `diff_method` keyword argument to the `integral` API, allowing the user to specify which differentiation method should be used when calculating differential element magnitudes throughout the integration domain.
- Implemented `DifferentiationMethod` types:
    - `FiniteDifference` for finite-difference approximation.
    - `AutoEnzyme` for using Enzyme.jl automatic differentiation (AD) via a package extension.
- Added `diff_method` as an optional third argument to the `jacobian` and `differential` API.
- Adds standardized support for integrating over `Tetrahedron` volumes.
- Generalizes integrand functions to support any `f::Any` with a method defined for `f(::Point)`.
- Refactored specialization methods by implementing an internal `_ParametricGeometry <: Meshes.Geometry` to define geometries with custom parametric functions, standardizing support for `BezierCurve`, `Line`, `Plane`, `Ray`, `Tetrahedron`, and `Triangle`.
- Significant performance improvements:
  - Achieved an 80x improvement when integrating over `BezierCurve`.
  - Achieved an up-to-4x improvement when integrating using `HAdaptiveCubature`.

### Deprecated

- Deprecated manual specification of `GaussKronrod` rules for surfaces, i.e. geometries where `Meshes.paramdim(geometry) == 2`. A warning is now generated recommending users switch to `HAdaptiveCubature`.

### Fixed

- Refactored the unit test system.
  - Standardized `combinations.jl` tests by constructing a `TestableGeometry` package and passing it to a `@test` generation function to provide more thorough and standardized test coverage.
  - Reorganized `@testsnippet`s to exist in same source file as relevant tests.
  - Removed `:extended` tag from `Tetrahedron` now that performance is significantly improved.

## [0.15.2] - 2024-10-25

MeshIntegrals.jl is now owned by the JuliaGeometry organization!

### Added

- Added a benchmarking suite using AirspeedVelocity.jl.
- Implemented more unit tests with analytical solutions.

### Changed

- Tagged unit tests for `Meshes.Box` (4D) and `Tetrahedron` as `:extended`, removing them from automatic CI testing due to lengthy compute times.
