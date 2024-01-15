# LineIntegrals.jl

The goal of this project is to provide a simple interface for computing line integrals
using geometric representations from [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl).

Verified to work with
- Meshes geometries with Unitful coordinate types, e.g. `Point(1.0u"m", 2.0u"m")`
- Functions that map to Real-valued scalars and vectors
- Functions that map to Real-valued Unitful scalars and vectors

Implements `QuadGK.quadgk` methods for
- `quadgk(f, ::Meshes.Segment)`
- `quadgk(f, ::Meshes.BezierCurve)`

## Roadmap to General

Roadmap:
- [ ] Docstrings available for all exports
- [ ] Expand README documentation to include usage examples, logo
- [ ] Request registration in General
- [ ] Implement Documenter documentation
- [ ] Complete implementation of `SurfaceTrajectory` concept

Planned tests
- `f: Point -> Complex`
- Integrate with DynamicQuantities.jl
