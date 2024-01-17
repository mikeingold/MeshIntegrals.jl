# LineIntegrals.jl

This package implements methods for computing line integrals using the adaptive
Gauss-Kronrod quadrature solver from [**QuadGK.jl**](https://github.com/JuliaMath/QuadGK.jl)
and the geometric 1-Dim polytopes representations from [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl).

Verified to work with
- Meshes.jl geometries with Unitful coordinate types, e.g. `Point(1.0u"m", 2.0u"m")`
- Functions that map to Real-valued scalars and vectors
- Functions that map to Real-valued Unitful scalars and vectors

Implements `QuadGK.quadgk` methods for
- `quadgk(f, ::Meshes.Point...) `
- `quadgk(f, ::Meshes.Segment)`
- `quadgk(f, ::Meshes.Ring)`
- `quadgk(f, ::Meshes.Rope)`
- `quadgk(f, ::Meshes.BezierCurve)`

## Roadmap to General

Roadmap:
- [ ] Implement `quadgk` for `Vector{<:Geometry}`
- [ ] Docstrings available for all exports
- [ ] Expand README documentation to include usage examples, logo
- [ ] Decide whether registration in General is appropriate (consult the brain trust)
- [ ] Implement Documenter documentation
- [ ] Continue implementation of `SurfaceTrajectory` concept

Planned tests
- `f: Point -> Complex`
- Also composes with DynamicQuantities.jl?
