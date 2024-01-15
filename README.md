# LineIntegrals.jl

The goal of this project is to provide a simple interface for computing line integrals
using geometric representations from [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl).

Verified to work with
- Meshes geometries with Unitful coordinate types, e.g. `Point(1.0u"m", 2.0u"m")`
- Functions that map to Real-valued scalars and vectors
- Functions that map to Real-valued Unitful scalars and vectors

## Roadmap to General

Planned tests to integrate
- `f: Point -> Complex`
- Integrate with DynamicQuantities.jl

Roadmap to release in General:
- Docstrings available for all exports
- Request registration in General

Later plans:
- Implement Documenter
- Add logo
- Add a method for passing through `kwargs` to `quadgk`
- Consider adding `quadgk(f,::Geometry)` for full control and outputs
- Complete implementation of `SurfaceTrajectory` stuff
