# LineIntegrals.jl

The goal of this project is to provide a simple interface for computing line integrals
using geometric representations from [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl).

## Plans

Roadmap to release in General:
- [ ] Solidify API: e.g. `integrate(f, line)`
- [ ] Docstrings available for all exports
- [ ] Implement integration over a `BezierCurve`
- [ ] Test integral results over a Unitful geometry, verify dimensionality

Later plans:
- [ ] Implement Documenter for docs
- [ ] Complete implementation of `SurfaceTrajectory` stuff
