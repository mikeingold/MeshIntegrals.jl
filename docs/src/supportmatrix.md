# Support Matrix

This library aims to enable users to calculate the value of integrals over all [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl)
geometry types using an array of numerical integration rules and techniques. However, some
combinations of geometry types and integration rules are ill-suited, and some others are simply
not yet implemented. The following Support Matrix captures the current state of support for
all geometry/rule combinations. Entries with a green check mark are fully supported and pass
unit tests designed to check for accuracy.

In general, Gauss-Kronrod integration rules are recommended (and the default) for geometries
with one parametric dimension, e.g.: `Segment`, `BezierCurve`, and `Rope`. For geometries with
more than one parametric dimension, e.g. surfaces and volumes, H-Adaptive Cubature rules are
recommended (and the default).

While it is possible to apply nested Gauss-Kronrod rules to numerically integrate geometries
with more than one parametric dimension, this produces results that are strictly inferior to
using an equivalent H-Adaptive Cubature rule, so support for this usage is not recommended.

| Symbol | Support Level |
|--------|---------|
| ✅ | Supported |
| 🎗️ | Planned to support in the future |
| ⚠️ | Deprecated |
| 🛑 | Not supported |

| `Meshes.Geometry` | Gauss-Legendre | Gauss-Kronrod | H-Adaptive Cubature |
|----------|----------------|---------------|---------------------|
| `Ball` in `𝔼{2}` | ✅ | ⚠️ | ✅ |
| `Ball` in `𝔼{3}` | ✅ | 🛑 | ✅ |
| `BezierCurve` | ✅ | ✅ | ✅ |
| `Box` in `𝔼{1}` | ✅ | ✅ | ✅ |
| `Box` in `𝔼{2}` | ✅ | ⚠️ | ✅ |
| `Box` in `𝔼{≥3}` | ✅ | 🛑 | ✅ |
| `Circle` | ✅ | ✅ | ✅ |
| `Cone` | ✅ | 🛑 | ✅ |
| `ConeSurface` | ✅ | ⚠️ | ✅ |
| `Cylinder` | ✅ | 🛑 | ✅ |
| `CylinderSurface` | ✅ | ⚠️ | ✅ |
| `Disk` | ✅ | ⚠️ | ✅ |
| `Ellipsoid` | ✅ | ✅ | ✅ |
| `Frustum` | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) |
| `FrustumSurface` | ✅ | ⚠️ | ✅ |
| `Hexahedron` | ✅ | ✅ | ✅ |
| `Line` | ✅ | ✅ | ✅ |
| `ParaboloidSurface` | ✅ | ⚠️ | ✅ |
| `ParametrizedCurve` | ✅ | ✅ | ✅ |
| `Plane` | ✅ | ✅ | ✅ |
| `Polyarea` | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) |
| `Pyramid` | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) |
| `Quadrangle` | ✅ | ⚠️ | ✅ |
| `Ray` | ✅ | ✅ | ✅ |
| `Ring` | ✅ | ✅ | ✅ |
| `Rope` | ✅ | ✅ | ✅ |
| `Segment` | ✅ | ✅ | ✅ |
| `SimpleMesh` | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/27) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/27) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/27) |
| `Sphere` in `𝔼{2}` | ✅ | ✅ | ✅ |
| `Sphere` in `𝔼{3}` | ✅ | ⚠️ | ✅ |
| `Tetrahedron` | ✅ | ⚠️ | ✅ |
| `Triangle` | ✅ | ✅ | ✅ |
| `Torus` | ✅ | ⚠️ | ✅ |
| `Wedge` | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) | [🎗️](https://github.com/JuliaGeometry/MeshIntegrals.jl/issues/28) |
