# Support Matrix

While this library aims to support all possible integration rules and **Meshes.jl**
geometry types, some combinations are ill-suited and some others are simplu not yet
implemented. The following Support Matrix aims to capture the current development state of
all geometry/rule combinations. Entries with a green check mark are fully supported
and have passing unit tests that provide some confidence they produce accurate results.

In general, Gauss-Kronrod integration rules are recommended (and the default) for geometries
with one parametric dimension, e.g.: `Segment`, `BezierCurve`, and `Rope`. Gauss-Kronrod
rules can also be applied to some geometries with more dimensions by nesting multiple
integration solves, but this is inefficient. These Gauss-Kronrod rules are supported (but
not recommended) for surface-like geometries, but not for volume-like geometries. For
geometries with more than one parametric dimension, e.g. surfaces and volumes, H-Adaptive
Cubature integration rules are recommended (and the default).

| Symbol | Support Level |
|--------|---------|
| ✅ | Supported, passes unit tests |
| 🎗️ | Planned to support in the future |
| 🛑 | Not supported |

| `Meshes.Geometry` | Gauss-Legendre | Gauss-Kronrod | H-Adaptive Cubature |
|----------|----------------|---------------|---------------------|
| `Ball` in `𝔼{2}` | ✅ | ✅ | ✅ |
| `Ball` in `𝔼{3}` | ✅ | 🛑 | ✅ |
| `BezierCurve` | ✅ | ✅ | ✅ |
| `Box` in `𝔼{1}` | ✅ | ✅ | ✅ |
| `Box` in `𝔼{2}` | ✅ | ✅ | ✅ |
| `Box` in `𝔼{3}` | ✅ | 🛑 | ✅ |
| `Circle` | ✅ | ✅ | ✅ |
| `Cone` | ✅ | ✅ | ✅ |
| `ConeSurface` | ✅ | ✅ | ✅ |
| `Cylinder` | ✅ | 🛑 | ✅ |
| `CylinderSurface` | ✅ | ✅ | ✅ |
| `Disk` | ✅ | ✅ | ✅ |
| `Frustum` | [🎗️](https://github.com/mikeingold/MeshIntegrals.jl/issues/57) | [🎗️](https://github.com/mikeingold/MeshIntegrals.jl/issues/57) | [🎗️](https://github.com/mikeingold/MeshIntegrals.jl/issues/57) |
| `FrustumSurface` | ✅ | ✅ | ✅ |
| `Line` | ✅ | ✅ | ✅ |
| `ParaboloidSurface` | ✅ | ✅ | ✅ |
| `Plane` | ✅ | ✅ | ✅ |
| `Quadrangle` | ✅ | ✅ | ✅ |
| `Ray` | ✅ | ✅ | ✅ |
| `Ring` | ✅ | ✅ | ✅ |
| `Rope` | ✅ | ✅ | ✅ |
| `Segment` | ✅ | ✅ | ✅ |
| `SimpleMesh` | [🎗️](https://github.com/mikeingold/MeshIntegrals.jl/issues/27) | [🎗️](https://github.com/mikeingold/MeshIntegrals.jl/issues/27) | [🎗️](https://github.com/mikeingold/MeshIntegrals.jl/issues/27) |
| `Sphere` in `𝔼{2}` | ✅ | ✅ | ✅ |
| `Sphere` in `𝔼{3}` | ✅ | ✅ | ✅ |
| `Tetrahedron` in `𝔼{3}` | [🎗️](https://github.com/mikeingold/MeshIntegrals.jl/issues/40) | ✅ | [🎗️](https://github.com/mikeingold/MeshIntegrals.jl/issues/40) |
| `Triangle` | ✅ | ✅ | ✅ |
| `Torus` | ✅ | ✅ | ✅ |
