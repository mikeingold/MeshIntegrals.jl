# LineIntegrals.jl

The goal of this project is to provide a simple interface for computing line integrals
using geometric representations from [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl).

## Roadmap to General

Status of Integration Methods:
| Function Output Type | Geometry | Implemented | Has Test |
| `Float64` | `Meshes.Segment` | :white_check_mark: | :white_check_mark: |
| `Float64` | `Meshes.BezierCurve` | :white_check_mark: | :white_check_mark: |
| `Float64` | `Vector{<:Meshes.Geometry}` | :white_check_mark: | :white_check_mark: |
| `Vector{Float64}` | `Meshes.Segment` | :white_check_mark: | :white_check_mark: |
| `Vector{Float64}` | `Meshes.BezierCurve` | :white_check_mark: | :white_check_mark: |
| `Vector{Float64}` | `Vector{<:Meshes.Geometry}` | :white_check_mark: | :white_check_mark: |

Roadmap to release in General:
- [ ] Solidify API: e.g. `integrate(f, line)`
- [ ] Docstrings available for all exports
- [ ] Test integral results over a Unitful geometry, verify dimensionality

Later plans:
- [ ] Implement Documenter for docs
- [ ] Complete implementation of `SurfaceTrajectory` stuff
- [ ] Add logo