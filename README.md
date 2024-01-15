# LineIntegrals.jl

The goal of this project is to provide a simple interface for computing line integrals
using geometric representations from [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl).

## Roadmap to General

Status of Integration Methods:

| Function Output Type | Geometry | Implemented | Has Test |
|:---:|:---:|:---:|:---:|
| `Float64`, `Vector{Float64}` | `Meshes.Segment` | :white_check_mark: | :white_check_mark: |
| `Float64`, `Vector{Float64}` | `Meshes.Ring` | :white_check_mark: | :white_check_mark: |
| `Float64`, `Vector{Float64}` | `Meshes.Rope` | :white_check_mark: | :white_check_mark: |
| `Float64`, `Vector{Float64}` | `Meshes.BezierCurve` | :white_check_mark: | :white_check_mark: |
| `Float64`, `Vector{Float64}` | `Vector{<:Meshes.Geometry}` | :white_check_mark: | :white_check_mark: |

Roadmap to release in General:
- Docstrings available for all exports
- Request registration in General

Later plans:
- Implement Documenter for docs
- Add logo
- Test with DynamicQuantities.jl
- Add a method for passing through `kwargs` to `quadgk`
- Consider adding `quadgk(f,::Geometry)` for full control and outputs
- Complete implementation of `SurfaceTrajectory` stuff
