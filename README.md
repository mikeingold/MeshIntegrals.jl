# LineIntegrals.jl

This package implements methods for computing line integrals using the adaptive
Gauss-Kronrod quadrature solver from [**QuadGK.jl**](https://github.com/JuliaMath/QuadGK.jl)
and the geometric 1-Dim polytopes representations from [**Meshes.jl**](https://github.com/JuliaGeometry/Meshes.jl).

Verified to work with
- Meshes.jl geometries with **Unitful.jl** coordinate types, e.g. `Point(1.0u"m", 2.0u"m")`
- Meshes.jl geometries with **DynamicQuantities.jl** coordinate types, e.g. `Point(1.0u"m", 2.0u"m")`
- Any `f(::Meshes.Point)` that maps to a value type that QuadGK can integrate, including:
    - Real or complex-valued scalars
    - Real or complex-valued vectors
    - Dimensionful scalars or vectors from Unitful.jl
    - Dimensionful scalars or vectors from DynamicQuantities.jl

Implements `QuadGK.quadgk` methods for
- `quadgk(f, ::Meshes.Point...) `
- `quadgk(f, ::Meshes.Segment)`
- `quadgk(f, ::Meshes.Ring)`
- `quadgk(f, ::Meshes.Rope)`
- `quadgk(f, ::Meshes.BezierCurve)`

## Example Usage

```julia
using Meshes
using QuadGK
using LineIntegrals

# Construct a path that approximates a unit circle on the xy-plane
#   embedded in 3D space using a Bezier curve
unit_circle = BezierCurve(
    [Point(cos(t), sin(t), 0.0) for t in range(0,2pi,length=361)]
)

# Real function
fr(x,y,z) = abs(x + y)
fr(p) = fr(p.coords...)
quadgk(fr, unit_circle)
    # (5.551055333711397, 1.1102230246251565e-16)

# Complex function
#   Treat the xy-plane as the complex plane: x + iy
fc(z::Complex) = 1/z
fc(p) = fc(complex(p.coords[1],p.coords[2]))
quadgk(fc, unit_circle)
    # (-0.017331713663560157 + 0.0im, 4.585229184817946e-12)
```