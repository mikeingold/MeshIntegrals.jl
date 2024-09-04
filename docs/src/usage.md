# Example Usage

## Integrating along a Bezier curve

```julia
using Meshes
using MeshIntegrals

# Define a unit circle on the xy-plane
origin = Point(0,0,0)
ẑ = Vec(0,0,1)
xy_plane = Plane(origin,ẑ)
unit_circle_xy = Circle(xy_plane, 1.0)

# Approximate unit_circle_xy with a high-order Bezier curve
unit_circle_bz = BezierCurve(
    [Point(cos(t), sin(t), 0.0) for t in range(0,2pi,length=361)]
)

# A Real-valued function
f(x, y, z) = abs(x + y)
f(p) = f(to(p)...)

integral(f, unit_circle_xy, GaussKronrod())
    # 0.000170 seconds (5.00 k allocations: 213.531 KiB)
    # ans == 5.656854249525293 m^2

integral(f, unit_circle_bz, GaussKronrod())
    # 0.017122 seconds (18.93 k allocations: 78.402 MiB)
    # ans = 5.551055333711397 m^2
```