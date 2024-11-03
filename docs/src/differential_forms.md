# Differential Forms (Work in Progress)

**MeshIntegrals.jl** uses differential forms to perform numerical integration. Every supported `Meshes.Geometry` type is defined as having a parametric function that maps from a local coordinate system to every point on the geometry. For example,
```julia
a = Meshes.Point(0, 0)
b = Meshes.Point(2, 4)
segment = Meshes.Segment(a, b)
```
defines a line segment beginning at point `a` and ending at point `b`. As a geometry with one parametric dimension (i.e. `paramdim(segment) == 1`), it can be treated as a function with one argument that returns a corresponding point on the segment.
```julia
segment(0) == a
segment(0.5) == Point(1, 2)
segment(1) == b
```

The general solution for integrating such a geometry in 3D space can look something like the following, where $t$ is a parametric coordinate used to generate points in the domain.

**TODO: update all of the above for a 2D geometry (Sphere?) to make it more interesting/relevant. Something simple enough to follow but non-trivial.**

Using differential forms, the general solution for integrating a geometry with two parametric dimensions, $t_1$ and $t_2$, is
```math
\iint f(x, y, z) ~ \text{d}t_1 \wedge \text{d}t_2
```

Since `Meshes.Geometry`s parametric functions have arguments that are defined on the domain $[0,1]$, this is equivalent to
```math
\int_0^1 \int_0^1 f(x, y, z) ~ \text{d}t_1 \wedge \text{d}t_2
```

For every point in the integration domain where the integrand function is evaluated, the differential element $\text{d}t_1 \wedge \text{d}t_2$ is calculated using [the Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) of the parametric function. For a two-dimensional surface this Jacobian consists of two vectors, each pointing in the direction that the parametric function's output point will move by changing each input argument. The differential element, then, is simply the magnitude of the exterior product of these vectors.