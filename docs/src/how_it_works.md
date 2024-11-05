# How it Works (By Example) (Work in Progress)

Let $f$ be a function of position ($\bar{r}$)
```julia
function f(r̄::Meshes.Point)
    x, y, z = to(r̄)
    ...
end
```

Let the integration domain be the space (a ball) enclosed by a sphere centered on the origin with a radius of 5 meters.
```julia
center = Meshes.Point(0u"m", 0u"m", 0u"m")
radius = 5.0
ball = Meshes.Ball(center, radius)
```

This integral is often expressed abstractly as simply the following, where the triple integral signs and $\text{d}V$ indicate that the integration domain is some three-dimensional volume.
```math
\iiint f(r̄) ~ \text{d}V
```

## Parametric Functions

Every supported `Meshes.Geometry` type is defined as having a parametric function that maps from a local parametric coordinate system to every point on the geometry. Curve-like geometries will have a single parametric dimension, surfaces will have two dimensions, and volumes will have three dimensions; this can be checked for a particular geometry via `Meshes.paramdim(geometry)`.

For consistency across geometry types, these parametric functions are defined to always take coordinates inside a normalized range $[0,1]$. In the example case of `ball`, Meshes.jl defines a parametric function mapped in normalized spherical coordinates $(t_\rho, ~t_\theta, ~t_\phi)$. We find, then:
```julia
Meshes.paramdim(ball) == 3    # a volume

ball(tρ, tθ, tφ)    # for args in range [0, 1], maps to a corresponding Meshes.Point

ball(0, tθ, tφ) == center
```



## Differential Forms

Using differential forms, the general solution for integrating a geometry with three parametric dimensions ($t_1$, $t_2$, and $t_3$) is
```math
\iiint f(x, y, z) ~ \left\| \text{d}t_1 \wedge \text{d}t_2 \wedge \text{d}t_3 \right\|
```

Since `Meshes.Geometry`s parametric functions have arguments that are defined on the domain $[0,1]$, this is equivalent to
```math
\int_0^1 \int_0^1 \int_0^1 f(x, y, z) ~ \left\| \text{d}t_1 \wedge \text{d}t_2 \wedge \text{d}t_3 \right\|
```

For every point in the integration domain where the integrand function is evaluated, the differential element $\text{d}t_1 \wedge \text{d}t_2$ is calculated using [the Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) of the parametric function. For a two-dimensional surface this Jacobian consists of two vectors, each pointing in the direction that the parametric function's output point will move by changing each input argument. The differential element, then, is simply the magnitude of the exterior product of these vectors.