# [Specializations](@id specializations)

There are several notable exceptions to how Meshes.jl defines [parametric functions](@ref how-parametric).
- `Meshes.ConeSurface` is essentially a composite type and has a parametric function that only maps the conical portion of the geometry, so the `Meshes.Disk` base element has to be integrated separately.
- `Meshes.CylinderSurface` is essentially a composite type and has a parametric function that only maps the cylindrical portion of the geometry, so the `Meshes.Disk` element has to be integrated separately.
- `Meshes.FrustumSurface` is essentially a composite type and has a parametric function that only maps the cylindrical portion of the geometry, so the top and bottom `Meshes.Disk` elements have to be integrated separately.
- `Meshes.Line` represents a line of infinite length that passes through two points, and it has a parametric function that is valid on the domain $(-\infty, \infty)$.
- `Meshes.Plane` represents a plane of infinite extent, and it has a parametric function that is valid on the domain $(-\infty, \infty)^2$.
- `Meshes.Ray` represents a line that begins at a point and extends in a particular direction with infinite length, and it has a parametric function that is valid on the domain $[0, \infty)$.
- `Meshes.Ring` is a composite type that lacks a parametric function, but can be decomposed into `Meshes.Segment`s and then integrated by parts.
- `Meshes.Rope` is a composite type that lacks a parametric function, but can be decomposed into `Meshes.Segment`s and then integrated by parts.
- `Meshes.Triangle` has a parametric function that takes coordinates on a 2D barycentric coordinate system. So, for `(::Meshes.Triangle)(t1, t2)`, the coordinates must obey: $t_1, t_2 \in [0,1]$ where $t_1 + t_2 \le 1$.
- `Meshes.Tetrahedron` has a parametric function that takes coordinates on a 3D barycentric coordinate system. So, for `(::Meshes.Tetrahedron)(t1, t2)`, the coordinates must obey: $t_1, t_2, t_3 \in [0,1]$ where $t_1 + t_2 + t_3 \le 1$.

## Triangle

For a specified `Meshes.Triangle` surface with area $A$, let $u$ and $v$ be Barycentric coordinates that span the surface.
```math
\int_\triangle f(\bar{r}) \, \text{d}A
    = \iint_\triangle f\left( \bar{r}(u,v) \right) \, \left( \text{d}u \wedge \text{d}v \right)
```

Since the geometric transformation from the originally-arbitrary domain to a Barycentric domain is linear, the magnitude of the surface element $\text{d}u \wedge \text{d}v$ is constant throughout the integration domain. This constant will be equal to twice the magnitude of $A$.
```math
\int_\triangle f(\bar{r}) \, \text{d}A
    = 2A \int_0^1 \int_0^{1-v} f\left( \bar{r}(u,v) \right) \, \text{d}u \, \text{d}v
```

Since the integral domain is a right-triangle in the Barycentric domain, a nested application of Gauss-Kronrod quadrature rules is capable of computing the result, albeit inefficiently. However, many numerical integration methods that require rectangular bounds can not be directly applied.

In order to enable integration methods that operate over rectangular bounds, two coordinate system transformations are applied: the first maps from Barycentric coordinates $(u, v)$ to polar coordinates $(r, \phi)$, and the second is a non-linear map from polar coordinates to a new curvilinear basis $(R, \phi)$.

For the first transformation, let $u = r~\cos\phi$ and $v = r~\sin\phi$ where $\text{d}u~\text{d}v = r~\text{d}r~\text{d}\phi$. The Barycentric triangle's hypotenuse boundary line is described by the function $v(u) = 1 - u$. Substituting in the previous definitions leads to a new boundary line function in polar coordinate space $r(\phi) = 1 / (\sin\phi + \cos\phi)$.
```math
\int_0^1 \int_0^{1-v} f\left( \bar{r}(u,v) \right) \, \text{d}u \, \text{d}v =
    \int_0^{\pi/2} \int_0^{1/(\sin\phi+\cos\phi)} f\left( \bar{r}(r,\phi) \right) \, r \, \text{d}r \, \text{d}\phi
```

These integral boundaries remain non-rectangular, so an additional transformation will be applied to a curvilinear $(R, \phi)$ space that normalizes all of the hypotenuse boundary line points to $R=1$. To achieve this, a function $R(r,\phi)$ is required such that $R(r_0, \phi) = 1$ where $r_0 = 1 / (\sin\phi + \cos\phi)$

To achieve this, let $R(r, \phi) = r~(\sin\phi + \cos\phi)$. Now, substituting some terms leads to
```math
\int_0^{\pi/2} \int_0^{1/(\sin\phi+\cos\phi)} f\left( \bar{r}(r,\phi) \right) \, r \, \text{d}r \, \text{d}\phi
    = \int_0^{\pi/2} \int_0^{r_0} f\left( \bar{r}(r,\phi) \right) \, \left(\frac{R}{\sin\phi + \cos\phi}\right) \, \text{d}r \, \text{d}\phi
```

Since $\text{d}R/\text{d}r = \sin\phi + \cos\phi$, a change of integral domain leads to
```math
\int_0^{\pi/2} \int_0^{r_0} f\left( \bar{r}(r,\phi) \right) \, \left(\frac{R}{\sin\phi + \cos\phi}\right) \, \text{d}r \, \text{d}\phi
    = \int_0^{\pi/2} \int_0^1 f\left( \bar{r}(R,\phi) \right) \, \left(\frac{R}{\left(\sin\phi + \cos\phi\right)^2}\right) \, \text{d}R \, \text{d}\phi
```

The second term in this new integrand function serves as a correction factor that corrects for the impact of the non-linear domain transformation. Since all of the integration bounds are now constants, specialized integration methods can be defined for triangles that performs these domain transformations and then solve the new rectangular integration problem using a wider range of solver options.
