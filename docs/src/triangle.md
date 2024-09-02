# Integrating a Triangle

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

This non-rectangular Barycentric domain prevents a direct application of most numerical integration methods. It can be directly integrated, albeit inefficiently, using nested Gauss-Kronrod quadrature rules. Alternatively, additional transformation could be applied to map this domain onto a rectangular domain.

**WORK IN PROGRESS:** continued derivation to detail this barycentric-rectangular domain transformation
