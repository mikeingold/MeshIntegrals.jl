# Surface Integral Over a Triangle

**Note: This explanation reflects my original method but is no longer valid.
An updated derivation is planned. **

A linear transformation can be applied that maps any triangle onto a Barycentric
coordinate system. A linear correction factor is then applied to correct for the
domain transformation, where the area of the original triangle is $A$ and the area
of the Barycentric triangle is $1/2$.
```math
\int_\triangle f(\bar{r}) \text{d}A
    = \frac{A}{1/2} \int_0^1 \int_0^{1-v} f(u,v) \text{d}v \text{d}u
```

This Barycentric integral can be directly estimated using a nested application of
the h-adaptive Gauss-Kronrod quadrature rules from QuadGK.jl (`quadgk_surface`).
Alternatively, with some additional modification it can be transformed onto a
fixed rectangular domain for integration via cubature rules or other nested
quadrature rules.

Let $g$ be a wrapper function for $f$ that is only non-zero at valid Barycentric
coordinates.
```math
g(u,v,) =
    \begin{cases}
        f(u,v) & \text{if } 0 \le u+v \le 1 \\
        0 & \text{otherwise}
    \end{cases}
```

Then
```math
\int_0^1 \int_0^{1-v} f(u,v) \text{d}v \text{d}u
    = \int_0^1 \int_0^1 g(u,v) \text{d}v \text{d}u
```

A domain transformation can be applied to map the Barycentric coordinate domains
from $u,v \in [0,1]$ to $s,t \in [-1,1]$, enabling the application of Gauss-Legendre
quadrature rules.
```math
s(u) = 2u - 1 \\
u(s) = \frac{s+1}{2}
\text{d}u = \frac{1}{2}\text{d}s
```
```math
t(v) = 2v - 1 \\
v(t) = \frac{t+1}{2}
\text{d}v = \frac{1}{2}\text{d}t
```

Leading to
```math
\int_0^1 \int_0^1 g(u,v) \text{d}v \text{d}u
    = \frac{1}{4} \int_{-1}^1 \int_{-1}^1 g(\frac{s+1}{2},\frac{t+1}{2}) \text{d}t \text{d}s
```

Gauss-Legendre nodes ($x \in [-1,1]$) and weights ($w$) for a rule of order $N$
can be efficiently calculated using the FastGaussQuadrature.jl package.
```math
\int_{-1}^1 \int_{-1}^1 g(\frac{s+1}{2},\frac{t+1}{2}) \text{d}t \text{d}s
    \approx \sum_{i=1}^N \sum_{j=1}_N w_i w_j f(x_i,x_j)
```

This approximation can be rolled back up the stack of equations, leading to an
expression that numerically approximates the original integral problem as
```math
\int_\triangle f(\bar{r}) \text{d}A
    = \frac{A}{4(1/2)} \sum_{i=1}^N \sum_{j=1}_N w_i w_j f(\frac{x_i+1}{2}, \frac{x_j+1}{2})
```
