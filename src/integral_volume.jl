################################################################################
#                               Gauss-Legendre
################################################################################

function integral(
    f,
    ball::Meshes.Ball{3,T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,3,T)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws, ws)
    xxs = Iterators.product(xs, xs, xs)

    # Domain transformation: xi,xj,xk [-1,1] ↦ s,t,u [0,1]
    s(xi) = 0.5xi + 0.5
    t(xj) = 0.5xj + 0.5
    u(xk) = 0.5xk + 0.5

    # Integrate the ball in parametric (s,t,u)-space [0,1]³
    integrand(s,t,u) = s^2 * sinpi(t) * f(ball(s,t,u))
    g(((wi,wj,wk), (xi,xj,xk))) = wi * wj * wk * integrand(s(xi),t(xj),u(xk))
    R = ball.radius
    return (1/8) * 2π^2 * R^3 .* sum(g, zip(wws,xxs))
end

function integral(
    f,
    box::Meshes.Box{3,T},
    settings::GaussLegendre
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,3,T)

    # Get Gauss-Legendre nodes and weights for a 3D region [-1,1]³
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws, ws)
    xxs = Iterators.product(xs, xs, xs)

    # Domain transformation: xi,xj,xk [-1,1] ↦ s,t,u [0,1]
    s(xi) = 0.5xi + 0.5
    t(xj) = 0.5xj + 0.5
    u(xk) = 0.5xk + 0.5
    point(xi,xj,xk) = box(s(xi), t(xj), u(xk))

    # Calculate weight-node product
    g(((wi,wj,wk), (xi,xj,xk))) = wi * wj * wk * f(point(xi,xj,xk))

    # Calculate 3D Gauss-Legendre integral of f over parametric coordinates [-1,1]³
    # Apply a linear domain-correction factor [-1,1]³ ↦ volume(box)
    return (1/8) * volume(box) .* sum(g, zip(wws,xxs))
end


################################################################################
#                             GaussKronrod
################################################################################

function integral(
    f,
    ball::Meshes.Ball{3,T},
    settings::GaussKronrod
) where {T}
    error("Integrating a Ball{3,T} with GaussKronrod not supported.")
end

function integral(
    f,
    box::Meshes.Box{3,T},
    settings::GaussKronrod
) where {T}
    error("Integrating a Ball{3,T} with GaussKronrod not supported.")
end


################################################################################
#                               HCubature
################################################################################

function integral(
    f,
    ball::Meshes.Ball{3,T},
    settings::HAdaptiveCubature
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,3,T)

    # Integrate the ball in parametric (s,t,u)-space [0,1]³
    integrand(s,t,u) = s^2 * sinpi(t) * f(ball(s,t,u))
    integrand(stu) = integrand(stu[1],stu[2],stu[3])
    intval = hcubature(stu -> integrand(stu), [0,0,0], [1,1,1], settings.kwargs...)[1]

    R = ball.radius
    return 2π^2 * R^3 .* intval
end

function integral(
    f,
    box::Meshes.Box{3,T},
    settings::HAdaptiveCubature
) where {T}
    # Validate the provided integrand function
    _validate_integrand(f,3,T)

    # Integrate the box in parametric (s,t,u)-space
    integrand(stu) = f(box(stu[1],stu[2],stu[3]))
    intval = hcubature(integrand, [0,0,0], [1,1,1], settings.kwargs...)[1]

    # Apply a linear domain-correction factor [0,1]³ ↦ volume(box)
    return volume(box) .* intval
end