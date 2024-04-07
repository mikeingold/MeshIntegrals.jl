################################################################################
#                       Generalized 3D Methods
################################################################################

function _integral_3d(
    f,
    geometry3d::G,
    settings::GaussLegendre
) where {Dim, T, G<:Meshes.Geometry{Dim,T}}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = _gausslegendre(T, settings.n)
    wws = Iterators.product(ws, ws, ws)
    xxs = Iterators.product(xs, xs, xs)

    # Domain transformation: x [-1,1] ↦ s,t,u [0,1]
    s(x) = T(1/2) * x + T(1/2)
    t(x) = T(1/2) * x + T(1/2)
    u(x) = T(1/2) * x + T(1/2)

    point(stu) = geometry3d(stu[1], stu[2], stu[3])

    function integrand(((wi,wj,wk), (xi,xj,xk)))
        stu = [s(xi),t(xj),u(xk)]
        wi * wj * wk * f(point(stu)) * differential(geometry3d, stu)
    end

    return T(1/8) .* sum(integrand, zip(wws,xxs))
end

function _integral_3d(
    f,
    geometry3d::G,
    settings::HAdaptiveCubature
) where {Dim, T, G<:Meshes.Geometry{Dim,T}}
    integrand(ts) = differential(geometry3d, ts) * f(geometry3d(ts...))
    return hcubature(integrand, zeros(T,3), ones(T,3); settings.kwargs...)[1]
end


################################################################################
#                         Unsupported Placeholders
################################################################################

function integral(
    f::F,
    ball::Meshes.Ball{3,T},
    settings::GaussKronrod
) where {F<:Function, T}
    error("Integrating a Ball{3,T} with GaussKronrod not supported.")
end

function integral(
    f::F,
    box::Meshes.Box{3,T},
    settings::GaussKronrod
) where {F<:Function, T}
    error("Integrating a Box{3,T} with GaussKronrod not supported.")
end

function integral(
    f::F,
    box::Meshes.Cylinder{T},
    settings::GaussKronrod
) where {F<:Function, T}
    error("Integrating a Cylinder{T} with GaussKronrod not supported.")
end