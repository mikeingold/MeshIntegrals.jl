################################################################################
#                               Gauss-Legendre
################################################################################

function _integral_3d(f, geometry3d, settings::GaussLegendre)
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = gausslegendre(settings.n)
    wws = Iterators.product(ws, ws, ws)
    xxs = Iterators.product(xs, xs, xs)

    # Domain transformation: x [-1,1] ↦ s,t,u [0,1]
    s(x) = 0.5x + 0.5
    t(x) = 0.5x + 0.5
    u(x) = 0.5x + 0.5

    point(stu) = geometry3d(stu[1], stu[2], stu[3])

    function paramfactor(stu)
        J = jacobian(geometry3d, stu)
        return abs((J[1] × J[2]) ⋅ J[3])
    end

    function integrand(((wi,wj,wk), (xi,xj,xk)))
        stu = [s(xi),t(xj),u(xk)]
        wi * wj * wk * f(point(stu)) * paramfactor(stu)
    end

    return (1/8) .* sum(integrand, zip(wws,xxs))
end


################################################################################
#                             GaussKronrod
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
    error("Integrating a Ball{3,T} with GaussKronrod not supported.")
end


################################################################################
#                               HCubature
################################################################################

# Generalized method
function _integral_3d(f, geometry3d, settings::HAdaptiveCubature)
    function paramfactor(ts)
        J = jacobian(geometry3d, ts)
        return abs((J[1] × J[2]) ⋅ J[3])
    end

    integrand(ts) = paramfactor(ts) * f(geometry3d(ts...))
    return hcubature(integrand, zeros(3), ones(3); settings.kwargs...)[1]
end
