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
    t(x) = T(1/2) * x + T(1/2)

    function integrand(((wi,wj,wk), (xi,xj,xk)))
        ts = t.([xi, xj, xk])
        wi * wj * wk * f(geometry3d(ts...)) * differential(geometry3d, ts)
    end

    return T(1/8) .* sum(integrand, zip(wws,xxs))
end

function _integral_3d(
    f,
    geometry3d::G,
    settings::HAdaptiveCubature
) where {Dim, T, G<:Meshes.Geometry{Dim,T}}
    integrand(ts) = f(geometry3d(ts...)) * differential(geometry3d, ts)
    return HCubature.hcubature(integrand, zeros(T,3), ones(T,3); settings.kwargs...)[1]
end


################################################################################
#                  Specialized Methods for Tetrahedron
################################################################################

function integral(
    f::F,
    tetrahedron::Meshes.Tetrahedron,
    settings::GaussLegendre
) where {F<:Function, T}
    error("Integrating a Tetrahedron with GaussLegendre not supported.")
end

function integral(
    f::F,
    tetrahedron::Meshes.Tetrahedron{3,T},
    settings::GaussKronrod
) where {F<:Function, T}
    # Validate the provided integrand function
    _validate_integrand(f,3,T)

    inner∫₂(v,w) = QuadGK.quadgk(u -> f(tetrahedron(u,v,w)), T(0), T(1-v-w); settings.kwargs...)[1]
    inner∫₁(w) = QuadGK.quadgk(v -> inner∫₂(v,w), T(0), T(1-w); settings.kwargs...)[1]
    outer∫ = QuadGK.quadgk(w -> inner∫₁(w), T(0), T(1); settings.kwargs...)[1]

    # Apply barycentric domain correction (volume: 1/6 → actual)
    return 6 * volume(tetrahedron) * outer∫
end

function integral(
    f::F,
    tetrahedron::Meshes.Tetrahedron,
    settings::HAdaptiveCubature
) where {F<:Function, T}
    error("Integrating a Tetrahedron with HAdaptiveCubature not supported.")
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
    box::Meshes.Cylinder,
    settings::GaussKronrod
) where {F<:Function}
    error("Integrating a Cylinder with GaussKronrod not supported.")
end
