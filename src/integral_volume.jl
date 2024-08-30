################################################################################
#                       Generalized 3D Methods
################################################################################

function _integral_3d(
    f,
    geometry3d,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]^2
    xs, ws = _gausslegendre(FP, settings.n)
    wws = Iterators.product(ws, ws, ws)
    xxs = Iterators.product(xs, xs, xs)

    # Domain transformation: x [-1,1] ↦ s,t,u [0,1]
    t(x) = FP(1/2) * x + FP(1/2)

    function integrand(((wi,wj,wk), (xi,xj,xk)))
        ts = t.([xi, xj, xk])
        wi * wj * wk * f(geometry3d(ts...)) * differential(geometry3d, ts)
    end

    return FP(1/8) .* sum(integrand, zip(wws,xxs))
end

# Integrating volumes with GaussKronrod not supported by default
function _integral_3d(
    f,
    geometry,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    error("Integrating this volume type with GaussKronrod not supported.")
end


################################################################################
#                  Specialized Methods for Tetrahedron
################################################################################

function integral(
    f::F,
    tetrahedron::Meshes.Tetrahedron,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    error("Integrating a Tetrahedron with GaussLegendre not supported.")
end

function integral(
    f::F,
    tetrahedron::Meshes.Tetrahedron,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    inner∫₂(v,w) = QuadGK.quadgk(u -> f(tetrahedron(u,v,w)), FP(0), FP(1-v-w); settings.kwargs...)[1]
    inner∫₁(w) = QuadGK.quadgk(v -> inner∫₂(v,w), FP(0), FP(1-w); settings.kwargs...)[1]
    outer∫ = QuadGK.quadgk(w -> inner∫₁(w), FP(0), FP(1); settings.kwargs...)[1]

    # Apply barycentric domain correction (volume: 1/6 → actual)
    return 6 * volume(tetrahedron) * outer∫
end

function integral(
    f::F,
    tetrahedron::Meshes.Tetrahedron,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64,
) where {F<:Function, T<:AbstractFloat}
    error("Integrating a Tetrahedron with HAdaptiveCubature not supported.")
end
