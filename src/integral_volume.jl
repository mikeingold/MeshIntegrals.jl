################################################################################
#                       Generalized 3D Methods
################################################################################

function _integral_3d(
    FP::Type{T} = Float64,
    f,
    geometry3d::G,
    settings::GaussLegendre
) where {T<:AbstractFloat, G<:Meshes.Geometry}
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

function _integral_3d(
    FP::Type{T} = Float64,
    f,
    geometry3d::G,
    settings::HAdaptiveCubature
) where {T<:AbstractFloat, G<:Meshes.Geometry}
    integrand(ts) = f(geometry3d(ts...)) * differential(geometry3d, ts)
    return HCubature.hcubature(integrand, zeros(FP,3), ones(FP,3); settings.kwargs...)[1]
end


################################################################################
#                  Specialized Methods for Tetrahedron
################################################################################

function integral(
    FP::Type{T} = Float64,
    f::F,
    tetrahedron::Meshes.Tetrahedron,
    settings::GaussLegendre
) where {T<:AbstractFloat, F<:Function}
    error("Integrating a Tetrahedron with GaussLegendre not supported.")
end

function integral(
    FP::Type{T} = Float64,
    f::F,
    tetrahedron::Meshes.Tetrahedron,
    settings::GaussKronrod
) where {T<:AbstractFloat, F<:Function}
    inner∫₂(v,w) = QuadGK.quadgk(u -> f(tetrahedron(u,v,w)), FP(0), FP(1-v-w); settings.kwargs...)[1]
    inner∫₁(w) = QuadGK.quadgk(v -> inner∫₂(v,w), FP(0), FP(1-w); settings.kwargs...)[1]
    outer∫ = QuadGK.quadgk(w -> inner∫₁(w), FP(0), FP(1); settings.kwargs...)[1]

    # Apply barycentric domain correction (volume: 1/6 → actual)
    return 6 * volume(tetrahedron) * outer∫
end

function integral(
    FP::Type{T} = Float64,
    f::F,
    tetrahedron::Meshes.Tetrahedron,
    settings::HAdaptiveCubature
) where {T<:AbstractFloat, F<:Function}
    error("Integrating a Tetrahedron with HAdaptiveCubature not supported.")
end


################################################################################
#                         Unsupported Placeholders
################################################################################

function integral(
    FP::Type{T} = Float64,
    f::F,
    ball::Meshes.Ball{Meshes.𝔼{3},CRS,ℒ},
    settings::GaussKronrod
) where {T<:AbstractFloat, F<:Function, CRS, ℒ}
    error("Integrating a Ball in 𝔼{3} with GaussKronrod not supported.")
end

function integral(
    FP::Type{T} = Float64,
    f::F,
    box::Meshes.Box{Meshes.𝔼{3},CRS,ℒ},
    settings::GaussKronrod
) where {T<:AbstractFloat, F<:Function, CRS, ℒ}
    error("Integrating a Box in 𝔼{3} with GaussKronrod not supported.")
end

function integral(
    FP::Type{T} = Float64,
    f::F,
    box::Meshes.Cylinder,
    settings::GaussKronrod
) where {T<:AbstractFloat, F<:Function}
    error("Integrating a Cylinder with GaussKronrod not supported.")
end
