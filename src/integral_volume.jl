################################################################################
#                       Generalized 3D Methods
################################################################################

function _integral_3d(
    f,
    geometry,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    N = Meshes.paramdim(geometry)

    # Get Gauss-Legendre nodes and weights for a region [-1,1]^N
    xs, ws = _gausslegendre(FP, settings.n)
    weights = Iterators.product(ntuple(Returns(ws), N)...)
    nodes = Iterators.product(ntuple(Returns(xs), N)...)

    # Domain transformation: x [-1,1] ↦ u [0,1]
    t(x) = FP(1//2) * x + FP(1//2)

    function integrand((weights, nodes))
        ts = t.(nodes)
        prod(weights) * f(geometry(ts...)) * differential(geometry, ts)
    end

    return FP(1//(2^N)) .* sum(integrand, zip(weights, nodes))
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
