################################################################################
#                  Specialized Methods for Tetrahedron
#
# Why Specialized?
#   The Tetrahedron geometry is a volumetric simplex whose parametric function
#   in Meshes.jl uses barycentric coordinates on a domain {u,v,w} with coordinates
#   that are non-negative and bound by the surface $u + v + w ≤ 1$. This requires
#   a multi-step domain transformation whose derivation is detailed in the package
#   documentation.
################################################################################

function integral(
        f::F,
        tetrahedron::Meshes.Tetrahedron,
        rule::GaussLegendre;
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    _error_unsupported_combination("Tetrahedron", "GaussLegendre")
end

function integral(
        f::F,
        tetrahedron::Meshes.Tetrahedron,
        rule::GaussKronrod;
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    nil = zero(FP)
    ∫uvw(u, v, w) = f(tetrahedron(u, v, w))
    ∫vw(v, w) = QuadGK.quadgk(u -> ∫uvw(u, v, w), nil, FP(1 - v - w); rule.kwargs...)[1]
    ∫w(w) = QuadGK.quadgk(v -> ∫vw(v, w), nil, FP(1 - w); rule.kwargs...)[1]
    ∫ = QuadGK.quadgk(∫w, nil, one(FP); rule.kwargs...)[1]

    # Apply barycentric domain correction (volume: 1/6 → actual)
    return 6 * Meshes.volume(tetrahedron) * ∫
end

function integral(
        f::F,
        tetrahedron::Meshes.Tetrahedron,
        rule::HAdaptiveCubature;
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    _error_unsupported_combination("Tetrahedron", "HAdaptiveCubature")
end
