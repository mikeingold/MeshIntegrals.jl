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
        rule::GaussLegendre,
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    error("Integrating a Tetrahedron with GaussLegendre not supported.")
end

function integral(
        f::F,
        tetrahedron::Meshes.Tetrahedron,
        rule::GaussKronrod,
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    function inner∫₂(v, w)
        QuadGK.quadgk(u -> f(tetrahedron(u, v, w)), FP(0), FP(1 - v - w); rule.kwargs...)[1]
    end
    inner∫₁(w) = QuadGK.quadgk(v -> inner∫₂(v, w), FP(0), FP(1 - w); rule.kwargs...)[1]
    outer∫ = QuadGK.quadgk(w -> inner∫₁(w), FP(0), FP(1); rule.kwargs...)[1]

    # Apply barycentric domain correction (volume: 1/6 → actual)
    return 6 * volume(tetrahedron) * outer∫
end

function integral(
        f::F,
        tetrahedron::Meshes.Tetrahedron,
        rule::HAdaptiveCubature,
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    error("Integrating a Tetrahedron with HAdaptiveCubature not supported.")
end
