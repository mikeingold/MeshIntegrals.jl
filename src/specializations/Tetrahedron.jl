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
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    _error_unsupported_combination("Tetrahedron", "GaussLegendre")
end

function integral(
        f::F,
        tetrahedron::Meshes.Tetrahedron,
        rule::GaussKronrod;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    _guarantee_analytical(Meshes.Tetrahedron, diff_method)

    o = zero(FP)
    ∫uvw(u, v, w) = f(tetrahedron(u, v, w))
    ∫vw(v, w) = QuadGK.quadgk(u -> ∫uvw(u, v, w), o, FP(1 - v - w); rule.kwargs...)[1]
    ∫w(w) = QuadGK.quadgk(v -> ∫vw(v, w), o, FP(1 - w); rule.kwargs...)[1]
    ∫ = QuadGK.quadgk(∫w, o, one(FP); rule.kwargs...)[1]

    # Apply barycentric domain correction (volume: 1/6 → actual)
    return 6 * Meshes.volume(tetrahedron) * ∫
end

function integral(
        f::F,
        tetrahedron::Meshes.Tetrahedron,
        rule::HAdaptiveCubature;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    _error_unsupported_combination("Tetrahedron", "HAdaptiveCubature")
end

################################################################################
#                               jacobian
################################################################################

_has_analytical(::Type{T}) where {T <: Meshes.Tetrahedron} = true
