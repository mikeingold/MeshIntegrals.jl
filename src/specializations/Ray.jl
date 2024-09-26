################################################################################
#                   Specialized Methods for Ray
#
# Why Specialized?
#   The Ray geometry is a special case, representing a ray, originating at a point
#   and extending an infinite length in a particular direction. This requires
#   a domain transformation mapping from the typical parametric region [0,1] to
#   an infinite one (-∞,∞).
################################################################################

function integral(
        f::F,
        ray::Meshes.Ray,
        rule::GaussLegendre;
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, rule.n)

    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Meshes.Ray(ray.p, Meshes.unormalize(ray.v))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ [0,∞)
    t₁(x) = FP(1 / 2) * x + FP(1 / 2)
    t₁′(x) = FP(1 / 2)
    t₂(x) = x / (1 - x^2)
    t₂′(x) = (1 + x^2) / (1 - x^2)^2
    t = t₂ ∘ t₁
    t′(x) = t₂′(t₁(x)) * t₁′(x)

    # Integrate f along the Ray
    domainunits = _units(ray(0))
    integrand(x) = f(ray(t(x))) * t′(x) * domainunits
    return sum(w .* integrand(x) for (w, x) in zip(ws, xs))
end

function integral(
        f::F,
        ray::Meshes.Ray,
        rule::GaussKronrod;
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Meshes.Ray(ray.p, Meshes.unormalize(ray.v))

    # Integrate f along the Ray
    domainunits = _units(ray(0))
    integrand(t) = f(ray(t)) * domainunits
    return QuadGK.quadgk(integrand, zero(FP), FP(Inf); rule.kwargs...)[1]
end

function integral(
        f::F,
        ray::Meshes.Ray,
        rule::HAdaptiveCubature;
        FP::Type{T} = Float64
) where {F <: Function, T <: AbstractFloat}
    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Meshes.Ray(ray.p, Meshes.unormalize(ray.v))

    # Domain transformation: x ∈ [0,1] ↦ t ∈ [0,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f along the Ray
    domainunits = _units(ray(0))
    integrand(x::AbstractVector) = f(ray(t(x[1]))) * t′(x[1]) * domainunits

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = zero(FP)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, zero(FP), one(FP); rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end
