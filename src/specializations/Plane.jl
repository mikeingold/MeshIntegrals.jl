################################################################################
#                      Specialized Methods for Plane
################################################################################

function integral(
    f::F,
    plane::Meshes.Plane,
    rule::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]²
    xs, ws = _gausslegendre(FP, rule.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Normalize the Plane's orthogonal vectors
    plane = Plane(plane.p, Meshes.unormalize(plane.u), Meshes.unormalize(plane.v))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f over the Plane
    domainunits = _units(plane(0,0))
    integrand(((wi,wj), (xi,xj))) = wi * wj * f(plane(t(xi), t(xj))) * t′(xi) * t′(xj) * domainunits^2
    return sum(integrand, zip(wws,xxs))
end

function integral(
    f::F,
    plane::Meshes.Plane,
    rule::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Plane's orthogonal vectors
    plane = Plane(plane.p, Meshes.unormalize(plane.u), Meshes.unormalize(plane.v))

    # Integrate f over the Plane
    domainunits = _units(plane(0,0))
    inner∫(v) = QuadGK.quadgk(u -> f(plane(u,v)) * domainunits, FP(-Inf), FP(Inf); rule.kwargs...)[1]
    return QuadGK.quadgk(v -> inner∫(v) * domainunits, FP(-Inf), FP(Inf); rule.kwargs...)[1]
end

function integral(
    f::F,
    plane::Meshes.Plane,
    rule::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Plane's orthogonal vectors
    plane = Plane(plane.p, Meshes.unormalize(plane.u), Meshes.unormalize(plane.v))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f over the Plane
    domainunits = _units(plane(0,0))
    integrand(x::AbstractVector) = f(plane(t(x[1]), t(x[2]))) * t′(x[1]) * t′(x[2]) * domainunits^2

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = FP[0.5, 0.5]
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, FP[-1,-1], FP[1,1]; rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end
