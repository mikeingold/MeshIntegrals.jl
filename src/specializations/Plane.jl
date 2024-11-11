############################################################################################
#                               Specialized Methods for Plane
#
# Why Specialized?
#   The Plane geometry is a special case, representing a planar surface with infinite extent
#   along two basis vectors. This requires a pair of domain transformations mapping from the
#   typical parametric region [0,1]² to an infinite one (-∞,∞)².
############################################################################################

function integral(
    f::Function,
    plane::Meshes.Plane,
    rule::IntegrationRule;
    kwargs...
)
    paramfunction(t1, t2) = _parametric(plane, t1, t2)
    param_plane = _ParametricGeometry(paramfunction, 2)
    return _integral(f, param_plane, rule; kwargs...)
end

############################################################################################
#                                       Parametric
############################################################################################

function _parametric(plane::Meshes.Plane, t1, t2)
    f1(t) = t / (1 - t^2)
    f2(t) = 2t - 1
    f(t) = f1 ∘ f2
    return plane(f(t1), f(t2))
end


#=

function integral(
        f::F,
        plane::Meshes.Plane,
        rule::GaussLegendre;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    _guarantee_analytical(Meshes.Plane, diff_method)

    # Get Gauss-Legendre nodes and weights for a 2D region [-1,1]²
    xs, ws = _gausslegendre(FP, rule.n)
    wws = Iterators.product(ws, ws)
    xxs = Iterators.product(xs, xs)

    # Normalize the Plane's orthogonal vectors
    uu = Meshes.unormalize(plane.u)
    uv = Meshes.unormalize(plane.v)
    uplane = Meshes.Plane(plane.p, uu, uv)

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f over the Plane
    domainunits = _units(uplane(0, 0))
    function integrand(((wi, wj), (xi, xj)))
        wi * wj * f(uplane(t(xi), t(xj))) * t′(xi) * t′(xj) * domainunits^2
    end
    return sum(integrand, zip(wws, xxs))
end

function integral(
        f::F,
        plane::Meshes.Plane,
        rule::GaussKronrod;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    _guarantee_analytical(Meshes.Plane, diff_method)

    # Normalize the Plane's orthogonal vectors
    uu = Meshes.unormalize(plane.u)
    uv = Meshes.unormalize(plane.v)
    uplane = Meshes.Plane(plane.p, uu, uv)

    # Integrate f over the Plane
    domainunits = _units(uplane(0, 0))^2
    integrand(u, v) = f(uplane(u, v)) * domainunits
    inner∫(v) = QuadGK.quadgk(u -> integrand(u, v), FP(-Inf), FP(Inf); rule.kwargs...)[1]
    return QuadGK.quadgk(inner∫, FP(-Inf), FP(Inf); rule.kwargs...)[1]
end

function integral(
        f::F,
        plane::Meshes.Plane,
        rule::HAdaptiveCubature;
        diff_method::DM = Analytical(),
        FP::Type{T} = Float64
) where {F <: Function, DM <: DifferentiationMethod, T <: AbstractFloat}
    _guarantee_analytical(Meshes.Plane, diff_method)

    # Normalize the Plane's orthogonal vectors
    uu = Meshes.unormalize(plane.u)
    uv = Meshes.unormalize(plane.v)
    uplane = Meshes.Plane(plane.p, uu, uv)

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ (-∞,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2

    # Integrate f over the Plane
    domainunits = _units(uplane(0, 0))
    function integrand(x::AbstractVector)
        f(uplane(t(x[1]), t(x[2]))) * t′(x[1]) * t′(x[2]) * domainunits^2
    end

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = zeros(FP, 2)
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, -ones(FP, 2), ones(FP, 2); rule.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end

################################################################################
#                               jacobian
################################################################################

_has_analytical(::Type{T}) where {T <: Meshes.Plane} = true
=#