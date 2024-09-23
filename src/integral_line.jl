################################################################################
#                          Generalized 1D Methods
################################################################################

function _integral_1d(
    f,
    geometry,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    integrand(t) = f(geometry(t)) * differential(geometry, [t])
    return QuadGK.quadgk(integrand, FP(0), FP(1); settings.kwargs...)[1]
end


################################################################################
#                   Specialized Methods for Ray
################################################################################

function integral(
    f::F,
    ray::Meshes.Ray,
    settings::GaussLegendre,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = _gausslegendre(FP, settings.n)

    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, Meshes.unormalize(ray.v))

    # Domain transformation: x ∈ [-1,1] ↦ t ∈ [0,∞)
    t₁(x) = FP(1/2) * x + FP(1/2)
    t₁′(x) = FP(1/2)
    t₂(x) = x / (1 - x^2)
    t₂′(x) = (1 + x^2) / (1 - x^2)^2
    t = t₂ ∘ t₁
    t′(x) = t₂′(t₁(x)) * t₁′(x)

    # Integrate f along the Ray
    domainunits = _units(ray(0))
    integrand(x) = f(ray(t(x))) * t′(x) * domainunits
    return sum(w .* integrand(x) for (w,x) in zip(ws, xs))
end

function integral(
    f::F,
    ray::Meshes.Ray,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, Meshes.unormalize(ray.v))

    # Integrate f along the Ray
    domainunits = _units(ray(0))
    return QuadGK.quadgk(t -> f(ray(t)) * domainunits, FP(0), FP(Inf); settings.kwargs...)[1]
end

function integral(
    f::F,
    ray::Meshes.Ray,
    settings::HAdaptiveCubature,
    FP::Type{T} = Float64
) where {F<:Function, T<:AbstractFloat}
    # Normalize the Ray s.t. ray(t) is distance t from origin point
    ray = Ray(ray.p, Meshes.unormalize(ray.v))

    # Domain transformation: x ∈ [0,1] ↦ t ∈ [0,∞)
    t(x) = x / (1 - x^2)
    t′(x) = (1 + x^2) / (1 - x^2)^2
    
    # Integrate f along the Ray
    domainunits = _units(ray(0))
    integrand(x::AbstractVector) = f(ray(t(x[1]))) * t′(x[1]) * domainunits

    # HCubature doesn't support functions that output Unitful Quantity types
    # Establish the units that are output by f
    testpoint_parametriccoord = FP[0.5]
    integrandunits = Unitful.unit.(integrand(testpoint_parametriccoord))
    # Create a wrapper that returns only the value component in those units
    uintegrand(uv) = Unitful.ustrip.(integrandunits, integrand(uv))
    # Integrate only the unitless values
    value = HCubature.hcubature(uintegrand, FP[0], FP[1]; settings.kwargs...)[1]

    # Reapply units
    return value .* integrandunits
end

################################################################################
#                    Specialized Methods for Ring, Rope
################################################################################

function integral(
    f::F,
    ring::Meshes.Ring,
    settings::HAdaptiveCubature
) where {F<:Function}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> _integral(f, segment, settings), segments(ring))
end

function integral(
    f::F,
    ring::Meshes.Ring,
    settings::HAdaptiveCubature,
    FP::Type{T}
) where {F<:Function, T<:AbstractFloat}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> _integral(f, segment, settings, FP), segments(ring))
end

function integral(
    f::F,
    ring::Meshes.Ring,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> integral(f, segment, settings), segments(ring))
end

function integral(
    f::F,
    ring::Meshes.Ring,
    settings::I,
    FP::Type{T}
) where {F<:Function, I<:IntegrationAlgorithm, T<:AbstractFloat}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> integral(f, segment, settings, FP), segments(ring))
end


function integral(
    f::F,
    rope::Meshes.Rope,
    settings::HAdaptiveCubature
) where {F<:Function}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> _integral(f, segment, settings), segments(rope))
end

function integral(
    f::F,
    rope::Meshes.Rope,
    settings::HAdaptiveCubature,
    FP::Type{T}
) where {F<:Function, T<:AbstractFloat}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> _integral(f, segment, settings, FP), segments(rope))
end

function integral(
    f::F,
    rope::Meshes.Rope,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> integral(f, segment, settings), segments(rope))
end

function integral(
    f::F,
    rope::Meshes.Rope,
    settings::I,
    FP::Type{T}
) where {F<:Function, I<:IntegrationAlgorithm, T<:AbstractFloat}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> integral(f, segment, settings, FP), segments(rope))
end
