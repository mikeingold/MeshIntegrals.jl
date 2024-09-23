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
