################################################################################
#                      Specialized Methods for Ring
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
