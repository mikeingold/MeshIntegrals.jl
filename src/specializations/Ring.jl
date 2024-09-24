################################################################################
#                      Specialized Methods for Ring
################################################################################

function integral(
        f::F,
        ring::Meshes.Ring,
        rule::HAdaptiveCubature
) where {F <: Function}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> _integral(f, segment, rule), segments(ring))
end

function integral(
        f::F,
        ring::Meshes.Ring,
        rule::HAdaptiveCubature,
        FP::Type{T}
) where {F <: Function, T <: AbstractFloat}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> _integral(f, segment, rule, FP), segments(ring))
end

function integral(
        f::F,
        ring::Meshes.Ring,
        rule::I
) where {F <: Function, I <: IntegrationRule}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> integral(f, segment, rule), segments(ring))
end

function integral(
        f::F,
        ring::Meshes.Ring,
        rule::I,
        FP::Type{T}
) where {F <: Function, I <: IntegrationRule, T <: AbstractFloat}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> integral(f, segment, rule, FP), segments(ring))
end
