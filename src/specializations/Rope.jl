################################################################################
#                      Specialized Methods for Rope
################################################################################

function integral(
        f::F, rope::Meshes.Rope, settings::HAdaptiveCubature) where {F <: Function}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> _integral(f, segment, settings), segments(rope))
end

function integral(
        f::F,
        rope::Meshes.Rope,
        settings::HAdaptiveCubature,
        FP::Type{T}
) where {F <: Function, T <: AbstractFloat}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> _integral(f, segment, settings, FP), segments(rope))
end

function integral(f::F, rope::Meshes.Rope,
        settings::I) where {F <: Function, I <: IntegrationAlgorithm}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> integral(f, segment, settings), segments(rope))
end

function integral(
        f::F,
        rope::Meshes.Rope,
        settings::I,
        FP::Type{T}
) where {F <: Function, I <: IntegrationAlgorithm, T <: AbstractFloat}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> integral(f, segment, settings, FP), segments(rope))
end
