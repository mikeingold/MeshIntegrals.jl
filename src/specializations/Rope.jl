################################################################################
#                      Specialized Methods for Rope
#
# Why Specialized?
#   The Rope geometry defines a polytope whose length spans segments between
#   consecutive points. Meshes.jl does not define a parametric function for
#   Rope's, but they can be decomposed into their constituent Segment's,
#   integrated separately, and then summed.
################################################################################

function integral(
        f::F,
        rope::Meshes.Rope,
        rule::HAdaptiveCubature
) where {F <: Function}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> _integral(f, segment, rule), segments(rope))
end

function integral(
        f::F,
        rope::Meshes.Rope,
        rule::HAdaptiveCubature,
        FP::Type{T}
) where {F <: Function, T <: AbstractFloat}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> _integral(f, segment, rule, FP), segments(rope))
end

function integral(
        f::F,
        rope::Meshes.Rope,
        rule::I
) where {F <: Function, I <: IntegrationRule}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> integral(f, segment, rule), segments(rope))
end

function integral(
        f::F,
        rope::Meshes.Rope,
        rule::I,
        FP::Type{T}
) where {F <: Function, I <: IntegrationRule, T <: AbstractFloat}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> integral(f, segment, rule, FP), segments(rope))
end
