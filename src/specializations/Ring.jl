################################################################################
#                      Specialized Methods for Ring
#
# Why Specialized?
#   The Ring geometry defines a polytope whose length spans segments between
#   consecutive points that form a closed path. Meshes.jl does not define a
#   parametric function for Ring's, but they can be decomposed into their
#   constituent Segment's, integrated separately, and then summed.
################################################################################

function integral(
        f::F,
        ring::Meshes.Ring,
        rule::I;
        kwargs...
) where {F <: Function, I <: IntegrationRule}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> _integral(f, segment, rule; kwargs...), segments(ring))
end
