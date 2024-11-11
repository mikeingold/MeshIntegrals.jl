############################################################################################
#                    Specialized Methods for Triangle
#
# Why Specialized?
#   The Triangle geometry is a surface simplex whose parametric function in Meshes.jl uses
#   barycentric coordinates on a domain {u,v} with coordinates that are non-negative and
#   bound by the surface $u + v ≤ 1$. A transformation is used to map this volume to a
#   rectangular domain [0,1]^3. 
############################################################################################

function integral(
    f::F,
    triangle::Meshes.Triangle,
    rule::IntegrationRule;
    kwargs...
) where {F <: Function}
    paramfunction(t1, t2) = _parametric(triangle, t1, t2)
    tri = _ParametricGeometry(paramfunction, 2)
    return _integral(f, tri, rule; kwargs...)
end

################################################################################
#                              Parametric
################################################################################

function _parametric(triangle::Meshes.Triangle, t1, t2)
    @assert (0 <= t1 <= 1) && (0 <= t2 <= 1)

    # Form a horizontal line segment at t2
    a = triangle(0, t2)
    b = triangle(_constrain(1 - t2), t2)
    segment = Meshes.Segment(a, b)

    return segment(t1)
end
