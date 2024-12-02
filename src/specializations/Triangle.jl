################################################################################
#                    Specialized Methods for Triangle
#
# Why Specialized?
#   The Triangle geometry is a surface simplex whose parametric function in
#   Meshes.jl uses barycentric coordinates on a domain {u,v} with coordinates
#   that are non-negative and bound by the surface $u + v ≤ 1$. This requires a
#   multi-step domain transformation whose derivation is detailed in the package
#   documentation.
################################################################################

function integral(
    f,
    triangle::Meshes.Triangle,
    rule::IntegrationRule;
    kwargs...
)
    # Generate a _ParametricGeometry whose parametric function domain spans [0,1]²
    param_triangle = _ParametricGeometry(_parametric(triangle), 2)

    # Integrate the _ParametricGeometry using the standard methods
    return _integral(f, param_triangle, rule; kwargs...)
end

################################################################################
#                              Parametric
################################################################################

function _parametric(triangle::Meshes.Triangle)
    function f(t1, t2)
        if any(Iterators.map(n -> (n < 0) || (n > 1), (t1, t2)))
            msg = "triangle(t1, t2) is not defined for (t1, t2) outside [0, 1]²."
            throw(DomainError((t1, t2), msg))
        end

        # Form a line segment between sides
        a = triangle(0, t2)
        b = triangle(t2, 0)
        segment = Meshes.Segment(a, b)

        return segment(t1)
    end
    return f
end
