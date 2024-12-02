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
    param_triangle = _ParametricGeometry(_parametric(triangle), Meshes.paramdim(triangle))

    # Integrate the _ParametricGeometry using the standard methods
    return _integral(f, param_triangle, rule; kwargs...)
end

################################################################################
#                              Parametric
################################################################################

# Map argument domain from [0, 1]² to Barycentric domain for (::Triangle)(t1, t2)
function _parametric(triangle::Meshes.Triangle)
    function f(t1, t2)
        if any(Iterators.map(n -> (n < 0) || (n > 1), (t1, t2)))
            msg = "triangle(t1, t2) is not defined for (t1, t2) outside [0, 1]²."
            throw(DomainError((t1, t2), msg))
        end

        t1t2 = t1 * t2
        return triangle(t1t2, t2 - t1t2)
    end
    return f
end
