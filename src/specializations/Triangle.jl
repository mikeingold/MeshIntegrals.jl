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
    function f(t₁, t₂)
        if any(Iterators.map(n -> (n < 0) || (n > 1), (t₁, t₂)))
            msg = "triangle(t₁, t₂) is not defined for (t₁, t₂) outside [0, 1]²."
            throw(DomainError((t₁, t₂), msg))
        end

        """
        # Algorithm:
        - Form a barycentric triangle bounded by the points [0, 0], [1, 0], and [0, 1].
        - Use t₂ to take a line segment cross-section of the triangle between points
            [0, t₂] and [t₂, 0].
        - Use t₁ to select a point along this line segment, i.e. ā + t₁(b̄ - ā).
        """
        u₁ = t₁ * t₂
        u₂ = t₂ - (t₁ * t₂)
        return triangle(u₁, u₂)
    end
    return f
end
