################################################################################
#                  Specialized Methods for Tetrahedron
#
# Why Specialized?
#   The Tetrahedron geometry is a volumetric simplex whose parametric function
#   in Meshes.jl uses barycentric coordinates on a domain {u,v,w} with coordinates
#   that are non-negative and bound by the surface $u + v + w ≤ 1$. This requires
#   a multi-step domain transformation whose derivation is detailed in the package
#   documentation.
################################################################################

function integral(
        f,
        tetrahedron::Meshes.Tetrahedron,
        rule::IntegrationRule;
        kwargs...
)
    # Generate a _ParametricGeometry whose parametric function domain spans [0,1]³
    tetra = _ParametricGeometry(_parametric(tetrahedron), Meshes.paramdim(tetrahedron))

    # Integrate the _ParametricGeometry using the standard methods
    return _integral(f, tetra, rule; kwargs...)
end

################################################################################
#                               Parametric
################################################################################

# Map argument domain from [0, 1]³ to Barycentric domain for (::Tetrahedron)(t1, t2, t3)
function _parametric(tetrahedron::Meshes.Tetrahedron)
    function f(t₁, t₂, t₃)
        if any(Iterators.map(n -> (n < 0) || (n > 1), (t₁, t₂, t₃)))
            msg = "tetrahedron(t₁, t₂, t₃) is not defined for (t₁, t₂, t₃) outside [0, 1]³."
            throw(DomainError((t₁, t₂, t₃), msg))
        end

        """
        # Algorithm:
        - Form a barycentric tetrahedron bounded by the points [0, 0, 0], [1, 0, 0],
            [0, 1, 0], and [0, 0, 1].
        - Use t₃ to take a triangular cross-section of the tetrahedron at points
            [t₃, 0, 0], [0, t₃, 0], and [0, 0, t₃].
        - Use t₂ to take a line segment cross-section of the triangle between
            points [t₂t₃, 0, t₃ - t₂t₃] and [0, t₂t₃, t₃ - t₂t₃].
        - Use t₁ to select a point along this line segment, i.e. ā + t₁(b̄ - ā).
        """
        u₁ = (t₂ * t₃) - (t₁ * t₂ * t₃)
        u₂ = t₁ * t₂ * t₃
        u₃ = t₃ - (t₂ * t₃)
        return tetrahedron(u₁, u₂, u₃)
    end
    return f
end
