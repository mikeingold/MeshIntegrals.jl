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
    function f(t1, t2, t3)
        if any(Iterators.map(n -> (n < 0) || (n > 1), (t1, t2, t3)))
            msg = "tetrahedron(t1, t2, t3) is not defined for t outside [0, 1]."
            throw(DomainError((t1, t2, t3), msg))
        end

        # Take a triangular cross-section at t3
        a = tetrahedron(t3, 0, 0)
        b = tetrahedron(0, t3, 0)
        c = tetrahedron(0, 0, t3)
        cross_section = _parametric(Meshes.Triangle(a, b, c))

        return cross_section(t1, t2)
    end
    return f
end
