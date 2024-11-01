################################################################################
#                           Misc. Internal Tools
################################################################################

# Calculate Gauss-Legendre nodes/weights and convert to type T
function _gausslegendre(T, n)
    xs, ws = FastGaussQuadrature.gausslegendre(n)
    return T.(xs), T.(ws)
end

# Common error message structure
function _error_unsupported_combination(geometry, rule)
    msg = "Integrating a $geometry using a $rule rule not supported."
    throw(ArgumentError(msg))
end

function _guarantee_analytical(G, diff_method)
    if diff_method != Analytical()
        ArgumentError("Geometry type $G requires kwarg diff_method = Analytical()")
    end
end

################################################################################
#                        CliffordNumbers and Units
################################################################################

# Meshes.Vec -> ::CliffordNumber.KVector
function _kvector(v::Meshes.Vec{Dim, T}) where {Dim, T}
    ucoords = Iterators.map(Unitful.ustrip, v.coords)
    return CliffordNumbers.KVector{1, VGA(Dim)}(ucoords...)
end

# Extract the length units used by the CRS of a Geometry
_units(::Geometry{M, CRS}) where {M, CRS} = Unitful.unit(CoordRefSystems.lentype(CRS))
