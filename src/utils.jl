################################################################################
#                           Misc. Internal Tools
################################################################################

# Calculate Gauss-Legendre nodes/weights and convert to type T
function _gausslegendre(T, n)
    xs, ws = FastGaussQuadrature.gausslegendre(n)
    return T.(xs), T.(ws)
end

# Extract the length units used by the CRS of a Geometry
function _units(g::Meshes.Geometry{M, CRS}) where {M, CRS}
    return Unitful.unit(CoordRefSystems.lentype(CRS))
end

################################################################################
#                        CliffordNumbers Interface
################################################################################

# Meshes.Vec -> ::CliffordNumber.KVector
function _kvector(v::Meshes.Vec{Dim, T}) where {Dim, T}
    ucoords = Iterators.map(Unitful.ustrip, v.coords)
    return CliffordNumbers.KVector{1, VGA(Dim)}(ucoords...)
end

# Meshes.Vec -> (units, ::CliffordNumber.KVector)
function _ukvector(v::Meshes.Vec{Dim, T}) where {Dim, T}
    units = Unitful.unit(T)
    ucoords = Iterators.map(x -> Unitful.ustrip(units, x), v.coords)
    kvec = CliffordNumbers.KVector{1, VGA(Dim)}(ucoords...)
    return (units, kvec)
end

################################################################################
#                           Error Conditions
################################################################################

@inline function _error_unsupported_gk()
    msg = "Integrating this geometry type with GaussKronrod not supported."
    throw(ArgumentError(msg))
end

@inline function _error_unsupported_combination(geometry, rule)
    msg = "Integrating a $geometry with $rule not supported."
    throw(ArgumentError(msg))
end