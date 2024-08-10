################################################################################
#                              Line Integral
################################################################################

"""
    lineintegral(f, geometry)
    lineintegral(f, geometry, algorithm)
    lineintegral(f, geometry, algorithm, FP)

Numerically integrate a given function `f(::Point)` along a line-like `geometry`
using a particular `integration algorithm` with floating point precision of type
`FP`.

Algorithm types available:
- GaussKronrod (default)
- GaussLegendre
- HAdaptiveCubature
"""
function lineintegral(
    f::F,
    geometry::G,
) where {F<:Function, G<:Meshes.Geometry}
    Dim = Meshes.paramdim(geometry)

    if Dim == 1
        return integral(f, geometry, GaussKronrod())
    else
        error("Performing a line integral on a geometry with $Dim parametric dimensions not supported.")
    end
end

function lineintegral(
    f::F,
    geometry::G,
    settings::I
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationAlgorithm}
    Dim = Meshes.paramdim(geometry)

    if Dim == 1
        return integral(f, geometry, settings)
    else
        error("Performing a line integral on a geometry with $Dim parametric dimensions not supported.")
    end
end

function lineintegral(
    f::F,
    geometry::G,
    settings::I,
    FP::Type{T}
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationAlgorithm, T<:AbstractFloat}
    Dim = Meshes.paramdim(geometry)

    if Dim == 1
        return integral(f, geometry, settings, FP)
    else
        error("Performing a line integral on a geometry with $Dim parametric dimensions not supported.")
    end
end


################################################################################
#                              Surface Integral
################################################################################

"""
    surfaceintegral(f, geometry)
    surfaceintegral(f, geometry, algorithm)
    surfaceintegral(f, geometry, algorithm, FP)

Numerically integrate a given function `f(::Point)` along a surface `geometry`
using a particular `integration algorithm` with floating point precision of type
`FP`.

Algorithm types available:
- GaussKronrod
- GaussLegendre
- HAdaptiveCubature (default)
"""
function surfaceintegral(
    f::F,
    geometry::G,
) where {F<:Function, G<:Meshes.Geometry}
    Dim = Meshes.paramdim(geometry)

    if Dim == 2
        return integral(f, geometry, GaussKronrod())
    else
        error("Performing a surface integral on a geometry with $Dim parametric dimensions not supported.")
    end
end

function surfaceintegral(
    f::F,
    geometry::G,
    settings::I
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationAlgorithm}
    Dim = Meshes.paramdim(geometry)

    if Dim == 2
        return integral(f, geometry, settings)
    else
        error("Performing a surface integral on a geometry with $Dim parametric dimensions not supported.")
    end
end

function surfaceintegral(
    f::F,
    geometry::G,
    settings::I,
    FP::Type{T}
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationAlgorithm, T<:AbstractFloat}
    Dim = Meshes.paramdim(geometry)

    if Dim == 2
        return integral(f, geometry, settings, FP)
    else
        error("Performing a surface integral on a geometry with $Dim parametric dimensions not supported.")
    end
end


################################################################################
#                              Volume Integral
################################################################################

"""
    volumeintegral(f, geometry)
    volumeintegral(f, geometry, algorithm)
    volumeintegral(f, geometry, algorithm, FP)

Numerically integrate a given function `f(::Point)` throughout a volumetric
`geometry` using a particular `integration algorithm` with floating point
precision of type `FP`.

Algorithm types available:
- GaussKronrod
- GaussLegendre
- HAdaptiveCubature (default)
"""
function volumeintegral(
    f::F,
    geometry::G,
) where {F<:Function, G<:Meshes.Geometry}
    Dim = Meshes.paramdim(geometry)

    if Dim == 3
        return integral(f, geometry, GaussKronrod())
    else
        error("Performing a volume integral on a geometry with $Dim parametric dimensions not supported.")
    end
end

function volumeintegral(
    f::F,
    geometry::G,
    settings::I
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationAlgorithm}
    Dim = Meshes.paramdim(geometry)

    if Dim == 3
        return integral(f, geometry, settings)
    else
        error("Performing a volume integral on a geometry with $Dim parametric dimensions not supported.")
    end
end

function volumeintegral(
    f::F,
    geometry::G,
    settings::I,
    FP::Type{T}
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationAlgorithm, T<:AbstractFloat}
    Dim = Meshes.paramdim(geometry)

    if Dim == 3
        return integral(f, geometry, settings, FP)
    else
        error("Performing a volume integral on a geometry with $Dim parametric dimensions not supported.")
    end
end
