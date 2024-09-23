################################################################################
#                              Line Integral
################################################################################

"""
    lineintegral(f, geometry)
    lineintegral(f, geometry, rule)
    lineintegral(f, geometry, rule, FP)

Numerically integrate a given function `f(::Point)` along a line-like `geometry`
using a particular numerical `integration rule` with floating point precision of
type `FP`.

Rule types available:
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
    rule::I
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationRule}
    Dim = Meshes.paramdim(geometry)

    if Dim == 1
        return integral(f, geometry, rule)
    else
        error("Performing a line integral on a geometry with $Dim parametric dimensions not supported.")
    end
end

function lineintegral(
    f::F,
    geometry::G,
    rule::I,
    FP::Type{T}
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationRule, T<:AbstractFloat}
    Dim = Meshes.paramdim(geometry)

    if Dim == 1
        return integral(f, geometry, rule, FP)
    else
        error("Performing a line integral on a geometry with $Dim parametric dimensions not supported.")
    end
end


################################################################################
#                              Surface Integral
################################################################################

"""
    surfaceintegral(f, geometry)
    surfaceintegral(f, geometry, rule)
    surfaceintegral(f, geometry, rule, FP)

Numerically integrate a given function `f(::Point)` along a surface `geometry`
using a particular numerical `integration rule` with floating point precision of
type `FP`.

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
        return integral(f, geometry, HAdaptiveCubature())
    else
        error("Performing a surface integral on a geometry with $Dim parametric dimensions not supported.")
    end
end

function surfaceintegral(
    f::F,
    geometry::G,
    rule::I
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationRule}
    Dim = Meshes.paramdim(geometry)

    if Dim == 2
        return integral(f, geometry, rule)
    else
        error("Performing a surface integral on a geometry with $Dim parametric dimensions not supported.")
    end
end

function surfaceintegral(
    f::F,
    geometry::G,
    rule::I,
    FP::Type{T}
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationRule, T<:AbstractFloat}
    Dim = Meshes.paramdim(geometry)

    if Dim == 2
        return integral(f, geometry, rule, FP)
    else
        error("Performing a surface integral on a geometry with $Dim parametric dimensions not supported.")
    end
end


################################################################################
#                              Volume Integral
################################################################################

"""
    volumeintegral(f, geometry)
    volumeintegral(f, geometry, rule)
    volumeintegral(f, geometry, rule, FP)

Numerically integrate a given function `f(::Point)` throughout a volumetric
`geometry` using a particular numerical `integration rule` with floating point
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
        return integral(f, geometry, HAdaptiveCubature())
    else
        error("Performing a volume integral on a geometry with $Dim parametric dimensions not supported.")
    end
end

function volumeintegral(
    f::F,
    geometry::G,
    rule::I
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationRule}
    Dim = Meshes.paramdim(geometry)

    if Dim == 3
        return integral(f, geometry, rule)
    else
        error("Performing a volume integral on a geometry with $Dim parametric dimensions not supported.")
    end
end

function volumeintegral(
    f::F,
    geometry::G,
    rule::I,
    FP::Type{T}
) where {F<:Function, G<:Meshes.Geometry, I<:IntegrationRule, T<:AbstractFloat}
    Dim = Meshes.paramdim(geometry)

    if Dim == 3
        return integral(f, geometry, rule, FP)
    else
        error("Performing a volume integral on a geometry with $Dim parametric dimensions not supported.")
    end
end
