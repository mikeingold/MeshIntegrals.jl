################################################################################
#                              Line Integral
################################################################################

"""
    lineintegral(f, geometry[, rule]; FP=Float64)

Numerically integrate a given function `f(::Point)` along a line-like `geometry`
using a particular numerical integration `rule` with floating point precision of
type `FP`.

Rule types available:
- GaussKronrod (default)
- GaussLegendre
- HAdaptiveCubature
"""
function lineintegral(
        f::Function,
        geometry::Geometry,
        rule::IntegrationRule = GaussKronrod();
        kwargs...
)
    N = Meshes.paramdim(geometry)

    if N == 1
        return integral(f, geometry, rule; kwargs...)
    else
        throw(ArgumentError("Performing a line integral on a geometry \
                            with $N parametric dimensions not supported."))
    end
end

################################################################################
#                              Surface Integral
################################################################################

"""
    surfaceintegral(f, geometry[, rule]; FP=Float64)

Numerically integrate a given function `f(::Point)` along a surface `geometry`
using a particular numerical integration `rule` with floating point precision of
type `FP`.

Algorithm types available:
- GaussKronrod
- GaussLegendre
- HAdaptiveCubature (default)
"""
function surfaceintegral(
        f::Function,
        geometry::Geometry,
        rule::IntegrationRule = HAdaptiveCubature();
        kwargs...
)
    N = Meshes.paramdim(geometry)

    if N == 2
        return integral(f, geometry, rule; kwargs...)
    else
        throw(ArgumentError("Performing a surface integral on a geometry \
                            with $N parametric dimensions not supported."))
    end
end

################################################################################
#                              Volume Integral
################################################################################

"""
    volumeintegral(f, geometry[, rule]; FP=Float64)

Numerically integrate a given function `f(::Point)` throughout a volumetric
`geometry` using a particular numerical integration `rule` with floating point
precision of type `FP`.

Algorithm types available:
- GaussKronrod
- GaussLegendre
- HAdaptiveCubature (default)
"""
function volumeintegral(
        f::Function,
        geometry::Geometry,
        rule::IntegrationRule = HAdaptiveCubature();
        kwargs...
)
    N = Meshes.paramdim(geometry)

    if N == 3
        return integral(f, geometry, rule; kwargs...)
    else
        throw(ArgumentError("Performing a volume integral on a geometry \
                            with $N parametric dimensions not supported."))
    end
end
