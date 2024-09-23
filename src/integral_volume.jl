################################################################################
#                       Generalized 3D Methods
################################################################################

# Integrating volumes with GaussKronrod not supported by default
function _integral_3d(
    f,
    geometry,
    settings::GaussKronrod,
    FP::Type{T} = Float64
) where {T<:AbstractFloat}
    error("Integrating this volume type with GaussKronrod not supported.")
end
