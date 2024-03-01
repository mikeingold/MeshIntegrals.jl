################################################################################
#                            Common Methods
################################################################################

function integral(
    f::F,
    ring::Meshes.Ring,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Ring into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(f, segment, settings), segments(ring))
end

function integral(
    f::F,
    rope::Meshes.Rope,
    settings::I
) where {F<:Function, I<:IntegrationAlgorithm}
    # Convert the Rope into Segments, sum the integrals of those 
    return sum(segment -> lineintegral(f, segment, settings), segments(rope))
end

function lineintegral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::I;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, Dim, T, V, I<:IntegrationAlgorithm}
    return integral(f, curve, settings; alg=alg)
end


################################################################################
#                            Gauss-Legendre
################################################################################

# Generalized method
function _integral_1d(f, geometry, settings::GaussLegendre)
    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x) = geometry(t(x))

    function paramfactor(x)
        J = jacobian(geometry,[t(x)])
        return norm(J[1])
    end

    # Integrate f along the geometry and apply a domain-correction factor for [-1,1] ↦ [0, 1]
    integrand((w,x)) = w * f(point(x)) * paramfactor(x)
    return 0.5 * sum(integrand, zip(ws, xs))
end

"""
    integral(f, curve::Meshes.BezierCurve, ::GaussLegendre; alg=Meshes.Horner())

Like [`integral`](@ref) but integrates along the domain defined a `curve`. By
default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function integral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::GaussLegendre;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function, Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(settings.n)

    # Change of variables: x [-1,1] ↦ t [0,1]
    t(x) = 0.5x + 0.5
    point(x) = curve(t(x), alg)

    # Integrate f along the line and apply a domain-correction factor for [-1,1] ↦ [0, length]
    return 0.5 * length(curve) * sum(w .* f(point(x)) for (w,x) in zip(ws, xs))
end

function integral(
    f::F,
    line::Meshes.Box{1,T},
    settings::GaussLegendre
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Box{1,T} is definitionally embedded in 1D-space
    _validate_integrand(f,1,T)

    return _integral_1d(f, line, settings)
end

function integral(
    f::F,
    circle::Meshes.Circle{T},
    settings::GaussLegendre
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Circle is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    return _integral_1d(f, circle, settings)
end

function integral(
    f::F,
    line::Meshes.Line{Dim,T},
    settings::GaussLegendre
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Compute Gauss-Legendre nodes/weights for x in interval [-1,1]
    xs, ws = gausslegendre(settings.n)

    # Get domain-corrected parametric locator
    len = length(Segment(line.a,line.b))
    point(t) = line(t/len)

    # Change of variables: t [-Inf,Inf] ↦ x [-1,1]
    integrand(x) = f(point((x/(1-x^2)))) * (1+x^2)/((1-x^2)^2)

    # Integrate f along the line
    return sum(w .* integrand(x) for (w,x) in zip(ws, xs))
end

function integral(
    f::F,
    segment::Meshes.Segment{Dim,T},
    settings::GaussLegendre
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    return _integral_1d(f, segment, settings)
end

function integral(
    f::F,
    circle::Meshes.Sphere{2,T},
    settings::GaussLegendre
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Sphere{2,T} is simply a circle in 2D-space
    _validate_integrand(f,2,T)

    return _integral_1d(f, circle, settings)
end


################################################################################
#                                Gauss-Kronrod
################################################################################

"""
    integral(f, curve::BezierCurve, ::GaussKronrod; alg=Horner(), kws...)

Like [`integral`](@ref) but integrates along the domain defined a `curve`. By
default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function integral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::GaussKronrod;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function,Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(curve)
    point(t) = curve(t, alg)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end

function integral(
    f::F,
    line::Meshes.Box{1,T},
    settings::GaussKronrod
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Box is definitionally embedded in 1D-space
    _validate_integrand(f,1,T)

    len = length(line)
    point(t) = line(t)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end

function integral(
    f::F,
    circle::Meshes.Circle{T},
    settings::GaussKronrod
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Circle is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    len = length(circle)
    point(t) = circle(t)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end

function integral(
    f::F,
    line::Meshes.Line{Dim,T},
    settings::GaussKronrod
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Get domain-corrected parametric locator
    len = length(Segment(line.a,line.b))
    point(t) = line(t/len)

    # Lines are infinite-length passing through defined points a and b
    return QuadGK.quadgk(t -> f(point(t)), -Inf, Inf; settings.kwargs...)[1]
end

function integral(
    f::F,
    segment::Meshes.Segment{Dim,T},
    settings::GaussKronrod
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(segment)
    point(t) = segment(t)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end

function integral(
    f::F,
    circle::Meshes.Sphere{2,T},
    settings::GaussKronrod
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Sphere{2,T} is simply a circle in 2D-space
    _validate_integrand(f,2,T)

    len = length(circle)
    point(t) = circle(t)
    return QuadGK.quadgk(t -> len * f(point(t)), 0, 1; settings.kwargs...)[1]
end


################################################################################
#                                HCubature
################################################################################

"""
    integral(f, curve::BezierCurve, ::HAdaptiveCubature; alg=Horner(), kws...)

Like [`integral`](@ref) but integrates along the domain defined a `curve`. By
default this uses Horner's method to improve performance when parameterizing
the `curve` at the expense of a small loss of precision. Additional accuracy
can be obtained by specifying the use of DeCasteljau's algorithm instead with
`alg=Meshes.DeCasteljau()` but can come at a steep cost in memory allocations,
especially for curves with a large number of control points.
"""
function integral(
    f::F,
    curve::Meshes.BezierCurve{Dim,T,V},
    settings::HAdaptiveCubature;
    alg::Meshes.BezierEvalMethod=Meshes.Horner()
) where {F<:Function,Dim, T, V}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(curve)
    point(t) = curve(t, alg)
    return hcubature(t -> len * f(point(t[1])), [0], [1]; settings.kwargs...)[1]
end

function integral(
    f::F,
    line::Meshes.Box{1,T},
    settings::HAdaptiveCubature
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Box is definitionally embedded in 1D-space
    _validate_integrand(f,1,T)

    len = length(line)
    point(t) = line(t)
    return hcubature(t -> len * f(point(t[1])), [0], [1]; settings.kwargs...)[1]
end

function integral(
    f::F,
    circle::Meshes.Circle{T},
    settings::HAdaptiveCubature
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Circle is definitionally embedded in 3D-space
    _validate_integrand(f,3,T)

    len = length(circle)
    point(t) = circle(t)
    return hcubature(t -> len * f(point(t[1])), [0], [1]; settings.kwargs...)[1]
end

function integral(
    f::F,
    line::Meshes.Line{Dim,T},
    settings::HAdaptiveCubature
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    # Get domain-corrected parametric locator
    len = length(Segment(line.a,line.b))
    point(t) = line(t/len)
    
    # Change of variables: t [-Inf,Inf] ↦ x [-1,1]
    integrand(x) = f(point((x/(1-x^2)))) * (1+x^2)/((1-x^2)^2)
    integrand(x::AbstractVector) = integrand(x[1])

    # Lines are infinite-length passing through defined points a and b
    return hcubature(integrand, [-1], [1]; settings.kwargs...)[1]
end

function integral(
    f::F,
    segment::Meshes.Segment{Dim,T},
    settings::HAdaptiveCubature
) where {F<:Function, Dim, T}
    # Validate the provided integrand function
    _validate_integrand(f,Dim,T)

    len = length(segment)
    point(t) = segment(t)
    return hcubature(t -> len * f(point(t[1])), [0], [1]; settings.kwargs...)[1]
end

function integral(
    f::F,
    circle::Meshes.Sphere{2,T},
    settings::HAdaptiveCubature
) where {F<:Function, T}
    # Validate the provided integrand function
    # A Sphere{2,T} is simply a circle in 2D-space
    _validate_integrand(f,2,T)

    len = length(circle)
    point(t) = circle(t)
    return hcubature(t -> len * f(point(t[1])), [0], [1]; settings.kwargs...)[1]
end
