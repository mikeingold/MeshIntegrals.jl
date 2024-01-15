################################################################################
#                              HELPER FUNCTIONS
################################################################################

# Validate that f is a f(::Point)
function _validate_integrand_point(f,Dim,T)
	# Verify that the provided function has the correct method available
	if !hasmethod(f, (Point{Dim,T},))
		error("The provided Function f must have a method f(::Point{$Dim,$T})")
	end

	return nothing
end

# Validate that f is a f(::Point,::Vec)
function _validate_integrand_pointvec(f,Dim,T)
	# Verify that the provided function has the correct method available
	if !hasmethod(f, (Point{Dim,T}, Vec{Dim,T}))
		error("The provided Function f must have a method f(::Point{$Dim,$T}, ::Vec{$Dim,$T})")
	end

	return nothing
end

################################################################################
#                         INTEGRALS OF f(position)
################################################################################

# Integrate f(::Point{Dim,T}) over an arbitrary geometry construct
function integrate(f::F, path::Vector{<:Meshes.Geometry{Dim,T}}; kwargs...) where {F<:Function,Dim,T}
	# Validate the provided integrand function
	_validate_integrand_point(f,Dim,T)

	return sum(line -> integrate(f, line; kwargs...), path)
end

# Integrate f(::Point{Dim,T}) over a Segment
function integrate(f::F, segment::Meshes.Segment{Dim,T}; kwargs...) where {F<:Function,Dim,T}
	# Validate the provided integrand function
	_validate_integrand_point(f,Dim,T)
	
    return length(segment) * quadgk(t -> f(segment(t)), 0, 1; kwargs...)[1]
end

# Integrate f(::Point{Dim,T}) over a Rope (an open Chain)
function integrate(f::F, rope::Meshes.Rope{Dim,T}; kwargs...) where {F<:Function,Dim,T}
	# Validate the provided integrand function
	_validate_integrand_point(f,Dim,T)
	
    sum(segment -> integrate(f, segment; kwargs...), segments(rope))
end

# Integrate f(::Point{Dim,T}) over a Ring (a closed Chain)
function integrate(f::F, ring::Meshes.Ring{Dim,T}; kwargs...) where {F<:Function,Dim,T}
	# Validate the provided integrand function
	_validate_integrand_point(f,Dim,T)
	
    sum(segment -> integrate(f, segment; kwargs...), segments(ring))
end

# Integrate f(::Point{Dim,T}) over a BezierCurve
function integrate(f::F, curve::Meshes.BezierCurve{Dim,T,V}; kwargs...) where {F<:Function,Dim,T,V}
	# Validate the provided integrand function
	_validate_integrand_point(f,Dim,T)
	
    return length(curve) * quadgk(t -> f(curve(t)), 0, 1; kwargs...)[1]
end

################################################################################
#                       QuadGK.quadgk Methods
################################################################################

# Implement QuadGK.quadgk over a Meshes.Segment
function QuadGK.quadgk(f::F, segment::Meshes.Segment{Dim,T}; kwargs...) where {F<:Function,Dim,T}
	# Validate the provided integrand function
	_validate_integrand_point(f,Dim,T)

	len = length(segment)
	return quadgk(t -> len * f(segment(t)), 0, 1; kwargs...)
end

# Implement QuadGK.quadgk over a Meshes.BezierCurve
function QuadGK.quadgk(f::F, curve::Meshes.BezierCurve{Dim,T,V}; kwargs...) where {F<:Function,Dim,T,V}
	# Validate the provided integrand function
	_validate_integrand_point(f,Dim,T)

	len = length(segment)
    return quadgk(t -> len * f(curve(t)), 0, 1; kwargs...)
end

################################################################################
#                       INTEGRALS OF f(position, normal)
################################################################################

# Integrate f(::Point{Dim,T}, ::Vec{Dim,T}) over a SurfaceSegment
function integrate(f::F, ss::SurfacePathSegment{Dim,T}) where {F<:Function,Dim,T}
	# Validate the provided integrand function
	_validate_integrand_pointvec(f,Dim,T)
	
    return length(ss.segment) * quadgk(t -> f(ss.segment(t),ss.normal), 0, 1)[1]
end

# Integrate f(::Point{Dim,T}, ::Vec{Dim,T}) over a SurfaceTrajectory
function integrate(f::F, traj::SurfaceTrajectory{Dim,T}) where {F<:Function,Dim,T}
	# Validate the provided integrand function
	_validate_integrand_pointvec(f,Dim,T)
	
    return sum(ss -> integrate(ss,f), traj.path)
end
