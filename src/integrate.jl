################################################################################
#                              HELPER FUNCTIONS
################################################################################

# Validate that f is a f(::Point)
function _validate_integrand_point(f,Dim,T)
	# Verify that the provided function has the correct method available
	if !hasmethod(f, (Point{Dim,T},))
		error("The provided Function f has no method f(::Point{Dim,T})")
	end

	return nothing
end

# Validate that f is a f(::Point,::Vec)
function _validate_integrand_pointvec(f,Dim,T)
	# Verify that the provided function has the correct method available
	if !hasmethod(f, (Point{Dim,T}, Vec{Dim,T}))
		error("The provided Function f has no method f(::Point{Dim,T})")
	end

	return nothing
end

################################################################################
#                         INTEGRALS OF f(position)
################################################################################

# Integrate f(::Point{Dim,T}) over a Segment
function integrate(f::F, segment::Meshes.Segment{Dim,T}) where {F<:Function,Dim,T}
	# Validate the provided integrand function
	_validate_integrand_point(f,Dim,T)
	
    return length(segment) * quadgk(t -> f(segment(t)), 0, 1)[1]
end

# Integrate f(::Point{Dim,T}) over a BezierCurve
function integrate(f::F, curve::Meshes.BezierCurve{Dim,T,V}) where {F<:Function,Dim,T}
	# Validate the provided integrand function
	_validate_integrand_point(f)
	
    return length(curve) * quadgk(t -> f(curve(t)), 0, 1)[1]
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
function integrate(f::F, traj::SurfaceTrajectory)
	# Validate the provided integrand function
	_validate_integrand_pointvec(f,Dim,T)
	
    return sum(ss -> integrate(ss,f), traj.path)
end
