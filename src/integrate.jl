# INTEGRALS OF f(position)

# Integrate f(::Point{Dim,T}) over a Segment
function integrate(segment::Meshes.Segment{Dim,T}, f) where {Dim,T}
	# Verify that the provided function has the correct method available
	if !hasmethod(f, (Point{Dim,T},))
		error("The provided Function f has no method f(::Point{Dim,T})")
	end
	
    return length(segment) * quadgk(t -> f(segment(t)), 0, 1)[1]
end

# Integrate f(::Point{Dim,T}) over a BezierCurve
function integrate(curve::Meshes.BezierCurve{Dim,T,V}, f) where {Dim,T,V}
	# Verify that the provided function has the correct method available
	if !hasmethod(f, (Point{Dim,T},))
		error("The provided Function f has no method f(::Point{Dim,T})")
	end
	
    return length(curve) * quadgk(t -> f(curve(t)), 0, 1)[1]
end

# INTEGRALs OF f(position, normal)

# Integrate f(::Point{Dim,T}, ::Vec{Dim,T}) over a SurfaceSegment
function integrate(ss::SurfaceSegment{Dim,T}, f) where {Dim,T}
	# Verify that the provided function has the correct method available
	if !hasmethod(f, (Point{Dim,T}, Vec{Dim,T}))
		error("The provided Function f has no method f(::Point{Dim,T})")
	end
	
    return length(ss.segment) * quadgk(t -> f(ss.segment(t),ss.normal), 0, 1)[1]
end

# Integrate f(::Point{Dim,T}, ::Vec{Dim,T}) over a SurfaceTrajectory
function integrate(traj::SurfaceTrajectory, f)
	# Verify that the provided function has the correct method available
	if !hasmethod(f, (Point{Dim,T}, Vec{Dim,T}))
		error("The provided Function f has no method f(::Point{Dim,T})")
	end
	
    return sum(ss -> integrate(ss,f), traj.path)
end
