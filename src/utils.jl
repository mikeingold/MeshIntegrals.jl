# Get the derivative vector of a Bezier Curve at parametric t
#   Ref https://en.wikipedia.org/wiki/B%C3%A9zier_curve#Derivative
function derivative(bz::Meshes.BezierCurve{Dim,T,V}, t) where {Dim,T,V}
	P = bz.controls
	N = length(bz.controls)
	
	# Binomial coefficient
	bc(n,i) = factorial(n) / (factorial(i)*factorial(n-i))
	
	# Bernstein basis polynomials
	b(i,n) = t -> bc(n,i) * t^i * (1-t)^(n-i)
	
	# Term to be summed
	sigma(i) = b(i,N-1)(t) .* (P[i+2] - P[i+1])
	
	# Derivative of curve
	return N .* sum(sigma, 0:(N-1))
end
