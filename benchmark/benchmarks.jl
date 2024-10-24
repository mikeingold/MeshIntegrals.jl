using BenchmarkTools
using Meshes
using MeshIntegrals
using Unitful

const SUITE = BenchmarkGroup()

SUITE["Integrals"] = BenchmarkGroup()

SUITE["Integrals"]["Meshes.Segment"] = let s = BenchmarkGroup()
    # Connect a line segment from the origin to an arbitrary point on the unit sphere
    φ, θ = (7pi / 6, pi / 3)  # Arbitrary spherical angles
    pt_a = Point(0.0u"m", 0.0u"m", 0.0u"m")
    pt_b = Point(sin(θ) * cos(φ) * u"m", sin(θ) * sin(φ) * u"m", cos(θ) * u"m")
    segment = Segment(pt_a, pt_b)

    f(p) = norm(to(p))
    fv(p) = fill(f(p), 3)

    s["Scalar-valued GaussLegendre"] = @benchmarkable integral($f, $segment, GaussLegendre(100)) evals=1000
    s["Scalar-valued GaussKronrod"] = @benchmarkable integral($f, $segment, GaussKronrod()) evals=1000
    s["Scalar-valued HAdaptiveCubature"] = @benchmarkable integral($f, $segment, HAdaptiveCubature()) evals=1000

    s["Vector-valued GaussLegendre"] = @benchmarkable integral($fv, $segment, GaussLegendre(100)) evals=1000
    s["Vector-valued GaussKronrod"] = @benchmarkable integral($fv, $segment, GaussKronrod()) evals=1000
    s["Vector-valued HAdaptiveCubature"] = @benchmarkable integral($fv, $segment, HAdaptiveCubature()) evals=1000

    s
end
