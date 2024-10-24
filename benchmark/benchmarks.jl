using BenchmarkTools
using LinearAlgebra: norm
using Meshes
using MeshIntegrals
using Unitful

const SUITE = BenchmarkGroup()
SUITE["Integrals"] = BenchmarkGroup()

segment = Segment(Point(0, 0, 0), Point(1, 1, 1))
f(p) = norm(to(p))
fv(p) = fill(f(p), 3)

SUITE["Integrals"]["Meshes.Segment"] = BenchmarkGroup()
let s = SUITE["Integrals"]["Meshes.Segment"]
    s["Scalar GaussLegendre"] = @benchmarkable integral($f, $segment, GaussLegendre(100)) evals=1000
    s["Scalar GaussKronrod"] = @benchmarkable integral($f, $segment, GaussKronrod()) evals=1000
    s["Scalar HAdaptiveCubature"] = @benchmarkable integral($f, $segment, HAdaptiveCubature()) evals=1000

    s["Vector GaussLegendre"] = @benchmarkable integral($fv, $segment, GaussLegendre(100)) evals=1000
    s["Vector GaussKronrod"] = @benchmarkable integral($fv, $segment, GaussKronrod()) evals=1000
    s["Vector HAdaptiveCubature"] = @benchmarkable integral($fv, $segment, HAdaptiveCubature()) evals=1000
end
