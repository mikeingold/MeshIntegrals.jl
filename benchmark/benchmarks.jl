using BenchmarkTools
using LinearAlgebra
using Meshes
using MeshIntegrals
using Unitful

const SUITE = BenchmarkGroup()
SUITE["Integrals"] = BenchmarkGroup()

segment = Segment(Point(0, 0, 0), Point(1, 1, 1))
f(p) = norm(to(p))
fv(p) = fill(f(p), 3)

integrands = (
    (name = "Scalar", f = f),
    (name = "Vector", f = fv)
)
rules = (
    (name = "GaussLegendre", rule = GaussLegendre(100)),
    (name = "GaussKronrod", rule = GaussKronrod()),
    (name = "HAdaptiveCubature", rule = HAdaptiveCubature())
)

SUITE["Integrals"]["Meshes.Segment"] = let s = BenchmarkGroup()
    for (int, rule) in Iterators.product(integrands, rules)
        name = "$(int.name) $(rule.name)"
        s[name] = @benchmarkable integral($int.f, $segment, $rule.rule) evals=1000
    end
    s
end

#tune!(SUITE)
#run(SUITE, verbose=true)
