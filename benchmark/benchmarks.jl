using BenchmarkTools
using LinearAlgebra
using Meshes
using MeshIntegrals
using Unitful

const SUITE = BenchmarkGroup()

integrands = (
    (name = "Scalar", f = p -> norm(to(p))),
    (name = "Vector", f = p -> fill(norm(to(p)), 3))
)
rules = (
    (name = "GaussLegendre", rule = GaussLegendre(100)),
    (name = "GaussKronrod", rule = GaussKronrod()),
    (name = "HAdaptiveCubature", rule = HAdaptiveCubature())
)
geometries = (
    (name = "Meshes.Segment", item = Segment(Point(0, 0, 0), Point(1, 1, 1))),
    (name = "Meshes.Sphere", item = Sphere(Point(0, 0, 0), 1.0))
)

SUITE["Integrals"] = let s = BenchmarkGroup()
    for (int, rule, g) in Iterators.product(integrands, rules, geometries)
        n1 = geometry.name
        n2 = "$(int.name) $(rule.name)"
        s[n1, n2] = @benchmarkable integral($int.f, $geometry.item, $rule.rule) evals=1000
    end
    s
end

#tune!(SUITE)
#run(SUITE, verbose=true)
