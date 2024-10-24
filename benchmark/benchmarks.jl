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
    (name = "GaussLegendre", rule = GaussLegendre(100), N = 100),
    (name = "GaussKronrod", rule = GaussKronrod(), N = 100),
    (name = "HAdaptiveCubature", rule = HAdaptiveCubature(), N = 500)
)
geometries = (
    (name = "Meshes.Segment", item = Segment(Point(0, 0, 0), Point(1, 1, 1))),
    (name = "Meshes.Sphere", item = Sphere(Point(0, 0, 0), 1.0))
)

SUITE["Integrals"] = let s = BenchmarkGroup()
    for (int, rule, geometry) in Iterators.product(integrands, rules, geometries)
        n1, n2 = geometry.name,"$(int.name) $(rule.name)"
        N = rule.N
        s[n1, n2] = @benchmarkable integral($int.f, $geometry.item, $rule.rule) evals=N
    end
    s
end

#tune!(SUITE)
#run(SUITE, verbose=true)
