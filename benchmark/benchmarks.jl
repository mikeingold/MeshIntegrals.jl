using BenchmarkTools
using LinearAlgebra
using Meshes
using MeshIntegrals

const SUITE = BenchmarkGroup()

############################################################################################
#                                      Integrals
############################################################################################

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
        n1, n2, N = geometry.name, "$(int.name) $(rule.name)", rule.N
        s[n1][n2] = @benchmarkable integral($int.f, $geometry.item, $rule.rule)
    end
    s
end

############################################################################################
#                                    Specializations
############################################################################################

# TODO after merge of PR #131
#   spec.f = p -> norm(to(p))
#   un-comment s["Tetrahedron"]

spec = (
    f = p -> 1.0,
    g = (
        bezier = BezierCurve([Point(t, sin(t), 0) for t in range(-pi, pi, length = 361)]),
        triangle = Triangle(Point(1, 0, 0), Point(0, 1, 0), Point(0, 0, 1)),
        tetrahedron = let
            a = Point(0, 3, 0)
            b = Point(-7, 0, 0)
            c = Point(8, 0, 0)
            ẑ = Vec(0, 0, 1)
            Tetrahedron(a, b, c, a + ẑ)
        end
    ),
    rule = HAdaptiveCubature()
)

SUITE["Specializations"] = let s = BenchmarkGroup()
    s["BezierCurve"] = @benchmarkable integral($spec.f, $spec.g.bezier, $spec.rule)
    s["Triangle"] = @benchmarkable integral($spec.f, $spec.g.triangle, $spec.rule)
    # s["Tetrahedron"] = @benchmarkable integral($spec.f, $spec.g.tetrahedron, $spec.rule)
    s
end

############################################################################################
#                                      Differentials
############################################################################################

sphere = Sphere(Point(0, 0, 0), 1.0)
differential = MeshIntegrals.differential

SUITE["Differentials"] = let s = BenchmarkGroup()
    s["Jacobian"] = @benchmarkable jacobian($sphere, $(0.5, 0.5)) evals=10
    s["Differential"] = @benchmarkable differential($sphere, $(0.5, 0.5)) evals=10
    s
end

#tune!(SUITE)
#run(SUITE, verbose=true)
