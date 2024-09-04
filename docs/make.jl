using Documenter
using MeshIntegrals

makedocs(
    sitename="MeshIntegrals.jl",
    pages = [
        "Home" => [
            "About" => "index.md",
            "Support Matrix" => "supportmatrix.md",
            "Integration Algorithms" => "algorithms.md"
        ],
        "Derivations" => [
            "Integrating a Triangle" => "triangle.md"
        ],
        "Public API" => "api.md"
    ]
)

deploydocs(repo = "github.com/mikeingold/MeshIntegrals.jl.git")
