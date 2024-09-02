using Documenter
using MeshIntegrals

makedocs(
    sitename="MeshIntegrals.jl",
    pages = [
        "index.md",
        "Derivations" => [
            "Integrating a Triangle" => "triangle.md"
        ]
    ]
)

deploydocs(repo = "github.com/mikeingold/MeshIntegrals.jl.git")
