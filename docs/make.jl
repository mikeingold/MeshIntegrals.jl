using Documenter
using MeshIntegrals

makedocs(
  sitename="MeshIntegrals.jl",
  pages=[
    "Home" => ["About" => "index.md", "Support Matrix" => "supportmatrix.md", "Example Usage" => "usage.md"],
    "Derivations" => ["Integrating a Triangle" => "triangle.md"],
    "Public API" => "api.md"
  ]
)

deploydocs(repo="github.com/mikeingold/MeshIntegrals.jl.git", devbranch="main", push_preview=true)
