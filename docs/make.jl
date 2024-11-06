using Documenter
using MeshIntegrals

# Dynamically set all files in subdirectories of the source directory to include all files in these subdirectories.
# This way they don't need to be listed explicitly.
SPECIALIZATIONS_FILES = joinpath.(Ref("specializations"),
    readdir(joinpath(dirname(@__DIR__), "src", "specializations")))

makedocs(
    sitename = "MeshIntegrals.jl",
    pages = [
        "Home" => [
            "About" => "index.md",
            "Support Matrix" => "supportmatrix.md",
            "Example Usage" => "usage.md"
        ],
        "Developer Notes" => [
            "How it Works" => "how_it_works.md",
            "Integrating a Triangle" => "triangle.md"
        ],
        "Public API" => "api.md"
    ]
)

deploydocs(repo = "github.com/JuliaGeometry/MeshIntegrals.jl.git",
    devbranch = "main",
    push_preview = true)
