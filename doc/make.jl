push!(LOAD_PATH,"../src/")

using Documenter, ProteinES

const pages = [
    "Home" => "index.md",
    "Manual" => [
        "man/introduction.md",
        "man/references.md"
    ],
    "ProteinES module" => [
        "base/born.md",
        "base/constants.md",
        "base/model.md",
        "base/potentials.md",
        "base/quadrature.md",
        "base/util.md"
    ],
    "ProteinES.IO module" => [
        "io/input.md",
        "io/output.md"
    ]
]

makedocs(
    modules   = [ProteinES, ProteinES.IO],
    clean     = true,
    doctest   = true,
    linkcheck = true,
    checkdocs = :all,
    pages     = pages,
    format    = :html,
    sitename  = "ProteinES.jl",
    repo      = "../..{path}"
)
