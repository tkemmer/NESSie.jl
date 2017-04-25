push!(LOAD_PATH,"../src/")

using Documenter, ProteinES

const pages = [
    "Home" => "index.md",
    "Manual" => [
        "man/introduction.md",
        "man/references.md"
    ],
    "Data formats" => [
        "formats/input.md",
        "formats/output.md"
    ],
    "Library" => [
        "lib/constants.md",
        "lib/electrostatics.md",
        "lib/models.md",
        "lib/quadrature.md",
        "lib/solvers.md",
        "lib/util.md"
    ],
    "Internals" => [
        "intern/constants.md",
        "intern/input.md",
        "intern/quadrature.md",
        "intern/util.md"
    ]
]

makedocs(
    modules   = [ProteinES, ProteinES.BEM, ProteinES.Born, ProteinES.IO],
    clean     = true,
    doctest   = true,
    linkcheck = true,
    checkdocs = :all,
    pages     = pages,
    format    = :html,
    sitename  = "ProteinES.jl",
    repo      = "../..{path}"
)
