push!(LOAD_PATH,"../src/")

using Documenter, NESSie

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
    modules   = [NESSie, NESSie.BEM, NESSie.Born, NESSie.Format],
    clean     = true,
    doctest   = true,
    linkcheck = true,
    checkdocs = :all,
    pages     = pages,
    format    = :html,
    sitename  = "NESSie.jl",
    repo      = "https://github.com/tkemmer/nessie.jl/tree/master{path}"
)
