push!(LOAD_PATH,"../src/")

using Documenter, ProteinES

const pages = [
    "Home" => "index.md",
    "Manual" => [
        "man/introduction.md",
        "man/references.md"
    ],
    "Data formats" => [
        "io/input.md",
        "io/output.md"
    ],
    "Library" => [
        "base/constants.md",
        "base/model.md",
        "base/potentials.md",
        "base/quadrature.md",
        "base/util.md"
    ],
    "Internals" => [
        "intern/constants.md",
        "intern/input.md",
        "intern/util.md"
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
