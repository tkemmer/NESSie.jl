push!(LOAD_PATH,"../src/")

using Documenter, ProteinES

const pages = [
    "Home" => "index.md",
    "Manual" => [
        "man/introduction.md",
        "man/references.md"
    ],
    "ProteinES" => [
        "base/constants.md",
        "base/model.md",
        "base/util.md"
    ]
]

makedocs(
    modules   = [ProteinES],
    clean     = true,
    doctest   = true,
    linkcheck = true,
    checkdocs = :all,
    pages     = pages,
    format    = :html,
    sitename  = "ProteinES.jl",
    repo      = "../..{path}"
)
