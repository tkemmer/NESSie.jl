using Documenter, NESSie

const pages = [
    "Home" => "index.md",
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
        "intern/base.md",
        "intern/bem.md",
        "intern/format.md",
        "intern/radon.md",
        "intern/rjasanow.md",
        "intern/testmodel.md"
    ]
]

makedocs(
    modules   = [NESSie, NESSie.BEM, NESSie.Format, NESSie.TestModel],
    clean     = true,
    doctest   = true,
    linkcheck = true,
    checkdocs = :all,
    pages     = pages,
    format    = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename  = "NESSie.jl",
    repo      = "https://github.com/tkemmer/NESSie.jl/blob/master{path}"
)

deploydocs(;
    repo = "github.com/tkemmer/NESSie.jl.git",
    devbranch = "master"
)
