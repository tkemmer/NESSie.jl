push!(LOAD_PATH,"../src/")

using ProteinES
using FactCheck

# All tests. Will be used if no command line arguments are given.
const tests = [
    "base/constants.jl",
    "base/model.jl",
    "base/potentials.jl",
    "base/quadrature.jl",
    "base/util.jl",
    "born.jl",
    "bem/local.jl",
    "bem/nonlocal.jl",
    "io/hmo.jl",
    "io/mcsf.jl",
    "io/msms.jl",
    "io/off.jl",
    "io/pqr.jl",
    "io/skel.jl",
    "io/vtk.jl",
    "io/xml3d.jl",
    "radon.jl",
    "rjasanow.jl"
]
const testtypes = (Float64, Float32)

# Check command line arguments
ok = true
for arg in ARGS
    arg in tests || "$(arg).jl" in tests || begin
        println("\n\e[1;31mERROR: Invalid test file $arg\e[0m\n")
        ok = false
    end
end
ok || exit(1)

# Run tests
include("testutils.jl")
for t in (length(ARGS) > 0 ? ARGS : tests)
    facts(t) do
        include(endswith(t, ".jl") ? t : "$(t).jl")
    end
    println()
end
FactCheck.exitstatus()
