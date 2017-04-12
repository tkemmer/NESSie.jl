using ProteinES
using FactCheck

# All tests. Will be used if no command line arguments are given.
const tests = [
    "base/born", "base/constants", "base/model", "base/potentials", "base/quadrature", "base/util",
    "bem/local", "bem/nonlocal",
    "io/hmo", "io/mcsf", "io/off", "io/pqr", "io/skel", "io/vtk", "io/xml3d",
    "radon", "rjasanow"]
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
