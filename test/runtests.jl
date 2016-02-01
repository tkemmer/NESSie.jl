using ProteinES
using FactCheck

# All tests. Will be used if no command line arguments are given.
const tests = [
    "base/types.jl", "base/util.jl",
    "io/hmo", "io/matlab", "io/off", "io/pqr", "io/xml3d",
    "local/bem",
    "nonlocal/bem", "nonlocal/fem",
    "radon", "rjasanow"]
const testtypes = (Float64, Float32)

for t in (length(ARGS) > 0 ? ARGS : tests)
    facts(t) do
        include(endswith(t, ".jl") ? t : "$(t).jl")
    end
    println()
end

FactCheck.exitstatus()
