using ProteinES
using FactCheck

# All tests. Will be used if no command line arguments are given.
const tests = [
    "base/born.jl", "base/constants.jl", "base/model.jl", "base/potentials.jl", "base/quadrature.jl", "base/util.jl",
    "bem/local.jl", "bem/nonlocal.jl",
    "io/hmo", "io/mcsf", "io/off", "io/pqr", "io/skel", "io/xml3d",
    "radon", "rjasanow"]
const testtypes = (Float64, Float32)

for t in (length(ARGS) > 0 ? ARGS : tests)
    facts(t) do
        include(endswith(t, ".jl") ? t : "$(t).jl")
    end
    println()
end

FactCheck.exitstatus()
