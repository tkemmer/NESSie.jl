using NonLocalES
using FactCheck

# All tests. Will be used if no command line arguments are given.
const tests = ["input/hmo", "input/matlab", "input/off", "input/pqr", "types", "util", "radon", "rjasanow", "bem"]
const testtypes = (Float64, Float32)

for t in (length(ARGS) > 0 ? ARGS : tests)
    facts(t) do
        include(endswith(t, ".jl") ? t : "$(t).jl")
    end
    println()
end

FactCheck.exitstatus()
