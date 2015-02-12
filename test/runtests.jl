using NonLocalBEM
using Base.Test

# All tests. Will be used if no command line arguments are given.
const tests = ["hmo", "radon", "rjasanow", "util", "nonlocal"]

println("Running tests:")
for t in (length(ARGS) > 0 ? ARGS : tests)
    println(" * $(t)")
    include("$(t).jl")
end
