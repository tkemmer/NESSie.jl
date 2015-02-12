using NonLocalBEM
using Base.Test

const tests = ["hmo", "radon", "rjasanow", "util", "nonlocal"]

println("Running tests:")
for t in tests
	println(" * $(t)")
	include("$(t).jl")
end
