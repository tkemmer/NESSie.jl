using NonLocalBEM
using Base.Test

const tests = ["hmo", "nonlocal"]

println("Running tests:")
for t in tests
	println(" * $(t)")
	include("$(t).jl")
end
