include("../src/nonlocal.jl")
using Base.Test
using NonLocalBEM

const testtypes = (Float64, Float32)

const tempfile = mktemp()
const testfiles = ((tempfile[1], tempfile[2], (0,0,0)),
				   ("./hmo.hmo", open("./test/hmo.hmo"), (2,2,2)))

try
	for (fname, fh, len) in testfiles
		try
			for dtype in testtypes
				nodes, elements, charges = readhmo(fh, dtype=dtype)
				@test isa(nodes, Vector{Vector{dtype}})
				@test isa(elements, Vector{Element{dtype}})
				@test isa(charges, Vector{Charge{dtype}})
				@test length(nodes) == len[1]
				@test length(elements) == len[2]
				@test length(charges) == len[3]
				seek(fh,0)
			end
		finally
			close(fh)
		end
	end
finally
	rm(testfiles[1][1])
end
