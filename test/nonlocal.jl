# check singularpot
for dtype in (Float64, Float32)
	# empty lists
	@test singularpot(Element{dtype}[], Charge{dtype}[]) == (dtype[], dtype[])
	# empty charge list
	elements = [Element(map(dtype, [0, 0, 0]), map(dtype, [0, 0, 3]), map(dtype, [0, 3, 0])),
				Element(map(dtype, [3, 0, 0]), map(dtype, [0, 4, 0]), map(dtype, [0, 0, 5]))]
	map(props!, elements)
	(umol, qmol) = singularpot(elements, Charge{dtype}[])
	@test isa(umol, Vector{dtype}) && isa(qmol, Vector{dtype})
	@test umol == qmol == [0, 0]
	# single charge (trivial options)
	opt = Option(one(dtype), zeros(dtype, 3)...)
	charges = [Charge(dtype, 0, 0, 0, √2)]
	(umol, qmol) = singularpot(elements, charges, opt)
	@test isa(umol, Vector{dtype}) && isa(qmol, Vector{dtype})
	@test_approx_eq umol [1, .6]
	@test_approx_eq qmol [0, -162/25/√769]
	# single charge (non-trivial options)
	opt = Option(convert(dtype, 2), zeros(dtype, 3)...)
	(umol, qmol) = singularpot(elements, charges, opt)
	@test isa(umol, Vector{dtype}) && isa(qmol, Vector{dtype})
	@test_approx_eq umol [.5, .3]
	@test_approx_eq qmol [0, -162/50/√769]
	# multiple charges (non-trivial options)
	push!(charges, Charge(dtype, 1, 1, 1, -√5))
	(umol, qmol) = singularpot(elements, charges, opt)
	@test isa(umol, Vector{dtype}) && isa(qmol, Vector{dtype})
	@test_approx_eq umol [2 \ (1-√5), -1.2]
	@test_approx_eq qmol [2 \ √5, 1593/50/√769]
end

# TODO cauchy
