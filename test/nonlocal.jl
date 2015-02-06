# check eye!
let m = zeros(Int, 3, 3)
	eye!(m);    @test m == [1 0 0; 0 1 0; 0 0 1]
	eye!(m, 2); @test m == [2 0 0; 0 2 0; 0 0 2]
	m = zeros(Int, 2, 3)
	eye!(m);    @test m == [1 0 0; 0 1 0]
	eye!(m, 2); @test m == [2 0 0; 0 2 0]
	m = zeros(Int, 3, 2)
	eye!(m);    @test m == [1 0; 0 1; 0 0]
	eye!(m, 2); @test m == [2 0; 0 2; 0 0]
end

# check compute_props! and isdegenerate
@test_throws ErrorException isdegenerate([1., 1.], [1., 1., 1.], [1., 1., 1.])
@test_throws ErrorException isdegenerate([1., 1., 1.], [1., 1.], [1., 1., 1.])
@test_throws ErrorException isdegenerate([1., 1., 1.], [1., 1., 1.], [1., 1.])
for dtype in (Float64, Float32)
	# degenerate triangles
	elem = Element(map(dtype, [0, 0, 0]), map(dtype, [0, 0, 0]), map(dtype, [1, 0, 0]))
	@test_throws ErrorException compute_props!(elem)
	elem = Element(map(dtype, [0, 0, 0]), map(dtype, [0, 0, 1]), map(dtype, [0, 0, 0]))
	@test_throws ErrorException compute_props!(elem)
	elem = Element(map(dtype, [1, 0, 0]), map(dtype, [0, 0, 0]), map(dtype, [0, 0, 0]))
	@test_throws ErrorException compute_props!(elem)
	elem = Element(map(dtype, [0, 1, 0]), map(dtype, [0, 2, 0]), map(dtype, [0, 3, 0]))
	@test_throws ErrorException compute_props!(elem)
	# simple 2D triangle
	elem = Element(map(dtype, [0, 0, 0]), map(dtype, [0, 0, 3]), map(dtype, [0, 3, 0]))
	compute_props!(elem)
	@test_approx_eq elem.center [0, 1, 1]
	@test_approx_eq elem.normal [-1, 0, 0]
	@test_approx_eq elem.distorig 0.
	@test_approx_eq elem.area 4.5
	# simple 3D triangle
	elem = Element(map(dtype, [3, 0, 0]), map(dtype, [0, 4, 0]), map(dtype, [0, 0, 5]))
	compute_props!(elem)
	@test_approx_eq elem.center [1., 4/3, 5/3]
	@test_approx_eq elem.normal (√769 \ [20, 15, 12])
	@test_approx_eq elem.distorig -60/√769
	@test_approx_eq elem.area √769/2
end

# check compute_singularpot
for dtype in (Float64, Float32)
	# empty lists
	@test compute_singularpot(Element{dtype}[], Charge{dtype}[]) == (dtype[], dtype[])
	# empty charge list
	elements = [Element(map(dtype, [0, 0, 0]), map(dtype, [0, 0, 3]), map(dtype, [0, 3, 0])),
				Element(map(dtype, [3, 0, 0]), map(dtype, [0, 4, 0]), map(dtype, [0, 0, 5]))]
	map(compute_props!, elements)
	(umol, qmol) = compute_singularpot(elements, Charge{dtype}[])
	@test isa(umol, Vector{dtype}) && isa(qmol, Vector{dtype})
	@test umol == qmol == [0, 0]
	# single charge (trivial options)
	opt = Option(one(dtype), zero(dtype), zero(dtype), zero(dtype))
	charges = [Charge(dtype, 0, 0, 0, √2)]
	(umol, qmol) = compute_singularpot(elements, charges, opt)
	@test isa(umol, Vector{dtype}) && isa(qmol, Vector{dtype})
	@test_approx_eq umol [1, .6]
	@test_approx_eq qmol [0, -162/25/√769]
	# single charge (non-trivial options)
	opt = Option(convert(dtype, 2), zero(dtype), zero(dtype), zero(dtype))
	(umol, qmol) = compute_singularpot(elements, charges, opt)
	@test isa(umol, Vector{dtype}) && isa(qmol, Vector{dtype})
	@test_approx_eq umol [.5, .3]
	@test_approx_eq qmol [0, -162/50/√769]
	# multiple charges (non-trivial options)
	push!(charges, Charge(dtype, 1, 1, 1, -√5))
	(umol, qmol) = compute_singularpot(elements, charges, opt)
	@test isa(umol, Vector{dtype}) && isa(qmol, Vector{dtype})
	@test_approx_eq umol [2 \ (1-√5), -1.2]
	@test_approx_eq qmol [2 \ √5, 1593/50/√769]
end

