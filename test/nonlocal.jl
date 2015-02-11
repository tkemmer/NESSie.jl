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

# check props! and isdegenerate
@test_throws ErrorException isdegenerate([1., 1.], [1., 1., 1.], [1., 1., 1.])
@test_throws ErrorException isdegenerate([1., 1., 1.], [1., 1.], [1., 1., 1.])
@test_throws ErrorException isdegenerate([1., 1., 1.], [1., 1., 1.], [1., 1.])
for dtype in (Float64, Float32)
	# degenerate triangles
	elem = Element(map(dtype, [0, 0, 0]), map(dtype, [0, 0, 0]), map(dtype, [1, 0, 0]))
	@test_throws ErrorException props!(elem)
	elem = Element(map(dtype, [0, 0, 0]), map(dtype, [0, 0, 1]), map(dtype, [0, 0, 0]))
	@test_throws ErrorException props!(elem)
	elem = Element(map(dtype, [1, 0, 0]), map(dtype, [0, 0, 0]), map(dtype, [0, 0, 0]))
	@test_throws ErrorException props!(elem)
	elem = Element(map(dtype, [0, 1, 0]), map(dtype, [0, 2, 0]), map(dtype, [0, 3, 0]))
	@test_throws ErrorException props!(elem)
	# simple 2D triangle
	elem = Element(map(dtype, [0, 0, 0]), map(dtype, [0, 0, 3]), map(dtype, [0, 3, 0]))
	props!(elem)
	@test_approx_eq elem.center [0, 1, 1]
	@test_approx_eq elem.normal [-1, 0, 0]
	@test_approx_eq elem.distorig 0.
	@test_approx_eq elem.area 4.5
	# simple 3D triangle
	elem = Element(map(dtype, [3, 0, 0]), map(dtype, [0, 4, 0]), map(dtype, [0, 0, 5]))
	props!(elem)
	@test_approx_eq elem.center [1., 4/3, 5/3]
	@test_approx_eq elem.normal (√769 \ [20, 15, 12])
	@test_approx_eq elem.distorig -60/√769
	@test_approx_eq elem.area √769/2
end

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

# check regularyukawapot
for dtype in (Float64, Float32)
	# x == ξ
	x = ones(dtype, 3)
	opt = Option(zeros(dtype, 4)..., convert(dtype, 7))
	ret = regularyukawapot(x, x, opt)
	@test isa(ret, dtype)
	@test_approx_eq ret -7
	# ξ not in origin (no cancellation)
	ξ = -ones(dtype, 3)
	ret = regularyukawapot(x, ξ, opt)
	@test isa(ret, dtype)
	@test_approx_eq ret -0.28867513458631466433435706025773931468027671112178301167202
	# ξ in origin (no cancellation)
	ξ = zeros(dtype, 3)
	ret = regularyukawapot(x, ξ, opt)
	@test isa(ret, dtype)
	@test_approx_eq ret -0.57734713663526745894988668695023080532035174827300561882957
	# ξ in origin (potential cancellation)
	ret = regularyukawapot(.001x, ξ, opt)
	@test isa(ret, dtype)
	@test_approx_eq ret -6.95773573664079611429119468507639364064769481076591007615189
	# ξ in origin (potential cancellation 2)
	ret = regularyukawapot(.0001x, ξ, opt)
	@test isa(ret, dtype)
	@test_approx_eq ret -6.99575819000175052904190394524100762244862345273988759731283
	# ξ in origin (potential cancellation 3)
	ret = regularyukawapot(.00001x, ξ, opt)
	@test isa(ret, dtype)
	@test_approx_eq ret -6.99957566470162580591945945320355718542304458045960164651420
end

# check regularyukawapot_dn
for dtype in (Float64, Float32)
	# x == ξ
	x = ones(dtype, 3)
	opt = Option(zeros(dtype, 4)..., convert(dtype, 7))
	ret = regularyukawapot_dn(x, x, x, opt)
	@test isa(ret, dtype)
	@test ret == 0
	# ξ not in origin (no cancellation)
	ξ = -ones(dtype, 3)
	n = map(dtype, [1, 0, 0])
	ret = regularyukawapot_dn(x, ξ, n, opt)
	@test isa(ret, dtype)
	@test_approx_eq ret 0.048112522396707305228639137148607498838982054839859863103104
	# ξ in origin (no cancellation)
	ξ = zeros(dtype, 3)
	ret = regularyukawapot_dn(x, ξ, n, opt)
	@test isa(ret, dtype)
	@test_approx_eq ret 0.192436385477375021033020985138032408973537816345889754008150
	# ξ in origin (potential cancellation)
	ret = regularyukawapot_dn(.001x, ξ, n, opt)
	@test isa(ret, dtype)
	@test_approx_eq ret 14.03126641709760402264474861719654017225786809145550365963283
	# ξ in origin (potential cancellation 2)
	ret = regularyukawapot_dn(.0001x, ξ, n, opt)
	@test isa(ret, dtype)
	@test_approx_eq ret 14.13365345844970855427401235115544109979343227826364222026794
	# ξ in origin (potential cancellation 3)
	ret = regularyukawapot_dn(.00001x, ξ, n, opt)
	@test isa(ret, dtype)
	@test_approx_eq ret 14.14393831379399210175378534648974379367907886588210581085651
end

# TODO regularyukawacoll

# TODO cauchy
