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
