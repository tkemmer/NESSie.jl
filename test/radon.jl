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

# TODO laplacepot
# TODO laplacepot_dn
# TODO radoncoll!

