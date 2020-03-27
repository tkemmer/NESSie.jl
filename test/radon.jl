using NESSie: yukawa
using NESSie.Radon
using NESSie.Radon: radoncoll!

@testset "regularyukawapot" begin
    for T in testtypes
        # x -> ξ
        x = ones(T, 3)
        yuk = T(7)
        ret = Radon.regularyukawapot(x, x, yuk)
        @test typeof(ret) == T
        @test ret == -7
        # ξ not in origin (no cancellation)
        ξ = -ones(T, 3)
        ret = Radon.regularyukawapot(x, ξ, yuk)
        @test typeof(ret) == T
        @test ret ≈ T(-0.28867513458631466433435706025773931468027671112178301167202)
        # ξ in origin (no cancellation)
        ξ = zeros(T, 3)
        ret = Radon.regularyukawapot(x, ξ, yuk)
        @test typeof(ret) == T
        @test ret ≈ T(-0.57734713663526745894988668695023080532035174827300561882957)
        # ξ in origin (potential cancellation)
        ret = Radon.regularyukawapot(T(.001) * x, ξ, yuk)
        @test typeof(ret) == T
        @test ret ≈ T(-6.95773573664079611429119468507639364064769481076591007615189)
        # ξ in origin (potential cancellation 2)
        ret = Radon.regularyukawapot(T(.0001) * x, ξ, yuk)
        @test typeof(ret) == T
        @test ret ≈ T(-6.99575819000175052904190394524100762244862345273988759731283)
        # ξ in origin (potential cancellation 3)
        ret = Radon.regularyukawapot(T(.00001) * x, ξ, yuk)
        @test typeof(ret) == T
        @test ret ≈ T(-6.99957566470162580591945945320355718542304458045960164651420)
    end
end

@testset "∂ₙregularyukawapot" begin
    for T in testtypes
        # x -> ξ
        x = ones(T, 3)
        n = map(T, [1, 0, 0])
        yuk = T(7)
        ret = Radon.∂ₙregularyukawapot(x, x, yuk, n)
        @test typeof(ret) == T
        @test ret ≈ T(14.145081595145832)
        # ξ not in origin (no cancellation)
        ξ = -ones(T, 3)
        ret = Radon.∂ₙregularyukawapot(x, ξ, yuk, n)
        @test typeof(ret) == T
        @test ret ≈ T(0.048112522396707305228639137148607498838982054839859863103104)
        # ξ in origin (no cancellation)
        ξ = zeros(T, 3)
        ret = Radon.∂ₙregularyukawapot(x, ξ, yuk, n)
        @test typeof(ret) == T
        @test ret ≈ T(0.192436385477375021033020985138032408973537816345889754008150)
        # ξ in origin (potential cancellation)
        ret = Radon.∂ₙregularyukawapot(T(.001) * x, ξ, yuk, n)
        @test typeof(ret) == T
        @test ret ≈ T(14.03126641709760402264474861719654017225786809145550365963283)
        # ξ in origin (potential cancellation 2)
        ret = Radon.∂ₙregularyukawapot(T(.0001) * x, ξ, yuk, n)
        @test typeof(ret) == T
        @test ret ≈ T(14.13365345844970855427401235115544109979343227826364222026794)
        # ξ in origin (potential cancellation 3)
        ret = Radon.∂ₙregularyukawapot(T(.00001) * x, ξ, yuk, n)
        @test typeof(ret) == T
        @test ret ≈ T(14.14393831379399210175378534648974379367907886588210581085651)
    end
end

@testset "radoncoll!" begin
    for T in testtypes
        elem = [Triangle(zeros(T, 3), T[2, 0, 0], T[0, 2, 0])]
        # compute triangle area ∫dx
        res = zeros(T, 1, 1)
        radoncoll!(res, elem, [zeros(T, 3)], (_, __, ___, ____) -> one(T), zero(T))
        @test res[1] ≈ 2one(T)
        # compute prism volume
        res = zeros(T, 2, 1)
        radoncoll!(res, elem, [T[0, 0, 4], T[0, 0, 6]], (x, ξ, _, __) -> ξ[3] - x[3], zero(T))
        @test res[1] ≈ 8one(T)
        @test res[2] ≈ 12one(T)
        # compute tetrahedron volume
        res = zeros(T, 1, 1)
        radoncoll!(res, elem, [T[0, 0, 2]], (x, ξ, _, __) -> ξ[3] - x[1] - x[2], zero(T))
        @test res[1] ≈ 4one(T)/3
        # a slightly more elaborate tetrahedron volume
        elem = [Triangle(zeros(T, 3), T[3, 0, 0], T[0, 2, 0])]
        res = zeros(T, 1, 1)
        radoncoll!(res, elem, [T[0, 0, 6]], (x, ξ, _, __) -> ξ[3] - 2x[1] - 3x[2], zero(T))
        @test res[1] ≈ 6one(T)
        # ∫∫6x²-40y dA
        elem = [Triangle(T[0, 3, 0], T[1, 1, 0], T[5, 3, 0])]
        res = zeros(T, 1, 1)
        radoncoll!(res, elem, [zeros(T, 3)], (x, _, __, ___) -> 6*x[1]^2-40x[2], zero(T))
        @test res[1] ≈ T(-935/3)
    end
end
