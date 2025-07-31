@testitem "Radon" begin
    include("testsetup.jl")

    using NESSie.Radon
    using NESSie.Radon: _regularyukawapot

    @testset "_regularyukawapot (SingleLayer)" begin
        for T in testtypes
            # x -> ξ
            x = ones(T, 3)
            yuk = T(7)
            ret = _regularyukawapot(SingleLayer, x, x, yuk)
            @test typeof(ret) == T
            @test ret == -7
            # ξ not in origin (no cancellation)
            ξ = -ones(T, 3)
            ret = _regularyukawapot(SingleLayer, x, ξ, yuk)
            @test typeof(ret) == T
            @test ret ≈ T(-0.28867513458631466433435706025773931468027671112178301167202)
            # ξ in origin (no cancellation)
            ξ = zeros(T, 3)
            ret = _regularyukawapot(SingleLayer, x, ξ, yuk)
            @test typeof(ret) == T
            @test ret ≈ T(-0.57734713663526745894988668695023080532035174827300561882957)
            # ξ in origin (potential cancellation)
            ret = _regularyukawapot(SingleLayer, T(.001) * x, ξ, yuk)
            @test typeof(ret) == T
            @test ret ≈ T(-6.95773573664079611429119468507639364064769481076591007615189)
            # ξ in origin (potential cancellation 2)
            ret = _regularyukawapot(SingleLayer, T(.0001) * x, ξ, yuk)
            @test typeof(ret) == T
            @test ret ≈ (T == Float64 ?
                -6.99575819000175052904190394524100762244862345273988759731283 :
                -yuk)
            # ξ in origin (potential cancellation 3)
            ret = _regularyukawapot(SingleLayer, T(.00001) * x, ξ, yuk)
            @test typeof(ret) == T
            @test ret ≈ (T == Float64 ?
                -6.99957566470162580591945945320355718542304458045960164651420 :
                -yuk)
        end
    end

    @testset "_regularyukawapot (DoubleLayer)" begin
        for T in testtypes
            # x -> ξ
            x = ones(T, 3)
            n = map(T, [1, 0, 0])
            yuk = T(7)
            ret = _regularyukawapot(DoubleLayer, x, x, yuk, n)
            @test typeof(ret) == T
            @test ret ≈ T(14.145081595145832)
            # ξ not in origin (no cancellation)
            ξ = -ones(T, 3)
            ret = _regularyukawapot(DoubleLayer, x, ξ, yuk, n)
            @test typeof(ret) == T
            @test ret ≈ T(0.048112522396707305228639137148607498838982054839859863103104)
            # ξ in origin (no cancellation)
            ξ = zeros(T, 3)
            ret = _regularyukawapot(DoubleLayer, x, ξ, yuk, n)
            @test typeof(ret) == T
            @test ret ≈ T(0.192436385477375021033020985138032408973537816345889754008150)
            # ξ in origin (potential cancellation)
            ret = _regularyukawapot(DoubleLayer, T(.001) * x, ξ, yuk, n)
            @test typeof(ret) == T
            @test ret ≈ T(14.03126641709760402264474861719654017225786809145550365963283)
            # ξ in origin (potential cancellation 2)
            ret = _regularyukawapot(DoubleLayer, T(.0001) * x, ξ, yuk, n)
            @test typeof(ret) == T
            @test ret ≈ (T == Float64 ?
                14.13365345844970855427401235115544109979343227826364222026794 :
                yuk^2 / 2 / Float32(√3))
            # ξ in origin (potential cancellation 3)
            ret = _regularyukawapot(DoubleLayer, T(.00001) * x, ξ, yuk, n)
            @test typeof(ret) == T
            @test ret ≈ (T == Float64 ?
                14.14393831379399210175378534648974379367907886588210581085651 :
                yuk^2 / 2 / Float32(√3))
        end
    end
end
