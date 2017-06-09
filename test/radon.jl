using NESSie: yukawa
using NESSie.Radon

context("regularyukawapot") do
    for T in testtypes
        # x --> ξ
        x = ones(T, 3)
        yuk = T(7)
        ret = Radon.regularyukawapot(x, x, yuk)
        @fact typeof(ret) --> T
        @fact ret --> -7
        # ξ not in origin (no cancellation)
        ξ = -ones(T, 3)
        ret = Radon.regularyukawapot(x, ξ, yuk)
        @fact typeof(ret) --> T
        @fact ret --> roughly(T(-0.28867513458631466433435706025773931468027671112178301167202))
        # ξ in origin (no cancellation)
        ξ = zeros(T, 3)
        ret = Radon.regularyukawapot(x, ξ, yuk)
        @fact typeof(ret) --> T
        @fact ret --> roughly(T(-0.57734713663526745894988668695023080532035174827300561882957))
        # ξ in origin (potential cancellation)
        ret = Radon.regularyukawapot(T(.001) * x, ξ, yuk)
        @fact typeof(ret) --> T
        @fact ret --> roughly(T(-6.95773573664079611429119468507639364064769481076591007615189))
        # ξ in origin (potential cancellation 2)
        ret = Radon.regularyukawapot(T(.0001) * x, ξ, yuk)
        @fact typeof(ret) --> T
        @fact ret --> roughly(T(-6.99575819000175052904190394524100762244862345273988759731283))
        # ξ in origin (potential cancellation 3)
        ret = Radon.regularyukawapot(T(.00001) * x, ξ, yuk)
        @fact typeof(ret) --> T
        @fact ret --> roughly(T(-6.99957566470162580591945945320355718542304458045960164651420))
    end
end

context("∂ₙregularyukawapot") do
    for T in testtypes
        # x --> ξ
        x = ones(T, 3)
        n = map(T, [1, 0, 0])
        yuk = T(7)
        ret = Radon.∂ₙregularyukawapot(x, x, n, yuk)
        @fact typeof(ret) --> T
        @fact ret --> 0
        # ξ not in origin (no cancellation)
        ξ = -ones(T, 3)
        ret = Radon.∂ₙregularyukawapot(x, ξ, n, yuk)
        @fact typeof(ret) --> T
        @fact ret --> roughly(T(0.048112522396707305228639137148607498838982054839859863103104))
        # ξ in origin (no cancellation)
        ξ = zeros(T, 3)
        ret = Radon.∂ₙregularyukawapot(x, ξ, n, yuk)
        @fact typeof(ret) --> T
        @fact ret --> roughly(T(0.192436385477375021033020985138032408973537816345889754008150))
        # ξ in origin (potential cancellation)
        ret = Radon.∂ₙregularyukawapot(T(.001) * x, ξ, n, yuk)
        @fact typeof(ret) --> T
        @fact ret --> roughly(T(14.03126641709760402264474861719654017225786809145550365963283))
        # ξ in origin (potential cancellation 2)
        ret = Radon.∂ₙregularyukawapot(T(.0001) * x, ξ, n, yuk)
        @fact typeof(ret) --> T
        @fact ret --> roughly(T(14.13365345844970855427401235115544109979343227826364222026794))
        # ξ in origin (potential cancellation 3)
        ret = Radon.∂ₙregularyukawapot(T(.00001) * x, ξ, n, yuk)
        @fact typeof(ret) --> T
        @fact ret --> roughly(T(14.14393831379399210175378534648974379367907886588210581085651))
    end
end

@pending radoncoll! --> :nothing
