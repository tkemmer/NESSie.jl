using ProteinES.Radon

context("laplacepot") do
    for T in testtypes
        # x --> ξ
        x = ones(T, 3)
        ret = Radon.laplacepot(x, x)
        @fact isa(ret, T) --> true
        @fact ret --> 0
        # ξ not in origin (moderate vecnorm)
        ξ = -ones(T, 3)
        ret = Radon.Radon.laplacepot(x, ξ)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(0.288675134594812882254574390250978727823800875635063438009301))
        # ξ in origin (moderate vecnorm)
        ξ = zeros(T, 3)
        ret = Radon.Radon.laplacepot(x, ξ)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(0.577350269189625764509148780501957455647601751270126876018602))
        @pending alternating_series --> :nothing
    end
end

context("∂ₙlaplacepot") do
    for T in testtypes
        # x --> ξ
        x = ones(T, 3)
        n = map(T, [1, 0, 0])
        ret = Radon.Radon.∂ₙlaplacepot(x, x, n)
        @fact isa(ret, T) --> true
        @fact ret --> 0
        # ξ not in origin (moderate vecnorm)
        ξ = -ones(T, 3)
        ret = Radon.Radon.∂ₙlaplacepot(x, ξ, n)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(-0.04811252243246881370909573170849645463730014593917723966821))
        # ξ in origin (moderate vecnorm)
        ξ = zeros(T, 3)
        ret = Radon.Radon.∂ₙlaplacepot(x, ξ, n)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(-0.19245008972987525483638292683398581854920058375670895867286))
        @pending alternating_series --> :nothing
    end
end

context("regularyukawapot") do
    for T in testtypes
        # x --> ξ
        x = ones(T, 3)
        opt = Option(zeros(T, 4)..., T(7))
        ret = Radon.regularyukawapot(x, x, opt)
        @fact isa(ret, T) --> true
        @fact ret --> -7
        # ξ not in origin (no cancellation)
        ξ = -ones(T, 3)
        ret = Radon.regularyukawapot(x, ξ, opt)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(-0.28867513458631466433435706025773931468027671112178301167202))
        # ξ in origin (no cancellation)
        ξ = zeros(T, 3)
        ret = Radon.regularyukawapot(x, ξ, opt)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(-0.57734713663526745894988668695023080532035174827300561882957))
        # ξ in origin (potential cancellation)
        ret = Radon.regularyukawapot(.001x, ξ, opt)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(-6.95773573664079611429119468507639364064769481076591007615189))
        # ξ in origin (potential cancellation 2)
        ret = Radon.regularyukawapot(.0001x, ξ, opt)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(-6.99575819000175052904190394524100762244862345273988759731283))
        # ξ in origin (potential cancellation 3)
        ret = Radon.regularyukawapot(.00001x, ξ, opt)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(-6.99957566470162580591945945320355718542304458045960164651420))
    end
end

context("∂ₙregularyukawapot") do
    for T in testtypes
        # x --> ξ
        x = ones(T, 3)
        n = map(T, [1, 0, 0])
        opt = Option(zeros(T, 4)..., T(7))
        ret = Radon.∂ₙregularyukawapot(x, x, n, opt)
        @fact isa(ret, T) --> true
        @fact ret --> 0
        # ξ not in origin (no cancellation)
        ξ = -ones(T, 3)
        ret = Radon.∂ₙregularyukawapot(x, ξ, n, opt)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(0.048112522396707305228639137148607498838982054839859863103104))
        # ξ in origin (no cancellation)
        ξ = zeros(T, 3)
        ret = Radon.∂ₙregularyukawapot(x, ξ, n, opt)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(0.192436385477375021033020985138032408973537816345889754008150))
        # ξ in origin (potential cancellation)
        ret = Radon.∂ₙregularyukawapot(.001x, ξ, n, opt)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(14.03126641709760402264474861719654017225786809145550365963283))
        # ξ in origin (potential cancellation 2)
        ret = Radon.∂ₙregularyukawapot(.0001x, ξ, n, opt)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(14.13365345844970855427401235115544109979343227826364222026794))
        # ξ in origin (potential cancellation 3)
        ret = Radon.∂ₙregularyukawapot(.00001x, ξ, n, opt)
        @fact isa(ret, T) --> true
        @fact ret --> roughly(T(14.14393831379399210175378534648974379367907886588210581085651))
    end
end

@pending radoncoll! --> :nothing
