using NESSie.Rjasanow
using NESSie.Rjasanow: laplacepot

@test_skip laplacecoll!

@testset "laplacepot" begin
    for T in testtypes
        # 1. Single layer potentials
        # 1.1 ξ is a triangle vertex
        elem = Triangle(T[0, 0, 0], T[0, 1, 0], T[0, 0, 1])
        ξ = elem.v1     # h = √2, φ1 = -π/4, φ2=π/4
        res = laplacepot(SingleLayer, ξ, elem, zero(T))
        @test typeof(res) == T
        @test res ≈ T(4π * 0.099189377627951192068415700486459476011275034113019436460)
        ξ = elem.v2     # h = 1, φ1 = -π/4, φ2 = 0
        res = laplacepot(SingleLayer, ξ, elem, zero(T))
        @test typeof(res) == T
        @test res ≈ T(4π * 0.070137481542397515891964948616904961601109053078512490737)
        ξ = elem.v3     # h = 1, φ1 = 0, φ2 = π/4
        res = laplacepot(SingleLayer, ξ, elem, zero(T))
        @test typeof(res) == T
        @test res ≈ T(4π * 0.070137481542397515891964948616904961601109053078512490737)

        # 1.2 ξ lies in the same plane as the triangle
        ξ = T[0, -1, 0]
        res = laplacepot(SingleLayer, ξ, elem, zero(T))
        @test typeof(res) == T
        @test res ≈ T(4π * 0.029051896085553676176450751869554514410165981034506945723)
        ξ = T[0, -1, -1]
        res = laplacepot(SingleLayer, ξ, elem, zero(T))
        @test typeof(res) == T
        @test res ≈ T(4π * (0.110553267370012295141653085597608171948437903548411687443
                                - 2 * 0.044743379406364904901413439333009046862112643066826198636))

        # 1.3 ξ lies on a line through a vertex, perpendicular
        #     to the triangle plane. In the following code, ξ
        #     is already proceted onto the plane!
        ξ = elem.v1
        res = laplacepot(SingleLayer, ξ, elem, one(T))
        @test typeof(res) == T
        @test res ≈ T(4π * 0.034775298546314300210123433855429220291493953452986787115)
        ξ = elem.v2
        res = laplacepot(SingleLayer, ξ, elem, one(T))
        @test typeof(res) == T
        @test res ≈ T(4π * 0.03156675645153950833618064292205638527751593336549)
        ξ = elem.v3
        res = laplacepot(SingleLayer, ξ, elem, one(T))
        @test typeof(res) == T
        @test res ≈ T(4π * 0.03156675645153950833618064292205638527751593336549)

        # 1.4 ξ lies somewhere else in space. Again, ξ is
        #     already projected onto the plane!
        ξ = T[0, -1, 0]
        res = laplacepot(SingleLayer, ξ, elem, one(T))
        @test typeof(res) == T
        @test res ≈ T(4π * (0.05486307812003175845441110089494128884416545308420259271
                                     -0.03156675645153950833618064292205638527751593336549))
        ξ = T[0, -1, -1]
        res = laplacepot(SingleLayer, ξ, elem, one(T))
        @test typeof(res) == T
        @test res ≈ T(4π * (2 * 0.035315594609376722705763636889615906447125307241104235442
                                    - 2 * 0.0260303288327139057352736693314968773933071931295084))

        # 2. Double layer potentials
        # 1.1 ξ is a triangle vertex
        ξ = elem.v1
        res = laplacepot(DoubleLayer, ξ, elem, zero(T))
        @test typeof(res) == T
        @test res ≈ zero(T)
        ξ = elem.v2
        res = laplacepot(DoubleLayer, ξ, elem, zero(T))
        @test typeof(res) == T
        @test res ≈ zero(T)
        ξ = elem.v3
        res = laplacepot(DoubleLayer, ξ, elem, zero(T))
        @test typeof(res) == T
        @test res ≈ zero(T)

        # 1.2 ξ lies in the same plane as the triangle
        ξ = T[0, -1, 0]
        res = laplacepot(DoubleLayer, ξ, elem, zero(T))
        @test typeof(res) == T
        @test res ≈ zero(T)
        ξ = T[0, -1, -1]
        res = laplacepot(DoubleLayer, ξ, elem, zero(T))
        @test typeof(res) == T
        @test res ≈ zero(T)

        @test_skip dl_ξ_above_vertex
        @test_skip dl_ξ_somewhere_in_space
    end
end
