@testitem "Nonlocal BEM" begin
    include("../testsetup.jl")

    using NESSie.BEM, NESSie.TestModel

    @testset "Comparison to Born test model" begin
        for T in testtypes, method in [:blas, :gmres]
            ion   = bornion("Ca", T)
            model = Model(ion)
            bem   = solve(NonlocalES, model; method = method)
            Ξ     = LinRange(T[0, 0, 0], T[0, 0, 10], 100)
            ξΩ    = T[0, 0, ion.radius - 0.2]
            ξΓ    = T[0, 0, ion.radius]
            ξΣ    = T[0, 0, ion.radius + 0.2]

            res = molpotential(Ξ, bem)
            @test res isa Vector{T}
            @test res == molpotential(Ξ, ion)

            res = rfpotential(ξΩ, bem)
            @test res isa T
            @test res ≈ rfpotential(:Ω, NonlocalES, ξΩ, ion) rtol=0.05
            @test res == rfpotential(:Ω, ξΩ, bem)
            @test res == only(rfpotential(:Ω, [ξΩ], bem))

            res = rfpotential(ξΓ, bem)
            @test res isa T
            @test res ≈ rfpotential(:Γ, NonlocalES, ξΓ, ion) rtol=0.05
            @test res == rfpotential(:Γ, ξΓ, bem)
            @test res == only(rfpotential(:Γ, [ξΓ], bem))

            res = rfpotential(ξΣ, bem)
            @test res isa T
            @test res ≈ rfpotential(:Σ, NonlocalES, ξΣ, ion) rtol=0.05
            @test res == rfpotential(:Σ, ξΣ, bem)
            @test res == only(rfpotential(:Σ, [ξΣ], bem))

            res = rfpotential(Ξ, bem)
            @test res isa Vector{T}
            @test res ≈ rfpotential(NonlocalES, Ξ, ion) rtol=0.05

            res = rfpotential((ξ for ξ in Ξ), bem)
            @test res isa Vector{T}
            @test res ≈ rfpotential(NonlocalES, Ξ, ion) rtol=0.05

            res = espotential(ξΩ, bem)
            @test res isa T
            @test res ≈ espotential(:Ω, NonlocalES, ξΩ, ion) rtol=0.05
            @test res == espotential(:Ω, ξΩ, bem)
            @test res == only(espotential(:Ω, [ξΩ], bem))

            res = espotential(ξΓ, bem)
            @test res isa T
            @test res ≈ espotential(:Γ, NonlocalES, ξΓ, ion) rtol=0.05
            @test res == espotential(:Γ, ξΓ, bem)
            @test res == only(espotential(:Γ, [ξΓ], bem))

            res = espotential(ξΣ, bem)
            @test res isa T
            @test res ≈ espotential(:Σ, NonlocalES, ξΣ, ion) rtol=0.05
            @test res == espotential(:Σ, ξΣ, bem)
            @test res == only(espotential(:Σ, [ξΣ], bem))

            res = espotential(Ξ, bem)
            @test res isa Vector{T}
            @test res ≈ espotential(NonlocalES, Ξ, ion) rtol=0.05
            @test res ≈ rfpotential(Ξ, bem) .+ molpotential(Ξ, bem)

            res = espotential((ξ for ξ in Ξ), bem)
            @test res isa Vector{T}
            @test res ≈ espotential(NonlocalES, Ξ, ion) rtol=0.05
            @test res ≈ rfpotential(Ξ, bem) .+ molpotential(Ξ, bem)

            res = rfenergy(bem)
            @test res isa T
            @test res ≈ rfenergy(NonlocalES, ion) rtol=0.05

            @test_throws ErrorException rfpotential(:Foo, ξΩ, bem)
            @test_throws ErrorException rfpotential(:Foo, Ξ, bem)
            @test_throws ErrorException espotential(:Foo, ξΩ, bem)
            @test_throws ErrorException espotential(:Foo, Ξ, bem)
        end
    end
end
