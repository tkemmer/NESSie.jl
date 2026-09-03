using NESSie.TestModel

@testset "XieSphere" begin
    using NESSie.TestModel: scalemodel

    for T in testtypes
        charges = [
            Charge(T[0, 0, 0, 1]...),
            Charge(T[1, 0, 0, 2]...),
            Charge(T[1, 1, 0, 3]...),
            Charge(T[0, 1, 0, 4]...)
        ]

        model = XieSphere(2one(T), charges)
        @test typeof(model) == XieSphere{T}
        @test model.radius == 2one(T)
        for (q1, q2) in zip(model.charges, scalemodel(charges, model.radius))
            @test q1.pos == q2.pos
            @test q1.val == q2.val
        end
        @test model.params == defaultopt(T)
    end
end

@testset "scalemodel" begin
    using NESSie.TestModel: scalemodel

    for T in testtypes
        charges = [
            Charge(T[0, 0, 0, 1]...),
            Charge(T[1, 0, 0, 2]...),
            Charge(T[1, 1, 0, 3]...),
            Charge(T[0, 1, 0, 4]...)
        ]

        # single charge only
        scaled = scalemodel([first(charges)], T(√(.5)/.8))
        @test typeof(scaled) == Vector{Charge{T}}
        @test scaled[1].pos ≈ T[0, 0, 0]
        @test scaled[1].val == one(T)
        scaled = scalemodel([Charge(T[1, 1, 1, 2]...)], T(√(.5)/.8))
        @test typeof(scaled) == Vector{Charge{T}}
        @test scaled[1].pos ≈ T[0, 0, 0]
        @test scaled[1].val == 2one(T)

        # translation only
        scaled = scalemodel(charges, T(√(.5)/.8))
        @test typeof(scaled) == Vector{Charge{T}}
        @test scaled[1].pos ≈ T[-.5, -.5, 0]
        @test scaled[1].val == one(T)
        @test scaled[2].pos ≈ T[.5, -.5, 0]
        @test scaled[2].val == 2one(T)
        @test scaled[3].pos ≈ T[.5, .5, 0]
        @test scaled[3].val == 3one(T)
        @test scaled[4].pos ≈ T[-.5, .5, 0]
        @test scaled[4].val == 4one(T)

        # translation and scaling
        scaled = scalemodel(charges, T(1.25))
        @test typeof(scaled) == Vector{Charge{T}}
        @test scaled[1].pos ≈ T[-√.5, -√.5, 0]
        @test scaled[1].val == one(T)
        @test scaled[2].pos ≈ T[√.5, -√.5, 0]
        @test scaled[2].val == 2one(T)
        @test scaled[3].pos ≈ T[√.5, √.5, 0]
        @test scaled[3].val == 3one(T)
        @test scaled[4].pos ≈ T[-√.5, √.5, 0]
        @test scaled[4].val == 4one(T)
    end
end

@testset "Model" begin
    for T in testtypes
        charges = [
            Charge(T[0, 0, 0, 1]...),
            Charge(T[1, 0, 0, 2]...),
            Charge(T[1, 1, 0, 3]...),
            Charge(T[0, 1, 0, 4]...)
        ]
        xie = XieSphere(2one(T), charges)

        let model = Model(xie)
            @test model isa Model{T, Triangle{T}}
            @test !isempty(model.nodes)
            @test !isempty(model.elements)
            @test model.charges == xie.charges
            @test model.params == defaultopt(T)
        end
    end
end

@testset "legendre" begin
    using NESSie.TestModel: legendre

    for T in testtypes
        # 0-length
        p = legendre(0, one(T))
        @test isa(p, Function)
        @test_throws BoundsError p(-1)
        @test_throws BoundsError p(0)
        @test_throws BoundsError p(1)

        # one element
        p = legendre(1, one(T))
        @test isa(p, Function)
        @test typeof(p(0)) == T
        @test p(0) == one(T)
        @test_throws BoundsError p(-1)
        @test_throws BoundsError p(1)

        # more elements
        p = legendre(5, one(T)/2)
        @test isa(p, Function)
        @test typeof(p(0)) == T
        @test p(0) == one(T)
        @test typeof(p(1)) == T
        @test p(1) == one(T)/2
        @test typeof(p(2)) == T
        @test p(2) ≈ T(-0.125)
        @test typeof(p(3)) == T
        @test p(3) ≈ T(-0.4375)
        @test typeof(p(4)) == T
        @test p(4) ≈ T(-0.2890625)
        @test_throws BoundsError p(-1)
        @test_throws BoundsError p(5)
    end
end

@testset "spherical_besseli" begin
    using NESSie.TestModel: spherical_besseli

    for T in testtypes
        # 0-length
        i = spherical_besseli(-2, one(T))
        @test isa(i, Function)
        @test_throws BoundsError i(-2)
        @test_throws BoundsError i(-1)
        @test_throws BoundsError i(0)

        # some elements
        i = spherical_besseli(0, one(T)/2)
        @test isa(i, Function)
        @test typeof(i(-1)) == T
        @test i(-1) ≈ T(√π * 1.2723896474148502)
        @test typeof(i(0)) == T
        @test i(0) ≈ T(√π * 0.5879930867904171)
        @test_throws BoundsError i(-2)
        @test_throws BoundsError i(1)
    end
end

@testset "spherical_besselk" begin
    using NESSie.TestModel: spherical_besselk

    for T in testtypes
        # 0-length
        k = spherical_besselk(-2, one(T))
        @test isa(k, Function)
        @test_throws BoundsError k(-2)
        @test_throws BoundsError k(-1)
        @test_throws BoundsError k(0)

        # some elements
        k = spherical_besselk(0, one(T)/2)
        @test isa(k, Function)
        @test typeof(k(-1)) == T
        @test k(-1) ≈ T(√π * 1.0750476034999203)
        @test typeof(k(0)) == T
        @test k(0) ≈ T(√π * 1.0750476034999203)
        @test_throws BoundsError k(-2)
        @test_throws BoundsError k(1)
    end
end

@testset "potentials and energies (interface variant identity)" begin
    for T in testtypes
        charges = [Charge(T[0, 0, 0, 1]...)]
        xie_sphere = XieSphere(T(2), charges)
        xie_local = LocalXieModel(xie_sphere, 10)
        xie_nl1   = NonlocalXieModel1(xie_sphere, 10)
        xie_nl2   = NonlocalXieModel2(xie_sphere, 10)

        r = T(2)
        Ξ = LinRange(T[0, 0, 0], T[0, 0, r + 0.5], 20)
        ξΩ   = T[0, 0, r - 0.2]
        ξΓ   = T[0, 0, r]
        ξΣ   = T[0, 0, r + 0.2]

        for (label, xie) in [("Local", xie_local), ("Nonlocal1", xie_nl1), ("Nonlocal2", xie_nl2)]
            # molpotential
            res = molpotential(Ξ, xie)
            @test res isa Vector{T}
            @test length(res) == length(Ξ)

            res = molpotential(reshape(Ξ, 1, length(Ξ)), xie)
            @test res isa Matrix{T}
            @test res == transpose(molpotential(Ξ, xie))

            res = molpotential((ξ for ξ in Ξ), xie)
            @test res isa Vector{T}
            @test res == molpotential(Ξ, xie)

            # rfpotential with explicit domain
            res = rfpotential(:Ω, ξΩ, xie)
            @test res isa T

            res = rfpotential(:Ω, [ξΩ], xie)
            @test res isa Vector{T}
            @test only(res) == rfpotential(:Ω, ξΩ, xie)

            res = rfpotential(:Ω, reshape([ξΩ], 1, 1), xie)
            @test res isa Matrix{T}
            @test only(res) == rfpotential(:Ω, ξΩ, xie)

            res = rfpotential(:Ω, (ξ for ξ in [ξΩ]), xie)
            @test res isa Vector{T}
            @test only(res) == rfpotential(:Ω, ξΩ, xie)

            res = rfpotential(:Γ, ξΓ, xie)
            @test res isa T
            @test res == rfpotential(:Ω, ξΓ, xie)

            res = rfpotential(:Γ, [ξΓ], xie)
            @test only(res) == rfpotential(:Ω, ξΓ, xie)

            res = rfpotential(:Σ, ξΣ, xie)
            @test res isa T

            res = rfpotential(:Σ, [ξΣ], xie)
            @test res isa Vector{T}
            @test only(res) == rfpotential(:Σ, ξΣ, xie)

            res = rfpotential(:Σ, reshape([ξΣ], 1, 1), xie)
            @test res isa Matrix{T}
            @test only(res) == rfpotential(:Σ, ξΣ, xie)

            res = rfpotential(:Σ, (ξ for ξ in [ξΣ]), xie)
            @test res isa Vector{T}
            @test only(res) == rfpotential(:Σ, ξΣ, xie)

            # rfpotential with auto-location
            res = rfpotential(ξΩ, xie)
            @test res isa T
            @test res == rfpotential(:Ω, ξΩ, xie)

            res = rfpotential(ξΓ, xie)
            @test res isa T
            @test res == rfpotential(:Γ, ξΓ, xie)

            res = rfpotential(ξΣ, xie)
            @test res isa T
            @test res == rfpotential(:Σ, ξΣ, xie)

            res = rfpotential(Ξ, xie)
            @test res isa Vector{T}
            @test length(res) == length(Ξ)

            res = rfpotential(reshape(Ξ, 1, length(Ξ)), xie)
            @test res isa Matrix{T}
            @test res == transpose(rfpotential(Ξ, xie))

            res = rfpotential((ξ for ξ in Ξ), xie)
            @test res isa Vector{T}
            @test res == rfpotential(Ξ, xie)

            # espotential with explicit domain
            res = espotential(:Ω, ξΩ, xie)
            @test res isa T

            res = espotential(:Ω, [ξΩ], xie)
            @test res isa Vector{T}
            @test only(res) == espotential(:Ω, ξΩ, xie)

            res = espotential(:Γ, ξΓ, xie)
            @test res isa T
            @test res == espotential(:Ω, ξΓ, xie)

            res = espotential(:Γ, [ξΓ], xie)
            @test only(res) == espotential(:Ω, ξΓ, xie)

            res = espotential(:Σ, ξΣ, xie)
            @test res isa T

            res = espotential(:Σ, [ξΣ], xie)
            @test res isa Vector{T}
            @test only(res) == espotential(:Σ, ξΣ, xie)

            res = espotential(:Σ, reshape([ξΣ], 1, 1), xie)
            @test res isa Matrix{T}
            @test only(res) == espotential(:Σ, ξΣ, xie)

            res = espotential(:Σ, (ξ for ξ in [ξΣ]), xie)
            @test res isa Vector{T}
            @test only(res) == espotential(:Σ, ξΣ, xie)

            # espotential with auto-location
            res = espotential(ξΩ, xie)
            @test res isa T
            @test res == espotential(:Ω, ξΩ, xie)

            res = espotential(ξΓ, xie)
            @test res isa T
            @test res == espotential(:Γ, ξΓ, xie)

            res = espotential(ξΣ, xie)
            @test res isa T
            @test res == espotential(:Σ, ξΣ, xie)

            res = espotential(Ξ, xie)
            @test res isa Vector{T}
            @test length(res) == length(Ξ)

            res = espotential(reshape(Ξ, 1, length(Ξ)), xie)
            @test res isa Matrix{T}
            @test res == transpose(espotential(Ξ, xie))

            res = espotential((ξ for ξ in Ξ), xie)
            @test res isa Vector{T}
            @test res == espotential(Ξ, xie)

            # espotential = rfpotential + molpotential identity
            res = espotential(Ξ, xie)
            @test res ≈ rfpotential(Ξ, xie) .+ molpotential(Ξ, xie)

            res = espotential(reshape(Ξ, 1, length(Ξ)), xie)
            @test res ≈ rfpotential(reshape(Ξ, 1, length(Ξ)), xie) .+ molpotential(reshape(Ξ, 1, length(Ξ)), xie)

            res = espotential((ξ for ξ in Ξ), xie)
            @test res ≈ rfpotential((ξ for ξ in Ξ), xie) .+ molpotential((ξ for ξ in Ξ), xie)

            # rfenergy
            res = rfenergy(xie)
            @test res isa T

            # error handling
            @test_throws ErrorException rfpotential(:Foo, ξΩ, xie)
            @test_throws ErrorException rfpotential(:Foo, Ξ, xie)
            @test_throws ErrorException espotential(:Foo, ξΩ, xie)
            @test_throws ErrorException espotential(:Foo, Ξ, xie)
        end
    end
end
