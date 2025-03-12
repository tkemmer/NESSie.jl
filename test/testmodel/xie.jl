@testitem "Test model: Xie" begin
    include("../testsetup.jl")

    using NESSie.TestModel

    @testset "XieModel" begin
        using NESSie.TestModel: scalemodel

        for T in testtypes
            charges = [
                Charge(T[0, 0, 0, 1]...),
                Charge(T[1, 0, 0, 2]...),
                Charge(T[1, 1, 0, 3]...),
                Charge(T[0, 1, 0, 4]...)
            ]

            model = XieModel(2one(T), charges)
            @test typeof(model) == XieModel{T}
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
            xie = XieModel(2one(T), charges)

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

    @test_skip NonlocalXieModel1
    @test_skip coefficients
    @test_skip φΣ
    @test_skip φΩ
    @test_skip espotential
end
