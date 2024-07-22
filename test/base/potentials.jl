@testitem "Potentials" begin
    include("../testsetup.jl")

    @testset "φmol and ∂ₙφmol" begin
        using NESSie: φmol, ∂ₙφmol

        # TODO check list variants of the functions
        for T in testtypes
            # empty lists
            model = Model{T, Triangle{T}}()
            @test φmol(model) == T[]
            # empty charge list
            elements = [Triangle(T[0, 0, 0], T[0, 0, 3], T[0, 3, 0]),
                        Triangle(T[3, 0, 0], T[0, 4, 0], T[0, 0, 5])]
            model = Model(Vector{T}[], elements)
            umol = φmol(model)
            @test typeof(umol) == Vector{T}
            @test umol == zeros(T, 2)
            qmol = ∂ₙφmol(model)
            @test typeof(qmol) == Vector{T}
            @test qmol == zeros(T, 2)
            # single charge
            charges = [Charge(T[0, 0, 0, √2]...)]
            model = Model(Vector{T}[], elements, charges)
            umol = φmol(model)
            @test typeof(umol) == Vector{T}
            @test umol ≈ T[1, .6]
            qmol = ∂ₙφmol(model)
            @test typeof(qmol) == Vector{T}
            @test qmol ≈ T[0, -162/25/√769]
            # multiple charges
            push!(model.charges, Charge(T[1, 1, 1, -√5]...))
            umol = φmol(model)
            @test typeof(umol) == Vector{T}
            @test umol ≈ T[1 - √5, -2.4]
            qmol = ∂ₙφmol(model)
            @test typeof(qmol) == Vector{T}
            @test qmol ≈ T[√5, 1593/25/√769]
            # singular case
            umol = φmol(zeros(T, 3), [Charge(T[0, 0, 0, √2]...)])
            @test typeof(umol) == T
            @test umol ≈ T(√2 / 1e-10)
            umol = φmol(zeros(T, 3), [Charge(T[0, 0, 0, √2]...)], tolerance=T(1e-12))
            @test typeof(umol) == T
            @test umol ≈ T(√2 / 1e-12)
        end
    end

    @test_skip ∇φmol
end
