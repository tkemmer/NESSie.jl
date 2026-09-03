using NESSie.TestModel
using NESSie: potprefactor, euclidean, ec, ε0

const ion_names = [
    "Li", "li", "Na", "na", "K", "k", "Rb", "rb", "Cs", "cs",
    "Mg", "mg", "Ca", "ca", "Sr", "sr", "Ba", "ba"
]

@testset "defaultopt" begin
    for T in testtypes
        @test typeof(defaultopt(BornIon{T})) == Option{T}
    end
end

@testset "bornion" begin
    for T in testtypes, name in ion_names
        @test typeof(bornion(name, T)) == BornIon{T}
    end
end

@testset "Model" begin
    for T in testtypes
        let model = Model(bornion(first(ion_names), T))
            @test model isa Model{T, Triangle{T}}
            @test !isempty(model.nodes)
            @test !isempty(model.elements)
            @test !isempty(model.charges)
            @test model.params == defaultopt(BornIon{T})
        end
    end
end

@testset "potentials and energies (analytical values)" begin
    for T in testtypes
        ion = bornion("Ca", T)
        q    = ion.charge.val
        r_ion = ion.radius
        opt  = ion.params
        εΣ   = opt.εΣ
        ε∞   = opt.ε∞
        λ    = opt.λ
        ξc   = T[0, 0, 0]
        pp   = potprefactor(T)
        tol  = T(1e-10)
        rtol = T == Float64 ? 1.45e-8 : 3.45f-4

        # Observation points
        ξΩ = T[0, 0, r_ion - 0.2]
        ξΣ = T[0, 0, r_ion + 0.2]
        d  = euclidean(ξΣ, ξc)

        # ---- Analytic pre-computations ----
        # molpotential(ξ) = potprefactor * q / |ξ - ξc|
        φ_mol_Σ = pp * q / max(d, tol)
        φ_mol_Ω = pp * q / max(euclidean(ξΩ, ξc), tol)

        # _rfpotential_Ω (LocalES) = potprefactor * q / r * (1/εΣ - 1)
        φ_rf_Ω_local = pp * q / r_ion * (1/εΣ - T(1))

        # _rfpotential_Ω (NonlocalES)
        ν = sqrt(εΣ / ε∞) * r_ion / λ
        φ_rf_Ω_nonlocal = pp * q / r_ion / εΣ *
            (T(1) - εΣ + (εΣ - ε∞)/ε∞ * sinh(ν)/ν * exp(-ν))

        # _rfpotential_Σ (LocalES) = molpotential * (1/εΣ - 1)
        φ_rf_Σ_local = φ_mol_Σ * (1/εΣ - T(1))

        # _rfpotential_Σ (NonlocalES) = molpotential * (c - 1)
        c = (T(1) + (εΣ - ε∞)/ε∞ * sinh(ν)/ν * exp(-ν * d / r_ion)) / εΣ
        φ_rf_Σ_nonlocal = φ_mol_Σ * (c - T(1))

        # _espotential_Σ (LocalES) = molpotential / εΣ
        φ_es_Σ_local = φ_mol_Σ / εΣ

        # _espotential_Σ (NonlocalES) = molpotential * c
        φ_es_Σ_nonlocal = φ_mol_Σ * c

        # ---- rfpotential :Ω ----
        @test rfpotential(:Ω, LocalES, ξΩ, ion) ≈ φ_rf_Ω_local rtol=rtol
        @test rfpotential(:Ω, NonlocalES, ξΩ, ion) ≈ φ_rf_Ω_nonlocal rtol=rtol

        # ---- rfpotential :Γ (maps to :Ω) ----
        @test rfpotential(:Γ, LocalES, ξΩ, ion) ≈ φ_rf_Ω_local rtol=rtol
        @test rfpotential(:Γ, NonlocalES, ξΩ, ion) ≈ φ_rf_Ω_nonlocal rtol=rtol

        # ---- rfpotential :Σ ----
        @test rfpotential(:Σ, LocalES, ξΣ, ion) ≈ φ_rf_Σ_local rtol=rtol
        @test rfpotential(:Σ, NonlocalES, ξΣ, ion) ≈ φ_rf_Σ_nonlocal rtol=rtol

        # ---- rfpotential auto-location ----
        @test rfpotential(LocalES, ξΩ, ion) ≈ φ_rf_Ω_local rtol=rtol
        @test rfpotential(NonlocalES, ξΩ, ion) ≈ φ_rf_Ω_nonlocal rtol=rtol
        @test rfpotential(LocalES, ξΣ, ion) ≈ φ_rf_Σ_local rtol=rtol
        @test rfpotential(NonlocalES, ξΣ, ion) ≈ φ_rf_Σ_nonlocal rtol=rtol

        # ---- molpotential ----
        @test molpotential(ξΣ, ion) ≈ φ_mol_Σ rtol=rtol
        @test molpotential(ξΩ, ion) ≈ φ_mol_Ω rtol=rtol

        # ---- espotential :Ω ----
        @test espotential(:Ω, LocalES, ξΩ, ion) ≈ φ_rf_Ω_local + φ_mol_Ω rtol=rtol
        @test espotential(:Ω, NonlocalES, ξΩ, ion) ≈ φ_rf_Ω_nonlocal + φ_mol_Ω rtol=rtol

        # ---- espotential :Γ (maps to :Ω) ----
        @test espotential(:Γ, LocalES, ξΩ, ion) ≈ φ_rf_Ω_local + φ_mol_Ω rtol=rtol
        @test espotential(:Γ, NonlocalES, ξΩ, ion) ≈ φ_rf_Ω_nonlocal + φ_mol_Ω rtol=rtol

        # ---- espotential :Σ ----
        @test espotential(:Σ, LocalES, ξΣ, ion) ≈ φ_es_Σ_local rtol=rtol
        @test espotential(:Σ, NonlocalES, ξΣ, ion) ≈ φ_es_Σ_nonlocal rtol=rtol

        # ---- espotential auto-location ----
        @test espotential(LocalES, ξΩ, ion) ≈ φ_rf_Ω_local + φ_mol_Ω rtol=rtol
        @test espotential(NonlocalES, ξΩ, ion) ≈ φ_rf_Ω_nonlocal + φ_mol_Ω rtol=rtol
        @test espotential(LocalES, ξΣ, ion) ≈ φ_es_Σ_local rtol=rtol
        @test espotential(NonlocalES, ξΣ, ion) ≈ φ_es_Σ_nonlocal rtol=rtol

        # ---- rfenergy ----
        expected_E_local = q * φ_rf_Ω_local * ec * T(6.022140857e10) / T(2)
        expected_E_nonlocal = q * φ_rf_Ω_nonlocal * ec * T(6.022140857e10) / T(2)
        @test rfenergy(LocalES, ion) ≈ expected_E_local rtol=rtol
        @test rfenergy(NonlocalES, ion) ≈ expected_E_nonlocal rtol=rtol
    end
end

@testset "potentials and energies (interface variant identity)" begin
    for T in testtypes, lt in (LocalES, NonlocalES)
        ion   = bornion("Ca", T)
        Ξ     = LinRange(T[0, 0, 0], T[0, 0, 2], 20)
        ξΩ    = T[0, 0, ion.radius - 0.2]
        ξΓ    = T[0, 0, ion.radius]
        ξΣ    = T[0, 0, ion.radius + 0.2]

        # molpotential
        res = molpotential(Ξ, ion)
        @test res isa Vector{T}
        @test length(res) == length(Ξ)

        res = molpotential(reshape(Ξ, 1, length(Ξ)), ion)
        @test res isa Matrix{T}
        @test res == transpose(molpotential(Ξ, ion))

        res = molpotential((ξ for ξ in Ξ), ion)
        @test res isa Vector{T}
        @test res == molpotential(Ξ, ion)

        # rfpotential with explicit domain
        res = rfpotential(:Ω, lt, ξΩ, ion)
        @test res isa T

        res = rfpotential(:Ω, lt, [ξΩ], ion)
        @test res isa Vector{T}
        @test only(res) == rfpotential(:Ω, lt, ξΩ, ion)

        res = rfpotential(:Ω, lt, reshape([ξΩ], 1, 1), ion)
        @test res isa Matrix{T}
        @test only(res) == rfpotential(:Ω, lt, ξΩ, ion)

        res = rfpotential(:Ω, lt, (ξ for ξ in [ξΩ]), ion)
        @test res isa Vector{T}
        @test only(res) == rfpotential(:Ω, lt, ξΩ, ion)

        res = rfpotential(:Γ, lt, ξΓ, ion)
        @test res isa T
        @test res == rfpotential(:Ω, lt, ξΓ, ion)

        res = rfpotential(:Γ, lt, [ξΓ], ion)
        @test only(res) == rfpotential(:Ω, lt, ξΓ, ion)

        res = rfpotential(:Σ, lt, ξΣ, ion)
        @test res isa T

        res = rfpotential(:Σ, lt, [ξΣ], ion)
        @test res isa Vector{T}
        @test only(res) == rfpotential(:Σ, lt, ξΣ, ion)

        res = rfpotential(:Σ, lt, reshape([ξΣ], 1, 1), ion)
        @test res isa Matrix{T}
        @test only(res) == rfpotential(:Σ, lt, ξΣ, ion)

        res = rfpotential(:Σ, lt, (ξ for ξ in [ξΣ]), ion)
        @test res isa Vector{T}
        @test only(res) == rfpotential(:Σ, lt, ξΣ, ion)

        # rfpotential with auto-location
        res = rfpotential(lt, ξΩ, ion)
        @test res isa T
        @test res == rfpotential(:Ω, lt, ξΩ, ion)

        res = rfpotential(lt, ξΓ, ion)
        @test res isa T
        @test res == rfpotential(:Γ, lt, ξΓ, ion)

        res = rfpotential(lt, ξΣ, ion)
        @test res isa T
        @test res == rfpotential(:Σ, lt, ξΣ, ion)

        res = rfpotential(lt, Ξ, ion)
        @test res isa Vector{T}
        @test length(res) == length(Ξ)

        res = rfpotential(lt, reshape(Ξ, 1, length(Ξ)), ion)
        @test res isa Matrix{T}
        @test res == transpose(rfpotential(lt, Ξ, ion))

        res = rfpotential(lt, (ξ for ξ in Ξ), ion)
        @test res isa Vector{T}
        @test res == rfpotential(lt, Ξ, ion)

        # espotential with explicit domain
        res = espotential(:Ω, lt, ξΩ, ion)
        @test res isa T

        res = espotential(:Ω, lt, [ξΩ], ion)
        @test res isa Vector{T}
        @test only(res) == espotential(:Ω, lt, ξΩ, ion)

        res = espotential(:Γ, lt, ξΓ, ion)
        @test res isa T
        @test res == espotential(:Ω, lt, ξΓ, ion)

        res = espotential(:Γ, lt, [ξΓ], ion)
        @test only(res) == espotential(:Ω, lt, ξΓ, ion)

        res = espotential(:Σ, lt, ξΣ, ion)
        @test res isa T

        res = espotential(:Σ, lt, [ξΣ], ion)
        @test res isa Vector{T}
        @test only(res) == espotential(:Σ, lt, ξΣ, ion)

        res = espotential(:Σ, lt, reshape([ξΣ], 1, 1), ion)
        @test res isa Matrix{T}
        @test only(res) == espotential(:Σ, lt, ξΣ, ion)

        res = espotential(:Σ, lt, (ξ for ξ in [ξΣ]), ion)
        @test res isa Vector{T}
        @test only(res) == espotential(:Σ, lt, ξΣ, ion)

        # espotential with auto-location
        res = espotential(lt, ξΩ, ion)
        @test res isa T
        @test res == espotential(:Ω, lt, ξΩ, ion)

        res = espotential(lt, ξΓ, ion)
        @test res isa T
        @test res == espotential(:Γ, lt, ξΓ, ion)

        res = espotential(lt, ξΣ, ion)
        @test res isa T
        @test res == espotential(:Σ, lt, ξΣ, ion)

        res = espotential(lt, Ξ, ion)
        @test res isa Vector{T}
        @test length(res) == length(Ξ)

        res = espotential(lt, reshape(Ξ, 1, length(Ξ)), ion)
        @test res isa Matrix{T}
        @test res == transpose(espotential(lt, Ξ, ion))

        res = espotential(lt, (ξ for ξ in Ξ), ion)
        @test res isa Vector{T}
        @test res == espotential(lt, Ξ, ion)

        # espotential = rfpotential + molpotential identity
        res = espotential(lt, Ξ, ion)
        @test res ≈ rfpotential(lt, Ξ, ion) .+ molpotential(Ξ, ion)

        res = espotential(lt, reshape(Ξ, 1, length(Ξ)), ion)
        @test res ≈ rfpotential(lt, reshape(Ξ, 1, length(Ξ)), ion) .+ molpotential(reshape(Ξ, 1, length(Ξ)), ion)

        res = espotential(lt, (ξ for ξ in Ξ), ion)
        @test res ≈ rfpotential(lt, (ξ for ξ in Ξ), ion) .+ molpotential((ξ for ξ in Ξ), ion)

        # rfenergy
        res = rfenergy(lt, ion)
        @test res isa T

        # error handling
        @test_throws ErrorException rfpotential(:Foo, lt, ξΩ, ion)
        @test_throws ErrorException rfpotential(:Foo, lt, Ξ, ion)
        @test_throws ErrorException espotential(:Foo, lt, ξΩ, ion)
        @test_throws ErrorException espotential(:Foo, lt, Ξ, ion)
    end
end
