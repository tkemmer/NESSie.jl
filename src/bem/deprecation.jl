# removed in v1.4
@deprecate solve_implicit(lt, model) solve(lt, model; method=:gmres)

# removed in v1.5
@deprecate NESSie.φΓ(ξorΞ, bem::BEMResult) espotential(:Γ, ξorΞ, bem) false
@deprecate NESSie.φΩ(ξorΞ, bem::BEMResult) espotential(:Ω, ξorΞ, bem) false
@deprecate NESSie.φΣ(ξorΞ, bem::BEMResult) espotential(:Σ, ξorΞ, bem) false
