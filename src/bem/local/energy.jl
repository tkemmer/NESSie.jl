#=
    Computes the local reaction field energy W* as

        W* = ∫φ*⋅ρ dΩ

    where φ* is the reaction field and ρ is the corresponding charge distribution.

    @param bem
        Result of `solvelocal`
    @param LaplaceMod
        Module to be used for Laplace potential; Valid values: Radon, Rjasanow
    @return T ([W*] = kJ/mol)
=#
function rfenergy{T}(bem::LocalBEMResult{T}, LaplaceMod::Module=Rjasanow)
    qposs = [charge.pos for charge in bem.model.charges]
    qvals = [charge.val for charge in bem.model.charges]

    # reaction field energy per point charge
    wstar = zeros(T, length(bem.model.charges))

    # W* = -[W ⋅ u](ξ)
    LaplaceMod.laplacecoll!(DoubleLayer, wstar, bem.model.elements, qposs, bem.u)
    scale!(wstar, -1)

    # W* += [V ⋅ q](ξ)
    LaplaceMod.laplacecoll!(SingleLayer, wstar, bem.model.elements, qposs, bem.q)

    # Apply ρ, integrate over Ω and apply remaining prefactors (in order)
    # ▶ 4π        for Vtilde, W
    # (in potprefactor:)
    # ▶ 4π⋅ε0     for u and q
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    #
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 6.022e23  for Avogadro constant Nₐ; [Nₐ] = 1/mol
    # ▶ 1e-3      for the conversion 1/J → 1/kJ
    # ▶ 0.5       since we have counted all interactions twice
    wstar ⋅ qvals / 4π * potprefactor(T) * ec * 6.022140857e10 / 2
end
