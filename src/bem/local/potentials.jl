# =========================================================================================
"""
    φΩ{T, B <: BEMResult{T}}(
        Ξ         ::Vector{Vector{T}},
        bem       ::B,
        LaplaceMod::Module=Rjasanow
    )

Computes the local or nonlocal interior electrostatic potential ``φ\_Ω`` for the given set
of observation points `Ξ`.

!!! warning
    This function does not verify whether all points in `Ξ` are located in ``Ω``!

# Arguments
 * `LaplaceMod` Module to be used for Laplace potential; Valid values: `Radon, Rjasanow`

# Unit
``V = \\frac{C}{F}``

# Return type
`Vector{T}`
"""
function φΩ{T}(Ξ::Vector{Vector{T}}, bem::LocalBEMResult{T}, LaplaceMod::Module=Rjasanow)
    # result vector
    φ = zeros(T, length(Ξ))

    # φ  = -[W ⋅ u](ξ)
    LaplaceMod.laplacecoll!(DoubleLayer, φ, bem.model.elements, Ξ, bem.u)
    scale!(φ, -1)

    # φ += [Vtilde ⋅ q](ξ)
    LaplaceMod.laplacecoll!(SingleLayer, φ, bem.model.elements, Ξ, bem.q)

    # φ *= 1/4π
    # (W and Vtilde were premultiplied by 4π! 4π⋅ε0 from u and q still to be applied)
    scale!(φ, 1 / 4π)

    # φ += 1/εΩ ⋅ φ*mol(ξ)
    # (φ*mol was premultiplied by 4π⋅ε0⋅εΩ; 4π⋅ε0 remain to be applied)
    axpy!(1/bem.model.params.εΩ, φmol(Ξ, bem.model.charges), φ)

    # Apply remaining prefactors:
    # ▶ 4π⋅ε0     for u, q, and umol
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    scale!(φ, potprefactor(T))

    φ
end


# =========================================================================================
"""
    φΣ{T, B <: BEMResult{T}}(
        Ξ         ::Vector{Vector{T}},
        bem       ::B,
        LaplaceMod::Module=Rjasanow
    )

Computes the local or nonlocal exterior electrostatic potential ``φ\_Σ`` for the given set
of observation points `Ξ`.

!!! warning
    This function does not verify whether all points in `Ξ` are located in ``Σ``!

# Arguments
 * `LaplaceMod` Module to be used for Laplace potential; Valid values: `Radon, Rjasanow`

# Unit
``V = \\frac{C}{F}``

# Return type
`Vector{T}`
"""
function φΣ{T}(Ξ::Vector{Vector{T}}, bem::LocalBEMResult{T}, LaplaceMod::Module=Rjasanow)
    # result vector
    φ = zeros(T, length(Ξ))
    buf = Array{T}(length(bem.model.elements))

    # φ  = -εΩ/εΣ ⋅ [Vtilde ⋅ (q + qmol)](ξ)
    copy!(buf, bem.q)
    axpy!(1, bem.qmol, buf)
    LaplaceMod.laplacecoll!(SingleLayer, φ, bem.model.elements, Ξ, buf)
    scale!(φ, -bem.model.params.εΩ/bem.model.params.εΣ)

    # φ += [W ⋅ (u + umol)](ξ)
    copy!(buf, bem.u)
    axpy!(1, bem.umol, buf)
    LaplaceMod.laplacecoll!(DoubleLayer, φ, bem.model.elements, Ξ, buf)

    # Apply remaining prefactors:
    # ▶ 4π        for Vtilde, W
    # ▶ 4π⋅ε0     for u, q, umol, and qmol
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    scale!(φ, potprefactor(T) / 4π)

    φ
end
