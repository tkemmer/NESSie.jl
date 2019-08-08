# =========================================================================================
# Documented in bem/local/potentials.jl
function φΩ(
        Ξ         ::Vector{Vector{T}},
        bem       ::NonlocalBEMResult{T}
    ) where T
    # result vector
    φ = zeros(T, length(Ξ))

    # φ  = -[K ⋅ u](ξ)
    Rjasanow.laplacecoll!(DoubleLayer, φ, bem.model.elements, Ξ, bem.u)
    scale!(φ, -1)

    # φ += [V ⋅ q](ξ)
    Rjasanow.laplacecoll!(SingleLayer, φ, bem.model.elements, Ξ, bem.q)

    # φ *= 1/4π
    # (K and V were premultiplied by 4π! 4π⋅ε0 from u and q still to be applied)
    scale!(φ, 1/4π)

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
# Documented in bem/local/potentials.jl
function φΣ(
        Ξ         ::Vector{Vector{T}},
        bem       ::NonlocalBEMResult{T}
    ) where T
    # result vector
    φ = zeros(T, length(Ξ))

    # convenience aliases
    εΩ  = bem.model.params.εΩ
    εΣ  = bem.model.params.εΣ
    ε∞  = bem.model.params.ε∞
    yuk = yukawa(bem.model.params)
    elements = bem.model.elements

    buf = Array{T}(length(elements))

    # φ  = -V[εΩ/ε∞ ⋅ (q + qmol)](ξ)
    copy!(buf, bem.q)
    axpy!(1, bem.qmol, buf)
    Rjasanow.laplacecoll!(SingleLayer, φ, elements, Ξ, buf)
    scale!(φ, -εΩ/ε∞)

    # φ += (Vʸ-V)[εΩ(1/εΣ - 1/ε∞) ⋅ (q + qmol)](ξ)
    scale!(buf, εΩ * (1/εΣ - 1/ε∞))
    Radon.regularyukawacoll!(SingleLayer, φ, elements, Ξ, yuk, buf)

    # φ += K[u + umol](ξ)
    copy!(buf, bem.u)
    axpy!(1, bem.umol, buf)
    Rjasanow.laplacecoll!(DoubleLayer, φ, elements, Ξ, buf)

    # φ += (Kʸ-K)[u + (1-εΩ/εΣ) ⋅ umol - ε∞\εΣ ⋅ w](ξ)
    copy!(buf, bem.u)
    axpy!(1-εΩ/εΣ, bem.umol, buf)
    axpy!(-ε∞/εΣ, bem.w, buf)
    Radon.regularyukawacoll!(DoubleLayer, φ, elements, Ξ, yuk, buf)

    # Apply remaining prefactors:
    # ▶ 4π        for V, K, (Vʸ-V), and (Kʸ-K)
    # ▶ 4π⋅ε0     for u, q, w, umol, and qmol
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    scale!(φ, potprefactor(T) / 4π)

    φ
end
