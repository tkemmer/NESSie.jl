#=
    Computes the interior nonlocal electrostatic potential φΩ for the given set of points Ξ. Note that the function
    does not verify whether the given points actually lie inside of the protein.

    @param Ξ
        List of observation points
    @param bem
        Result of `solvelocal`
    @param LaplaceMod
        Module to be used for Laplace potential; Valid values: Radon, Rjasanow
    @return Vector{T}  ([φΩ] = V = C/F)
=#
function φΩ{T}(Ξ::Vector{Vector{T}}, bem::NonlocalBEMResult{T}, LaplaceMod::Module=Rjasanow)
    # result vector
    φ = zeros(T, length(Ξ))

    # φ  = -[K ⋅ u](ξ)
    LaplaceMod.laplacecoll!(DoubleLayer, φ, bem.model.elements, Ξ, bem.u)
    scale!(φ, -1)

    # φ += [V ⋅ q](ξ)
    LaplaceMod.laplacecoll!(SingleLayer, φ, bem.model.elements, Ξ, bem.q)

    # φ *= 2/4π  (K and V were premultiplied by 4π! 4π⋅ε0 from u and q still to be applied)
    scale!(φ, 1/2π)

    # φ += 1/εΩ ⋅ φ*mol(ξ)  (φ*mol was premultiplied by 4π⋅ε0⋅εΩ; 4π⋅ε0 remain to be applied)
    axpy!(1/bem.opt.εΩ, φmol(Ξ, bem.model.charges), φ)

    # Apply remaining prefactors:
    # ▶ 4π⋅ε0     for u, q, and umol
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    scale!(φ, potprefactor(T))

    φ
end

#=
    Computes the exterior nonlocal electrostatic potential φΣ for the given set of points Ξ. Note that the function
    does not verify whether the given points actually lie in the surrounding space.

    @param Ξ
        List of observation points
    @param bem
        Result of `solvelocal`
    @param LaplaceMod
        Module to be used for Laplace potential; Valid values: Radon, Rjasanow
    @return Vector{T}  ([φΣ] = V = C/F)
=#
function φΣ{T}(Ξ::Vector{Vector{T}}, bem::NonlocalBEMResult{T}, LaplaceMod::Module=Rjasanow)
    # result vector
    φ = zeros(T, length(Ξ))

    # convenience aliases
    const εΩ  = bem.opt.εΩ
    const εΣ  = bem.opt.εΣ
    const ε∞  = bem.opt.ε∞
    const yuk = yukawa(bem.opt)
    const elements = bem.model.elements

    buf = Array(T, length(elements))

    # φ  = -V[εΩ/ε∞ ⋅ (q + qmol)](ξ)
    copy!(buf, bem.q)
    axpy!(1, bem.qmol, buf)
    LaplaceMod.laplacecoll!(SingleLayer, φ, elements, Ξ, buf)
    scale!(φ, -εΩ/ε∞)

    # φ += (Vʸ-V)[εΩ(1/εΣ - 1/ε∞) ⋅ (q + qmol)](ξ)
    scale!(buf, εΩ * (1/εΣ - 1/ε∞))
    Radon.regularyukawacoll!(SingleLayer, φ, elements, Ξ, buf, yuk)

    # φ += K[u + umol](ξ)
    copy!(buf, bem.u)
    axpy!(1, bem.umol, buf)
    LaplaceMod.laplacecoll!(DoubleLayer, φ, elements, Ξ, buf)

    # φ += (Kʸ-K)[u + (1-εΩ/εΣ) ⋅ umol - ε∞\εΣ ⋅ w](ξ)
    copy!(buf, bem.u)
    axpy!(1-εΩ/εΣ, bem.umol, buf)
    axpy!(-ε∞/εΣ, bem.w, buf)
    Radon.regularyukawacoll!(DoubleLayer, φ, elements, Ξ, buf, yuk)

    # Apply remaining prefactors:
    # ▶ 4π        for V, K, (Vʸ-V), and (Kʸ-K)
    # ▶ 4π⋅ε0     for u, q, w, umol, and qmol
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    scale!(φ, potprefactor(T) / 4π)

    φ
end
