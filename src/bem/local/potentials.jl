#=
    Computes the interior local electrostatic potential φΩ for the given set of points Ξ. Note that the function does
    not verify whether the given points actually lie inside of the protein.

    @param Ξ
        List of observation points
    @param bem
        Result of `solvelocal`
    @param LaplaceMod
        Module to be used for Laplace potential; Valid values: Radon, Rjasanow
    @return Vector{T}  ([φΩ] = V = C/F)
=#
function φΩ{T}(Ξ::Vector{Vector{T}}, bem::LocalBEMResult{T}, LaplaceMod::Module=Rjasanow)
    # result vector
    φ = zeros(T, length(Ξ))

    # φ  = -[W ⋅ u](ξ)
    LaplaceMod.laplacecoll!(DoubleLayer, φ, bem.model.elements, Ξ, bem.u)
    scale!(φ, -1)

    # φ += [Vtilde ⋅ q](ξ)
    LaplaceMod.laplacecoll!(SingleLayer, φ, bem.model.elements, Ξ, bem.q)

    # φ *= 1/4π  (W and Vtilde were premultiplied by 4π! 4π⋅ε0 from u and q still to be applied)
    scale!(φ, T(1 / 4π))

    # φ += 1/εΩ ⋅ φ*mol(ξ)  (φ*mol was premultiplied by 4π⋅ε0⋅εΩ; 4π⋅ε0 remain to be applied)
    axpy!(1/bem.opt.εΩ, φmol(Ξ, bem.model.charges), φ)

    # Apply remaining prefactors:
    # ▶ 4π⋅ε0     for u, q, and umol
    # ▶ 1.69e-19  for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    scale!(φ, T(1.69e-9 / 4π / ε0))

    φ
end

#=
    Computes the exterior local electrostatic potential φΣ for the given set of points Ξ. Note that the function does
    not verify whether the given points actually lie in the surrounding space.

    @param Ξ
        List of observation points
    @param bem
        Result of `solvelocal`
    @param LaplaceMod
        Module to be used for Laplace potential; Valid values: Radon, Rjasanow
    @return Vector{T}  ([φΣ] = V = C/F)
=#
function φΣ{T}(Ξ::Vector{Vector{T}}, bem::LocalBEMResult{T}, LaplaceMod::Module=Rjasanow)
    # result vector
    φ = zeros(T, length(Ξ))
    buf = Array(T, length(bem.model.elements))

    # φ  = -εΩ/εΣ ⋅ [Vtilde ⋅ (q + qmol)](ξ)
    copy!(buf, bem.q)
    axpy!(1/bem.opt.εΩ, bem.qmol, buf)
    LaplaceMod.laplacecoll!(SingleLayer, φ, bem.model.elements, Ξ, buf)
    scale!(φ, -bem.opt.εΩ/bem.opt.εΣ)

    # φ += [W ⋅ (u + umol)](ξ)
    copy!(buf, bem.u)
    axpy!(1/bem.opt.εΩ, bem.umol, buf)
    LaplaceMod.laplacecoll!(DoubleLayer, φ, bem.model.elements, Ξ, buf)

    # Apply remaining prefactors:
    # ▶ 4π        for Vtilde, W
    # ▶ 4π⋅ε0     for u, q, umol, and qmol
    # ▶ 1.69e-19  for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    scale!(φ, T(1.69e-9 / 4π / 4π / ε0))

    φ
end
