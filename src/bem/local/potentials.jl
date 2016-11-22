#=
    TODO
=#
function φΩ{T}(Ξ::Vector{Vector{T}}, bem::LocalBEMResult{T}, LaplaceMod::Module=Rjasanow)
    # result vector
    φ = zeros(T, length(Ξ))

    # φ  = -[W ⋅ u](ξ)
    LaplaceMod.laplacecoll!(DoubleLayer, φ, bem.model.elements, Ξ, bem.u)
    scale!(φ, -1)

    # φ += [Vtilde ⋅ q](ξ)
    LaplaceMod.laplacecoll!(SingleLayer, φ, bem.model.elements, Ξ, bem.q)

    # TODO document coefficients (soon...)
    T(1.69e-9 / 4π / ε0) * (T(4π) \ φ + bem.opt.εΩ \ φmol(Ξ, bem.model.charges)) # [φΩ] = V = C/F
end

#=
    TODO
=#
function φΣ{T}(Ξ::Vector{Vector{T}}, bem::LocalBEMResult{T}, LaplaceMod::Module=Rjasanow)
    # result vector
    φ = zeros(T, length(Ξ))

    # φ  = -εΩ/εΣ ⋅ [Vtilde ⋅ (q + qmol)](ξ)
    LaplaceMod.laplacecoll!(SingleLayer, φ, bem.model.elements, Ξ, bem.q + bem.opt.εΩ \ bem.qmol)
    scale!(φ, -bem.opt.εΩ/bem.opt.εΣ)

    # φ += [W ⋅ (u + umol)](ξ)
    LaplaceMod.laplacecoll!(DoubleLayer, φ, bem.model.elements, Ξ, bem.u + bem.opt.εΩ \ bem.umol)
    scale!(φ, T(1.69e-9 / 4π / 4π / ε0)) # TODO doc!

    φ # [φΣ] = V = C/F
end
