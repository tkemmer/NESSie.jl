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

    T(1.69e-9 / 4π / ε0) * (4π \ φ + bem.opt.εΩ \ φmol(Ξ, bem.model.charges)) # [φΩ] = V = C/F
end
