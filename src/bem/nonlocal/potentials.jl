#=
    TODO
=#
function φΩ{T}(Ξ::Vector{Vector{T}}, bem::NonlocalBEMResult{T}, LaplaceMod::Module=Rjasanow)
    # result vector
    φ = zeros(T, length(Ξ))

    # φ  = -[K ⋅ u](ξ)
    LaplaceMod.laplacecoll!(DoubleLayer, φ, bem.model.elements, Ξ, bem.u)
    scale!(φ, -1)

    # φ += [V ⋅ q](ξ)
    LaplaceMod.laplacecoll!(SingleLayer, φ, bem.model.elements, Ξ, bem.q)

    # φ *= 2
    scale!(φ, 2)

    # TODO document coefficients (soon...)
    T(1.69e-9 / 4π / ε0) * (T(4π) \ φ + bem.opt.εΩ \ φmol(Ξ, bem.model.charges)) # [φΩ] = V = C/F
end
