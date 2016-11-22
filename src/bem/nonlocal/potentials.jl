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

#=
    TODO
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
    axpy!(1/εΩ, bem.qmol, buf)
    LaplaceMod.laplacecoll!(SingleLayer, φ, elements, Ξ, buf)
    scale!(φ, -εΩ/ε∞)

    # φ += (Vʸ-V)[εΩ(ε0/εΣ - 1/ε∞) ⋅ (q + qmol)](ξ)
    scale!(buf, εΩ * (ε0/εΣ - 1/ε∞))
    Radon.regularyukawacoll!(SingleLayer, φ, elements, Ξ, buf, yuk)

    # φ += K[u + umol](ξ)
    copy!(buf, bem.u)
    axpy!(1/εΩ, bem.umol, buf)
    LaplaceMod.laplacecoll!(DoubleLayer, φ, elements, Ξ, buf)

    # φ += (Kʸ-K)[u + (1-εΩ/εΣ) ⋅ umol - ε∞\εΣ ⋅ w](ξ)
    copy!(buf, bem.u)
    axpy!(1/εΩ-1/εΣ, bem.umol, buf)
    axpy!(-ε∞/εΣ, bem.w, buf)
    Radon.regularyukawacoll!(DoubleLayer, φ, elements, Ξ, buf, yuk)

    # apply coefficients
    scale!(φ, T(1.69e-9 / 4π / 4π / ε0))

    φ # [φΣ] = V = C/F
end
