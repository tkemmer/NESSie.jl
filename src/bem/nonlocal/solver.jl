# =========================================================================================
"""
    type NonlocalBEMResult{T} <: BEMResult{T}
        model::SurfaceModel{T}
        u    ::SubArray{T,1}   # [γ₀int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
        q    ::SubArray{T,1}   # [γ₁int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
        w    ::SubArray{T,1}   # [γ₀ext(Ψ)](ξ)     ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
        umol ::Vector{T}       # [γ₀int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
        qmol ::Vector{T}       # [γ₁int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
    end

Result data of the nonlocal solving process to be used for potential computation and
post-processing, with `Ξ` being the list of observation points, that is, the set of
triangle centroids.
"""
type NonlocalBEMResult{T} <: BEMResult{T}
    model::SurfaceModel{T}
    u::SubArray{T,1}
    q::SubArray{T,1}
    w::SubArray{T,1}
    umol::Vector{T}
    qmol::Vector{T}
end


# =========================================================================================
# Documented in bem/local/solver.jl
function solve{T}(
        ::Type{NonlocalES},
        model::SurfaceModel{T},
        LaplaceMod::Module=Rjasanow
    )
    # convenient access
    const elements = model.elements
    const εΩ       = model.params.εΩ
    const εΣ       = model.params.εΣ
    const ε∞       = model.params.ε∞
    const yuk      = yukawa(model.params)

    # create system matrix
    const numelem = length(elements)
    m = zeros(T, 3 * numelem, 3 * numelem)

    # convenient access to 9 blocks of the system matrix
    m11 = view(m,          1:numelem,           1:numelem )
    m12 = view(m,          1:numelem,   1+numelem:2numelem)
    m13 = view(m,          1:numelem,  1+2numelem:3numelem)
    m21 = view(m,  1+numelem:2numelem,          1:numelem )
    m22 = view(m,  1+numelem:2numelem,  1+numelem:2numelem)
    m23 = view(m,  1+numelem:2numelem, 1+2numelem:3numelem)
    m31 = view(m, 1+2numelem:3numelem,          1:numelem )
    m32 = view(m, 1+2numelem:3numelem,  1+numelem:2numelem)
    m33 = view(m, 1+2numelem:3numelem, 1+2numelem:3numelem)

    # initialize the system matrix;
    # since all other components of the system matrix will be premultiplied by 4π,
    # do the same for σ here
    pluseye!(m11, 4π * σ)
    pluseye!(m21, 4π * σ)
    pluseye!(m33, 4π * σ)

    # compute molecular potential for the point charges;
    # molecular potentials are initially premultiplied by 4π⋅ε0⋅εΩ
    const umol = εΩ \   φmol(model)
    const qmol = εΩ \ ∂ₙφmol(model)

    # create right hand side
    rhs = zeros(T, 3 * numelem)

    # convenient access to the first block of rhs
    β = view(rhs, 1:numelem)

    # initialize rhs;
    # again, we apply a prefactor of 4π to σ to match the other components of the vector
    copy!(β, umol)
    scale!(β, -4π * σ)

    # create list of observation points
    Ξ = [e.center for e in elements]

    #=
        generate and apply Kʸ-K
    =#
    buffer = Array{T}(numelem, numelem)
    Radon.regularyukawacoll!(DoubleLayer, buffer, elements, Ξ, yuk)

    # β += (1-εΩ/εΣ)(Kʸ-K)umol
    gemv!(1-εΩ/εΣ, buffer, umol, β)

    # m11 -= Kʸ-K
    axpy!(-1, buffer, m11)

    # m13 += ε∞/εΣ * (Kʸ-K)
    axpy!(ε∞/εΣ, buffer, m13)

    #=
        generate and apply Vʸ-V
    =#
    Radon.regularyukawacoll!(SingleLayer, buffer, elements, Ξ, yuk)

    # β += (εΩ/εΣ - εΩ/ε∞)(Vʸ-V)qmol
    gemv!(εΩ * (1/εΣ - 1/ε∞), buffer, qmol, β)

    # m12 += (εΩ/ε∞ - εΩ/εΣ)(Vʸ-V)
    axpy!(εΩ * (1/ε∞ - 1/εΣ), buffer, m12)

    #=
        generate and apply K
    =#
    LaplaceMod.laplacecoll!(DoubleLayer, buffer, elements, Ξ)

    # β += K⋅umol
    gemv!(one(T), buffer, umol, β)

    # m11 -= K
    axpy!(-1, buffer, m11)

    # m21 += K
    axpy!(1, buffer, m21)

    # m33 -= K
    axpy!(-1, buffer, m33)

    #=
        generate and apply V
    =#
    LaplaceMod.laplacecoll!(SingleLayer, buffer, elements, Ξ)

    # β -= εΩ/ε∞ * V * qmol
    gemv!(-εΩ/ε∞, buffer, qmol, β)

    # m12 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m12)

    # m22 -= V
    axpy!(-1, buffer, m22)

    # m32 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m32)

    # solve system
    cauchy = m\rhs
    NonlocalBEMResult(
        model,
        view(cauchy, 1:          numelem),
        view(cauchy, 1+numelem: 2numelem),
        view(cauchy, 1+2numelem:3numelem),
        umol,
        qmol
    )
end
