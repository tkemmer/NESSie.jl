# =========================================================================================
"""
    struct NonlocalBEMResult{T, E} <: BEMResult{T, E}
        model::Model{T, E}
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
struct NonlocalBEMResult{T, E} <: BEMResult{T, E}
    model::Model{T, E}
    u::SubArray{T,1}
    q::SubArray{T,1}
    w::SubArray{T,1}
    umol::Vector{T}
    qmol::Vector{T}
end


# =========================================================================================
# Documented in bem/local/solver.jl
function solve(
                  ::Type{NonlocalES},
        model     ::Model{T, Triangle{T}}
    ) where T
    # convenient access
    elements = model.elements
    εΩ       = model.params.εΩ
    εΣ       = model.params.εΣ
    ε∞       = model.params.ε∞
    yuk      = yukawa(model.params)

    # create system matrix
    numelem = length(elements)
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
    pluseye!(m11, T(4π * σ))
    pluseye!(m21, T(4π * σ))
    pluseye!(m33, T(4π * σ))

    # compute molecular potential for the point charges;
    # molecular potentials are initially premultiplied by 4π⋅ε0⋅εΩ
    umol = εΩ \   φmol(model)
    qmol = εΩ \ ∂ₙφmol(model)

    # create right hand side
    rhs = zeros(T, 3 * numelem)

    # convenient access to the first block of rhs
    β = view(rhs, 1:numelem)

    # initialize rhs;
    # again, we apply a prefactor of 4π to σ to match the other components of the vector
    copyto!(β, umol)
    rmul!(β, -T(4π * σ))

    # create list of observation points
    Ξ = [e.center for e in elements]

    #=
        generate and apply Kʸ-K
    =#
    buffer = Array{T}(undef, numelem, numelem)
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
    Rjasanow.laplacecoll!(DoubleLayer, buffer, elements, Ξ)

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
    Rjasanow.laplacecoll!(SingleLayer, buffer, elements, Ξ)

    # β -= εΩ/ε∞ * V * qmol
    gemv!(-εΩ/ε∞, buffer, qmol, β)

    # m12 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m12)

    # m22 -= V
    axpy!(-1, buffer, m22)

    # m32 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m32)

    # solve system
    cauchy = m \ rhs

    NonlocalBEMResult(
        model,
        view(cauchy, 1:          numelem),
        view(cauchy, 1+numelem: 2numelem),
        view(cauchy, 1+2numelem:3numelem),
        umol,
        qmol
    )
end

# TODO
function solve_implicit(
                  ::Type{NonlocalES},
        model     ::Model{T, Triangle{T}}
    ) where T
    # convenient access
    elements = model.elements
    εΩ       = model.params.εΩ
    εΣ       = model.params.εΣ
    ε∞       = model.params.ε∞
    yuk      = yukawa(model.params)
    numelem  = length(elements)

    # compute molecular potential for the point charges;
    # molecular potentials are initially premultiplied by 4π⋅ε0⋅εΩ
    umol = εΩ \   φmol(model)
    qmol = εΩ \ ∂ₙφmol(model)

    # observation points and elements
    Ξ      = [e.center for e in elements]
    iΞ     = collect(enumerate(Ξ))
    ielems = collect(enumerate(elements))


    #=
        right-hand side
    =#

    # -[1 - σ - Kʸ + εΩ/εΣ (Kʸ - K)]
    Ks = InteractionMatrix(T, iΞ, ielems,
        ((i, ξ), (j, elem)) -> -T(i == j) * 4π * (1-σ) +
            (1 - εΩ/εΣ) * Radon.regularyukawacoll(DoubleLayer, ξ, elem, yuk) +
            Rjasanow.laplacecoll(DoubleLayer, ξ, elem)
    )

    # -[εΩ/ε∞ Vʸ - εΩ/εΣ (Vʸ - V)]
    Vs = InteractionMatrix(T, Ξ, elements,
        (ξ, elem) -> εΩ * (1/εΣ - 1/ε∞) *
            Radon.regularyukawacoll(SingleLayer, ξ, elem, yuk) -
            εΩ/ε∞ * Rjasanow.laplacecoll(SingleLayer, ξ, elem)
    )

    # rhs = [β, 0 , 0]ᵀ
    rhs = BlockMatrix(3, 1,
        reshape(Ks * umol + Vs * qmol, (numelem, 1)),
        FixedValueArray(zero(T), numelem, 1),
        FixedValueArray(zero(T), numelem, 1)
    )


    #=
        system matrix
    =#

    # M₁₁ = 1 - σ - Kʸ
    M₁₁ = InteractionMatrix(T, iΞ, ielems,
        ((i, ξ), (j, elem)) -> T(i == j) * 4π * (1-σ) -
            Radon.regularyukawacoll(DoubleLayer, ξ, elem, yuk) -
            Rjasanow.laplacecoll(DoubleLayer, ξ, elem)
    )

    # M₁₂ = εΩ/ε∞ Vʸ - εΩ/εΣ (Vʸ - V)
    M₁₂ = InteractionMatrix(T, Ξ, elements,
        (ξ, elem) -> εΩ * (1/ε∞ - 1/εΣ) *
            Radon.regularyukawacoll(SingleLayer, ξ, elem, yuk) +
            εΩ / ε∞ * Rjasanow.laplacecoll(SingleLayer, ξ, elem)
    )

    # M₁₃ = ε∞/εΣ (Kʸ - K)
    M₁₃ = InteractionMatrix(T, Ξ, elements,
        (ξ, elem) -> ε∞/εΣ * Radon.regularyukawacoll(DoubleLayer, ξ, elem, yuk)
    )

    # M₂₁ = σ + K
    M₂₁ = InteractionMatrix(T, iΞ, ielems,
        ((i, ξ), (j, elem)) -> T(i == j) * 4π * σ  +
            Rjasanow.laplacecoll(DoubleLayer, ξ, elem)
    )

    # M₂₂ = -V
    M₂₂ = InteractionMatrix(T, Ξ, elements,
        (ξ, elem) -> -Rjasanow.laplacecoll(SingleLayer, ξ, elem)
    )

    # M₂₃ = M₃₁ = 0
    M₂₃ = FixedValueArray(zero(T), numelem, numelem)
    M₃₁ = FixedValueArray(zero(T), numelem, numelem)

    # M₃₂ = εΩ/ε∞ V
    M₃₂ = InteractionMatrix(T, Ξ, elements,
        (ξ, elem) -> εΩ/ε∞ * Rjasanow.laplacecoll(SingleLayer, ξ, elem)
    )

    # M₃₃ = 1 - σ - K
    M₃₃ = InteractionMatrix(T, iΞ, ielems,
        ((i, ξ), (j, elem)) -> T(i == j) * 4π * (1-σ) -
            Rjasanow.laplacecoll(DoubleLayer, ξ, elem)
    )

    M = BlockMatrix(3, 3,
        M₁₁, M₁₂, M₁₃,
        M₂₁, M₂₂, M₂₃,
        M₃₁, M₃₂, M₃₃
    )

    # solve system
    cauchy = M \ rhs

    NonlocalBEMResult(
        model,
        view(cauchy, 1:          numelem),
        view(cauchy, 1+numelem: 2numelem),
        view(cauchy, 1+2numelem:3numelem),
        umol,
        qmol
    )
end
