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

    # shortcuts
    params  = model.params
    elems   = model.elements
    ielems  = collect(enumerate(elems))
    Ξ       = [e.center for e in elems]
    iΞ      = collect(enumerate(Ξ))
    numelem = length(elems)

    # compute molecular potential for the point charges;
    # molecular potentials are initially premultiplied by 4π⋅ε0⋅εΩ
    umol = params.εΩ \   φmol(model)
    qmol = params.εΩ \ ∂ₙφmol(model)

    # rhs = [β, 0 , 0]ᵀ
    Ks = InteractionMatrix(iΞ, ielems, KSfun{T}(params))
    Vs = InteractionMatrix(Ξ,  elems,  VSfun{T}(params))

    rhs = BlockMatrix(3, 1,
        reshape(Ks * umol + Vs * qmol, (numelem, 1)),
        FixedValueArray(zero(T), numelem, 1),
        FixedValueArray(zero(T), numelem, 1)
    )

    # system matrix
    M₁₁ = InteractionMatrix(iΞ, ielems, M11fun{T}(params))
    M₁₂ = InteractionMatrix(Ξ,  elems,  M12fun{T}(params))
    M₁₃ = InteractionMatrix(Ξ,  elems,  M13fun{T}(params))
    M₂₁ = InteractionMatrix(iΞ, ielems, M21fun{T}())
    M₂₂ = InteractionMatrix(Ξ,  elems,  M22fun{T}())
    M₂₃ = FixedValueArray(zero(T), numelem, numelem)
    M₃₁ = FixedValueArray(zero(T), numelem, numelem)
    M₃₂ = InteractionMatrix(Ξ,  elems,  M32fun{T}(params))
    M₃₃ = InteractionMatrix(iΞ, ielems, M33fun{T}())

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

# Vs = -[εΩ/ε∞ Vʸ - εΩ/εΣ (Vʸ - V)]
struct VSfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T}
    opt::Option{T}
end
function (f::VSfun{T})(ξ::Vector{T}, elem::Triangle{T}) where T
    f.opt.εΩ * (1/f.opt.εΣ - 1/f.opt.ε∞) *
        Radon.regularyukawacoll(SingleLayer, ξ, elem, yukawa(f.opt)) -
        f.opt.εΩ/f.opt.ε∞ * Rjasanow.laplacecoll(SingleLayer, ξ, elem)
end

# Ks = -[1 - σ - Kʸ + εΩ/εΣ (Kʸ - K)]
struct KSfun{T} <: InteractionFunction{Tuple{Int, Vector{T}}, Tuple{Int, Triangle{T}}, T}
    opt::Option{T}
end
function (f::KSfun{T})(iξ::Tuple{Int, Vector{T}}, ielem::Tuple{Int, Triangle{T}}) where T
    -T(iξ[1] == ielem[1]) * 4π * (1-σ) + (1 - f.opt.εΩ/f.opt.εΣ) *
        Radon.regularyukawacoll(DoubleLayer, iξ[2], ielem[2], yukawa(f.opt)) +
        Rjasanow.laplacecoll(DoubleLayer, iξ[2], ielem[2])
end

# M₁₁ = 1 - σ - Kʸ
struct M11fun{T} <: InteractionFunction{Tuple{Int, Vector{T}}, Tuple{Int, Triangle{T}}, T}
    opt::Option{T}
end
function (f::M11fun{T})(iξ::Tuple{Int, Vector{T}}, ielem::Tuple{Int, Triangle{T}}) where T
    T(iξ[1] == ielem[1]) * 4π * (1-σ) -
        Radon.regularyukawacoll(DoubleLayer, iξ[2], ielem[2], yukawa(f.opt)) -
        Rjasanow.laplacecoll(DoubleLayer, iξ[2], ielem[2])
end

# M₁₂ = εΩ/ε∞ Vʸ -εΩ/εΣ (Vʸ - V)
struct M12fun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T}
    opt::Option{T}
end
function (f::M12fun{T})(ξ::Vector{T}, elem::Triangle{T}) where T
    f.opt.εΩ * (1/f.opt.ε∞ - 1/f.opt.εΣ) *
        Radon.regularyukawacoll(SingleLayer, ξ, elem, yukawa(f.opt)) +
        f.opt.εΩ / f.opt.ε∞ * Rjasanow.laplacecoll(SingleLayer, ξ, elem)
end

# M₁₃ = ε∞/εΣ (Kʸ - K)
struct M13fun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T}
    opt::Option{T}
end
function (f::M13fun{T})(ξ::Vector{T}, elem::Triangle{T}) where T
    f.opt.ε∞/f.opt.εΣ * Radon.regularyukawacoll(DoubleLayer, ξ, elem, yukawa(f.opt))
end

# M₂₁ = σ + K
struct M21fun{T} <: InteractionFunction{Tuple{Int, Vector{T}}, Tuple{Int, Triangle{T}}, T}
end
function (::M21fun{T})(iξ::Tuple{Int, Vector{T}}, ielem::Tuple{Int, Triangle{T}}) where T
    T(iξ[1] == ielem[1]) * 4π * σ  + Rjasanow.laplacecoll(DoubleLayer, iξ[2], ielem[2])
end

# M₂₂ = -V
struct M22fun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T} end
function (::M22fun{T})(ξ::Vector{T}, elem::Triangle{T}) where T
    -Rjasanow.laplacecoll(SingleLayer, ξ, elem)
end

# M₃₂ = εΩ/ε∞ V
struct M32fun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T}
    opt::Option{T}
end
function (f::M32fun{T})(ξ::Vector{T}, elem::Triangle{T}) where T
    f.opt.εΩ/f.opt.ε∞ * Rjasanow.laplacecoll(SingleLayer, ξ, elem)
end

# M₃₃ = 1 - σ - K
struct M33fun{T} <: InteractionFunction{Tuple{Int, Vector{T}}, Tuple{Int, Triangle{T}}, T}
end
function (::M33fun{T})(iξ::Tuple{Int, Vector{T}}, ielem::Tuple{Int, Triangle{T}}) where T
    T(iξ[1] == ielem[1]) * 4π * (1-σ) - Rjasanow.laplacecoll(DoubleLayer, iξ[2], ielem[2])
end
