# =========================================================================================
"""
    struct LocalBEMResult{T, E} <: BEMResult{T, E}
        model::Model{T, E}
        u    ::Vector{T}   # [γ₀int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
        q    ::Vector{T}   # [γ₁int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
        umol ::Vector{T}   # [γ₀int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
        qmol ::Vector{T}   # [γ₁int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
    end

Result data of the local solving process to be used for potential computation and
post-processing, with `Ξ` being the list of observation points, that is, the set of
triangle centroids.
"""
struct LocalBEMResult{T, E} <: BEMResult{T, E}
    """Surface model"""
    model::Model{T, E}
    """[γ₀int(φ*)](ξ) for all observation points ξ"""
    u::Vector{T}
    """[γ₁int(φ*)](ξ) for all observation points ξ"""
    q::Vector{T}
    """[γ₀int(φ*mol)](ξ) for all observation points ξ"""
    umol::Vector{T}
    """[γ₁int(φ*mol)](ξ) for all observation points ξ"""
    qmol::Vector{T}
end


# =========================================================================================
"""
    solve{T, L <: LocalityType}(
                  ::L,
        model     ::Model{T, Triangle{T}}
    )

Computes the full local or nonlocal cauchy data on the surface of the biomolecule.

# Return type
`LocalBEMResult{T, Triangle{T}}` or `NonlocalBEMResult{T, Triangle{T}}`
"""
function solve(
                  ::Type{LocalES},
        model     ::Model{T, Triangle{T}}
    ) where T
    # observation points ξ
    Ξ = [e.center for e in model.elements]

    # compute molecular potentials for the point charges;
    # molecular potentials are initially premultiplied by 4π⋅ε0⋅εΩ
    umol = model.params.εΩ \   φmol(model)
    qmol = model.params.εΩ \ ∂ₙφmol(model)

    # convenience aliases
    εΩ = model.params.εΩ
    εΣ = model.params.εΣ
    numelem = length(model.elements)

    # system matrix M
    m = zeros(T, numelem, numelem)
    k = Array{T}(undef, numelem, numelem)
    v = zeros(T, numelem, numelem)

    # right-hand sides bᵤ and b_q
    b = zeros(T, numelem)

    # Mᵤ = (1 + εΩ/εΣ) ⋅ σ;
    # since all other components of the system matrix will be premultiplied by 4π,
    # do the same for σ here
    pluseye!(m, (1 + εΩ/εΣ) * T(4π * σ))

    #=
        generate and apply V
    =#
    # M_q = V
    Rjasanow.laplacecoll!(SingleLayer, v, model.elements, Ξ)

    # bᵤ -= εΩ/εΣ ⋅ V ⋅ qmol
    gemv!(-εΩ/εΣ, v, qmol, b)

    #=
        generate and apply K
    =#
    Rjasanow.laplacecoll!(DoubleLayer, k, model.elements, Ξ)

    # Mᵤ += (εΩ/εΣ - 1) ⋅ K
    axpy!(εΩ/εΣ - 1, k, m)

    # bᵤ = (K - σ) ⋅ umol
    # again, we apply a prefactor of 4π to σ to match the other components of the vector
    pluseye!(k, -T(4π * σ))
    gemv!(one(T), k, umol, b)

    #=
        u = b / M, with
        b = (K - σ) ⋅ umol - εΩ/εΣ ⋅ V ⋅ qmol, and
        M = (1 + εΩ/εΣ) ⋅ σ + (εΩ/εΣ - 1) ⋅ K.
    =#
    u = m \ b

    # b_q = (σ + K) ⋅ u
    fill!(b, zero(T))
    pluseye!(k, T(8π * σ)) # meh...
    gemv!(one(T), k, u, b)

    #=
        q = V^{-1}⋅(σ + K)u
    =#
    q = v \ b

    LocalBEMResult(model, u, q, umol, qmol)
end
