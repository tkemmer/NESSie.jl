#=
    Result data of the local solving process to be used for potential computation and post-processing:
    ▶ u:    [γ0int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
    ▶ q:    [γ1int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
    ▶ umol: [γ0int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
    ▶ qmol: [γ1int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
    with Ξ being the list of observation points, that is, the set of triangle centroids.
=#
type LocalBEMResult{T} <: BEMResult{T}
    model::SurfaceModel{T}
    opt::Option{T}
    u::Vector{T}
    q::Vector{T}
    umol::Vector{T}
    qmol::Vector{T}
end

#=
    Computes the full cauchy data on the surface of the biomolecule.

    See `LocalBEMResult` for remarks on the present prefactors.

    @param model
            Surface model
    @param LaplaceMod
            Module to be used for Laplace potential; Valid values: Radon, Rjasanow
    @param opt
            Constants to be used
    @return LocalBEMResult{T}
=#
function solve{T}(
        ::Type{LocalES},
        model::SurfaceModel{T},
        LaplaceMod::Module=Rjasanow,
        opt::Option{T}=defaultopt(T)
    )
    # observation points ξ
    const Ξ = [e.center for e in model.elements]

    # compute molecular potentials for the point charges;
    # molecular potentials are initially premultiplied by 4π⋅ε0⋅εΩ
    const umol = opt.εΩ \   φmol(model)
    const qmol = opt.εΩ \ ∂ₙφmol(model)

    # convenience aliases
    const εΩ = opt.εΩ
    const εΣ = opt.εΣ
    const numelem = length(model.elements)

    # system matrix M
    m = zeros(T, numelem, numelem)
    k = Array(T, numelem, numelem)
    v = zeros(T, numelem, numelem)

    # right-hand sides bᵤ and b_q
    b = zeros(T, numelem)

    # Mᵤ = (1 + εΩ/εΣ) ⋅ σ;
    # since all other components of the system matrix will be premultiplied by 4π, do the same for σ here
    pluseye!(m, (1 + εΩ/εΣ) * 4π * σ)

    #=
        generate and apply V
    =#
    # M_q = V
    LaplaceMod.laplacecoll!(SingleLayer, v, model.elements, Ξ)

    # bᵤ -= εΩ/εΣ ⋅ V ⋅ qmol
    gemv!(-εΩ/εΣ, v, qmol, b)

    #=
        generate and apply K
    =#
    LaplaceMod.laplacecoll!(DoubleLayer, k, model.elements, Ξ)

    # Mᵤ += (εΩ/εΣ - 1) ⋅ K
    axpy!(εΩ/εΣ - 1, k, m)

    # bᵤ = (K - σ) ⋅ umol
    # again, we apply a prefactor of 4π to σ to match the other components of the vector
    pluseye!(k, -4π * σ)
    gemv!(one(T), k, umol, b)

    #=
        u = b / M, with
        b = (K - σ) ⋅ umol - εΩ/εΣ ⋅ V ⋅ qmol, and
        M = (1 + εΩ/εΣ) ⋅ σ + (εΩ/εΣ - 1) ⋅ K.
    =#
    u = m \ b

    # b_q = (σ + K) ⋅ u
    fill!(b, zero(T))
    pluseye!(k, 8π * σ) # meh...
    gemv!(one(T), k, u, b)

    #=
        q = V^{-1}⋅(σ + K)u
    =#
    q = v \ b

    LocalBEMResult(model, opt, u, q, umol, qmol)
end
