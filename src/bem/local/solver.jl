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
    Computes all vectors required for potential computation.

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

    const u = solve_u(model.elements, umol, qmol, Ξ, LaplaceMod, opt)
    const q = solve_q(model.elements, u, Ξ, LaplaceMod, opt)

    LocalBEMResult(model, opt, u, q, umol, qmol)
end

#=
    Computes u = b / M, with

        b = (K - σ) ⋅ umol - εΩ/εΣ ⋅ V ⋅ qmol, and

        M = (1 + εΩ/εΣ) ⋅ σ + (εΩ/εΣ - 1) ⋅ K.

    @param elements
            List of surface triangles
    @param umol
            Molecular potential on the surface
    @param qmol
            Normal derivative of the molecular potential on the surface
    @param Ξ
            List of observation points
    @param LaplaceMod
            Module to be used for Laplace potential; Valid values: Radon, Rjasanow
    @param opt
            Constants to be used
    @return Vector{T}
=#
function solve_u{T}(
        elements::Vector{Triangle{T}},
        umol::Vector{T},
        qmol::Vector{T},
        Ξ::Vector{Vector{T}},
        LaplaceMod::Module=Rjasanow,
        opt::Option{T}=defaultopt(T),
    )
    # convenient access to constants
    const εΩ = opt.εΩ
    const εΣ = opt.εΣ
    const numelem = length(elements)

    #=
        system matrix
    =#
    m = zeros(T, numelem, numelem)

    # M = (1 + εΩ/εΣ) ⋅ σ;
    # since all other components of the system matrix will be premultiplied by 4π, do the same for σ here
    pluseye!(m, (1 + εΩ/εΣ) * 4π * σ)

    # generate K
    buf = Array(T, numelem, numelem)
    LaplaceMod.laplacecoll!(DoubleLayer, buf, elements, Ξ)

    # M += (εΩ/εΣ - 1) ⋅ K
    axpy!(εΩ/εΣ - 1, buf, m)

    #=
        right-hand side
    =#
    b = zeros(T, numelem)

    # b = -σ ⋅ umol;
    # again, we apply a prefactor of 4π to σ to match the other components of the vector
    copy!(b, umol)
    scale!(b, -4π * σ)

    # b += K ⋅ umol
    gemv!(one(T), buf, umol, b)

    # generate V
    LaplaceMod.laplacecoll!(SingleLayer, buf, elements, Ξ)

    # b -= εΩ/εΣ ⋅ V ⋅ qmol
    gemv!(εΩ/εΣ, buf, qmol, b)

    # u = b / M
    m \ b
end

#=
    Computes q = V^{-1}⋅(σ + K)u.

    @param elements
            List of surface triangles
    @param u
            Vector u
    @param Ξ
            List of observation points
    @param LaplaceMod
            Module to be used for Laplace potential; Valid values: Radon, Rjasanow
    @param opt
            Constants to be used
    @return Vector{T}
=#
function solve_q{T}(
        elements::Vector{Triangle{T}},
        u::Vector{T},
        Ξ::Vector{Vector{T}},
        LaplaceMod::Module=Rjasanow,
        opt::Option{T}=defaultopt(T)
    )
    # constants
    const numelem = length(elements)

    # system matrix
    m = zeros(T, numelem, numelem)

    # m = (σ + K)
    LaplaceMod.laplacecoll!(DoubleLayer, m, elements, Ξ)
    pluseye!(m, 4π * σ)

    # m = V^{-1}⋅(σ + K)
    v = Array(T, numelem, numelem)
    LaplaceMod.laplacecoll!(SingleLayer, v, elements, Ξ)
    m = gemm(one(T), pinv(v), m)     # TODO check tolerance

    # q = m⋅u
    gemv(one(T), m, u)
end
