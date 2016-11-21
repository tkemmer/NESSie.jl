#=
    Result data of the local solving process to be used for potential computation and post-processing:
    ▶ u:    [γ0int(φ*)](ξ) ∀ ξ ∈ Ξ
    ▶ q:    [γ1int(φ*)](ξ) ∀ ξ ∈ Ξ
    ▶ umol: [γ0int(φ*mol)](ξ) ∀ ξ ∈ Ξ
    ▶ qmol: [γ1int(φ*mol)](ξ) ∀ ξ ∈ Ξ
    with Ξ being the list of observation points, that is, the set of triangle centroids.
=#
type LocalBEMResult{T} <: BEMResult{T}
    u::Vector{T}
    q::Vector{T}
    umol::Vector{T}
    qmol::Vector{T}
end

#=
    Computes all vectors required for potential computation.

    Note that the result is premultiplied by 4π!

    @param model
            Surface model
    @param LaplaceMod
            Module to be used for Laplace potential; Valid values: Radon, Rjasanow
    @param opt
            Constants to be used
    @return LocalBEMResult{T}
=#
function solvelocal{T}(
        model::SurfaceModel{T},
        LaplaceMod::Module=Rjasanow,
        opt::Option{T}=defaultopt(T)
    )
    # observation points ξ
    const Ξ = [e.center for e in model.elements]

    # compute molecular potentials for the point charges
    const umol = opt.εΩ \   φmol(model)
    const qmol = opt.εΩ \ ∂ₙφmol(model)

    const u = solve_u(model.elements, umol, qmol, Ξ, LaplaceMod, opt)
    const q = solve_q(model.elements, u, Ξ, LaplaceMod, opt)

    LocalBEMResult(u, q, umol, qmol)
end

#=
    Computes u.

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

    # m = (1 + εΩ/εΣ) ⋅ σ
    pluseye!(m, (1 + εΩ/εΣ) * 4π * σ)

    # generate K
    buf = Array(T, numelem, numelem)
    LaplaceMod.laplacecoll!(DoubleLayer, buf, elements, Ξ)

    # m += (εΩ/εΣ - 1) ⋅ K
    axpy!(εΩ/εΣ - 1, buf, m)

    #=
        right-hand side
    =#
    rhs = zeros(T, numelem)

    # rhs = -σ ⋅ umol
    copy!(rhs, umol)
    scale!(rhs, -4π * σ)

    # rhs += K ⋅ umol
    gemv!(one(T), buf, umol, rhs)

    # generate V
    LaplaceMod.laplacecoll!(SingleLayer, buf, elements, Ξ)

    # rhs -= εΩ/εΣ ⋅ V ⋅ qmol
    gemv!(εΩ/εΣ, buf, qmol, rhs)

    # solve system
    m\rhs
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
