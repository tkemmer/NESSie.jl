
#=
    TODO
=#
type LocalBEMResult{T} <: BEMResult{T}
    u::Vector{T}
    q::Vector{T}
    umol::Vector{T}
    qmol::Vector{T}
end

#=
    TODO
=#
function solvelocal{T}(
        model::SurfaceModel{T},
        LaplaceMod::Module=Rjasanow,
        opt::Option{T}=defaultopt(T)
    )
    # convenience aliases
    const elements = model.elements
    const charges  = model.charges

    # observation points ξ
    const ξlist = [e.center for e in elements]

    # compute molecular potentials for the point charges
    const umol = opt.εΩ \   φmol(elements, charges)
    const qmol = opt.εΩ \ ∂ₙφmol(elements, charges)

    const u = solve_u(elements, umol, qmol, ξlist, LaplaceMod, opt)
    const q = solve_q(elements, u, ξlist, LaplaceMod, opt)

    LocalBEMResult(u, q, umol, qmol)
end


#=
    TODO
=#
function solve_u{T}(
        elements::Vector{Triangle{T}},
        umol::Vector{T},
        qmol::Vector{T},
        ξlist::Vector{Vector{T}},
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
    LaplaceMod.laplacecoll!(DoubleLayer, buf, elements, ξlist)

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
    LaplaceMod.laplacecoll!(SingleLayer, buf, elements, ξlist)

    # rhs -= εΩ/εΣ ⋅ V ⋅ qmol
    gemv!(εΩ/εΣ, buf, qmol, rhs)

    # solve system
    m\rhs
end


#=
    q = V^{-1}⋅(σ + K)u
=#
function solve_q{T}(
        elements::Vector{Triangle{T}},
        u::Vector{T},
        ξlist::Vector{Vector{T}},
        LaplaceMod::Module=Rjasanow,
        opt::Option{T}=defaultopt(T)
    )
    # constants
    const numelem = length(elements)

    # system matrix
    m = zeros(T, numelem, numelem)

    # m = (σ + K)
    LaplaceMod.laplacecoll!(DoubleLayer, m, elements, ξlist)
    pluseye!(m, 4π * σ)

    # m = V^{-1}⋅(σ + K)
    v = Array(T, numelem, numelem)
    LaplaceMod.laplacecoll!(SingleLayer, v, elements, ξlist)
    m = gemm(one(T), pinv(v), m)     # TODO check tolerance

    # q = m⋅u
    gemv(one(T), m, u)
end
