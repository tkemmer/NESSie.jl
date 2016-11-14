#=
    TODO
=#
function solve_u{T}(
        elements::Vector{Triangle{T}},
        charges::Vector{Charge{T}},
        LaplaceMod::Module=Rjasanow,
        opt::Option{T}=defaultopt(T)
    )
    # convenient access to constants
    const εΩ = opt.εΩ
    const εΣ = opt.εΣ
    const numelem = length(elements)

    # compute molecular potentials for the point charges
    const umol = εΩ \ φmol(elements, charges)
    const qmol = εΩ \ ∂ₙφmol(elements, charges)

    # observation points
    const ξlist = [e.center for e in elements]

    #=
        system matrix
    =#
    m = zeros(T, numelem, numelem)

    # m = (1 + εΩ/εΣ) ⋅ σ
    eye!(m, (1 + εΩ/εΣ) * 4π * σ)

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
