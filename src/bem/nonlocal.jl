#=
    Computes the full cauchy data on the surface of the biomolecule.

    Note that the result is premultiplied by 4π!

    @param elements
                List of surface elements
    @param charges
                List of charges in the biomolecule
    @param LaplaceMod
                Module to be used for Laplace potential; Valid values: Radon, Rjasanow
    @param opt
                Constants to be used
    @return Vector{T}
=#
function cauchy{T}(elements::Vector{Triangle{T}}, charges::Vector{Charge{T}}, LaplaceMod::Module=Rjasanow, opt::Option{T}=defaultopt(T))
    # convient access to constants
    const εΩ  = opt.εΩ
    const εΣ  = opt.εΣ
    const ε∞  = opt.ε∞
    const yuk = yukawa(opt)

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

    # initialize the system matrix
    eye!(m11, 2π)
    eye!(m21, 2π)
    eye!(m33, 2π)

    # compute molecular potential for the point charges
    umol = εΩ \ φmol(elements, charges)
    qmol = εΩ \ ∂ₙφmol(elements, charges)

    # create right hand side
    rhs = zeros(T, 3 * numelem)

    # convenient access to the first block of rhs
    β = view(rhs, 1:numelem)

    # initialize rhs
    copy!(β, umol)
    scale!(β, -2π)

    # create list of observation points
    ξlist = [e.center for e in elements]

    #=
        generate and apply Kʸ-K
    =#
    buffer = Array(T, numelem, numelem)
    Radon.regularyukawacoll!(DoubleLayer, buffer, elements, ξlist, yuk)

    # β += (1-εΩ/εΣ)(Kʸ-K)umol
    gemv!(1-εΩ/εΣ, buffer, umol, β)

    # m11 -= Kʸ-K
    axpy!(-1, buffer, m11)

    # m13 += ε∞/εΣ * (Kʸ-K)
    axpy!(ε∞/εΣ, buffer, m13)

    #=
        generate and apply Vʸ-V
    =#
    Radon.regularyukawacoll!(SingleLayer, buffer, elements, ξlist, yuk)

    # β += (εΩ/εΣ - εΩ/ε∞)(Vʸ-V)qmol
    gemv!(εΩ * (1/εΣ - 1/ε∞), buffer, qmol, β)

    # m12 += (εΩ/ε∞ - εΩ/εΣ)(Vʸ-V)
    axpy!(εΩ * (1/ε∞ - 1/εΣ), buffer, m12)

    #=
        generate and apply K
    =#
    LaplaceMod.laplacecoll!(DoubleLayer, buffer, elements, ξlist)

    # β += K
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
    LaplaceMod.laplacecoll!(SingleLayer, buffer, elements, ξlist)

    # β -= εΩ/ε∞ * V * qmol
    gemv!(-εΩ/ε∞, buffer, qmol, β)

    # m12 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m12)

    # m22 -= V
    axpy!(-1, buffer, m22)

    # m32 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m32)

    # solve system
    m\rhs
end

function φΩ{T}(nodes::Vector{Vector{T}}, elements::Vector{Triangle{T}}, charges::Vector{Charge{T}}, LaplaceMod::Module=Rjasanow, opt::Option{T}=defaultopt(T))
    c = cauchy(elements, charges, LaplaceMod, opt)
    numelem = length(elements)
    φ = zeros(T, length(nodes))

    # convenient access
    γ0intφstar = c[1:numelem]
    γ1intφstar = c[1+numelem:2numelem]

    LaplaceMod.laplacecoll!(DoubleLayer, φ, elements, nodes, γ0intφstar)
    scale!(φ, -1)
    LaplaceMod.laplacecoll!(SingleLayer, φ, elements, nodes, γ1intφstar)
    scale!(φ, 2)
    T(1.69e-9 / 4π / ε0) * (4π \ φ + opt.εΩ \ φmol(nodes, charges))
end
