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
    const εΩ = opt.εΩ
    const εΣ = opt.εΣ
    const ε∞ = opt.ε∞

    # create system matrix
    numelem = length(elements)
    m = zeros(T, 3 * numelem, 3 * numelem)

    # convenient access to 9 blocks of the system matrix
    m11 = sub(m,          1:numelem,           1:numelem )
    m12 = sub(m,          1:numelem,   1+numelem:2numelem)
    m13 = sub(m,          1:numelem,  1+2numelem:3numelem)
    m21 = sub(m,  1+numelem:2numelem,          1:numelem )
    m22 = sub(m,  1+numelem:2numelem,  1+numelem:2numelem)
    m23 = sub(m,  1+numelem:2numelem, 1+2numelem:3numelem)
    m31 = sub(m, 1+2numelem:3numelem,          1:numelem )
    m32 = sub(m, 1+2numelem:3numelem,  1+numelem:2numelem)
    m33 = sub(m, 1+2numelem:3numelem, 1+2numelem:3numelem)

    # initialize the system matrix
    eye!(m11, 2π)
    eye!(m21, 2π)
    eye!(m33, 2π)

    # compute molecular potential for the point charges
    umol = opt.εΩ \ φmol(elements, charges)
    qmol = opt.εΩ \ ∂ₙφmol(elements, charges)

    # create right hand side
    rhs = zeros(T, 3 * numelem)

    # convenient access to the first block of rhs
    β = sub(rhs, 1:numelem)

    # initialize rhs
    copy!(β, umol)
    scale!(β, -2π)

    # create list of observation points
    ξlist = [e.center for e in elements]

    #=
        generate and apply Kʸ-K
    =#
    buffer = Array(T, numelem, numelem)
    Radon.regularyukawacoll!(DoubleLayer, buffer, elements, ξlist, opt)

    # β += (1-εΩ/εΣ)(Kʸ-K)umol
    gemv!(1-εΩ/εΣ, buffer, umol, β)

    # m11 -= Kʸ-K
    axpy!(-one(T), buffer, m11)

    # m13 += ε∞/εΣ * (Kʸ-K)
    axpy!(ε∞/εΣ, buffer, m13)

    #=
        generate and apply Vʸ-V
    =#
    Radon.regularyukawacoll!(SingleLayer, buffer, elements, ξlist, opt)

    # β += (εΩ/εΣ - εΩ/ε∞)(Vʸ-V)qmol
    gemv!(εΩ * (1/εΣ - 1/ε∞), buffer, qmol, β)

    # m12 += (εΩ/ε∞ - εΩ/εΣ)(Vʸ-V)
    axpy!(εΩ * (1/ε∞ - 1/εΣ), buffer, m12)

    #=
        generate and apply K
    =#
    LaplaceMod.laplacecoll!(DoubleLayer, buffer, elements, ξlist, opt)

    # β += K
    gemv!(one(T), buffer, umol, β)

    # m11 -= K
    axpy!(-one(T), buffer, m11)

    # m21 += K
    axpy!(one(T), buffer, m21)

    # m33 -= K
    axpy!(-one(T), buffer, m33)

    #=
        generate and apply V
    =#
    LaplaceMod.laplacecoll!(SingleLayer, buffer, elements, ξlist, opt)

    # β -= εΩ/ε∞ * V * qmol
    gemv!(-εΩ/ε∞, buffer, qmol, β)

    # m12 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m12)

    # m22 -= V
    axpy!(-one(T), buffer, m22)

    # m32 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m32)

    # solve system
    m\rhs
end