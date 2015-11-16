module NonLocalBEM

import Base: cos, sign
import Base.LinAlg.BLAS: gemv!, axpy!
import JSON: json

export
    # types.jl
    Element,
    Triangle,
    Tetrahedron,
    Charge,
    Option,
    SingleLayer,
    DoubleLayer,

    # util.jl
    props!,
    indexmap,
    unpack,
    vertexnormals,
    xml3d_mesh,
    eye!,
    isdegenerate,
    cos,
    cathetus,
    sign,
    distance,

    # input file readers
    readhmo,
    readmatlab,
    readoff,

    # this file
    defaultopt,
    singularpot,
    cauchy

# Include module files
include("types.jl")
include("util.jl")
include("input/hmo.jl")
include("input/matlab.jl")
include("input/off.jl")

# Default options
const defaultopt64 = Option(2., 78., 1.8, 20.)
const defaultopt32 = Option(2f0, 78f0, 1.8f0, 20f0)
defaultopt(::Type{Float64}) = defaultopt64
defaultopt(::Type{Float32}) = defaultopt32

# Include submodule files
include("radon.jl")
include("rjasanow.jl")

#=
    Computes the molecular potential (and the normal derivative) of the given system of
    point charges in a structureless medium.

    Note that the results are premultiplied by 4π!

    @param elements
                List of elements in the system
    @param charges
                List of charges in the system
    @param opt
                Constants to be used, including the dielectric constant of the solute
    @return (Vector{T}, Vector{T})
=#
function singularpot{T}(elements::Vector{Triangle{T}}, charges::Vector{Charge{T}}, opt::Option{T}=defaultopt(T))
    umol = T[]; qmol = T[]
    for elem in elements
        push!(umol, zero(T))
        push!(qmol, zero(T))
        for charge in charges
            r = elem.center - charge.pos
            rnorm = vecnorm(r)

            umol[end] += charge.val / rnorm
            qmol[end] -= charge.val * (r ⋅ elem.normal) / rnorm^3
        end
    end
    (opt.εΩ \ umol, opt.εΩ \ qmol)
end

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
    umol, qmol = singularpot(elements, charges, opt)

    # create right hand side
    rhs = zeros(T, 3 * numelem)

    # convenient access to the first block of rhs
    β = sub(rhs, 1:numelem)

    # initialize rhs
    copy!(β, umol)
    scale!(β, -2π)

    #=
        generate and apply Kʸ-K
    =#
    buffer = Array(T, numelem, numelem)
    Radon.regularyukawacoll!(DoubleLayer, buffer, elements, opt)

    # β += (1-εΩ/εΣ)(Kʸ-K)umol
    gemv!(1-εΩ/εΣ, buffer, umol, β)

    # m11 -= Kʸ-K
    axpy!(-1., buffer, m11)

    # m13 += ε∞/εΣ * (Kʸ-K)
    axpy!(ε∞/εΣ, buffer, m13)

    #=
        generate and apply Vʸ-V
    =#
    Radon.regularyukawacoll!(SingleLayer, buffer, elements, opt)

    # β += (εΩ/εΣ - εΩ/ε∞)(Vʸ-V)qmol
    gemv!(εΩ * (1/εΣ - 1/ε∞), buffer, qmol, β)

    # m12 += (εΩ/ε∞ - εΩ/εΣ)(Vʸ-V)
    axpy!(εΩ * (1/ε∞ - 1/εΣ), buffer, m12)

    #=
        generate and apply K
    =#
    LaplaceMod.laplacecoll!(DoubleLayer, buffer, elements, opt)

    # β += K
    gemv!(1., buffer, umol, β)

    # m11 -= K
    axpy!(-1., buffer, m11)

    # m21 += K
    axpy!(1., buffer, m21)

    # m33 -= K
    axpy!(-1., buffer, m33)

    #=
        generate and apply V
    =#
    LaplaceMod.laplacecoll!(SingleLayer, buffer, elements, opt)

    # β -= εΩ/ε∞ * V * qmol
    gemv!(-εΩ/ε∞, buffer, qmol, β)

    # m12 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m12)

    # m22 -= V
    axpy!(-1., buffer, m22)

    # m32 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m32)

    # solve system
    m\rhs
end

end # module
