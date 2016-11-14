#=
    Global constants
=#
const ε0 = 1 / (4π * 1e-7 * 299792458^2) # vacuum permittivity [F/m]

#=
    Model parameters
=#
immutable Option{T <: AbstractFloat}
    εΩ::T       # dielectric constant of the solute []
    εΣ::T       # dielectric constant of the solvent []
    ε∞::T       # large-scale (bulk) solvent response []
    λ::T        # correlation length scale [Å]
end

# Default options
const defaultopt64 = Option(2., 78., 1.8, 20.)
const defaultopt32 = Option(2f0, 78f0, 1.8f0, 20f0)
defaultopt(::Type{Float64}) = defaultopt64
defaultopt(::Type{Float32}) = defaultopt32

# exponent for fundamental solution of yukawa operator -1/Λ = -1/(λ√(ε∞/εΣ))
yukawa{T}(opt::Option{T}) = √(opt.εΣ/opt.ε∞)/opt.λ
