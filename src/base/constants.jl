#=
    Global constants

    References:
    [1] O. Steinbach, Numerische Näherungsverfahren für elliptische Randwertprobleme -
        Finite Elemente und Randelemente (in German). Advances in Numerical Matheamtics.
        Teubner Verlag/GWV Fachverlage GmbH, Wiesbaden, 2003.
=#
const ε0 = 1 / (4π * 1e-7 * 299792458^2)    # vacuum permittivity [ε0] = F/m
const σ  = 0.5                              # σ(ξ) = lim_{ε→0} 1/4πε² ∫dΓᵣ = 1/2 for almost all ξ ∈ Γ [1]

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
for T in [:Float64, :Float32]
    varname = Symbol("defaultopt_", T)
    @eval begin
        const $(varname) = Option($(T)[2, 78, 1.8, 20]...)
        defaultopt(::Type{$(T)}) = $(varname)
    end
end

# exponent for fundamental solution of yukawa operator -1/Λ = -1/(λ√(ε∞/εΣ))
yukawa{T}(opt::Option{T}) = √(opt.εΣ/opt.ε∞)/opt.λ
