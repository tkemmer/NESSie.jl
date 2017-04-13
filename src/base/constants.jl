# =========================================================================================
"""
Vacuum permittivity

# Unit
``\\frac{F}{m}``
"""
const ε0 = 1 / (4π * 1e-7 * 299792458^2)


# =========================================================================================
"""
Geometric quantity ``σ(ξ)``

```math
σ(ξ) = \\lim_{ε→0} \\frac{1}{4πε²} ∫_{r ∈ Ω:|r-ξ|=ε} dΓᵣ = \\frac{1}{2}
```
for almost all ``ξ ∈ Γ`` (cf. [[Ste03]](@ref Bibliography)).
"""
const σ  = 0.5


# =========================================================================================
"""
``10^{10}`` times the elementary charge (for ``Å → m`` conversion)

# Unit
``C``
"""
const ec = 1.602176e-9


# =========================================================================================
@doc """
    potprefactor(T::Type{Float64})
    potprefactor(T::Type{Float32})

Common prefactor for all potentials ``φ\_Ω`` and ``φ\_Σ``:

```math
\\frac{1.602 ⋅ 10^{-19}}{10^{-10} ⋅ 4π  ⋅ ε₀} ≈ 1.145 ⋅ 4π
```

# Return type
`T`
""" potprefactor
for T in [:Float64, :Float32]
    varname = Symbol("potprefactor_", T)
    @eval begin
        const $(varname) = $(T)(ec / 4π / ε0)
        potprefactor(::Type{$(T)}) = $(varname)
    end
end


# =========================================================================================
"""
    immutable Option{T <: AbstractFloat}
        εΩ::T       # dielectric constant of the solute
        εΣ::T       # dielectric constant of the solvent
        ε∞::T       # large-scale (bulk) solvent response
        λ ::T       # correlation length scale [λ] = Å
    end

System parameters
"""
immutable Option{T <: AbstractFloat}
    "dielectric constant of the solute"
    εΩ::T
    "dielectric constant of the solvent"
    εΣ::T
    "large-scale (bulk) solvent response"
    ε∞::T
    "correlation length scale [λ] = Å"
    λ::T
end


# =========================================================================================
@doc """
    defaultopt(T::Type{Float64})
    defaultopt(T::Type{Float32})

Default system parameters

# Return type
[`Option{T}`](@ref ProteinES.Option)

# Default values
 * ``ε\_Ω = 2``
 * ``ε\_Σ = 78``
 * ``ε\_∞ = 1.8``
 * ``λ = 20 Å``
""" defaultopt
for T in [:Float64, :Float32]
    varname = Symbol("defaultopt_", T)
    @eval begin
        const $(varname) = Option($(T)[2, 78, 1.8, 20]...)
        defaultopt(::Type{$(T)}) = $(varname)
    end
end


# =========================================================================================
"""
    yukawa{T}(opt::Option{T})

Exponent ``-Λ^{-1}`` for the fundamental solution of the yukawa operator

```math
Λ := λ\\sqrt{\\frac{ε\_∞}{ε\_Σ}}
```

# Return type
`T`
"""
yukawa{T}(opt::Option{T}) = √(opt.εΣ/opt.ε∞)/opt.λ
