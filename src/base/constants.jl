# =========================================================================================
@doc raw"""
    _etol(::Float64)
    _etol(::Type{Float64})
    _etol(::Float32)
    _etol(::Type{Float32})

Common tolerance values internally used by NESSie.

# Return type
`T`
""" _etol
const _etol_f64 = 1.45e-8
const _etol_f32 = 3.45f-4

@inline _etol(::Float64) = _etol_f64
@inline _etol(::Type{Float64}) = _etol_f64
@inline _etol(::Float32) = _etol_f32
@inline _etol(::Type{Float32}) = _etol_f32


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
@doc raw"""
    potprefactor(T::Type{Float64} = Float64)
    potprefactor(T::Type{Float32})

Common prefactor for all potentials ``φ_Ω`` and ``φ_Σ``:

```math
\frac{1.602 ⋅ 10^{-19}}{10^{-10} ⋅ 4π  ⋅ ε₀} ≈ 1.145 ⋅ 4π
```

# Return type
`T`
""" potprefactor
for T in [:Float64, :Float32]
    varname = Symbol("potprefactor_", T)
    @eval begin
        const $(varname) = $(T)(ec / 4π / ε0)
        @inline potprefactor(::Type{$(T)}) = $(varname)
    end
end
@inline potprefactor() = potprefactor(Float64)


# =========================================================================================
"""
    mutable struct Option{T <: AbstractFloat}
        εΩ::T       # dielectric constant of the solute
        εΣ::T       # dielectric constant of the solvent
        ε∞::T       # large-scale (bulk) solvent response
        λ ::T       # correlation length scale [λ] = Å
    end

System parameters
"""
@auto_hash_equals mutable struct Option{T <: AbstractFloat}
    "dielectric constant of the solute"
    εΩ::T
    "dielectric constant of the solvent"
    εΣ::T
    "large-scale (bulk) solvent response"
    ε∞::T
    "correlation length scale [λ] = Å"
    λ::T
end

@inline function Base.show(io::IO, ::MIME"text/plain", opt::Option)
    show(io, opt)
end

@inline function Base.show(io::IO, opt::Option)
    print(io,
        "$(typeof(opt))",
        "(εΩ = $(opt.εΩ), εΣ = $(opt.εΣ), ε∞ = $(opt.ε∞), λ = $(opt.λ))"
    )
end


# =========================================================================================
@doc raw"""
    defaultopt(T::Type{Float64} = Float64)
    defaultopt(T::Type{Float32})

Default system parameters for proteins in water.

# Return type
[`Option{T}`](@ref Option)

# Default values
 * ``ε_Ω = 2``
 * ``ε_Σ = 78``
 * ``ε_∞ = 1.8``
 * ``λ = 20 Å``
""" defaultopt
for T in [:Float64, :Float32]
    varname = Symbol("defaultopt_", T)
    @eval begin
        const $(varname) = Option($(T)[2, 78, 1.8, 20]...)
        @inline defaultopt(::Type{$(T)}) = $(varname)
    end
end
@inline defaultopt() = defaultopt(Float64)


# =========================================================================================
"""
    yukawa(opt::Option{T})

Exponent ``1/Λ`` for the fundamental solution of the yukawa operator

```math
Λ := λ\\sqrt{\\frac{ε_∞}{ε_Σ}}
```

# Return type
`T`
"""
@inline yukawa(opt::Option) = √(opt.εΣ/opt.ε∞)/opt.λ
