# =========================================================================================
"""
    rfenergy(::Type{<: LocalityType}, ::BornIon{T})

Computes the local or nonlocal reaction field energy W* as
```math
W^* = ∫φ^* ρ \\quad dΩ
```
where ``φ^*`` is the reaction field and ``ρ`` is the corresponding charge distribution.

# Unit
``\\frac{kJ}{mol}``

# Return type
`T`
"""
function NESSie.rfenergy(lt::Type{<: LocalityType}, ion::BornIon{T}; kwargs...) where T
    ξ = ion.charge.pos
    ζ = ion.charge.val

    ζ * rfpotential(:Ω, lt, ξ, ion; kwargs...) * T(ec * 6.022140857e10 / 2)
end


# =========================================================================================
"""
    espotential(::Type{<: LocalityType}, ξ::Vector{T}, ion::BornIon{T})
    espotential(::Type{<: LocalityType}, Ξ::AbstractVector{Vector{T}}, ion::BornIon{T})

Computes the local or nonlocal electrostatic potential(s) at the given observation point(s)
ξ (Ξ) for the given born ion. This function automatically locates the observation point(s).

The electrostatic potential is computed as the sum of the corresponding
[reaction field potential](@ref rfpotential) and the [molecular potential](@ref molpotential).

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`

# Alias
    espotential(domain::Symbol, ::Type{<: LocalityType}, ξ::Vector{T}, ::BornIon{T})
    espotential(domain::Symbol, ::Type{<: LocalityType}, Ξ::AbstractVector{Vector{T}}, ::BornIon{T})

Computes the electrostatic potential(s) for the given observation point(s) ξ (Ξ) and the
given domain `:Ω`, `:Σ`, or `:Γ`.
"""
@inline function NESSie.espotential(
    domain::Symbol,
    lt::Type{<: LocalityType},
    ξorΞ::Union{Vector{T}, <: AbstractVector{Vector{T}}, <: Base.Generator},
    ion::BornIon{T};
    kwargs...
) where T
    domain === :Ω && return _espotential_Ω(lt, ξorΞ, ion; kwargs...)
    domain === :Γ && return _espotential_Ω(lt, ξorΞ, ion; kwargs...)
    domain === :Σ && return _espotential_Σ(lt, ξorΞ, ion; kwargs...)
    error("unknown domain $domain")
end

@inline function NESSie.espotential(
    lt::Type{<: LocalityType},
    ξ::Vector{T},
    ion::BornIon{T};
    kwargs...
) where T
    euclidean(ξ, ion.charge.pos) <= ion.radius ?
        espotential(:Ω, lt, ξ, ion; kwargs...) :
        espotential(:Σ, lt, ξ, ion; kwargs...)
end

@inline function NESSie.espotential(
    lt::Type{<: LocalityType},
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    ion::BornIon{T};
    kwargs...
) where T
    collect(T, espotential(lt, ξ, ion; kwargs...) for ξ in Ξ)
end


# =========================================================================================
"""
    molpotential(ξ::Vector{T}, ion::BornIon{T})
    molpotential(Ξ::AbstractVector{Vector{T}}, ion::BornIon{T})

Computes the molecular potential(s) at the given observation point(s) ξ (Ξ) for the given
born ion.

# Supported keyword arguments
 - `tolerance::T = 1e-10` minimum distance assumed between any observation point and point
   charge. Closer distances are replaced by this value.

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
@inline function NESSie.molpotential(
    ξorΞ::Union{Vector{T}, <: AbstractVector{Vector{T}}, <: Base.Generator},
    ion::BornIon{T};
    kwargs...
) where T
    potprefactor(T) .* _molpotential(ξorΞ, [ion.charge]; kwargs...)
end


# =========================================================================================
"""
    rfpotential(::Type{<: LocalityType}, ξ::Vector{T}, ion::BornIon{T})
    rfpotential(::Type{<: LocalityType}, Ξ::AbstractVector{Vector{T}}, ion::BornIon{T})

Computes the local or nonlocal reaction field potential(s) at the given observation point(s)
ξ (Ξ) for the given born ion. This function automatically locates the observation point(s).

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`

# Alias
    rfpotential(domain::Symbol, ::Type{<: LocalityType}, ξ::Vector{T}, ::BornIon{T})
    rfpotential(domain::Symbol, ::Type{<: LocalityType}, Ξ::AbstractVector{Vector{T}}, ::BornIon{T})

Computes the reaction field potential(s) for the given observation point(s) ξ (Ξ) and the
given domain `:Ω`, `:Σ`, or `:Γ`.
"""
@inline function NESSie.rfpotential(
    domain::Symbol,
    lt::Type{<: LocalityType},
    ξorΞ::Union{Vector{T}, <: AbstractVector{Vector{T}}, <: Base.Generator},
    ion::BornIon{T};
    kwargs...
) where T
    domain === :Ω && return _rfpotential_Ω(lt, ξorΞ, ion)
    domain === :Γ && return _rfpotential_Ω(lt, ξorΞ, ion)
    domain === :Σ && return _rfpotential_Σ(lt, ξorΞ, ion; kwargs...)
    error("unknown domain $domain")
end

@inline function NESSie.rfpotential(
    lt::Type{<: LocalityType},
    ξ::Vector{T},
    ion::BornIon{T};
    kwargs...
) where T
    euclidean(ξ, ion.charge.pos) <= ion.radius ?
        rfpotential(:Ω, lt, ξ, ion; kwargs...) :
        rfpotential(:Σ, lt, ξ, ion; kwargs...)
end

@inline function NESSie.rfpotential(
    lt::Type{<: LocalityType},
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    ion::BornIon{T};
    kwargs...
) where T
    collect(T, rfpotential(lt, ξ, ion; kwargs...) for ξ in Ξ)
end


# =========================================================================================
"""
    _espotential_Ω(::Type{<: LocalityType}, ξ::Vector{T}, ion::BornIon{T})
    _espotential_Ω(::Type{<: LocalityType}, Ξ::AbstractVector{Vector{T}}, ion::BornIon{T})

Computes the local or nonlocal electrostatic potential for (an) observation point(s) ξ (Ξ)
inside the Born sphere.

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
@inline function _espotential_Ω(
    lt::Type{<: LocalityType},
    ξorΞ::Union{Vector{T}, <: AbstractVector{Vector{T}}, <: Base.Generator},
    ion::BornIon{T};
    kwargs...
) where T
    _rfpotential_Ω(lt, ξorΞ, ion) .+ molpotential(ξorΞ, ion; kwargs...)
end


# =========================================================================================
"""
    _rfpotential_Ω(::Type{<: LocalityType}, ξ::Vector{T}, ion::BornIon{T})
    _rfpotential_Ω(::Type{<: LocalityType}, Ξ::AbstractVector{Vector{T}}, ion::BornIon{T})

Computes the local or nonlocal reaction field potential for (an) observation point(s) ξ (Ξ)
inside the Born sphere.

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
function _rfpotential_Ω(
    ::Type{LocalES},
    ξ::Vector{T},
    ion::BornIon{T}
) where T
    potprefactor(T) * ion.charge.val / ion.radius * (1/ion.params.εΣ - 1)
end

function _rfpotential_Ω(
    ::Type{NonlocalES},
    ξ::Vector{T},
    ion::BornIon{T}
) where T
    opt = ion.params
    ν = sqrt(opt.εΣ/opt.ε∞) * ion.radius / opt.λ
    potprefactor(T) * ion.charge.val / ion.radius / opt.εΣ *
        (1 - opt.εΣ + (opt.εΣ - opt.ε∞)/opt.ε∞ * sinh(ν)/ν * exp(-ν))
end

@inline function _rfpotential_Ω(
    lt::Type{<: LocalityType},
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    ion::BornIon{T}
) where T
    collect(T, _rfpotential_Ω(lt, ξ, ion) for ξ in Ξ)
end


# =========================================================================================
"""
    _espotential_Σ(::Type{<: LocalityType}, ξ::Vector{T}, ion::BornIon{T})
    _espotential_Σ(::Type{<: LocalityType}, Ξ::AbstractVector{Vector{T}}, ion::BornIon{T})

Computes the local or nonlocal electrostatic potential for (an) observation point(s) ξ (Ξ)
outside the Born sphere.

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
@inline function _espotential_Σ(
    ::Type{LocalES},
    ξ::Vector{T},
    ion::BornIon{T};
    kwargs...
) where T
    molpotential(ξ, ion; kwargs...) / ion.params.εΣ
end

function _espotential_Σ(
    ::Type{NonlocalES},
    ξ::Vector{T},
    ion::BornIon{T};
    tolerance::T = T(1e-10)
) where T
    opt = ion.params
    ν = sqrt(opt.εΣ/opt.ε∞) * ion.radius / opt.λ
    r = max(tolerance, euclidean(ξ, ion.charge.pos))
    c = (1 + (opt.εΣ - opt.ε∞)/opt.ε∞ * sinh(ν)/ν * exp(-ν * r / ion.radius)) / opt.εΣ
    molpotential(ξ, ion; tolerance = tolerance) * c
end

@inline function _espotential_Σ(
    lt::Type{<: LocalityType},
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    ion::BornIon{T};
    kwargs...
) where T
    collect(T, _espotential_Σ(lt, ξ, ion; kwargs...) for ξ in Ξ)
end


# =========================================================================================
"""
    _rfpotential_Σ(::Type{<: LocalityType}, ξ::Vector{T}, ion::BornIon{T})
    _rfpotential_Σ(::Type{<: LocalityType}, Ξ::AbstractVector{Vector{T}}, ion::BornIon{T})

Computes the local or nonlocal reaction field potential for (an) observation point(s) ξ (Ξ)
outside the Born sphere.

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
function _rfpotential_Σ(
    ::Type{LocalES},
    ξ::Vector{T},
    ion::BornIon{T};
    kwargs...
) where T
    molpotential(ξ, ion; kwargs...) * (1/ion.params.εΣ - 1)
end

function _rfpotential_Σ(
    ::Type{NonlocalES},
    ξ::Vector{T},
    ion::BornIon{T};
    tolerance::T = T(1e-10)
) where T
    opt = ion.params
    ν = sqrt(opt.εΣ/opt.ε∞) * ion.radius / opt.λ
    r = max(tolerance, euclidean(ξ, ion.charge.pos))
    c = (1 + (opt.εΣ - opt.ε∞)/opt.ε∞ * sinh(ν)/ν * exp(-ν * r / ion.radius)) / opt.εΣ
    molpotential(ξ, ion; tolerance = tolerance) * (c - 1)
end

@inline function _rfpotential_Σ(
    lt::Type{<: LocalityType},
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    ion::BornIon{T};
    kwargs...
) where T
    collect(T, _rfpotential_Σ(lt, ξ, ion; kwargs...) for ξ in Ξ)
end
