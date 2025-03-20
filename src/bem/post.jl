# =========================================================================================
"""
    rfenergy(::BEMResult{T})

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
function NESSie.rfenergy(bem::R) where {T, R <: BEMResult{T}}
    qposs = [charge.pos for charge in bem.model.charges]
    qvals = [charge.val for charge in bem.model.charges]

    # reaction field energy per point charge
    wstar = zeros(T, length(bem.model.charges))

    # W* = -[W ⋅ u](ξ)
    Rjasanow.laplacecoll!(DoubleLayer, wstar, bem.model.elements, qposs, bem.u)
    rmul!(wstar, -1)

    # W* += [V ⋅ q](ξ)
    Rjasanow.laplacecoll!(SingleLayer, wstar, bem.model.elements, qposs, bem.q)

    # Apply ρ, integrate over Ω and apply remaining prefactors (in order)
    # ▶ 4π        for Vtilde, W
    # (in potprefactor:)
    # ▶ 4π⋅ε0     for u and q
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    #
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 6.022e23  for Avogadro constant Nₐ; [Nₐ] = 1/mol
    # ▶ 1e-3      for the conversion 1/J → 1/kJ
    # ▶ 0.5       since we have counted all interactions twice
    wstar ⋅ qvals / T(4π) * potprefactor(T) * T(ec * 6.022140857e10 / 2)
end


# =========================================================================================
"""
    espotential(ξ::Vector{T}, bem::BEMResult{T})
    espotential(Ξ::AbstractVector{Vector{T}}, bem::BEMResult{T})

Computes the local or nonlocal electrostatic potential(s) at the given observation point(s)
ξ (Ξ) for the given BEM result. This function tries to automatically locate the observation
point(s) using [`guess_domain`](@ref).

The electrostatic potential is computed as the sum of the corresponding
[reaction field potential](@ref rfpotential) and the [molecular potential](@ref molpotential).

# Supported keyword arguments
 - `surface_margin::T = 1e-6` see [`guess_domain`](@ref)
 - `tolerance::T = 1e-10` see [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`

# Alias
    espotential(domain::Symbol, ξ::Vector{T}, bem::BEMResult{T})
    espotential(domain::Symbol, Ξ::AbstractVector{T}, bem::BEMResult{T})

Computes the electrostatic potential(s) for the given observation point(s) ξ (Ξ) and the
given domain `:Ω`, `:Σ`, or `:Γ`.
"""
function NESSie.espotential(
    domain::Symbol,
    ξorΞ::Union{Vector{T}, <: AbstractVector{Vector{T}}, <: Base.Generator},
    bem::BEMResult{T};
    tolerance::T = T(1e-10)
) where T
    domain === :Ω && return _espotential_Ω(ξorΞ, bem; tolerance = tolerance)
    domain === :Σ && return _espotential_Σ(ξorΞ, bem)
    domain === :Γ && return _espotential_Γ(ξorΞ, bem; tolerance = tolerance)
    error("unknown domain $domain")
end

@inline function NESSie.espotential(
    ξ::Vector{T},
    bem::BEMResult{T};
    surface_margin::T = T(1e-6),
    kwargs...
) where T
    domain = guess_domain(ξ, bem.model; surface_margin = surface_margin)
    espotential(domain, ξ, bem; kwargs...)
end

function NESSie.espotential(
    Ξ::AbstractVector{Vector{T}},
    bem::BEMResult{T};
    surface_margin::T = T(1e-6),
    kwargs...
) where T
    domains = guess_domain.(Ξ, Ref(bem.model); surface_margin = surface_margin)
    unknown_domains = setdiff(domains, [:Ω, :Σ, :Γ])
    !isempty(unknown_domains) && error("unknown domains $unknown_domains")

    ret = Array{T}(undef, length(Ξ))
    view(ret, domains .== :Ω) .= espotential(:Ω, view(Ξ, domains .== :Ω), bem; kwargs...)
    view(ret, domains .== :Σ) .= espotential(:Σ, view(Ξ, domains .== :Σ), bem; kwargs...)
    view(ret, domains .== :Γ) .= espotential(:Γ, view(Ξ, domains .== :Γ), bem; kwargs...)
    ret
end

@inline function NESSie.espotential(
    Ξ::Base.Generator,
    bem::BEMResult{T};
    kwargs...
) where T
    espotential(collect(Vector{T}, Ξ), bem; kwargs...)
end


# =========================================================================================
"""
    molpotential(ξ::Vector{T}, bem::BEMResult{T})
    molpotential(Ξ::AbstractVector{Vector{T}}, bem::BEMResult{T})

Computes the molecular potential(s) at the given observation point(s) ξ (Ξ) for the given
BEM result.

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
    bem::BEMResult{T};
    kwargs...
) where T
    molpotential(ξorΞ, bem.model; kwargs...)
end


# =========================================================================================
"""
    rfpotential(ξ::Vector{T}, bem::BEMResult{T})
    rfpotential(Ξ::AbstractVector{Vector{T}}, bem::BEMResult{T})

Computes the local or nonlocal reaction field potential(s) at the given observation point(s)
ξ (Ξ) for the given BEM result. This function tries to automatically locate the observation
point(s) using [`guess_domain`](@ref).

# Supported keyword arguments
 - `surface_margin::T = 1e-6` see [`guess_domain`](@ref)
 - `tolerance::T = 1e-10` see [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`

# Alias
    rfpotential(domain::Symbol, ξ::Vector{T}, bem::BEMResult{T})
    rfpotential(domain::Symbol, Ξ::AbstractVector{T}, bem::BEMResult{T})

Computes the reaction field potential(s) for the given observation point(s) ξ (Ξ) and the
given domain `:Ω`, `:Σ`, or `:Γ`.
"""
@inline function NESSie.rfpotential(
    domain::Symbol,
    ξorΞ::Union{Vector{T}, <: AbstractVector{Vector{T}}, <: Base.Generator},
    bem::BEMResult{T};
    tolerance::T = T(1e-10)
) where T
    domain === :Ω && return _rfpotential_Ω(ξorΞ, bem)
    domain === :Σ && return _rfpotential_Σ(ξorΞ, bem; tolerance = tolerance)
    domain === :Γ && return _rfpotential_Γ(ξorΞ, bem)
    error("unknown domain $domain")
end

@inline function NESSie.rfpotential(
    ξ::Vector{T},
    bem::BEMResult{T};
    surface_margin::T = T(1e-6),
    kwargs...
) where T
    domain = guess_domain(ξ, bem.model; surface_margin = surface_margin)
    rfpotential(domain, ξ, bem; kwargs...)
end

function NESSie.rfpotential(
    Ξ::AbstractVector{Vector{T}},
    bem::BEMResult{T};
    surface_margin::T = T(1e-6),
    kwargs...
) where T
    domains = guess_domain.(Ξ, Ref(bem.model); surface_margin = surface_margin)
    unknown_domains = setdiff(domains, [:Ω, :Σ, :Γ])
    !isempty(unknown_domains) && error("unknown domains $unknown_domains")

    ret = Array{T}(undef, length(Ξ))
    view(ret, domains .== :Ω) .= rfpotential(:Ω, view(Ξ, domains .== :Ω), bem; kwargs...)
    view(ret, domains .== :Σ) .= rfpotential(:Σ, view(Ξ, domains .== :Σ), bem; kwargs...)
    view(ret, domains .== :Γ) .= rfpotential(:Γ, view(Ξ, domains .== :Γ), bem; kwargs...)
    ret
end

@inline function NESSie.rfpotential(
    Ξ::Base.Generator,
    bem::BEMResult{T};
    kwargs...
) where T
    rfpotential(collect(Vector{T}, Ξ), bem; kwargs...)
end


# =========================================================================================
"""
    _espotential_Γ(ξ::Vector{T}, bem::BEMResult{T})
    _espotential_Γ(Ξ::AbstractVector{T}, bem::BEMResult{T})

Computes the local or nonlocal electrostatic potential for (an) observation point(s) ξ (Ξ)
on the molecular surface.

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
@inline function _espotential_Γ(
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    bem::BEMResult{T};
    kwargs...
) where T
    _espotential_Γ.(Ξ, Ref(bem); kwargs...)
end

@inline function _espotential_Γ(ξ::Vector{T}, bem::BEMResult{T}; kwargs...) where T
    _rfpotential_Γ(ξ, bem) + molpot(ξ, bem.model; kwargs...)
end


# =========================================================================================
"""
    _rfpotential_Γ(ξ::Vector{T}, bem::BEMResult{T})
    _rfpotential_Γ(Ξ::AbstractVector{T}, bem::BEMResult{T})

Computes the local or nonlocal reaction field potential for (an) observation point(s) ξ (Ξ)
on the molecular surface.

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
@inline function _rfpotential_Γ(ξ::Vector{T}, bem::BEMResult{T}) where T
    bem.u[_closest_element_id(ξ, bem.model)] * potprefactor(T)
end

@inline function _rfpotential_Γ(
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    bem::BEMResult{T}
) where T
    _rfpotential_Γ.(Ξ, Ref(bem))
end


# =========================================================================================
"""
    _espotential_Ω(ξ::Vector{T}, bem::BEMResult{T})
    _espotential_Ω(Ξ::AbstractVector{T}, bem::BEMResult{T})

Computes the local or nonlocal electrostatic potential for (an) observation point(s) ξ (Ξ)
inside the molecule.

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
@inline function _espotential_Ω(
    Ξ::AbstractVector{Vector{T}},
    bem::BEMResult{T};
    kwargs...
) where T
    _rfpotential_Ω(Ξ, bem) .+ molpotential(Ξ, bem; kwargs...)
end

@inline function _espotential_Ω(Ξ::Base.Generator, bem::BEMResult{T}; kwargs...) where T
    _espotential_Ω(collect(Vector{T}, Ξ), bem; kwargs...)
end

@inline function _espotential_Ω(ξ::Vector{T}, bem::BEMResult{T}; kwargs...) where T
    only(_espotential_Ω([ξ], bem); kwargs...)
end


# =========================================================================================
"""
    _rfpotential_Ω(ξ::Vector{T}, bem::BEMResult{T})
    _rfpotential_Ω(Ξ::AbstractVector{T}, bem::BEMResult{T})

Computes the local or nonlocal reaction field potential for (an) observation point(s) ξ (Ξ)
inside the molecule.

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
function _rfpotential_Ω(Ξ::AbstractVector{Vector{T}}, bem:: BEMResult{T}) where T
    # result vector
    φ = zeros(T, length(Ξ))

    # φ  = -[W ⋅ u](ξ)
    Rjasanow.laplacecoll!(DoubleLayer, φ, bem.model.elements, Ξ, bem.u)
    rmul!(φ, -1)

    # φ += [Vtilde ⋅ q](ξ)
    Rjasanow.laplacecoll!(SingleLayer, φ, bem.model.elements, Ξ, bem.q)

    # φ *= 1/4π
    # (W and Vtilde were premultiplied by 4π! 4π⋅ε0 from u and q still to be applied)
    rmul!(φ, T(1 / 4π))

    # Apply remaining prefactors:
    # ▶ 4π⋅ε0     for u and q
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    rmul!(φ, potprefactor(T))
end

@inline function _rfpotential_Ω(Ξ::Base.Generator, bem:: BEMResult{T}) where T
    _rfpotential_Ω(collect(Vector{T}, Ξ), bem)
end

@inline function _rfpotential_Ω(ξ::Vector{T}, bem::BEMResult{T}) where T
    only(_rfpotential_Ω([ξ], bem))
end


# =========================================================================================
"""
    _espotential_Σ(ξ::Vector{T}, bem::BEMResult{T})
    _espotential_Σ(Ξ::AbstractVector{T}, bem::BEMResult{T})

Computes the local or nonlocal electrostatic potential for (an) observation point(s) ξ (Ξ)
in the solvent domain.

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
function _espotential_Σ(
    Ξ::AbstractVector{Vector{T}},
    bem::LocalBEMResult{T}
) where T
    # result vector
    φ = zeros(T, length(Ξ))
    buf = Array{T}(undef, length(bem.model.elements))

    # φ  = -εΩ/εΣ ⋅ [Vtilde ⋅ (q + qmol)](ξ)
    copyto!(buf, bem.q)
    _axpy!(1, bem.qmol, buf)
    Rjasanow.laplacecoll!(SingleLayer, φ, bem.model.elements, Ξ, buf)
    rmul!(φ, -bem.model.params.εΩ/bem.model.params.εΣ)

    # φ += [W ⋅ (u + umol)](ξ)
    copyto!(buf, bem.u)
    _axpy!(1, bem.umol, buf)
    Rjasanow.laplacecoll!(DoubleLayer, φ, bem.model.elements, Ξ, buf)

    # Apply remaining prefactors:
    # ▶ 4π        for Vtilde, W
    # ▶ 4π⋅ε0     for u, q, umol, and qmol
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    rmul!(φ, potprefactor(T) / T(4π))

    φ
end

function _espotential_Σ(
    Ξ::AbstractVector{Vector{T}},
    bem::NonlocalBEMResult{T}
) where T
    # result vector
    φ = zeros(T, length(Ξ))

    # convenience aliases
    εΩ  = bem.model.params.εΩ
    εΣ  = bem.model.params.εΣ
    ε∞  = bem.model.params.ε∞
    yuk = yukawa(bem.model.params)
    elements = bem.model.elements

    buf = Array{T}(undef, length(elements))

    # φ  = -V[εΩ/ε∞ ⋅ (q + qmol)](ξ)
    copyto!(buf, bem.q)
    _axpy!(1, bem.qmol, buf)
    Rjasanow.laplacecoll!(SingleLayer, φ, elements, Ξ, buf)
    rmul!(φ, -εΩ/ε∞)

    # φ += (Vʸ-V)[εΩ(1/εΣ - 1/ε∞) ⋅ (q + qmol)](ξ)
    rmul!(buf, εΩ * (1/εΣ - 1/ε∞))
    Radon.regularyukawacoll!(SingleLayer, φ, elements, Ξ, yuk, buf)

    # φ += K[u + umol](ξ)
    copyto!(buf, bem.u)
    _axpy!(1, bem.umol, buf)
    Rjasanow.laplacecoll!(DoubleLayer, φ, elements, Ξ, buf)

    # φ += (Kʸ-K)[u + (1-εΩ/εΣ) ⋅ umol - ε∞/εΣ ⋅ w](ξ)
    copyto!(buf, bem.u)
    _axpy!(1-εΩ/εΣ, bem.umol, buf)
    _axpy!(-ε∞/εΣ, bem.w, buf)
    Radon.regularyukawacoll!(DoubleLayer, φ, elements, Ξ, yuk, buf)

    # Apply remaining prefactors:
    # ▶ 4π        for V, K, (Vʸ-V), and (Kʸ-K)
    # ▶ 4π⋅ε0     for u, q, w, umol, and qmol
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    rmul!(φ, potprefactor(T) / T(4π))

    φ
end

@inline function _espotential_Σ(Ξ::Base.Generator, bem::BEMResult{T}) where T
    _espotential_Σ(collect(Vector{T}, Ξ), bem)
end

@inline function _espotential_Σ(ξ::Vector{T}, bem::BEMResult{T}) where T
    only(_espotential_Σ([ξ], bem))
end


# =========================================================================================
"""
    _rfpotential_Σ(ξ::Vector{T}, bem::BEMResult{T})
    _rfpotential_Σ(Ξ::AbstractVector{T}, bem::BEMResult{T})

Computes the local or nonlocal reaction field potential for (an) observation point(s) ξ (Ξ)
in the solvent domain.

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
@inline function _rfpotential_Σ(
    Ξ::AbstractVector{Vector{T}},
    bem::BEMResult{T};
    kwargs...
) where T
    _espotential_Σ(Ξ, bem) .- molpotential(Ξ, bem; kwargs...)
end

@inline function _rfpotential_Σ(Ξ::Base.Generator, bem:: BEMResult{T}; kwargs...) where T
    _rfpotential_Σ(collect(Vector{T}, Ξ), bem; kwargs...)
end

@inline function _rfpotential_Σ(ξ::Vector{T}, bem::BEMResult{T}; kwargs...) where T
    only(_rfpotential_Σ([ξ], bem; kwargs...))
end
