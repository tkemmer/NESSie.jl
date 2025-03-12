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

Computes the electrostatic potential at the given observation point ξ for the given BEM
result. This function tries to automatically locate the given observation point using
[`guess_domain`](@ref).

# Supported keyword arguments
See [`guess_domain`](@ref)

# Return type
`T`

# Alias
    espotential(Ξ::AbstractVector{Vector{T}}, bem::BEMResult{T})

Computes the electrostatic potentials for all observation points ``ξ \\in Ξ``.
"""
function NESSie.espotential(ξ::Vector{T}, bem::BEMResult{T}; kwargs...) where T
    loc = guess_domain(ξ, bem.model; kwargs...)
    loc === :Ω && return φΩ(ξ, bem)
    loc === :Σ && return φΣ(ξ, bem)
    loc === :Γ && return φΓ(ξ, bem)
    error("unknown location $loc")
end

function NESSie.espotential(Ξ::AbstractVector{Vector{T}}, bem::BEMResult{T}; kwargs...) where T
    locs = guess_domain.(Ξ, Ref(bem.model); kwargs...)
    unknown_locs = setdiff(locs, [:Ω, :Σ, :Γ])
    !isempty(unknown_locs) && error("unknown locations $unknown_locs")

    ret = Array{T}(undef, length(Ξ))
    view(ret, locs .== :Ω) .= φΩ(view(Ξ, locs .== :Ω), bem)
    view(ret, locs .== :Σ) .= φΣ(view(Ξ, locs .== :Σ), bem)
    view(ret, locs .== :Γ) .= φΓ(view(Ξ, locs .== :Γ), bem)
    ret
end

@inline function NESSie.espotential(Ξ::Base.Generator, bem::BEMResult{T}; kwargs...) where T
    espotential(collect(Vector{T}, Ξ), bem; kwargs...)
end


# =========================================================================================
"""
    function NESSie.φΓ(Ξ::AbstractVector{Vector{T}}, bem::BEMResult{T})

Computes the local or nonlocal surface potential ``φ_Γ`` for the given set of observation
points `Ξ`.

!!! warning
    This function does not verify whether all points in `Ξ` are located in ``Γ``!

# Unit
``V = \\frac{C}{F}``

# Return type
`Vector{T}`
"""
@inline function NESSie.φΓ(
    Ξ::Union{AbstractVector{Vector{T}}, <: Base.Generator},
    bem::BEMResult{T}
) where T
    φΓ.(Ξ, Ref(bem))
end

function NESSie.φΓ(ξ::Vector{T}, bem::BEMResult{T}) where T
    (
        bem.u[_closest_element_id(ξ, bem.model)] +
        φmol(ξ, bem.model.charges) / bem.model.params.εΩ
    ) / T(4π / ε0 * ec)
end


# =========================================================================================
"""
    φΩ(Ξ::Vector{Vector{T}}, bem::BEMResult{T})

Computes the local or nonlocal interior electrostatic potential ``φ_Ω`` for the given set
of observation points `Ξ`.

!!! warning
    This function does not verify whether all points in `Ξ` are located in ``Ω``!

# Unit
``V = \\frac{C}{F}``

# Return type
`Vector{T}`
"""
function NESSie.φΩ(
    Ξ::AbstractVector{Vector{T}},
    bem::R
) where {T, R <: BEMResult{T}}
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

    # φ += 1/εΩ ⋅ φ*mol(ξ)
    # (φ*mol was premultiplied by 4π⋅ε0⋅εΩ; 4π⋅ε0 remain to be applied)
    _axpy!(1/bem.model.params.εΩ, φmol(Ξ, bem.model.charges), φ)

    # Apply remaining prefactors:
    # ▶ 4π⋅ε0     for u, q, and umol
    # ▶ 1.602e-19 for elemental charge e; [e] = C
    # ▶ 1e10      for the conversion Å → m; [ε0] = F/m
    rmul!(φ, potprefactor(T))

    φ
end

@inline function NESSie.φΩ(Ξ::Base.Generator, bem::BEMResult{T}) where T
    φΩ(collect(Vector{T}, Ξ), bem)
end

@inline function NESSie.φΩ(ξ::Vector{T}, bem::BEMResult{T}) where T
    only(φΩ([ξ], bem))
end


# =========================================================================================
"""
    φΣ(Ξ::Vector{Vector{T}}, bem::BEMResult{T})

Computes the local or nonlocal exterior electrostatic potential ``φ_Σ`` for the given set
of observation points `Ξ`.

!!! warning
    This function does not verify whether all points in `Ξ` are located in ``Σ``!

# Unit
``V = \\frac{C}{F}``

# Return type
`Vector{T}`
"""
function NESSie.φΣ(
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

function NESSie.φΣ(
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

@inline function NESSie.φΣ(Ξ::Base.Generator, bem::BEMResult{T}) where T
    φΣ(collect(Vector{T}, Ξ), bem)
end

@inline function NESSie.φΣ(ξ::Vector{T}, bem::BEMResult{T}) where T
    only(φΣ([ξ], bem))
end
