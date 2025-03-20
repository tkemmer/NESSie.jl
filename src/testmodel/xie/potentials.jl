# =========================================================================================
"""
    rfenergy(::XieTestModel{T})

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
function NESSie.rfenergy(xie::XieTestModel{T}; kwargs...) where T
    Ξ = collect(Vector{T}, charge.pos for charge in xie.charges)
    ζ = collect(T, charge.val for charge in xie.charges)

    ζ ⋅ rfpotential(:Ω, Ξ, xie; kwargs...) * T(ec * 6.022140857e10 / 2)
end


# =========================================================================================
"""
    espotential(ξ::Vector{T}, xie::XieTestModel{T})
    espotential(Ξ::AbstractVector{Vector{T}}, xie::XieTestModel{T})

Computes the local or nonlocal electrostatic potential(s) at the given observation point(s)
ξ (Ξ) for the given Xie test model.

The electrostatic potential is computed as the sum of the corresponding
[reaction field potential](@ref rfpotential) and the [molecular potential](@ref molpotential).

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`

# Alias
    espotential(domain::Symbol, ξ::Vector{T}, xie::XieTestModel{T})
    espotential(domain::Symbol, Ξ::AbstractVector{T}, xie::XieTestModel{T})

Computes the electrostatic potential(s) for the given observation point(s) ξ (Ξ) and the
given domain `:Ω`, `:Σ`, or `:Γ`.
"""
@inline function NESSie.espotential(
    domain::Symbol,
    ξorΞ::Union{Vector{T}, <: AbstractVector{Vector{T}}, <: Base.Generator},
    xie::XieTestModel{T};
    kwargs...
) where T
    domain === :Ω && return _espotential_Ω(ξorΞ, xie; kwargs...)
    domain === :Γ && return _espotential_Ω(ξorΞ, xie; kwargs...)
    domain === :Σ && return _espotential_Σ(ξorΞ, xie; kwargs...)
    error("unknown domain $domain")
end

@inline function NESSie.espotential(
    ξ::Vector{T},
    xie::XieTestModel{T};
    kwargs...
) where T
    norm(ξ) <= xie.radius ?
        espotential(:Ω, ξ, xie; kwargs...) :
        espotential(:Σ, ξ, xie; kwargs...)
end

@inline function NESSie.espotential(
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    xie::XieTestModel{T};
    kwargs...
) where T
    espotential.(Ξ, Ref(xie); kwargs...)
end


# =========================================================================================
"""
    molpotential(ξ::Vector{T}, xie::XieSphere{T})
    molpotential(ξ::Vector{T}, xie::XieTestModel{T})
    molpotential(Ξ::AbstractVector{Vector{T}}, xie::XieSphere{T})
    molpotential(Ξ::AbstractVector{Vector{T}}, xie::XieTestModel{T})

Computes the molecular potential(s) at the given observation point(s) ξ (Ξ) for the given
Xie sphere or Xie test model.

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
    xie::Union{XieSphere{T}, <: XieTestModel{T}};
    kwargs...
) where T
    potprefactor(T) .* _molpotential(ξorΞ, xie.charges; kwargs...) ./ xie.params.εΩ
end


# =========================================================================================
"""
    rfpotential(ξ::Vector{T}, xie::XieTestModel{T})
    rfpotential(Ξ::AbstractVector{Vector{T}}, xie::XieTestModel{T})

Computes the local or nonlocal reaction field potential(s) at the given observation point(s)
ξ (Ξ) for the given Xie test model.

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`

# Alias
    rfpotential(domain::Symbol, ξ::Vector{T}, xie::XieTestModel{T})
    rfpotential(domain::Symbol, Ξ::AbstractVector{T}, xie::XieTestModel{T})

Computes the reaction field potential(s) for the given observation point(s) ξ (Ξ) and the
given domain `:Ω`, `:Σ`, or `:Γ`.
"""
@inline function NESSie.rfpotential(
    domain::Symbol,
    ξorΞ::Union{Vector{T}, <: AbstractVector{Vector{T}}, <: Base.Generator},
    xie::XieTestModel{T};
    kwargs...
) where T
    domain === :Ω && return _rfpotential_Ω(ξorΞ, xie)
    domain === :Γ && return _rfpotential_Ω(ξorΞ, xie)
    domain === :Σ && return _rfpotential_Σ(ξorΞ, xie; kwargs...)
    error("unknown domain $domain")
end

@inline function NESSie.rfpotential(
    ξ::Vector{T},
    xie::XieTestModel{T};
    kwargs...
) where T
    norm(ξ) <= xie.radius ?
        rfpotential(:Ω, ξ, xie; kwargs...) :
        rfpotential(:Σ, ξ, xie; kwargs...)
end

@inline function NESSie.rfpotential(
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    xie::XieTestModel{T};
    kwargs...
) where T
    rfpotential.(Ξ, Ref(xie); kwargs...)
end


# =========================================================================================
"""
    _espotential_Ω(ξ::Vector{T}, xie::XieTestModel{T})
    _espotential_Ω(Ξ::AbstractVector{T}, xie::XieTestModel{T})

Computes the local or nonlocal electrostatic potential for (an) observation point(s) ξ (Ξ)
inside the test model sphere.

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
@inline function _espotential_Ω(
    ξorΞ::Union{Vector{T}, <: AbstractVector{Vector{T}}, <: Base.Generator},
    xie::XieTestModel{T};
    kwargs...
) where T
    _rfpotential_Ω(ξorΞ, xie) .+ molpotential(ξorΞ, xie; kwargs...)
end


# =========================================================================================
"""
    _rfpotential_Ω(ξ::Vector{T}, xie::XieTestModel{T})
    _rfpotential_Ω(Ξ::AbstractVector{T}, xie::XieTestModel{T})

Computes the local or nonlocal reaction field potential for (an) observation point(s) ξ (Ξ)
inside the test model sphere.

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
function _rfpotential_Ω(
    ξ::Vector{T},
    xie::XieTestModel{T}
) where T
    a  = xie.radius
    λ  = xie.params.λ
    εΩ = xie.params.εΩ
    εΣ = xie.params.εΣ
    ε∞ = xie.params.ε∞
    M₃ = xie.M₃
    r  = _norm(ξ)

    φ = zero(T)
    for (qi, q) in enumerate(xie.charges)   # Eq. (17a)

        # if q is close to the origin, compute nonlocal Born potential
        # Note: this test model uses a different definition than our Born implementation
        # https://doi.org/10.4208/cicp.170811.211011s
        if _norm(q.pos) < 1e-10
            _term1 = (a * εΣ + λ * (εΩ - εΣ) * sinh(a / λ)) / (
                (a * √(ε∞ * εΣ) + λ * (ε∞ - εΣ)) * sinh(a / λ) + a * εΣ * cosh(a / λ)
            )
            _term2 = (εΩ - εΣ - (ε∞ - εΣ) * _term1) / a / εΣ
            φ += q.val / T(4π) / εΩ * _term2
            continue
        end

        # otherwise, use Eq. (18)/(38a)
        # if ξ is close to the origin, all terms for n > 0 become negligible
        if r < T(1e-10)
            φ += M₃[1, qi] * q.val
            continue
        end

        P = legendre(xie.len, _cos(ξ, q.pos, r))
        φj = zero(T)
        for n in 0:xie.len-1
            φj += M₃[n+1, qi] * r^n * P(n)
        end
        φ += φj * q.val
    end

    φ * T(ec/ε0)
end

@inline function _rfpotential_Ω(
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    xie::XieTestModel{T}
) where T
    _rfpotential_Ω.(Ξ, Ref(xie))
end


# =========================================================================================
"""
    _espotential_Σ(ξ::Vector{T}, xie::XieTestModel{T})
    _espotential_Σ(Ξ::AbstractVector{T}, xie::XieTestModel{T})

Computes the local or nonlocal electrostatic potential for (an) observation point(s) ξ (Ξ)
outside the test model sphere.

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
function _espotential_Σ(
    ξ::Vector{T},
    xie::XieTestModel{T};
    tolerance::T = T(1e-10)
) where T
    a  = xie.radius
    λ  = xie.params.λ
    εΩ = xie.params.εΩ
    εΣ = xie.params.εΣ
    ε∞ = xie.params.ε∞
    M₁ = xie.M₁
    M₂ = xie.M₂
    κ  = λ \ √(εΣ/ε∞)
    r  = _norm(ξ)

    kᵣ = spherical_besselk(xie.len, κ * r)

    φ = zero(T)
    for (qi, q) in enumerate(xie.charges)   # Eq. (17a)

        # if q is close to origin, compute nonlocal Born potential
        # Note: this test model uses a different definition than our Born implementation
        # https://doi.org/10.4208/cicp.170811.211011s
        if _norm(q.pos) < 1e-10
            _term1 = exp(κ * a) * (εΣ - ε∞)/εΩ * (a * εΣ + λ * (εΩ - εΣ) * sinh(a / λ))
            _term2 = (a * √(ε∞ * εΣ) + λ *(ε∞ - εΣ)) * sinh(a / λ) + a * εΣ * cosh(a / λ)
            _term3 = _term1 / _term2 * exp(-κ * r)
            φ += (1 + _term3) / εΣ * q.val / T(4π) / max(r, tolerance)
            continue
        end

        # otherwise, use Eq. (18)/(38b)
        P = legendre(xie.len, _cos(ξ, q.pos, r))
        φj = zero(T)
        for n in 0:xie.len-1
            φj += (ε∞ - εΣ)/ε∞ * M₂[n+1, qi] * kᵣ(n) * P(n) +
                  M₁[n+1, qi] / r^(n+1) * P(n)
        end
        φ += φj * q.val
    end

    φ * T(ec/ε0)
end

@inline function _espotential_Σ(
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    xie::XieTestModel{T};
    kwargs...
) where T
    _espotential_Σ.(Ξ, Ref(xie); kwargs...)
end


# =========================================================================================
"""
    _rfpotential_Σ(ξ::Vector{T}, xie::XieTestModel{T})
    _rfpotential_Σ(Ξ::AbstractVector{T}, xie::XieTestModel{T})

Computes the local or nonlocal reaction field potential for (an) observation point(s) ξ (Ξ)
outside the test model sphere.

# Supported keyword arguments
See [`molpotential`](@ref)

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
@inline function _rfpotential_Σ(
    ξorΞ::Union{Vector{T}, <: AbstractVector{Vector{T}}, <: Base.Generator},
    xie::XieTestModel{T};
    kwargs...
) where T
    _espotential_Σ(ξorΞ, xie; kwargs...) .- molpotential(ξorΞ, xie; kwargs...)
end
