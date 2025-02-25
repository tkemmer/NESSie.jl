# =========================================================================================
"""
    function φΩ(
        ξ    ::Vector{T},
        model::NonlocalXieModel1{T}
    )

Computes the interior nonlocal electrostatic potential ``φ_Ω`` for the given observation
point ``ξ``.

# Supported keyword arguments
 - `tolerance::T = T(1e-10)` minimum distance assumed between ξ and a point charge when
   computing molecular potentials (cf. [`φmol`](@ref)).

# Unit
``V = \\frac{C}{F}``

# Return type
`T`

## Alias
    φΩ(Ξ::Vector{Vector{T}}, model::NonlocalXieModel1{T})

Computes the potentials for all observation points ``ξ \\in Ξ``.

!!! warning
    This function does not verify whether ξ is located inside of the sphere!
"""
function NESSie.φΩ(ξ::Vector{T}, model::NonlocalXieModel1{T}; tolerance::T = T(1e-10)) where T
    a  = model.radius
    λ  = model.params.λ
    εΩ = model.params.εΩ
    εΣ = model.params.εΣ
    ε∞ = model.params.ε∞
    A₃ = model.A₃
    r  = _norm(ξ)

    φ = zero(T)
    for (qi, q) in enumerate(model.charges)   # Eq. (17a)

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

        # otherwise, use Eq. (18)
        # if ξ is close to the origin, all terms for n > 0 become negligible
        if r < 1e-10
            φ += A₃[1, qi] * q.val
            continue
        end

        P = legendre(model.len, _cos(ξ, q.pos, r))
        φj = zero(T)
        for n in 0:model.len-1
            φj += A₃[n+1, qi] * r^n * P(n)
        end
        φ += φj * q.val
    end

    (φ + φmol(ξ, model.charges; tolerance = tolerance) / T(4π) / εΩ) * T(ec/ε0)
end

@inline function NESSie.φΩ(
    Ξ::Union{AbstractVector{Vector{T}}, <: Base.Generator},
    model::NonlocalXieModel1{T}
) where T
    φΩ.(Ξ, Ref(model))
end


# =========================================================================================
"""
    function φΣ(
        ξ    ::Vector{T},
        model::NonlocalXieModel1{T}
    )

Computes the exterior nonlocal electrostatic potential ``φ_Σ`` for the given observation
point ``ξ``.

# Supported keyword arguments
 - `tolerance::T = T(1e-10)` minimum distance assumed between ξ and a point charge when
   computing molecular potentials (cf. [`φmol`](@ref)).

# Unit
``V = \\frac{C}{F}``

# Return type
`T`

## Alias
    φΣ(Ξ::Vector{Vector{T}}, model::NonlocalXieModel1{T})

Computes the potentials for all observation points ``ξ \\in Ξ``.

!!! warning
    This function does not verify whether ξ is located outside of the sphere!
"""
function NESSie.φΣ(ξ::Vector{T}, model::NonlocalXieModel1{T}; tolerance::T = T(1e-10)) where T
    a  = model.radius
    λ  = model.params.λ
    εΩ = model.params.εΩ
    εΣ = model.params.εΣ
    ε∞ = model.params.ε∞
    A₁ = model.A₁
    A₂ = model.A₂
    κ  = λ \ √(εΣ/ε∞)
    r  = _norm(ξ)

    kᵣ = spherical_besselk(model.len, κ * r)

    φ = zero(T)
    for (qi, q) in enumerate(model.charges)   # Eq. (17a)

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

        # otherwise, use Eq. (18)
        P = legendre(model.len, _cos(ξ, q.pos, r))
        φj = zero(T)
        for n in 0:model.len-1
            φj += (ε∞ - εΣ)/ε∞ * A₂[n+1, qi] * kᵣ(n) * P(n) +
                  A₁[n+1, qi] / r^(n+1) * P(n)
        end
        φ += φj * q.val
    end

    φ * T(ec/ε0)
end

@inline function NESSie.φΣ(
    Ξ::Union{AbstractVector{Vector{T}}, <: Base.Generator},
    model::NonlocalXieModel1{T}
) where T
    φΣ.(Ξ, Ref(model))
end
