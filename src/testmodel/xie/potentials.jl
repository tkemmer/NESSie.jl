# =========================================================================================
"""
    function φΩ{T}(
        ξ    ::Vector{T},
        model::NonlocalXieModel1{T}
    )

Computes the interior nonlocal electrostatic potential ``φ_Ω`` for the given observation
point ``ξ``.

# Unit
``V = \\frac{C}{F}``

# Return type
`T`

!!! warning
    This function does not verify whether ξ is located inside of the sphere!
"""
function φΩ(ξ::Vector{T}, model::NonlocalXieModel1{T}) where T
    const a  = model.radius
    const λ  = model.params.λ
    const εΩ = model.params.εΩ
    const εΣ = model.params.εΣ
    const ε∞ = model.params.ε∞
    const A₃ = model.A₃
    const r  = vecnorm(ξ)

    φ = zero(T)
    for (qi, q) in enumerate(model.charges)   # Eq. (17a)

        # if ξ is close to the origin, all terms for n > 0 become negligible
        if r < 1e-10
            φ += A₃[1, qi]
            continue
        end

        # if q is close to the origin, compute nonlocal Born potential
        # (see born/potentials.jl)
        if vecnorm(q.pos) < 1e-10
            ν = √(εΣ/ε∞) * a/λ
            φ += q.val * (1/r + 1/a/εΣ * (1 - εΣ + (εΣ - ε∞)/ε∞ * sinh(ν)/ν * exp(-ν)))
            continue
        end

        # otherwise, use Eq. (18)
        P = legendre(model.len, cos(ξ, q.pos, r))
        φj = zero(T)
        for n in 0:model.len-1
            φj += A₃[n+1, qi] * r^n * P(n)
        end
        φ += φj * q.val
    end

    (φ + φmol(ξ, model.charges) / 4π / εΩ) * T(ec/ε0)
end


# =========================================================================================
"""
    function φΣ{T}(
        ξ    ::Vector{T},
        model::NonlocalXieModel1{T}
    )

Computes the exterior nonlocal electrostatic potential ``φ_Σ`` for the given observation
point ``ξ``.

# Unit
``V = \\frac{C}{F}``

# Return type
`T`

!!! warning
    This function does not verify whether ξ is located outside of the sphere!
"""
function φΣ(ξ::Vector{T}, model::NonlocalXieModel1{T}) where T
    const a  = model.radius
    const λ  = model.params.λ
    const εΣ = model.params.εΣ
    const ε∞ = model.params.ε∞
    const A₁ = model.A₁
    const A₂ = model.A₂
    const κ  = λ \ √(εΣ/ε∞)
    const r  = vecnorm(ξ)

    kᵣ = spherical_besselk(model.len, κ * r)

    φ = zero(T)
    for (qi, q) in enumerate(model.charges)   # Eq. (17a)

        # if q is close to origin, compute nonlocal Born potential
        # (see born/potentials.jl)
        if vecnorm(q.pos) < 1e-10
            ν = √(εΣ/ε∞) * a/λ
            φ += q.val / εΣ / r * (1 + (εΣ - ε∞)/ε∞ * sinh(ν)/ν * exp(-ν * r/a))
            continue
        end

        # otherwise, use Eq. (18)
        P = legendre(model.len, cos(ξ, q.pos, r))
        φj = zero(T)
        for n in 0:model.len-1
            φj += (ε∞ - εΣ)/ε∞ * A₂[n+1, qi] * kᵣ(n) * P(n) +
                  A₁[n+1, qi] / r^(n+1) * P(n)
        end
        φ += φj * q.val
    end

    φ * T(ec/ε0)
end
