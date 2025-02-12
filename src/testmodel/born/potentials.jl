# =========================================================================================
"""
    φΩ(::Type{<: LocalityType}, ξ::Vector{T}, ion::BornIon{T})

Computes the interior local or nonlocal electrostatic potential ``φ_Ω`` for the given
observation point ``ξ``.

## Unit
``V = \\frac{C}{F}``

## Return type
`T`

## Alias
    φΩ(::Type{<: LocalityType}, Ξ::Vector{Vector{T}}, ion::BornIon{T})

Computes the potentials for all observation points ``ξ \\in Ξ``.

!!! warning
    This function does not verify whether ξ is located inside of the sphere!
"""
function NESSie.φΩ(
        ::Type{LocalES},
        ξ::Vector{T},
        ion::BornIon{T}
    ) where T
    potprefactor(T) * ion.charge.val *
        (1/euclidean(ion.charge.pos, ξ) + 1/ion.radius * (1/ion.params.εΣ - 1))
end

function NESSie.φΩ(
        ::Type{NonlocalES},
        ξ::Vector{T},
        ion::BornIon{T}
    ) where T
    opt = ion.params
    r = euclidean(ion.charge.pos, ξ)
    ν = sqrt(opt.εΣ/opt.ε∞) * ion.radius / opt.λ
    potprefactor(T) * ion.charge.val * (1/r + 1/ion.radius/opt.εΣ *
        (1 - opt.εΣ + (opt.εΣ - opt.ε∞)/opt.ε∞ * sinh(ν)/ν * exp(-ν)))
end

@inline function NESSie.φΩ(
    lt::Type{<: LocalityType},
    Ξ,
    ion::BornIon{T}
) where T
    φΩ.(lt, Ξ, Ref(ion))
end

# =========================================================================================
"""
    φΣ(::Type{<: LocalityType}, ξ::Vector{T}, ion::BornIon{T})

Computes the exterior local or nonlocal electrostatic potential ``φ_Σ`` for the given
observation point ``ξ``.

## Unit
``V = \\frac{C}{F}``

## Return type
`T`

## Alias
    φΣ(::Type{<: LocalityType}, Ξ::Vector{Vector{T}}, ion::BornIon{T})

Computes the potentials for all observation points ``ξ \\in Ξ``.

!!! warning
    This function does not verify whether ξ is located outside of the sphere!
"""
function NESSie.φΣ(
        ::Type{LocalES},
        ξ::Vector{T},
        ion::BornIon{T}
    ) where T
    potprefactor(T) * ion.charge.val / ion.params.εΣ / euclidean(ion.charge.pos, ξ)
end

function NESSie.φΣ(
        ::Type{NonlocalES},
        ξ::Vector{T},
        ion::BornIon{T}
    ) where T
    opt = ion.params
    r = euclidean(ion.charge.pos, ξ)
    ν = sqrt(opt.εΣ/opt.ε∞) * ion.radius / opt.λ
    potprefactor(T) * ion.charge.val / opt.εΣ /
        r * (1 + (opt.εΣ - opt.ε∞)/opt.ε∞ * sinh(ν)/ν * exp(-ν * r/ion.radius))
end

@inline function NESSie.φΣ(
    lt::Type{<: LocalityType},
    Ξ,
    ion::BornIon{T}
) where T
    φΣ.(lt, Ξ, Ref(ion))
end
