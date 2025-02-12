# =========================================================================================
"""
    φΩ(
           ::Type{<: LocalityType},
        ξ  ::Vector{T},
        ion::BornIon{T},
        opt::Option{T} = defaultopt(T)
    )

Computes the interior local or nonlocal electrostatic potential ``φ_Ω`` for the given
observation point ``ξ``.

## Unit
``V = \\frac{C}{F}``

## Return type
`T`

!!! warning
    This function does not verify whether ξ is located inside of the sphere!
"""
function NESSie.φΩ(
        ::Type{LocalES},
        ξ::Vector{T},
        ion::BornIon{T},
        opt::Option{T}=defaultopt(T)
    ) where T
    potprefactor(T) * ion.charge.val *
        (1/euclidean(ion.charge.pos, ξ) + 1/ion.radius * (1/opt.εΣ - 1))
end

function NESSie.φΩ(
        ::Type{NonlocalES},
        ξ::Vector{T},
        ion::BornIon{T},
        opt::Option{T}=defaultopt(T)
    ) where T
    r = euclidean(ion.charge.pos, ξ)
    ν = sqrt(opt.εΣ/opt.ε∞) * ion.radius / opt.λ
    potprefactor(T) * ion.charge.val * (1/r + 1/ion.radius/opt.εΣ *
        (1 - opt.εΣ + (opt.εΣ - opt.ε∞)/opt.ε∞ * sinh(ν)/ν * exp(-ν)))
end


# =========================================================================================
"""
    φΣ(
           ::Type{<: LocalityType},
        ξ  ::Vector{T},
        ion::BornIon{T},
        opt::Option{T} = defaultopt(T)
    )

Computes the exterior local or nonlocal electrostatic potential ``φ_Σ`` for the given
observation point ``ξ``.

## Unit
``V = \\frac{C}{F}``

## Return type
`T`

!!! warning
    This function does not verify whether ξ is located outside of the sphere!
"""
function NESSie.φΣ(
        ::Type{LocalES},
        ξ::Vector{T},
        ion::BornIon{T},
        opt::Option{T}=defaultopt(T)
    ) where T
    potprefactor(T) * ion.charge.val / opt.εΣ / euclidean(ion.charge.pos, ξ)
end

function NESSie.φΣ(
        ::Type{NonlocalES},
        ξ::Vector{T},
        ion::BornIon{T},
        opt::Option{T}=defaultopt(T)
    ) where T
    r = euclidean(ion.charge.pos, ξ)
    ν = sqrt(opt.εΣ/opt.ε∞) * ion.radius / opt.λ
    potprefactor(T) * ion.charge.val / opt.εΣ /
        r * (1 + (opt.εΣ - opt.ε∞)/opt.ε∞ * sinh(ν)/ν * exp(-ν * r/ion.radius))
end
