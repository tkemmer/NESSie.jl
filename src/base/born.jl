# =========================================================================================
"""
    type BornIon{T <: AbstractFloat}
        charge::Charge{T}  # point charge at the sphere's center
        radius::T          # sphere radius in Å
    end

Single Born ion, that is, a monoatomic ion represented as a spherically symmetric domain
with a single point charge located at its center (``ε\_Ω = 1``).

## Special constructors

    BornIon{T}(charge::T, radius::T)

Centers the sphere at ``(0, 0, 0)\^T``.
"""
type BornIon{T <: AbstractFloat}
    """Point charge at the sphere's center"""
    charge::Charge{T}
    """Sphere radius in Å"""
    radius::T
end
BornIon{T}(charge::T, radius::T) = BornIon(Charge(T[0, 0, 0], charge), radius)


# =========================================================================================
@doc """
    bornion(::Type{Float64}, name::String)
    bornion(::Type{Float32}, name::String)

Generator function for built-in Born ions:

| Name | Charge | Radius [[Åqv90]](@ref Bibliography) |
|------|-------:|------------------------------------:|
| Li   | +1     | 0.645                               |
| Na   | +1     | 1.005                               |
| K    | +1     | 1.365                               |
| Rb   | +1     | 1.505                               |
| Cs   | +1     | 1.715                               |
| Mg   | +2     | 0.615                               |
| Ca   | +2     | 1.015                               |
| Sr   | +2     | 1.195                               |
| Ba   | +2     | 1.385                               |

## Return type
`BornIon`
""" bornion
for T in [:Float64, :Float32]
    varname = Symbol("bornions_", T)
    @eval begin
        const $(varname) = Dict{String, BornIon{$(T)}}(
            "Li" => BornIon($(T)[1, 0.645]...),
            "Na" => BornIon($(T)[1, 1.005]...),
            "K"  => BornIon($(T)[1, 1.365]...),
            "Rb" => BornIon($(T)[1, 1.505]...),
            "Cs" => BornIon($(T)[1, 1.715]...),
            "Mg" => BornIon($(T)[2, 0.615]...),
            "Ca" => BornIon($(T)[2, 1.015]...),
            "Sr" => BornIon($(T)[2, 1.195]...),
            "Ba" => BornIon($(T)[2, 1.385]...)
        )
        bornion(::Type{$(T)}, name::String) = $(varname)[name]
    end
end


# =========================================================================================
"""
    φΩ{T, L <: LocalityType}(
           ::Type{L},
        ξ  ::Vector{T},
        ion::BornIon{T},
        opt::Option{T} = defaultopt(T)
    )

Computes the interior local or nonlocal electrostatic potential ``φ\_Ω`` for the given
observation point ``ξ``.

## Unit
``V = \\frac{C}{F}``

## Return type
`T`

!!! warning
    This function does not verify whether ξ is located inside of the sphere!
"""
function φΩ{T}(
        ::Type{LocalES},
        ξ::Vector{T},
        ion::BornIon{T},
        opt::Option{T}=defaultopt(T)
    )
    potprefactor(T) * ion.charge.val *
        (1/euclidean(ion.charge.pos, ξ) + 1/ion.radius * (1/opt.εΣ - 1))
end

function φΩ{T}(
        ::Type{NonlocalES},
        ξ::Vector{T},
        ion::BornIon{T},
        opt::Option{T}=defaultopt(T)
    )
    r = euclidean(ion.charge.pos, ξ)
    ν = sqrt(opt.εΣ/opt.ε∞) * ion.radius / opt.λ
    potprefactor(T) * ion.charge.val * (1/r + 1/ion.radius/opt.εΣ *
        (1 - opt.εΣ + (opt.εΣ - opt.ε∞)/opt.ε∞ * sinh(ν)/ν * exp(-ν)))
end


# =========================================================================================
"""
    φΣ{T, L <: LocalityType}(
           ::Type{L},
        ξ  ::Vector{T},
        ion::BornIon{T},
        opt::Option{T} = defaultopt(T)
    )

Computes the exterior local or nonlocal electrostatic potential ``φ\_Σ`` for the given
observation point ``ξ``.

## Unit
``V = \\frac{C}{F}``

## Return type
`T`

!!! warning
    This function does not verify whether ξ is located outside of the sphere!
"""
function φΣ{T}(
        ::Type{LocalES},
        ξ::Vector{T},
        ion::BornIon{T},
        opt::Option{T}=defaultopt(T)
    )
    potprefactor(T) * ion.charge.val / opt.εΣ / euclidean(ion.charge.pos, ξ)
end

function φΣ{T}(
        ::Type{NonlocalES},
        ξ::Vector{T},
        ion::BornIon{T},
        opt::Option{T}=defaultopt(T)
    )
    r = euclidean(ion.charge.pos, ξ)
    ν = sqrt(opt.εΣ/opt.ε∞) * ion.radius / opt.λ
    potprefactor(T) * ion.charge.val / opt.εΣ /
        r * (1 + (opt.εΣ - opt.ε∞)/opt.ε∞ * sinh(ν)/ν * exp(-ν * r/ion.radius))
end
