# =========================================================================================
"""
    mutable struct BornIon{T <: AbstractFloat}
        charge::Charge{T}                          # point charge at the sphere center
        radius::T                                  # sphere radius in Å
        params::Option{T} = defaultopt(BornIon{T}) # system constants
    end

Single Born ion, that is, a monoatomic ion represented as a spherically symmetric domain
with a single point charge located at its center (``ε_Ω = 1``).

## Special constructors

    BornIon(charge::T, radius::T, params::Option{T} = defaultopt(BornIon{T}))

Centers the sphere at ``(0, 0, 0)^T``.
"""
@auto_hash_equals mutable struct BornIon{T <: AbstractFloat}
    """Point charge at the sphere's center"""
    charge::Charge{T}
    """Sphere radius in Å"""
    radius::T
    """System constants"""
    params::Option{T}

    @inline function BornIon{T}(
        charge::Charge{T},
        radius::T,
        params::Option{T} = defaultopt(BornIon{T})
    ) where T
        new(charge, radius, params)
    end
end

@inline BornIon(
    charge::T,
    radius::T,
    params::Option{T} = defaultopt(BornIon{T})
) where T = BornIon{T}(Charge(T[0, 0, 0], charge), radius, params)


# =========================================================================================
@doc raw"""
    defaultopt(::Type{BornIon{T}})

Default system parameters for vacuum-like systems in water.

# Return type
[`Option{T}`](@ref)

# Default values
 * ``ε_Ω = 1``
 * ``ε_Σ = 78``
 * ``ε_∞ = 1.8``
 * ``λ = 20 Å``
""" defaultopt
for T in [:Float64, :Float32]
    varname = Symbol("_defaultopt_born_", T)
    @eval begin
        const $(varname) = Option($(T)[1, 78, 1.8, 20]...)
        @inline NESSie.defaultopt(::Type{BornIon{$(T)}}) = $(varname)
    end
end


# =========================================================================================
@doc """
    bornion(name::String, ::Type{Float64} = Float64)
    bornion(name::String, ::Type{Float32})

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
[`BornIon{T}`](@ref)
""" bornion
for T in [:Float64, :Float32]
    varname = Symbol("bornions_", T)
    @eval begin
        const $(varname) = Dict{String, BornIon{$(T)}}(
            "li" => BornIon($(T)[1, 0.645]...),
            "na" => BornIon($(T)[1, 1.005]...),
            "k"  => BornIon($(T)[1, 1.365]...),
            "rb" => BornIon($(T)[1, 1.505]...),
            "cs" => BornIon($(T)[1, 1.715]...),
            "mg" => BornIon($(T)[2, 0.615]...),
            "ca" => BornIon($(T)[2, 1.015]...),
            "sr" => BornIon($(T)[2, 1.195]...),
            "ba" => BornIon($(T)[2, 1.385]...)
        )
        @inline bornion(name::String, ::Type{$(T)}) = $(varname)[lowercase(name)]
    end
end
@inline bornion(name::String) = bornion(name, Float64)
