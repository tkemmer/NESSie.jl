# =========================================================================================
"""
    struct BornIon{T <: AbstractFloat}
        charge::Charge{T}  # point charge at the sphere's center
        radius::T          # sphere radius in Å
    end

Single Born ion, that is, a monoatomic ion represented as a spherically symmetric domain
with a single point charge located at its center (``ε\_Ω = 1``).

## Special constructors

    BornIon{T}(charge::T, radius::T)

Centers the sphere at ``(0, 0, 0)\^T``.
"""
struct BornIon{T <: AbstractFloat}
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
        bornion(::Type{$(T)}, name::String) = $(varname)[lowercase(name)]
    end
end
