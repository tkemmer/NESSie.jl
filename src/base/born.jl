#=
    Single Born ion, that is, a monoatomic ion represented as a spherically symmetric domain with a single point charge
    located at its center (εΩ = 1).
=#
type BornIon{T <: AbstractFloat}
    charge::Charge{T}
    radius::T # Å
end
BornIon{T}(charge::T, radius::T) = BornIon(Charge(T[0, 0, 0], charge), radius)

#=
    A selection of Born ions.

    Reference:
    [1]  J. Åqvist. Ion-water interaction potentials derived from free energy pertubation simulations.
         J. Phys. Chem., 94:8021, 1990.
=#
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

#=
    Computes the interior local electrostatic potential φΩ for the given observation point ξ. Note that the function
    does not verify whether the given points actually lie inside of the Born sphere.

    @param ξ
        Observation point; has to be located inside the Born sphere!
    @param ion
        A Born ion
    @param opt
        Constants to be used
    @return T
=#
function φΩ{T}(::Type{LocalES}, ξ::Vector{T}, ion::BornIon{T}, opt::Option{T}=defaultopt(T))
    T(1.69e-9 / 4π / ε0) * ion.charge.val * (1/euclidean(ion.charge.pos, ξ) + 1/ion.radius * (1/opt.εΣ - 1))
end

#=
    Computes the exterior local electrostatic potential φΣ for the given observation point ξ. Note that the function
    does not verify whether the given points actually lie in the surrounding space.

    @param ξ
        Observation point; has to be located outside the Born sphere!
    @param ion
        A Born ion
    @param opt
        Constants to be used
    @return T
=#
function φΣ{T}(::Type{LocalES}, ξ::Vector{T}, ion::BornIon{T}, opt::Option{T}=defaultopt(T))
    T(1.69e-9 / 4π / ε0) * ion.charge.val / opt.εΣ / euclidean(ion.charge.pos, ξ)
end

#=
    Computes the interior nonlocal electrostatic potential φΩ for the given observation point ξ. Note that the function
    does not verify whether the given points actually lie inside of the Born sphere.

    @param ξ
        Observation point; has to be located inside the Born sphere!
    @param ion
        A Born ion
    @param opt
        Constants to be used
    @return T
=#
function φΩ{T}(::Type{NonlocalES}, ξ::Vector{T}, ion::BornIon{T}, opt::Option{T}=defaultopt(T))
    r = euclidean(ion.charge.pos, ξ)
    ν = sqrt(opt.εΣ/opt.ε∞) * ion.radius / opt.λ
    T(1.69e-9 / 4π / ε0) * ion.charge.val * (1/r + 1/ion.radius/opt.εΣ * (1 - opt.εΣ + (opt.εΣ - opt.ε∞)/opt.ε∞ * sinh(ν)/ν * exp(-ν)))
end

#=
    Computes the exterior nonlocal electrostatic potential φΣ for the given observation point ξ. Note that the function
    does not verify whether the given points actually lie in the surrounding space.

    @param ξ
        Observation point; has to be located outside the Born sphere!
    @param ion
        A Born ion
    @param opt
        Constants to be used
    @return T
=#
function φΣ{T}(::Type{NonlocalES}, ξ::Vector{T}, ion::BornIon{T}, opt::Option{T}=defaultopt(T))
    r = euclidean(ion.charge.pos, ξ)
    ν = sqrt(opt.εΣ/opt.ε∞) * ion.radius / opt.λ
    T(1.69e-9 / 4π / ε0) * ion.charge.val / opt.εΣ / r * (1 + (opt.εΣ - opt.ε∞)/opt.ε∞ * sinh(ν)/ν * exp(-ν * r/ion.radius))
end
