# `opt` argument removed in v1.5
function NESSie.φΩ(
    lt::Type{<: LocalityType},
    ξ::Vector{T},
    ion::BornIon{T},
    opt::Option{T}
) where T
    Base.depwarn(
        "The `opt` parameter is deprecated and will be removed in a future " *
        "release! Instead, please set `ion.params` accordingly.",
        :φΩ
    )
    ion_copy = deepcopy(ion)
    ion_copy.params = opt
    φΩ(lt, ξ, ion_copy)
end

# `opt` argument removed in v1.5
function NESSie.φΣ(
    lt::Type{<: LocalityType},
    ξ::Vector{T},
    ion::BornIon{T},
    opt::Option{T}
) where T
    Base.depwarn(
        "The `opt` parameter is deprecated and will be removed in a future " *
        "release! Instead, please set `ion.params` accordingly.",
        :φΣ
    )
    ion_copy = deepcopy(ion)
    ion_copy.params = opt
    φΣ(lt, ξ, ion_copy)
end

# `XieModel` renamed to `XieSphere` in v1.5
Base.@deprecate_binding XieModel XieSphere

# removed in v1.5
@deprecate NESSie.φΩ(lt, ξorΞ, ion::BornIon) espotential(:Ω, ξorΞ, ion) false
@deprecate NESSie.φΣ(lt, ξorΞ, ion::BornIon) espotential(:Σ, ξorΞ, ion) false
@deprecate NESSie.φΩ(ξorΞ, xie::XieTestModel) espotential(:Ω, ξorΞ, xie) false
@deprecate NESSie.φΣ(ξorΞ, xie::XieTestModel) espotential(:Σ, ξorΞ, xie) false
