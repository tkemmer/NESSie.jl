# `opt` argument removed in v1.5
function NESSie.φΩ(
    lt::Type{<: LocalityType},
    ξ::Vector{T},
    ion::BornIon{T},
    opt::Option{T}
) where T
    Base.depwarn(
        "The `opt` parameter is deprecated and will be removed in a future" *
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
        "The `opt` parameter is deprecated and will be removed in a future" *
        "release! Instead, please set `ion.params` accordingly.",
        :φΣ
    )
    ion_copy = deepcopy(ion)
    ion_copy.params = opt
    φΣ(lt, ξ, ion_copy)
end
