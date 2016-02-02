#=
    Computes the molecular potential of the given system of point charges in a structureless medium.

    Note that the results are premultiplied by 4π * εΩ!

    @param ξ/ξlist
                Observation point(s)
    @param charges
                Point charges
    @return T or Vector{T}
=#
φmol{T}(ξlist::Vector{Vector{T}}, charges::Vector{Charge{T}}) = [φmol(ξ, charges) for ξ in ξlist]
φmol{T}(ξlist::Vector{Triangle{T}}, charges::Vector{Charge{T}}) = [φmol(ξ.center, charges) for ξ in ξlist]
φmol{T}(ξ::Vector{T}, charges::Vector{Charge{T}}) = sum([q.val / euclidean(ξ, q.pos) for q in charges])

#=
    Computes the molecular potential of the given system of point charges in a structureless medium.

    Note that the results are premultiplied by 4π * εΩ!

    @param ξ/ξlist
                Observation point(s)
    @param charges
                Point charges
    @return T or Vector{T}
=#
∂ₙφmol{T}(ξlist::Vector{Triangle{T}}, charges::Vector{Charge{T}}) = [∂ₙφmol(ξ, charges) for ξ in ξlist]
∂ₙφmol{T}(ξ::Triangle{T}, charges::Vector{Charge{T}}) = -sum([q.val * ((ξ.center - q.pos) ⋅ ξ.normal) / euclidean(ξ.center, q.pos)^3 for q in charges])
