#=
    Enum-like representation of potential types
=#
abstract PotentialType
type SingleLayer <: PotentialType end
type DoubleLayer <: PotentialType end

#=
    Computes the molecular potential of the given system of point charges in a structureless medium.

    φmol(r) = 1/(4π*ε0*εΩ) * Σ qᵢ/|rᵢ-r|

    Note that the results are premultiplied by 4π * ε0 * εΩ!

    @param ξ/ξlist
                Observation point(s)
    @param charges
                Point charges
    @return T or Vector{T}
=#
φmol{T}(ξlist::Vector{Vector{T}}, charges::Vector{Charge{T}}) = [φmol(ξ, charges) for ξ in ξlist]
φmol{T}(ξlist::Vector{Triangle{T}}, charges::Vector{Charge{T}}) = [φmol(ξ.center, charges) for ξ in ξlist]
φmol{T}(ξ::Vector{T}, charges::Vector{Charge{T}}) = sum([q.val / euclidean(ξ, q.pos) for q in charges])
# TODO devectorize!

#=
    Computes the molecular potential of the given system of point charges in a structureless medium.

    TODO φmol(r) =

    Note that the results are premultiplied by 4π * εΩ * ε0!

    @param ξ/ξlist
                Observation point(s)
    @param charges
                Point charges
    @return T or Vector{T}
=#
∂ₙφmol{T}(ξlist::Vector{Triangle{T}}, charges::Vector{Charge{T}}) = [∂ₙφmol(ξ, charges) for ξ in ξlist]
∂ₙφmol{T}(ξ::Triangle{T}, charges::Vector{Charge{T}}) = -sum([q.val * ddot(ξ.center, q.pos, ξ.normal) / euclidean(ξ.center, q.pos)^3 for q in charges])
# TODO devectorize!

#=
    TODO
=#
∇φmol{T}(ξlist::Vector{Vector{T}}, charges::Vector{Charge{T}}) = [∇φmol(ξ, charges) for ξ in ξlist]
∇φmol{T}(ξ::Vector{T}, charges::Vector{Charge{T}}) = -sum([q.val * (ξ - q.pos) / euclidean(ξ, q.pos)^3 for q in charges])
# TODO devectorize!
