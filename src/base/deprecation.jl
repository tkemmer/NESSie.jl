# removed in v1.5
@deprecate obspoints_line(u, v, n) LinRange(u, v, n)

# implementations removed in v1.5
export φΣ, φΩ, φΓ

function φΣ end
function φΩ end
function φΓ end

# deprecated since v1.5
export φmol, ∂ₙφmol, ∇φmol

@inline function φmol(
    ξorΞ::Union{Vector{T}, <: AbstractVector{Vector{T}}, <: Base.Generator},
    charges::AbstractVector{Charge{T}};
    kwargs...
) where T
    Base.depwarn(
        "The unscaled, model-less molecular potential function `φmol` is deprecated " *
        "and will be removed in a future release. Please use `molpotential` instead!",
        :φmol
    )
    _molpotential(ξorΞ, charges; kwargs...)
end

@inline function φmol(
    model::Model{T, Triangle{T}};
    kwargs...
) where T
    Base.depwarn(
        "The unscaled, model-less molecular potential function `φmol` is deprecated " *
        "and will be removed in a future release. Please use `molpotential` instead!",
        :φmol
    )
    _molpotential(model; kwargs...)
end

@inline function ∂ₙφmol(
    ξorΞ::Union{Vector{T}, <: AbstractVector{Triangle{T}}, <: Base.Generator},
    charges::AbstractVector{Charge{T}}
) where T
    Base.depwarn(
        "The unscaled, model-less molecular potential function `∂ₙφmol` is deprecated " *
        "and will be removed in a future release. No replacement is currently planned!",
        :∂ₙφmol
    )
    _molpotential_dn(ξorΞ, charges)
end

@inline function ∂ₙφmol(model::Model{T, Triangle{T}}) where T
    Base.depwarn(
        "The unscaled, model-less molecular potential function `∂ₙφmol` is deprecated " *
        "and will be removed in a future release. No replacement is currently planned!",
        :∂ₙφmol
    )
    _molpotential_dn(model.elements, model.charges)
end

@inline function ∇φmol(
    ξ::Vector{T},
    charges::AbstractVector{Charge{T}}
) where T
    Base.depwarn(
        "The unscaled, model-less molecular potential function `∇φmol` is deprecated " *
        "and will be removed in a future release. No replacement is currently planned!",
        :∇φmol
    )
    -sum(q.val .* (ξ .- q.pos) ./ euclidean(ξ, q.pos)^3 for q in charges; init = zeros(T, 3))
end

@inline function ∇φmol(
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    charges::AbstractVector{Charge{T}}
) where T
    Base.depwarn(
        "The unscaled, model-less molecular potential function `∇φmol` is deprecated " *
        "and will be removed in a future release. No replacement is currently planned!",
        :∇φmol
    )
    collect(Vector{T}, ∇φmol(ξ, charges) for ξ in Ξ)
end
