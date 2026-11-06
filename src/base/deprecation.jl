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

# deprecated in v1.6 (to be removed in v2.0)
export meshunion

function meshunion(
        model1::Model{T, Tetrahedron{T}},
        model2::Model{T, Tetrahedron{T}}
    ) where T
    Base.depwarn(
        "`meshunion` is deprecated and will be removed in the next major version (v2.0). " *
        "Please update your code to use alternative functionality before upgrading.",
        :meshunion
    )
    # find nodes that are to be replaced (tbr)
    obsolete = model1.nodes ∩ model2.nodes
    indices = something.(indexin(obsolete, model1.nodes)) # Vector{Union{Nothing, Int}} -> Vector{Int}
    tbr = Dict{Vector{T}, Int}(zip(obsolete, indices))

    elements = copy(model1.elements)
    for elem in model2.elements
        v1 = haskey(tbr, elem.v1) ? model1.nodes[tbr[elem.v1]] : elem.v1
        v2 = haskey(tbr, elem.v2) ? model1.nodes[tbr[elem.v2]] : elem.v2
        v3 = haskey(tbr, elem.v3) ? model1.nodes[tbr[elem.v3]] : elem.v3
        v4 = haskey(tbr, elem.v4) ? model1.nodes[tbr[elem.v4]] : elem.v4
        push!(elements, Tetrahedron(v1, v2, v3, v4, elem.domain))
    end

    Model(
        collect(Vector{T}, Set(model2.nodes) ∪ Set(model1.nodes)),
        elements,
        model1.charges ∪ model2.charges,
        model1.params
    )
end
