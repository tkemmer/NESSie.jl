# =========================================================================================
for T in [:PotentialType, :SingleLayer, :DoubleLayer]
    @eval @doc """
        abstract type PotentialType end
        struct SingleLayer <: PotentialType end
        struct DoubleLayer <: PotentialType end

    Enum-like representation of single and double layer potentials
    """ $T
end
abstract type PotentialType end
struct SingleLayer <: PotentialType end
struct DoubleLayer <: PotentialType end


# =========================================================================================
for T in [:LocalityType, :LocalES, :NonlocalES]
    @eval @doc """
        abstract type LocalityType end
        struct NonlocalES <: LocalityType end
        struct LocalES    <: LocalityType end

    Enum-like representation of locality assumption:
     * *Local electrostatics*:
       Complete independence of solvent molecules
     * *Nonlocal electrostatics*:
       Allow solvent molecule correlation effects (with area-of-effect radius λ)
    """ $T
end
abstract type LocalityType end
struct NonlocalES <: LocalityType end
struct LocalES    <: LocalityType end


# =========================================================================================
"""
    molpotential(ξ::Vector{T}, model::Model{T})
    molpotential(Ξ::AbstractVector{Vector{T}}, model::Model{T})

Computes the molecular potential(s) at the given observation point(s) ξ (Ξ) for the given
model.

# Supported keyword arguments
 - `tolerance::T = 1e-10` minimum distance assumed between any observation point and point
   charge. Closer distances are replaced by this value.

# Unit
``V = \\frac{C}{F}``

# Return type
`T` or `Vector{T}`
"""
@inline function molpotential(ξ::Vector{T}, model::Model{T}; kwargs...) where T
    _molpotential(ξ, model.charges; kwargs...) / model.params.εΩ * potprefactor(T)
end

@inline function molpotential(
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    model::Model{T};
    kwargs...
) where T
    collect(T, molpotential(ξ, model; kwargs...) for ξ in Ξ)
end


# =========================================================================================
"""
    _molpotential(ξ::Vector{T}, charges::AbstractVector{Charge{T}})
    _molpotential(Ξ::AbstractVector{Vector{T}}, charges::AbstractVector{Charge{T}})

Computes and returns the molecular potential(s) of the given system of point charges in a
structureless medium for the given observation point(s) ξ (Ξ):

```math
φ_{mol}(ξ) = \\frac{1}{4π ε_0 ε_Ω} \\sum_i \\frac{qᵢ}{|rᵢ-ξ|}
```

# Supported keyword arguments
 - `tolerance::T = 1e-10` minimum distance assumed between any observation point and point
   charge. Closer distances are replaced by this value.

## Return type
`T` or `Vector{T}`

## Aliases
    _molpotential(model::Model{T, Triangle{T}})

Computes the molecular potentials for the given surface model, using each triangle center
as observation point.

!!! note
    The return value is premultiplied by ``4π ⋅ ε₀ ⋅ ε_Ω``
"""
@inline function _molpotential(
    ξ::Vector{T},
    charges::AbstractVector{Charge{T}};
    tolerance::T=T(1e-10)
) where T
    sum(q.val / max(euclidean(ξ, q.pos), tolerance) for q in charges; init = zero(T))
end

@inline function _molpotential(
    Ξ::Union{<: AbstractVector{Vector{T}}, <: Base.Generator},
    charges::AbstractVector{Charge{T}};
    kwargs...
) where T
    collect(T, _molpotential(ξ, charges; kwargs...) for ξ in Ξ)
end

@inline function _molpotential(
    model::Model{T, Triangle{T}};
    kwargs...
) where T
    _molpotential((elem.center for elem in model.elements), model.charges; kwargs...)
end


# =========================================================================================
"""
    _molpotential_dn(ξ::Vector{T}, charges::AbstractVector{Charge{T}})
    _molpotential_dn(Ξ::AbstractVector{Vector{T}}, charges::AbstractVector{Charge{T}})

Computes and returns the normal derivative(s) of the molecular potential(s) of the given
system of point charges in a structureless medium for the given observation point(s) ξ (Ξ):

```math
∂ₙφ_{mol}(ξ) = -\\frac{1}{4π ε_0 ε_Ω} \\sum_i \\frac{qᵢ}{|rᵢ-ξ|³} (rᵢ-ξ) ⋅ n
```

# Supported keyword arguments
 - `tolerance::T = 1e-10` minimum distance cubed assumed between any observation point and
   point charge. Smaller values are replaced by this.

## Return type
`T` or `Vector{T}`

## Aliases
    _molpotential_dn(model::Model{T, Triangle{T}})

Computes the molecular potentials for the given surface model, using each triangle center
as observation point.

!!! note
    The return value is premultiplied by ``4π ⋅ ε₀ ⋅ ε_Ω``
"""
@inline function _molpotential_dn(
    ξ::Triangle{T},
    charges::AbstractVector{Charge{T}};
    tolerance::T = T(1e-10)
) where T
    -sum(
        q.val * ddot(ξ.center, q.pos, ξ.normal) / max(euclidean(ξ.center, q.pos)^3, tolerance)
        for q in charges;
        init = zero(T)
    )
end

@inline function _molpotential_dn(
    Ξ::Union{<: AbstractVector{Triangle{T}}, <: Base.Generator},
    charges::AbstractVector{Charge{T}};
    kwargs...
) where T
    collect(T, _molpotential_dn(ξ, charges; kwargs...) for ξ in Ξ)
end

@inline function _molpotential_dn(model::Model{T, Triangle{T}}; kwargs...) where T
    _molpotential_dn(model.elements, model.charges; kwargs...)
end


# =========================================================================================
# Function stubs for potentials in submodules
function espotential end
function rfpotential end
function rfenergy end
