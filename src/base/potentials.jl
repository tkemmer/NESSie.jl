# =========================================================================================
for T in [:PotentialType, :SingleLayer, :DoubleLayer]
    @eval @doc """
        abstract type PotentialType
        type SingleLayer <: PotentialType end
        type DoubleLayer <: PotentialType end

    Enum-like representation of single and double layer potentials
    """ $T
end
abstract type PotentialType end
type SingleLayer <: PotentialType end
type DoubleLayer <: PotentialType end


# =========================================================================================
for T in [:LocalityType, :LocalES, :NonlocalES]
    @eval @doc """
        abstract type PotentialType
        type SingleLayer <: PotentialType end
        type DoubleLayer <: PotentialType end

    Enum-like representation of locality assumption:
     * *Local electrostatics*:
       Complete independence of solvent molecules
     * *Nonlocal electrostatics*:
       Allow solvent molecule correlation effects (with area-of-effect radius λ)
    """ $T
end
abstract type LocalityType end
type NonlocalES <: LocalityType end
type LocalES <: LocalityType end


# =========================================================================================
"""
    φmol{T}(ξ::Vector{T}, charges::Vector{Charge{T}})

Computes and returns the molecular potential of the given system of point charges in a
structureless medium for the given observation point ξ:

```math
φ\_{mol}(ξ) = \\frac{1}{4π ε\_0 ε\_Ω} \\sum_i \\frac{qᵢ}{|rᵢ-ξ|}
```

!!! note
    The return value is premultiplied by ``4π ⋅ ε₀ ⋅ ε\_Ω``

## Return type
`T`

## Aliases
    φmol{T}(Ξ::Vector{Vector{T}}, charges::Vector{Charge{T}})

Computes the molecular potentials for a list of observation points.

    φmol{T}(model::SurfaceModel{T})

Computes the molecular potentials for the given surface model, using each triangle center
as observation point.
"""
function φmol{T}(ξ::Vector{T}, charges::Vector{Charge{T}})
    # TODO devectorize!
    sum([q.val / euclidean(ξ, q.pos) for q in charges])
end
φmol{T}(Ξ::Vector{Vector{T}}, charges::Vector{Charge{T}}) = [φmol(ξ, charges) for ξ in Ξ]
φmol{T}(model::SurfaceModel{T}) = [φmol(ξ.center, model.charges) for ξ in model.elements]


# =========================================================================================
"""
    ∂ₙφmol{T}(ξ::Triangle{T}, charges::Vector{Charge{T}})

Computes and returns the normal derivative of the given system's molecular potential in a
structureless medium, using the given triangle's center as observation point and the
triangle's normal as reference normal.

```math
∂ₙφ\_{mol}(ξ) = -\\frac{1}{4π ε\_0 ε\_Ω} \\sum_i \\frac{qᵢ}{|rᵢ-ξ|³} (rᵢ-ξ) ⋅ n
```

!!! note
    The return value is premultiplied by ``4π ⋅ ε₀ ⋅ ε\_Ω``

## Return type
`T`

## Aliases
    ∂ₙφmol{T}(model::SurfaceModel{T})

Computes the normal derivatives of the molecular potentials for the given surface model,
using each triangle center and normal as observation point.
"""
function ∂ₙφmol{T}(ξ::Triangle{T}, charges::Vector{Charge{T}})
    # TODO devectorize!
    - sum([q.val * ddot(ξ.center, q.pos, ξ.normal) /
        euclidean(ξ.center, q.pos)^3 for q in charges])
end
∂ₙφmol{T}(model::SurfaceModel{T}) = [∂ₙφmol(ξ, model.charges) for ξ in model.elements]


# =========================================================================================
"""
    ∇φmol{T}(ξ::Vector{T}, charges::Vector{Charge{T}})

Computes and returns the gradient of the given system's molecular potential in a
structureless medium for the given observation point ξ.

```math
∇φ\_{mol}(ξ) = -\\frac{1}{4π ε\_0 ε\_Ω} \\sum_i \\frac{qᵢ}{|rᵢ-ξ|³} (rᵢ-ξ)
```

!!! note
    The return value is premultiplied by ``4π ⋅ ε₀ ⋅ ε\_Ω``

## Return type
`Vector{T}`

## Aliases
    ∇φmol{T}(Ξ::Vector{Vector{T}})

Computes the molecular potential gradients for a list of observation points.
"""
function ∇φmol{T}(ξ::Vector{T}, charges::Vector{Charge{T}})
    # TODO devectorize!
    -sum([q.val * (ξ - q.pos) / euclidean(ξ, q.pos)^3 for q in charges])
end
∇φmol{T}(Ξ::Vector{Vector{T}}, charges::Vector{Charge{T}}) = [∇φmol(ξ, charges) for ξ in Ξ]
