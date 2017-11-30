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
    φmol{T}(
        ξ        ::Vector{T},
        charges  ::Vector{Charge{T}};
        # kwargs
        tolerance::T                 = T(1e-10)
    )

Computes and returns the molecular potential of the given system of point charges in a
structureless medium for the given observation point ξ:

```math
φ\_{mol}(ξ) = \\frac{1}{4π ε\_0 ε\_Ω} \\sum_i \\frac{qᵢ}{|rᵢ-ξ|}
```

If ``|rᵢ-ξ|`` is smaller than the given `tolerance`, the value is replaced by `tolerance`
for the affected charge.

!!! note
    The return value is premultiplied by ``4π ⋅ ε₀ ⋅ ε\_Ω``

## Return type
`T`

## Aliases
    φmol{T}(
        Ξ        ::Vector{Vector{T}},
        charges  ::Vector{Charge{T}};
        # kwargs
        tolerance::T                 = T(1e-10)
    )

Computes the molecular potentials for a list of observation points.

    φmol{T}(
        model    ::Model{T, Triangle{T}};
        # kwargs
        tolerance::T                     = T(1e-10)
    )

Computes the molecular potentials for the given surface model, using each triangle center
as observation point.
"""
function φmol(
        ξ        ::Vector{T},
        charges  ::Vector{Charge{T}};
        tolerance::T=T(1e-10)
    ) where T
    # devectorized version of
    # sum([q.val / euclidean(ξ, q.pos) for q in charges])
    ret = zero(T)
    for q in charges
        ret += q.val / max(euclidean(ξ, q.pos), tolerance)
    end
    ret
end

function φmol(
        Ξ        ::Vector{Vector{T}},
        charges  ::Vector{Charge{T}};
        tolerance::T=T(1e-10)
    ) where T
    [φmol(ξ, charges, tolerance=tolerance) for ξ in Ξ]
end

function φmol(
        model    ::Model{T, Triangle{T}};
        tolerance::T=T(1e-10)
    ) where T
    [φmol(ξ.center, model.charges, tolerance=tolerance) for ξ in model.elements]
end


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
    ∂ₙφmol{T}(model::Model{T, Triangle{T}})

Computes the normal derivatives of the molecular potentials for the given surface model,
using each triangle center and normal as observation point.
"""
function ∂ₙφmol(
        ξ      ::Triangle{T},
        charges::Vector{Charge{T}}
    ) where T
    ret = zero(T)
    for q in charges
        ret -= q.val * ddot(ξ.center, q.pos, ξ.normal) / euclidean(ξ.center, q.pos)^3
    end
    ret
end

function ∂ₙφmol(model::Model{T, Triangle{T}}) where T
    [∂ₙφmol(ξ, model.charges) for ξ in model.elements]
end


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
function ∇φmol(
        ξ      ::Vector{T},
        charges::Vector{Charge{T}}
    ) where T
    # TODO devectorize!
    -sum([q.val * (ξ - q.pos) / euclidean(ξ, q.pos)^3 for q in charges])
end

function ∇φmol(
        Ξ      ::Vector{Vector{T}},
        charges::Vector{Charge{T}}
    ) where T
    [∇φmol(ξ, charges) for ξ in Ξ]
end
