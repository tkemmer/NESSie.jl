module Radon

using ..NESSie
using ..NESSie: _etol, ddot
using ChunkSplitters
using Distances: euclidean
using LinearAlgebra

export regularyukawacoll, regularyukawacoll!


# =========================================================================================
"""
    _regularyukawapot(
              ::Type{SingleLayer},
        x     ::AbstractVector{T},
        ξ     ::AbstractVector{T},
        yukawa::T,
              ::AbstractVector{T}=T[]
    )

Computes the regular part of the Yukawa potential, that is, Yukawa minus Laplace:
```math
\\mathcal{G}^Y-\\mathcal{G}^L = \\frac{1}{4π}\\frac{e^{-\\frac{|x - ξ|}{Λ}} - 1}{|x-ξ|}
```

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `yukawa` [Exponent](@ref int-constants) of the Yukawa operator's fundamental solution

# Return type
`T`
"""
function _regularyukawapot(
    ::Type{SingleLayer},
    x::AbstractVector{T},
    ξ::AbstractVector{T},
    yukawa::T,
    ::AbstractVector{T}=T[]
) where T
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm::T = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= _etol(T) && return -yukawa

    scalednorm = yukawa * rnorm

    # guard against cancellation
    if scalednorm < T(.1)
        # use alternating series to approximate
        # e^(-c) - 1 = Σ((-c)^i / i!) for i=1 to ∞
        term = -scalednorm
        tolerance = _etol(T) * abs(term)
        tsum = zero(T)
        for i in 1:15
            abs(term) <= tolerance && break

            tsum += term
            term *= -scalednorm / (i+1)
        end
        return tsum / rnorm
    end

    # no danger of cancellation
    (exp(-scalednorm) - 1) / rnorm
end


# =========================================================================================
"""
    _regularyukawapot(
              ::Type{DoubleLayer},
        x     ::AbstractVector{T},
        ξ     ::AbstractVector{T},
        yukawa::T,
        normal::AbstractVector{T}
    )

Computes the normal derivative of the regular part of the Yukawa potential, that is, Yukawa
minus Laplace:
```math
\\frac{∂}{∂n} \\frac{1}{4π} \\frac{e^{-\\frac{|x - ξ|}{Λ}} - 1}{|x-ξ|}
= \\frac{1}{4π} \\frac{1 - (1 + Λ^{-1} |x - ξ|)e^{-\\frac{|x - ξ|}{Λ}}}{|x-ξ|²}
\\frac{(x - ξ) ⋅ n}{|x - ξ|}
```

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `yukawa` [Exponent](@ref int-constants) of the Yukawa operator's fundamental solution

# Return type
`T`
"""
function _regularyukawapot(
    ::Type{DoubleLayer},
    x::AbstractVector{T},
    ξ::AbstractVector{T},
    yukawa::T,
    normal::AbstractVector{T}
) where T
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm::T = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= _etol(T) && return yukawa^2 / 2 / T(√3)

    cosovernorm2 = ddot(x, ξ, normal) / rnorm^3
    scalednorm = yukawa * rnorm

    # guard against cancellation
    if scalednorm < T(.1)
        # use alternating series to approximate
        # 1 - (c+1)e^(-c) = Σ((-c)^i * (i-1) / i!) for i=2 to ∞
        term = scalednorm * scalednorm / 2
        tolerance = _etol(T) * abs(term)
        tsum = zero(T)
        for i in 2:16
            abs(term) <= tolerance && break

            tsum += term * (i-1)
            term *= -scalednorm / (i+1)
        end
        return tsum * cosovernorm2
    end

    # no danger of cancellation
    (1 - (1 + scalednorm) * exp(-scalednorm)) * cosovernorm2
end


# =========================================================================================
"""
    regularyukawacoll!(
                ::Type{<: PotentialType},
        dest    ::AbstractVector{T},
        elements::AbstractVector{Triangle{T}},
        Ξ       ::AbstractVector{Vector{T}},
        yukawa  ::T,
        fvals   ::AbstractVector{T}
    )

    regularyukawacoll!(
                ::Type{<: PotentialType},
        dest    ::AbstractMatrix{T},
        elements::AbstractVector{Triangle{T}},
        Ξ       ::AbstractVector{Vector{T}},
        yukawa  ::T
    )

Computes the regular part of the single or double layer Yukawa potential (that is, Yukawa
minus Laplace) using a seven-point Radon cubature [[Rad48]](@ref Bibliography) for a given
list of triangles and observation points `Ξ`.

The first version of this function uses a vector as destination `dest`, where each element
represents the dot product of the corresponding coefficient matrix row and the `fvals`
vector.

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `yukawa` [Exponent](@ref int-constants) of the Yukawa operator's fundamental solution

# Return type
`Nothing`
"""
function regularyukawacoll!(
    ptype::Type{<: PotentialType},
    dest::AbstractVector{T},
    elements::AbstractVector{Triangle{T}},
    Ξ::AbstractVector{Vector{T}},
    yukawa::T,
    fvals::AbstractVector{T}
) where T
    cubpts = collect(TriangleQuad{T}, TriangleQuad(elem) for elem in elements)
    tasks = map(index_chunks(Ξ; n = Threads.nthreads(), minsize = 100)) do idx
        Threads.@spawn _regularyukawacoll!(ptype, view(dest, idx), cubpts, view(Ξ, idx), yukawa, fvals)
    end
    wait.(tasks)
    nothing
end

@inline function _regularyukawacoll!(
    ptype::Type{<: PotentialType},
    dest::AbstractVector{T},
    elements::AbstractVector{TriangleQuad{T}},
    Ξ::AbstractVector{Vector{T}},
    yukawa::T,
    fvals::AbstractVector{T}
) where T
    #=== MEMORY-CRITICAL CODE! ===#
    @assert length(fvals) == length(elements)
    @assert length(dest) == length(Ξ)

    dest .+= _regularyukawacoll.(ptype, Ξ, (elements,), yukawa, (fvals,))
    nothing
end

@inline function _regularyukawacoll(
    ptype::Type{<: PotentialType},
    ξ::Vector{T},
    tquads::AbstractVector{TriangleQuad{T}},
    yukawa::T,
    fvals::AbstractVector{T}
) where T
    value = zero(T)
    for i in eachindex(tquads)
        value += regularyukawacoll(ptype, ξ, tquads[i], yukawa) * fvals[i]
    end
    value
end

function regularyukawacoll!(
    ptype::Type{<: PotentialType},
    dest::AbstractMatrix{T},
    elements::AbstractVector{Triangle{T}},
    Ξ::AbstractVector{Vector{T}},
    yukawa::T
) where T
    cubpts = collect(TriangleQuad{T}, TriangleQuad(elem) for elem in elements)
    tasks = map(index_chunks(Ξ; n = Threads.nthreads(), minsize = 100)) do idx
        Threads.@spawn _regularyukawacoll!(ptype, view(dest, idx, :), cubpts, view(Ξ, idx), yukawa)
    end
    wait.(tasks)
    nothing
end

function _regularyukawacoll!(
    ptype::Type{<: PotentialType},
    dest::AbstractMatrix{T},
    cubpts::AbstractVector{TriangleQuad{T}},
    Ξ::AbstractVector{Vector{T}},
    yukawa::T
) where T
    #=== MEMORY-CRITICAL CODE! ===#
    @assert size(dest) == (length(Ξ), length(cubpts))

    @inbounds for (eidx, cubpt) in enumerate(cubpts)
        for (oidx, ξ) in enumerate(Ξ)
            dest[oidx, eidx] = regularyukawacoll(ptype, ξ, cubpt, yukawa)
        end
    end
    nothing
end


# =========================================================================================
"""
    regularyukawacoll(
              ::Type{<: PotentialType},
        ξ     ::AbstractVector{T},
        elem  ::Union{Triangle{T}, TriangleQuad{T}},
        yukawa::T
    )

Computes the regular part of the single or double layer Yukawa potential (that is, Yukawa
minus Laplace) using a seven-point Radon cubature [[Rad48]](@ref Bibliography) for a given
triangle and observation point `ξ`.

!!! note
    The result is premultiplied by 4π.

!!! note
    If you need to call this method repeatedly for the same surface triangle, you should
    consider passing the `elem` parameter as `TriangleQuad` object to avoid expensive
    recomputations.

# Arguments
 * `yukawa` [Exponent](@ref int-constants) of the Yukawa operator's fundamental solution

# Return type
`T`
"""
function regularyukawacoll(
    ptype::Type{<: PotentialType},
    ξ::Vector{T},
    tquad::TriangleQuad{T},
    yukawa::T
) where T
    value = zero(T)
    for (i, qpt) in enumerate(eachcol(tquad.qpts))
        value += _regularyukawapot(ptype, qpt, ξ, yukawa, tquad.elem.normal) * tquad.weights[i]
    end
    value * 2 * tquad.elem.area
end

regularyukawacoll(
    ptype::Type{<: PotentialType},
    ξ::Vector{T},
    elem::Triangle{T},
    yukawa::T
) where T = regularyukawacoll(ptype, ξ, TriangleQuad(elem), yukawa)


# =========================================================================================
# Deprecation

# renamed in v1.5 (unexported but used by CuNESSie's tests)
@deprecate regularyukawapot(x, ξ, yuk) _regularyukawapot(SingleLayer, x, ξ, yuk) false
@deprecate regularyukawapot(x, ξ, yuk, normal) _regularyukawapot(SingleLayer, x, ξ, yuk, normal) false
@deprecate ∂ₙregularyukawapot(x, ξ, yuk, normal) _regularyukawapot(DoubleLayer, x, ξ, yuk, normal) false

end # module
