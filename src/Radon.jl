module Radon

using ..NESSie
using ..NESSie: _etol, ddot
using Distances: euclidean

export regularyukawacoll, regularyukawacoll!


# =========================================================================================
"""
    regularyukawapot{T}(
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
function regularyukawapot(
        x     ::AbstractVector{T},
        ξ     ::AbstractVector{T},
        yukawa::T,
              ::AbstractVector{T}=T[]
    ) where T
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= _etol(T) && return -yukawa

    scalednorm = yukawa * rnorm

    # guard against cancellation
    if scalednorm < .1
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
    ∂ₙregularyukawapot{T}(
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
function ∂ₙregularyukawapot(
        x     ::AbstractVector{T},
        ξ     ::AbstractVector{T},
        yukawa::T,
        normal::AbstractVector{T}
    ) where T
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= _etol(T) && return yukawa^2 / 2 / T(√3)

    cosovernorm2 = ddot(x, ξ, normal) / rnorm^3
    scalednorm = yukawa * rnorm

    # guard against cancellation
    if scalednorm < .1
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
    radoncoll!{T}(
            dest    ::AbstractVector{T},
            elements::AbstractVector{Triangle{T}},
            Ξ       ::AbstractVector{Vector{T}},
            solution::Function,
            yukawa  ::T,
            fvals   ::AbstractVector{T}
    )

    radoncoll!{T}(
            dest    ::AbstractMatrix{T},
            elements::AbstractVector{Triangle{T}},
            Ξ       ::AbstractVector{Vector{T}},
            solution::Function,
            yukawa  ::T
    )

Seven-point Radon cubature [[Rad48]](@ref Bibliography) for a given function and a list of
triangles and observation points `Ξ`. If `dest` is a vector, the function values f for each
surface triangle have to be specified, since each element of the vector represents the dot
product of the corresponding coefficient matrix row and the `fvals` vector.

If you intend computing single/double layer potentials with this function, you might want
to use the shorthand signature `regularyukawacoll!` instead.

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `solution` Fundamental solution; supported functions: `regularyukawapot`,
   `∂ₙregularyukawapot`
 * `yukawa` [Exponent](@ref int-constants) of the Yukawa operator's fundamental solution

# Return type
`Void`
"""
function radoncoll!(
        dest    ::AbstractVector{T},
        elements::AbstractVector{Triangle{T}},
        Ξ       ::AbstractVector{Vector{T}},
        solution::Function,
        yukawa  ::T,
        fvals   ::AbstractVector{T}
    ) where T
    #=== MEMORY-CRITICAL CODE! ===#
    @assert length(fvals) == length(elements)
    @assert length(dest) == length(Ξ)

    # pre-allocate memory for cubature points
    cubpts = quadraturepoints(elements)

    @inbounds for eidx in eachindex(cubpts)
        elem = cubpts[eidx].elem
        qpts = cubpts[eidx].qpts
        weig = cubpts[eidx].weights
        area = 2 * elem.area

        for oidx in eachindex(Ξ)
            ξ = Ξ[oidx]
            value = zero(T)
            for i in eachindex(weig)
                value += solution(view(qpts, :, i), ξ, yukawa, elem.normal) * weig[i]
            end
            dest[oidx] += value * area * fvals[eidx]
        end
    end
    nothing
end

function radoncoll!(
        dest    ::AbstractMatrix{T},
        elements::AbstractVector{Triangle{T}},
        Ξ       ::AbstractVector{Vector{T}},
        solution::Function,
        yukawa  ::T
    ) where T
    #=== MEMORY-CRITICAL CODE! ===#
    @assert size(dest) == (length(Ξ), length(elements))

    # pre-allocate memory for cubature points
    cubpts = quadraturepoints(elements)

    @inbounds for eidx in eachindex(cubpts)
        elem = cubpts[eidx].elem
        qpts = cubpts[eidx].qpts
        weig = cubpts[eidx].weights
        area = 2 * elem.area

        for oidx in eachindex(Ξ)
            ξ = Ξ[oidx]
            value = zero(T)
            for i in eachindex(weig)
                value += solution(view(qpts, :, i), ξ, yukawa, elem.normal) * weig[i]
            end
            dest[oidx, eidx] = value * area
        end
    end
    nothing
end


# =========================================================================================
"""
    radoncoll{T}(
        ξ       ::AbstractVector{T},
        tquad   ::TriangleQuad{T},
        yukawa  ::T,
        solution::Function
    )

Seven-point Radon cubature [[Rad48]](@ref Bibliography) for a given function and a pair of
surface element and observations point.

If you intend computing single/double layer potentials with this function, you might want
to use the shorthand signature `regularyukawacoll` instead.

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `solution` Fundamental solution; supported functions: `regularyukawapot`,
   `∂ₙregularyukawapot`
 * `yukawa` [Exponent](@ref int-constants) of the Yukawa operator's fundamental solution

# Return type
`Void`
"""
function radoncoll(
        ξ       ::AbstractVector{T},
        tquad   ::TriangleQuad{T},
        yukawa  ::T,
        solution::Function
    ) where T

    area = 2 * tquad.elem.area
    value = zero(T)
    for i in 1:length(tquad.weights)
        value += solution(view(tquad.qpts, :, i), ξ, yukawa, tquad.elem.normal)::T * tquad.weights[i]
    end
    value * area
end


# ========================================================================================
"""
    regularyukawacoll!{T, P <: PotentialType}(
                ::Type{P},
        dest    ::AbstractVector{T},
        elements::AbstractVector{Triangle{T}},
        Ξ       ::AbstractVector{Vector{T}},
        yukawa  ::T,
        fvals   ::AbstractVector{T}
    )

    regularyukawacoll!{T, P <: PotentialType}(
                ::Type{P},
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
`Void`
"""
@inline regularyukawacoll!(
            ::Type{SingleLayer},
    dest    ::AbstractVector{T},
    elements::AbstractVector{Triangle{T}},
    Ξ       ::AbstractVector{Vector{T}},
    yukawa  ::T,
    fvals   ::AbstractVector{T}
) where T = radoncoll!(dest, elements, Ξ, regularyukawapot, yukawa, fvals)

@inline regularyukawacoll!(
            ::Type{DoubleLayer},
    dest    ::AbstractVector{T},
    elements::AbstractVector{Triangle{T}},
    Ξ       ::AbstractVector{Vector{T}},
    yukawa  ::T,
    fvals   ::AbstractVector{T}
) where T = radoncoll!(dest, elements, Ξ, ∂ₙregularyukawapot, yukawa, fvals)

@inline regularyukawacoll!(
            ::Type{SingleLayer},
    dest    ::AbstractMatrix{T},
    elements::AbstractVector{Triangle{T}},
    Ξ       ::AbstractVector{Vector{T}},
    yukawa  ::T
) where T = radoncoll!(dest, elements, Ξ, regularyukawapot, yukawa)

@inline regularyukawacoll!(
            ::Type{DoubleLayer},
    dest    ::AbstractMatrix{T},
    elements::AbstractVector{Triangle{T}},
    Ξ       ::AbstractVector{Vector{T}},
    yukawa  ::T
) where T = radoncoll!(dest, elements, Ξ, ∂ₙregularyukawapot, yukawa)


# ========================================================================================
"""
    regularyukawacoll{T, P <: PotentialType}(
              ::Type{P},
        ξ     ::AbstractVector{T},
        elem  ::TriangleQuad{T},
        yukawa::T
    )

Computes the regular part of the single or double layer Yukawa potential (that is, Yukawa
minus Laplace) using a seven-point Radon cubature [[Rad48]](@ref Bibliography) for a given
triangle and observation point `ξ`.

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `yukawa` [Exponent](@ref int-constants) of the Yukawa operator's fundamental solution

# Return type
`Void`
"""
@inline regularyukawacoll(
          ::Type{SingleLayer},
    ξ     ::AbstractVector{T},
    elem  ::TriangleQuad{T},
    yukawa::T
) where T = radoncoll(ξ, elem, yukawa, regularyukawapot)

@inline regularyukawacoll(
          ::Type{DoubleLayer},
    ξ     ::AbstractVector{T},
    elem  ::TriangleQuad{T},
    yukawa::T
) where T = radoncoll(ξ, elem, yukawa, ∂ₙregularyukawapot)

end # module
