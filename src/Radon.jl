module Radon

using ..NESSie
using ..NESSie: ddot
using Distances: euclidean

export regularyukawacoll, regularyukawacoll!


# =========================================================================================
"""
    regularyukawapot{T}(
        x     ::DenseArray{T,1},
        ξ     ::Vector{T},
        yukawa::T,
              ::Vector{T}=T[]
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
        x     ::DenseArray{T,1},
        ξ     ::Vector{T},
        yukawa::T,
              ::Vector{T}=T[]
    ) where T
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= 1e-10 && return -yukawa

    scalednorm = yukawa * rnorm

    # guard against cancellation
    if scalednorm < .1
        # use alternating series to approximate
        # e^(-c) - 1 = Σ((-c)^i / i!) for i=1 to ∞
        term = -scalednorm
        tolerance = 1e-10 * abs(term)
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
        x     ::Vector{T},
        ξ     ::Vector{T},
        yukawa::T,
        normal::Vector{T}
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
        x     ::Vector{T},
        ξ     ::Vector{T},
        yukawa::T,
        normal::Vector{T}
    ) where T
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= 1e-10 && return zero(T)

    cosovernorm2 = ddot(x, ξ, normal) / rnorm^3
    scalednorm = yukawa * rnorm

    # guard against cancellation
    if scalednorm < .1
        # use alternating series to approximate
        # 1 - (c+1)e^(-c) = Σ((-c)^i * (i-1) / i!) for i=2 to ∞
        term = scalednorm * scalednorm / 2
        tolerance = 1e-10 * abs(term)
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
    setcubpts!{T}(
        dest::Vector{Vector{T}},
        qpts::QuadPts2D{T},
        elem::Triangle{T}
    )

Prepare cubature points for one surface element.

# Return type
`Void`
"""
function setcubpts!(dest::Vector{Vector{T}}, qpts::QuadPts2D{T}, elem::Triangle{T}) where T
    u = elem.v2 - elem.v1
    v = elem.v3 - elem.v1

    # devectorized version of
    # cubpts = [u * qpts.x[i] + v * qpts.y[i] + elem.v1 for i in 1:7]
    for i in 1:qpts.num, j in 1:3
        dest[i][j] = qpts.x[i] * u[j] + qpts.y[i] * v[j] + elem.v1[j]
    end
    nothing
end


# =========================================================================================
"""
    radoncoll!{T}(
            dest    ::DenseArray{T,1},
            elements::Vector{Triangle{T}},
            Ξ       ::Vector{Vector{T}},
            solution::Function,
            yukawa  ::T,
            fvals   ::DenseArray{T,1}
    )

    radoncoll!{T}(
            dest    ::DenseArray{T,2},
            elements::Vector{Triangle{T}},
            Ξ       ::Vector{Vector{T}},
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
        dest    ::DenseArray{T,1},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}},
        solution::Function,
        yukawa  ::T,
        fvals   ::DenseArray{T,1}
    ) where T
    #=== MEMORY-CRITICAL CODE! ===#
    @assert length(fvals) == length(elements)
    @assert length(dest) == length(Ξ)

    # pre-allocate memory for cubature points
    qpts = quadraturepoints(Triangle{T})
    cubpts = [zeros(T, 3) for _ in 1:qpts.num]

    @inbounds for eidx in 1:length(elements)
        elem = elements[eidx]
        area = 2 * elem.area
        setcubpts!(cubpts, qpts, elem)

        for oidx in 1:length(Ξ)
            ξ = Ξ[oidx]
            value = zero(T)
            for i in 1:qpts.num
                value += solution(cubpts[i], ξ, yukawa, elem.normal) * qpts.weight[i]
            end
            dest[oidx] += value * area * fvals[eidx]
        end
    end
    nothing
end

function radoncoll!(
        dest    ::DenseArray{T,2},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}},
        solution::Function,
        yukawa  ::T
    ) where T
    #=== MEMORY-CRITICAL CODE! ===#
    @assert size(dest) == (length(Ξ), length(elements))

    # pre-allocate memory for cubature points
    qpts = quadraturepoints(Triangle{T})
    cubpts = [zeros(T, 3) for _ in 1:qpts.num]

    @inbounds for eidx in 1:length(elements)
        elem = elements[eidx]
        area = 2 * elem.area
        setcubpts!(cubpts, qpts, elem)

        for oidx in 1:length(Ξ)
            ξ = Ξ[oidx]
            value = zero(T)
            for i in 1:qpts.num
                value += solution(cubpts[i], ξ, yukawa, elem.normal) * qpts.weight[i]
            end
            dest[oidx, eidx] = value * area
        end
    end
    nothing
end

# TODO
function radoncoll(
        ξ       ::Vector{T},
        elem    ::Triangle{T},
        yukawa  ::T,
        solution::Function
    ) where T

    qpts = quadraturepoints(Triangle{T})
    cubpts = [zeros(T, 3) for _ in 1:qpts.num]
    area = 2 * elem.area
    setcubpts!(cubpts, qpts, elem)

    value = zero(T)
    for i in 1:qpts.num
        value += solution(cubpts[i], ξ, yukawa, elem.normal)::T * qpts.weight[i]
    end
    value * area
end


# ========================================================================================
"""
    regularyukawacoll!{T, P <: PotentialType}(
                ::Type{P},
        dest    ::DenseArray{T,1},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}},
        yukawa  ::T,
        fvals   ::Union{DenseArray{T,1},SubArray{T,1}}
    )

    regularyukawacoll!{T, P <: PotentialType}(
                ::Type{P},
        dest    ::DenseArray{T,2},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}},
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
regularyukawacoll!(
            ::Type{SingleLayer},
    dest    ::DenseArray{T,1},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}},
    yukawa  ::T,
    fvals   ::Union{DenseArray{T,1},SubArray{T,1}}
) where T = radoncoll!(dest, elements, Ξ, regularyukawapot, yukawa, fvals)

regularyukawacoll!(
            ::Type{DoubleLayer},
    dest    ::DenseArray{T,1},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}},
    yukawa  ::T,
    fvals   ::Union{DenseArray{T,1},SubArray{T,1}}
) where T = radoncoll!(dest, elements, Ξ, ∂ₙregularyukawapot, yukawa, fvals)

regularyukawacoll!(
            ::Type{SingleLayer},
    dest    ::DenseArray{T,2},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}},
    yukawa  ::T
) where T = radoncoll!(dest, elements, Ξ, regularyukawapot, yukawa)

regularyukawacoll!(
            ::Type{DoubleLayer},
    dest    ::DenseArray{T,2},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}},
    yukawa  ::T
) where T = radoncoll!(dest, elements, Ξ, ∂ₙregularyukawapot, yukawa)

regularyukawacoll(
          ::Type{SingleLayer},
    ξ     ::Vector{T},
    elem  ::Triangle{T},
    yukawa::T
) where T = radoncoll(ξ, elem, yukawa, regularyukawapot)

regularyukawacoll(
          ::Type{DoubleLayer},
    ξ     ::Vector{T},
    elem  ::Triangle{T},
    yukawa::T
) where T = radoncoll(ξ, elem, yukawa, ∂ₙregularyukawapot)

end # module
