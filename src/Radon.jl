module Radon

using ..ProteinES
using ..ProteinES: ddot
using Distances: euclidean

export laplacecoll!, regularyukawacoll!


# =========================================================================================
"""
    laplacepot{T}(x::DenseArray{T,1}, ξ::Vector{T})

Computes the Laplace potential:
```math
\\frac{1}{|x-ξ|}
```

!!! note
    The result is premultiplied by 4π.

# Return type
`T`
"""
function laplacepot(x::DenseArray{T,1}, ξ::Vector{T}) where T
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= 1e-10 && return zero(T)

    # TODO check
    # guard against small rnorm
    if rnorm < .1
        # use alternating series to approximate
        # 1/c = Σ((-1)^i * (x-1)^i) for i = 0 to ∞
        term = one(T)
        tolerance = 1e-10
        tsum = zero(T)
        for i in 1:15
            abs(term) < tolerance && break

            tsum += term
            term *= -(rnorm - 1)
        end
        return tsum
    end

    1 / rnorm
end
laplacepot(x::DenseArray{T,1}, ξ::Vector{T}, ::Vector{T}, ::T) where T = laplacepot(x, ξ)


# =========================================================================================
"""
    ∂ₙlaplacepot{T}(x::DenseArray{T,1}, ξ::Vector{T}, normal::Vector{T})

Computes the normal derivative of the Laplace potential:
```math
\\frac{(x - ξ) ⋅ n}{|x-ξ|³}
```

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `normal` Normal unit vector at `x`

# Return type
`T`
"""
function ∂ₙlaplacepot(x::DenseArray{T,1}, ξ::Vector{T}, normal::Vector{T}) where T
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= 1e-10 && return zero(T)

    # TODO check
    # guard against small rnorm
    if rnorm < .1
        # use alternating series to approximate
        # -1/c^3 = -1/2 * Σ((-1)^i * (x-1)^i * (i+1) * (i+2)) for i = 0 to ∞
        term = T(2)
        tolerance = 1e-10
        tsum = zero(T)
        for i in 1:15
            abs(term) < tolerance && break

            tsum += term * (i+1) * (i+2)
            term *= -(rnorm - 1)
        end
        return -2 \ tsum * ddot(x, ξ, normal)
    end

    -1 / rnorm^3 * ddot(x, ξ, normal)
end

function ∂ₙlaplacepot{T}(x::DenseArray{T,1}, ξ::Vector{T}, normal::Vector{T}, ::T)
    ∂ₙlaplacepot(x, ξ, normal)
end


# =========================================================================================
"""
    regularyukawapot{T}(x::DenseArray{T,1}, ξ::Vector{T}, yukawa::T)

Computes the regular part of the Yukawa potential, that is, Yukawa minus Laplace:
```math
\\frac{e\^{Λ |x - ξ|} - 1}{|x-ξ|}
```

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `yukawa` [Exponent](@ref int-constants) of the Yukawa operator's fundamental solution

# Return type
`T`
"""
function regularyukawapot(x::DenseArray{T,1}, ξ::Vector{T}, yukawa::T) where T
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
        tsum = zero(T)     # DON'T EVER USE 0 HERE! Time: x2, Memory: x3
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

function regularyukawapot(x::DenseArray{T,1}, ξ::Vector{T}, ::Vector{T}, yukawa::T) where T
        regularyukawapot(x, ξ, yukawa)
end


# =========================================================================================
"""
    ∂ₙregularyukawapot{T}(
        x     ::Vector{T},
        ξ     ::Vector{T},
        normal::Vector{T},
        yukawa::T
    )

Computes the normal derivative of the regular part of the Yukawa potential, that is, Yukawa
minus Laplace:
```math
\\frac{∂}{∂n} \\frac{e\^{Λ  |x - ξ|} - 1}{|x-ξ|}
= \\frac{1 - (1 - Λ  |x - ξ|)e\^{Λ  |x - ξ|}}{|x-ξ|²} \\frac{(x - ξ) ⋅ n}{|x - ξ|}
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
        normal::Vector{T},
        yukawa::T
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
        tsum = zero(T)  # DON'T EVER USE 0 HERE!
        for i in 2:16
            abs(term #=* (i-1)=#) <= tolerance && continue

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
    const u = elem.v2 - elem.v1
    const v = elem.v3 - elem.v1

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
            fvals   ::DenseArray{T,1};
            # kwargs
            yukawa  ::T=zero(T)
    )

    radoncoll!{T}(
            dest    ::DenseArray{T,2},
            elements::Vector{Triangle{T}},
            Ξ       ::Vector{Vector{T}},
            solution::Function;
            # kwargs
            yukawa  ::T=zero(T)
    )

Seven-point Radon cubature [[Rad48]](@ref Bibliography) for a given function and a list of
triangles and observation points `Ξ`. If `dest` is a vector, the function values f for each
surface triangle have to be specified, since each element of the vector represents the dot
product of the corresponding coefficient matrix row and the `fvals` vector.

If you intend computing single/double layer potentials with this function, you might want
to use the shorthand signatures `laplacecoll!` and `regularyukawacoll!` instead.

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `solution` Fundamental solution; supported functions: `regularyukawapot`,
   `∂ₙregularyukawapot`, `laplacepot`, `∂ₙlaplacepot`
 * `yukawa` [Exponent](@ref int-constants) of the Yukawa operator's fundamental solution

# Return type
`Void`
"""
function radoncoll!(
        dest    ::DenseArray{T,1},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}},
        solution::Function,
        fvals   ::DenseArray{T,1};
        yukawa  ::T=zero(T)
    ) where T
    #=== MEMORY-CRITICAL CODE! ===#
    @assert length(fvals) == length(elements)
    @assert length(dest) == length(Ξ)

    # pre-allocate memory for cubature points
    const qpts = quadraturepoints(Triangle{T})
    cubpts = [zeros(T, 3) for _ in 1:qpts.num]

    @inbounds for (eidx, elem) in enumerate(elements)
        area = 2 * elem.area
        setcubpts!(cubpts, qpts, elem)

        for (oidx, ξ) in enumerate(Ξ)
            value = zero(T)
            for i in 1:qpts.num
                value += solution(cubpts[i], ξ, elem.normal, yukawa) * qpts.weight[i]
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
        solution::Function;
        yukawa  ::T=zero(T)
    ) where T
    #=== MEMORY-CRITICAL CODE! ===#
    @assert size(dest) == (length(Ξ), length(elements))

    # pre-allocate memory for cubature points
    const qpts = quadraturepoints(Triangle{T})
    cubpts = [zeros(T, 3) for _ in 1:qpts.num]

    @inbounds for (eidx, elem) in enumerate(elements)
        area = 2 * elem.area
        setcubpts!(cubpts, qpts, elem)

        for (oidx, ξ) in enumerate(Ξ)
            value = zero(T)
            for i in 1:qpts.num
                value += solution(cubpts[i], ξ, elem.normal, yukawa) * qpts.weight[i]
            end
            dest[oidx, eidx] = value * area
        end
    end
    nothing
end


# ========================================================================================
"""
    laplacecoll!{T, P <: PotentialType}(
                ::Type{P},
        dest    ::DenseArray{T,1},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}},
        fvals   ::Union{DenseArray{T,1},SubArray{T,1}}
    )

    laplacecoll!{T, P <: PotentialType}(
                ::Type{P},
        dest    ::DenseArray{T,2},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}}
    )

Computes the single or double layer Laplace potential using a seven-point Radon cubature
[[Rad48]](@ref Bibliography) for a given list of triangles and observation points `Ξ`.

The first version of this function uses a vector as destination `dest`, where each element
represents the dot product of the corresponding coefficient matrix row and the `fvals`
vector.

!!! note
    The result is premultiplied by 4π.

# Return type
`Void`
"""
laplacecoll!(
            ::Type{SingleLayer},
    dest    ::DenseArray{T,1},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}},
    fvals   ::Union{DenseArray{T,1},SubArray{T,1}}
) where T = radoncoll!(dest, elements, Ξ, laplacepot, fvals)

laplacecoll!(
            ::Type{DoubleLayer},
    dest    ::DenseArray{T,1},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}},
    fvals   ::Union{DenseArray{T,1},SubArray{T,1}}
) where T = radoncoll!(dest, elements, Ξ, ∂ₙlaplacepot, fvals)

laplacecoll!(
            ::Type{SingleLayer},
    dest    ::DenseArray{T,2},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}}
) where T = radoncoll!(dest, elements, Ξ, laplacepot)

laplacecoll!(
            ::Type{DoubleLayer},
    dest    ::DenseArray{T,2},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}}
) where T = radoncoll!(dest, elements, Ξ, ∂ₙlaplacepot)


# ========================================================================================
"""
    regularyukawacoll!{T, P <: PotentialType}(
                ::Type{P},
        dest    ::DenseArray{T,1},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}},
        fvals   ::Union{DenseArray{T,1},SubArray{T,1}}
    )

    regularyukawacoll!{T, P <: PotentialType}(
                ::Type{P},
        dest    ::DenseArray{T,2},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}}
    )

Computes the regular part of the single or double layer Yukawa potential (that is, Yukawa
minus Laplace) using a seven-point Radon cubature [[Rad48]](@ref Bibliography) for a given
list of triangles and observation points `Ξ`.

The first version of this function uses a vector as destination `dest`, where each element
represents the dot product of the corresponding coefficient matrix row and the `fvals`
vector.

!!! note
    The result is premultiplied by 4π.

# Return type
`Void`
"""
regularyukawacoll!(
            ::Type{SingleLayer},
    dest    ::DenseArray{T,1},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}},
    fvals   ::Union{DenseArray{T,1},SubArray{T,1}},
    yukawa  ::T
) where T = radoncoll!(dest, elements, Ξ, regularyukawapot, fvals, yukawa=yukawa)

regularyukawacoll!(
            ::Type{DoubleLayer},
    dest    ::DenseArray{T,1},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}},
    fvals   ::Union{DenseArray{T,1},SubArray{T,1}},
    yukawa  ::T
) where T = radoncoll!(dest, elements, Ξ, ∂ₙregularyukawapot, fvals, yukawa=yukawa)

regularyukawacoll!(
            ::Type{SingleLayer},
    dest    ::DenseArray{T,2},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}},
    yukawa  ::T
) where T = radoncoll!(dest, elements, Ξ, regularyukawapot, yukawa=yukawa)

regularyukawacoll!(
            ::Type{DoubleLayer},
    dest    ::DenseArray{T,2},
    elements::Vector{Triangle{T}},
    Ξ       ::Vector{Vector{T}},
    yukawa  ::T
) where T = radoncoll!(dest, elements, Ξ, ∂ₙregularyukawapot, yukawa=yukawa)

end # module
