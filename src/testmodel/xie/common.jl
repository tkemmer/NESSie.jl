# =========================================================================================
"""
    legendre{T <: AbstractFloat}(
        maxn::Int,
        x   ::T
    )

Precomputes the Legendre polynomials ``Pₙ(x)`` for ``n = 0, 1, ..., \\texttt{maxn}-1`` and
returns a generic function to access the values.

# Return type
`(generic function)`
```julia
(n::Int) -> Pₙ(x)   # return type: T
```

# Example

```jldoctest; setup = :(using NESSie.TestModel: legendre)
julia> p = legendre(3, 0.5);

julia> [p(n) for n in 0:2]    # [P₀(0.5), P₁(0.5), P₂(0.5)]
3-element Vector{Float64}:
  1.0
  0.5
 -0.125
```
"""
function legendre(maxn::Int, x::T) where T <: AbstractFloat
    maxn == 0 && return n::Int -> T[][n+1]
    maxn == 1 && return n::Int -> [one(T)][n+1]

    P = Array{T}(undef, maxn)
    P[1] = one(T) # P₀(x) = 0
    P[2] = x      # P₁(x) = x

    # Use Bonnet's recursion formula for remaining values:
    # Pₙ₊₁(x) = 1/(n+1) ⋅ [(2n + 1)⋅x⋅Pₙ(x) - n⋅Pₙ₋₁}(x)]
    @inbounds for n in 1:maxn-2
        P[n + 2] = (n + 1) \ ((2n + 1) * x * P[n + 1] - n * P[n])
    end
    n::Int -> P[n+1]
end


# =========================================================================================
"""
    spherical_besseli{T <: AbstractFloat}(
        maxn::Int,
        r   ::T
    )

Precomputes the modified spherical Bessel function of the first kind ``iₙ(r)`` for
``n=-1, 0, ..., \\texttt{maxn}`` and returns a generic function to access the values.
``iₙ(r)`` is defined as

```math
iₙ(r) = \\sqrt{\\frac{π}{2r}} I_{n+0.5}(r),
```

where ``I_ν`` is the modified Bessel function of the first kind
[[Xie16]](@ref Bibliography).

# Return type
`(generic function)`
```julia
(n::Int) -> iₙ(r)   # return type: T
```
"""
function spherical_besseli(maxn::Int, r::T) where T <: AbstractFloat
    i = √(π / 2r) * [besseli(n + one(T)/2, r) for n in -1:maxn]
    n::Int -> i[n+2]
end


# =========================================================================================
"""
    spherical_besselk{T <: AbstractFloat}(
        maxn::Int,
        r   ::T
    )

Precomputes the modified spherical Bessel function of the second kind ``kₙ(r)`` for
``n=-1, 0, ..., \\texttt{maxn}`` and returns a generic function to access the values.
``kₙ(r)`` is defined as

```math
kₙ(r) = \\sqrt{\\frac{π}{2r}} K_{n+0.5}(r),
```

where ``K_ν`` is the modified Bessel function of the second kind
[[Xie16]](@ref Bibliography).

# Return type
`(generic function)`
```julia
(n::Int) -> kₙ(r)   # return type: T
```
"""
function spherical_besselk(maxn::Int, r::T) where T <: AbstractFloat
    k = √(π / 2r) * [besselk(n + one(T)/2, r) for n in -1:maxn]
    n::Int -> k[n+2]
end
