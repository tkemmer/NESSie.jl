module Rjasanow

using ..NESSie
using ..NESSie: _etol, cos, cathetus, sign, distance

using LinearAlgebra: norm

export laplacecoll, laplacecoll!


# =========================================================================================
for T in [:ObservationPosition, :InPlane, :InSpace]
    @eval @doc """
        abstract type ObservationPosition end
        struct InPlane <: ObservationPosition end
        stryct InSpace <: ObservationPosition end

    Enum-like representation of the obseration point's position relative to the
    corresponding surface element
    """ $T
end
abstract type ObservationPosition end
struct InPlane <: ObservationPosition end
struct InSpace <: ObservationPosition end


# =========================================================================================
"""
    function laplacepot{T, P <: PotentialType}(
        ptype::Type{P},
        ξ    ::Vector{T},
        elem ::Triangle{T},
        dist ::T
    )

Computes the single or double layer Laplace potential of the given triangle for the given
observation point `ξ`. The latter needs to be projected onto the surface element plane
[[Rja90]](@ref Bibliography).

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `dist` Distance from the original `ξ` to the surface element plane

# Return type
`T`
"""
@inline function laplacepot(
        ptype::Type{P},
        ξ    ::Vector{T},
        elem ::Triangle{T},
        dist ::T
    ) where {T, P <: PotentialType}
    laplacepot(ptype, ξ, elem.v1, elem.v2, elem.normal, dist) +
    laplacepot(ptype, ξ, elem.v2, elem.v3, elem.normal, dist) +
    laplacepot(ptype, ξ, elem.v3, elem.v1, elem.normal, dist)
end


# =========================================================================================
"""
    laplacepot{T, P <: PotentialType}(
        ptype::Type{P},
        ξ::Vector{T},
        x1::Vector{T},
        x2::Vector{T},
        normal::Vector{T},
        dist::T
    )

Computes the single or double layer Laplace potential of the triangle defined by the
observation point `ξ` and two nodes `x1` and `x2` of the surface element. `x2` is required
to be `x1`'s next neighbor in counterclockwise direction. Also, `ξ` needs to be projected
onto the surface element plane  [[Rja90]](@ref Bibliography).

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `normal` Unit normal vector of the surface element
 * `dist` Distance from the original `ξ` to the surface element plane

# Return type
`T`
"""
function laplacepot(
        ptype ::Type{P},
        ξ     ::Vector{T},
        x1    ::Vector{T},
        x2    ::Vector{T},
        normal::Vector{T},
        dist ::T
    ) where {T, P <: PotentialType}
    # Construct triangle sides. Later on, we will use the height h of the triangle at point
    # ξ as well as the angles φ1 abd φ2 between h and the triangle sides extending from ξ.
    u1 = x1 - ξ
    u2 = x2 - ξ
    v  = x2 - x1

    # Compute side lengths
    u1norm = norm(u1)
    u2norm = norm(u2)
    vnorm  = norm(v)

    # Compute the sine of the angle φ1 (φ2) between u1 (u2) and the height, that is, the
    # cosine of the angle θ1 (θ2) between u1 (u2) and v if θ1 (θ2) is an acute or right
    # angle. Otherwise, compute the sine of -φ1 (-φ2), which is again the cosine of θ1
    # (θ2) in these cases.
    # Note that the equation (as given above) uses a polar coordinate system with ξ being
    # the pole and h giving the polar axis. The negative angles are needed whenever the
    # corresponding triangle side lies below the polar axis.
    sinφ1 = max(-one(T), min(one(T), cos(u1, v, u1norm, vnorm)))
    sinφ2 = max(-one(T), min(one(T), cos(u2, v, u2norm, vnorm)))

    # Compute the height of the triangle
    h = cathetus(u1norm, sinφ1)

    # Degenerate triangles
    if max(zero(T), h) < _etol(T) ||      # ξ on the line through v1 and v2 (|φ1| = |φ2| = π/2)
        one(T) - abs(sinφ1) < _etol(T) || # ξ on the line through v1 and v2 (φ1 = π/2)
        one(T) - abs(sinφ2) < _etol(T) || # ξ on the line through v1 and v2 (φ2 = π/2)
        abs(sinφ1 - sinφ2) < _etol(T)     # v1 and v2 on the line through ξ (φ1 = φ2 = 0)
        return zero(T)
    end

    # Since the observation point (projection) lies in the same plane as the surface
    # element, we can decide whether the result of this function is to be added or
    # subtracted from the whole surface triangle's Laplace potential by checking on which
    # side of v the observation point lies. This is equivalent to checking whether the
    # normal of the triangle here and the one of the surface element (which are obviously
    # (anti)parallel) are oriented alike.
    pot = abs(dist) < _etol(T) ?
        laplacepot(ptype, InPlane, sinφ1, sinφ2, h, dist) :
        laplacepot(ptype, InSpace, sinφ1, sinφ2, h, dist)
    sign(u1, u2, normal) * pot
end


# =========================================================================================
"""
    laplacepot{T, P <: PotentialType, O <: ObservationPosition}(
             ::Type{P},
             ::Type{O},
        sinφ1::T,
        sinφ2::T,
        h    ::T,
        d    ::T
    )

Computes the Laplace potential (or its normal derivative) of the triangle with the given
height `h` at the observation point `ξ` (projected onto the surface element plane) and the
sines of the angles ``φ₁`` and ``φ₂`` between `h` and the triangle sides extending from `ξ`
[[Rja90]](@ref Bibliography).

!!! note
    The result is premultiplied by 4π.

# Arguments
 * `d` Distance from the original `ξ` to the surface element plane

# Return type
`T`
"""
@inline function laplacepot(
             ::Type{SingleLayer},
             ::Type{InPlane},
        sinφ1::T,
        sinφ2::T,
        h    ::T,
        d    ::T
    ) where T
    #=
        h/8π * [ln((1 + sin(φ))/(1 - sin(φ)))]   from φ1 to φ2
        = h/8π * [ln((1 + sin(φ2))/(1 - sin(φ2))) - ln((1 + sin(φ1))/(1 - sin(φ1)))]
        = h/8π * ln((1 + sin(φ2))(1 - sin(φ1))/((1 - sin(φ2))(1 + sin(φ1))))

        The equation in the original source misses a factor of 0.5!
    =#
    2 \ h * log((1+sinφ2) * (1-sinφ1) / ((1-sinφ2) * (1+sinφ1)))
end

@inline function laplacepot(
             ::Type{SingleLayer},
             ::Type{InSpace},
        sinφ1::T,
        sinφ2::T,
        h    ::T,
        d    ::T
    ) where T
    #=
        h/8π * <1> + d/4π * <2>

        <1>: ln((√(1 - χ² sin²(φ)) + √(1 - χ²) sin(φ))
             / (√(1 - χ² sin²(φ)) - √(1 - χ²) sin(φ)))
             from φ1 to φ2
             = ln((√(1 - χ² sin²(φ2)) + √(1 - χ²) sin(φ2))
               * (√(1 - χ² sin²(φ1)) - √(1 - χ²) sin(φ1))
               / (√(1 - χ² sin²(φ2)) - √(1 - χ²) sin(φ2))
               * √(1 - χ² sin²(φ1)) + √(1 - χ²) sin(φ1))

        <2>: arcsin(χ sin(φ)) - φ   from φ1 to φ2
             = arcsin(χ sin(φ2)) - arcsin(sin(φ2)) - arcsin(χ sin(φ1)) + arcsin(sin(φ1))

        with χ = d / √(d² + h²).
    =#

    d  = abs(d)
    χ2 = d^2 / (d^2 + h^2)
    χ  = √χ2

    # h/8π * <1>
    result = 2 \ h * log(logterm(χ2, sinφ2) / logterm(χ2, sinφ1))
    # + d/4π * <2>
    result + d * (asin(χ * sinφ2) - asin(sinφ2) - asin(χ * sinφ1) + asin(sinφ1))
end

@inline function laplacepot(
             ::Type{DoubleLayer},
             ::Type{InPlane},
        sinφ1::T,
        sinφ2::T,
        h    ::T,
        d    ::T
    ) where T
    #=
        1/4π * ∫-1/|ξ-r'|^3 * (ξ-r')⋅n dr'
        = 1/4π * ∫-1/|ξ-r'|^3 * 0 dr'
        = 0
    =#
    zero(T)
end

@inline function laplacepot(
             ::Type{DoubleLayer},
             ::Type{InSpace},
        sinφ1::T,
        sinφ2::T,
        h    ::T,
        d    ::T
    ) where T
    χ  = abs(d) / √(d^2 + h^2)
    sign(d) * (asin(χ * sinφ1) - asin(sinφ1) - asin(χ * sinφ2) + asin(sinφ2))
end


# =========================================================================================
"""
    laplacecoll!{T, P <: PotentialType}(
        ptype   ::Type{P},
        dest    ::DenseArray{T,1},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}},
        fvals   ::Union{DenseArray{T,1},SubArray{T,1}}
    )

    laplacecoll!{T, P <: PotentialType}(
        ptype   ::Type{P},
        dest    ::DenseArray{T, 2},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}}
    )

Analytical solution for the single or double layer Laplace potential for a given list of
triangles and observation points `Ξ` [[Rja90]](@ref Bibliography).

The first version of this function uses a vector as destination `dest`, where each element
represents the dot product of the corresponding coefficient matrix row and the `fvals`
vector.

!!! note
    The result is premultiplied by 4π.

# Return type
`Void`
"""
function laplacecoll!(
        ptype   ::Type{P},
        dest    ::DenseArray{T,1},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}},
        fvals   ::Union{DenseArray{T,1},SubArray{T,1}}
    ) where {T, P <: PotentialType}
    @assert length(dest) == length(Ξ)
    @assert length(fvals) == length(elements)

    @inbounds for eidx in 1:length(elements)
        elem = elements[eidx]
        for oidx in 1:length(Ξ)
            ξ, dist = projectξ(Ξ[oidx], elem)
            dest[oidx] += laplacepot(ptype, ξ, elem, dist) * fvals[eidx]
        end
    end
    nothing
end

function laplacecoll!(
        ptype   ::Type{P},
        dest    ::DenseArray{T, 2},
        elements::Vector{Triangle{T}},
        Ξ       ::Vector{Vector{T}},
    ) where {T, P <: PotentialType}
    @assert size(dest) == (length(Ξ), length(elements))

    @inbounds for eidx in 1:length(elements)
        elem = elements[eidx]
        for oidx in 1:length(Ξ)
            ξ, dist = projectξ(Ξ[oidx], elem)
            dest[oidx, eidx] = laplacepot(ptype, ξ, elem, dist)
        end
    end
    nothing
end


# =========================================================================================
"""
    laplacecoll{T, P <: PotentialType}(
        ptype::Type{P},
        ξ    ::Vector{T},
        elem ::Vector{Triangle{T}}
    )

Analytical solution for the single or double layer Laplace potential for a given triangle
and observation point `ξ` [[Rja90]](@ref Bibliography).

!!! note
    The result is premultiplied by 4π.

# Return type
`Void`
"""
@inline function laplacecoll(
        ptype::Type{P},
        ξ    ::Vector{T},
        elem ::Triangle{T}
    ) where {T, P <: PotentialType}

    ξ_p, dist = projectξ(ξ, elem)
    laplacepot(ptype, ξ_p, elem, dist)
end


# =========================================================================================
"""
    logterm{T}(χ2::T, sinφ::T)

Utility function to compute
```math
\\frac
{\\sqrt{1 - χ² \\sin²(φ)} + \\sqrt{1 - χ² \\sin(φ)}}
{\\sqrt{1 - χ² \\sin²(φ)} - \\sqrt{1 - χ² \\sin(φ)}}
```

# Return type
`T`
"""
@inline function logterm(χ2::T, sinφ::T) where T
    term1 = √(1 - χ2 * sinφ^2)
    term2 = √(1 - χ2) * sinφ
    (term1 + term2) / (term1 - term2)
end


# =========================================================================================
"""
    projectξ{T}(ξ::Vector{T}, elem::Triangle{T})

Projects ξ onto the surface element plane. Returns the projection and the original distance
to the plane.

# Return type
`Tuple{Vector{T}, T}`
"""
function projectξ(ξ::Vector{T}, elem::Triangle{T}) where T
    dist = distance(ξ, elem)
    abs(dist) < _etol(T) && return (ξ, dist)
    # Devectorized version of ξ -= dist * elem.normal
    ([ξ[i] - dist * elem.normal[i] for i in 1:3], dist)
end

end # module
