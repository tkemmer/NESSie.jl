module Rjasanow

using ..ProteinES
import ..ProteinES: cos, cathetus, sign, distance

export laplacecoll!

#=
    Enum-like representation of the obseration point's position relative to the corresponding
    surface element.
=#
abstract ObservationPosition
type InPlane <: ObservationPosition end
type InSpace <: ObservationPosition end

#=
    Computes the single or double layer Laplace potential of the triangle defined by the
    observation point ξ and two nodes x1 and x2 of the surface element. x2 is required to
    be x1's next neighbor in counterclockwise direction. Also, ξ needs to be projected
    onto the surface element plane.

    Note that the result is premultiplied by 4π!

    @param ptype
        SingleLayer or DoubleLayer
    @param ξ
        Observation point (projection)
    @param x1
        One node of the surface element
    @param x2
        x1's next neighbor in counterclockwise direction
    @param normal
        Unit normal vector of the surface element
    @param dist
        Distance from the original ξ to the surface element plane
    @return T
=#
function laplacepot{T, P <: PotentialType}(ptype::Type{P}, ξ::Vector{T}, x1::Vector{T}, x2::Vector{T}, normal::Vector{T}, dist::T)
    # Construct triangle sides. Later on, we will use the height h of the triangle at point
    # ξ as well as the angles φ1 abd φ2 between h and the triangle sides extending from ξ.
    u1 = x1 - ξ
    u2 = x2 - ξ
    v  = x2 - x1

    # Compute side lengths
    u1norm = vecnorm(u1)
    u2norm = vecnorm(u2)
    vnorm  = vecnorm(v)

    # Compute the sine of the angle φ1 (φ2) between u1 (u2) and the height, that is, the
    # cosine of the angle θ1 (θ2) between u1 (u2) and v if θ1 (θ2) is an acute or right
    # angle. Otherwise, compute the sine of -φ1 (-φ2), which is again the cosine of θ1
    # (θ2) in these cases.
    # Note that the equation (as given above) uses a polar coordinate system with ξ being
    # the pole and h giving the polar axis. The negative angles are needed whenever the
    # corresponding triangle side lies below the polar axis.
    sinφ1 = max(-one(T), min(one(T), cos(u1, u1norm, v, vnorm)))
    sinφ2 = max(-one(T), min(one(T), cos(u2, u2norm, v, vnorm)))

    # Compute the height of the triangle
    h = cathetus(u1norm, sinφ1)

    # Degenerate triangles
    if max(zero(T), h) < 1e-10 ||         # ξ on the line through v1 and v2 (|φ1| = |φ2| = π/2)
        one(T) - abs(sinφ1) < 1e-10 ||    # ξ on the line through v1 and v2 (φ1 = π/2)
        one(T) - abs(sinφ2) < 1e-10 ||    # ξ on the line through v1 and v2 (φ2 = π/2)
        abs(sinφ1 - sinφ2) < 1e-10        # v1 and v2 on the line through ξ (φ1 = φ2 = 0)
        return zero(T)
    end

    # Check whether or not the original ξ lies in the surface element plane
    ξloc = abs(dist) < 1e-10 ? InPlane : InSpace

    # Since the observation point (projection) lies in the same plane as the surface element,
    # we can decide whether the result of this function is to be added or subtracted from the
    # whole surface triangle's Laplace potential by checking on which side of v the observation
    # point lies. This is equivalent to checking whether the normal of the triangle here and
    # the one of the surface element (which are obviously (anti)parallel) are oriented alike.
    sign(u1, u2, normal) * laplacepot(ptype, ξloc, sinφ1, sinφ2, h, dist)
end

#=
    Computes the Laplace potential of the triangle with the given height h at the observation
    point ξ and the sines of the angles φ1 and φ2 between h and the triangle sides extending
    from ξ.

    h/8π * [ln((1 + sin(φ))/(1 - sin(φ)))]   from φ1 to φ2
    = h/8π * [ln((1 + sin(φ2))/(1 - sin(φ2))) - ln((1 + sin(φ1))/(1 - sin(φ1)))]
    = h/8π * ln((1 + sin(φ2))(1 - sin(φ1))/((1 - sin(φ2))(1 + sin(φ1))))

    Note that the result is premultiplied by 4π! Also note that the equation in [1] misses a
    factor of 0.5.

    References:
    [1] S. Rjasanow. Vorkonditionierte iterative Auflösung von Randelementgleichungen für die
        Dirichlet-Aufgabe. Wissenschaftliche Schriftreihe der Technischen Universität Karl-
        Marx-Stadt, 7/1990.

    @param sinφ1
        Sine of the angle between h and one triangle side extending from ξ
    @param sinφ2
        Sine of the angle between h and the other triangle side extending from ξ
    @param h
        Height of the triangle
    @param d
        Distance from ξ to the plane the original surface element lies in (unused)
    @return T
=#
laplacepot{T}(::Type{SingleLayer}, ::Type{InPlane}, sinφ1::T, sinφ2::T, h::T, d::T) = T(.5) * h * log((one(T)+sinφ2) * (one(T)-sinφ1) / ((one(T)-sinφ2) * (one(T)+sinφ1)))

#=
    Computes the Laplace potential of the triangle with the given height h at the observation
    point ξ and the sines of the angles φ1 and φ2 between h and the triangle sides extending
    from ξ.

    This function projects ξ onto the plane the original surface element lies in before
    computing the potential.

    h/8π * <1> + d/4π * <2>

    <1>: ln((√(1 - χ² sin²(φ)) + √(1 - χ²) sin(φ)) / (√(1 - χ² sin²(φ)) - √(1 - χ²) sin(φ)))
           from φ1 to φ2
         = ln((√(1 - χ² sin²(φ2)) + √(1 - χ²) sin(φ2)) * (√(1 - χ² sin²(φ1)) - √(1 - χ²) sin(φ1))
           / (√(1 - χ² sin²(φ2)) - √(1 - χ²) sin(φ2)) * √(1 - χ² sin²(φ1)) + √(1 - χ²) sin(φ1))

    <2>: arcsin(χ sin(φ)) - φ   from φ1 to φ2
         = arcsin(χ sin(φ2)) - arcsin(sin(φ2)) - arcsin(χ sin(φ1)) + arcsin(sin(φ1))

    with χ = d / √(d² + h²).

    Note that the result is premultiplied by 4π!

    References:
    [1] S. Rjasanow. Vorkonditionierte iterative Auflösung von Randelementgleichungen für die
        Dirichlet-Aufgabe. Wissenschaftliche Schriftreihe der Technischen Universität Karl-
        Marx-Stadt, 7/1990.

    @param sinφ1
        Sine of the angle between h and one triangle side extending from ξ
    @param sinφ2
        Sine of the angle between h and the other triangle side extending from ξ
    @param h
        Height of the triangle
    @param d
        Distance from ξ to the plane the original surface element lies in
    @return T
=#
laplacepot{T}(::Type{SingleLayer}, ::Type{InSpace}, sinφ1::T, sinφ2::T, h::T, d::T) = begin
    d  = abs(d)
    χ  = d / √(d^2 + h^2)
    χ2 = χ^2

    # h/8π * <1>
    result = T(.5) * h * log(logterm(χ2, sinφ2) / logterm(χ2, sinφ1))
    # + d/4π * <2>
    result + d * (asin(χ * sinφ2) - asin(sinφ2) - asin(χ * sinφ1) + asin(sinφ1))
end

#=
    TODO
=#
laplacepot{T}(::Type{DoubleLayer}, ::Type{InPlane}, sinφ1::T, sinφ2::T, h::T, d::T) = zero(T)

#=
    TODO
=#
laplacepot{T}(::Type{DoubleLayer}, ::Type{InSpace}, sinφ1::T, sinφ2::T, h::T, d::T) = begin
    χ  = abs(d) / √(d^2 + h^2)
    sign(d) * (asin(χ * sinφ1) - asin(sinφ1) - asin(χ * sinφ2) + asin(sinφ2))
end

#=
    TODO
=#
laplacepot{T, P <: PotentialType}(ptype::Type{P}, ξ::Vector{T}, elem::Triangle{T}, dist::T) = begin
    laplacepot(ptype, ξ, elem.v1, elem.v2, elem.normal, dist) +
    laplacepot(ptype, ξ, elem.v2, elem.v3, elem.normal, dist) +
    laplacepot(ptype, ξ, elem.v3, elem.v1, elem.normal, dist)
end

#=
    Generates a single or double layer Laplace potential matrix/vector using equations by
    S. Rjasanow. If `dest` is a vector, the function values f for each surface triangle have
    to be specified, since each element of the vector represents the dot product of the
    corresponding coefficient matrix row and the `fvals` vector.

    Note that the result is premultiplied by 4π!

    References:
    [1] S. Rjasanow. Vorkonditionierte iterative Auflösung von Randelementgleichungen für die
        Dirichlet-Aufgabe. Wissenschaftliche Schriftreihe der Technischen Universität Karl-
        Marx-Stadt, 7/1990.

    @param ptype
        SingleLayer or DoubleLayer
    @param dest
        Destination matrix/vector
    @param elements
        Surface elements
    @param ξlist
        Observation points
    @param fvals
        Function values of the elements
=#
function laplacecoll!{T, P <: PotentialType}(ptype::Type{P}, dest::Union{DenseArray{T,1}, DenseArray{T, 2}}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, fvals::DenseArray{T,1}=T[])
    isvec  = isa(dest, DenseArray{T, 1})
    isvec && @assert length(dest) == length(ξlist)
    isvec && @assert length(fvals) == length(elements)
    isvec || @assert size(dest) == (length(ξlist), length(elements))

    @inbounds for (eidx, elem) in enumerate(elements), (oidx, ξ) in enumerate(ξlist)

        #TODO check whether zerodiag is necessary
        if !isvec && ptype == DoubleLayer && eidx == oidx
            dest[oidx, eidx] = zero(T)
            continue
        end

        dist = distance(ξ, elem)

        # Project ξ onto the surface element plane if necessary
        if abs(dist) >= 1e-10
            ξ = ξ - dist .* elem.normal
        end

        isvec ?
            dest[oidx] += laplacepot(ptype, ξ, elem, dist) * fvals[eidx] :
            dest[oidx, eidx] = laplacepot(ptype, ξ, elem, dist)
    end
    nothing
end
laplacecoll!{T, P <: PotentialType}(ptype::Type{P}, dest::DenseArray{T,1}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, fvals::DenseArray{T, 1}, ::Option{T}) = laplacecoll!(ptype, dest, elements, ξlist, fvals)
laplacecoll!{T, P <: PotentialType}(ptype::Type{P}, dest::DenseArray{T,2}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, ::Option{T}) = laplacecoll!(ptype, dest, elements, ξlist)


#=
    Helper function to compute
    √(1 - χ2 * sinφ^2) + √(1 - χ2) * sinφ / (√(1 - χ2 * sinφ^2) - √(1 - χ2) * sinφ)

    @param χ2
    @param sinφ
    @return T
=#
logterm{T}(χ2::T, sinφ::T) = begin
    term1 = √(one(T) - χ2 * sinφ^2)
    term2 = √(one(T) - χ2) * sinφ
    (term1 + term2) / (term1 - term2)
end

end # module
