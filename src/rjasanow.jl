module Rjasanow

import NonLocalBEM: Element, SingleLayer, DoubleLayer, cos, cathetus, sign, dist

export laplacecoll!

#=
    Enum-like representation of the obseration point's position relative to the
    corresponding surface element.
=#
abstract ObservationPosition
type InPlane <: ObservationPosition end
type InSpace <: ObservationPosition end

#=
    Computes a part of the Laplace potential of a surface element given the observation
    point ξ and two nodes x1 and x2 of the surface element as well as its normal vector.
    x2 is required to be x1's next neighbor in counterclockwise direction. Also, ξ is
    required to lie in the same plane as the surface element.

    h/4π * [ln((1 + sin φ)/(1 - sin φ))] from φ1 to φ2
    = h/4π * [ln((1 + sin φ2)/(1 - sin φ2)) - ln((1 + sin φ1)/(1 - sin φ1))]
    = h/4π * ln((1 + sin φ2)(1 - sin φ1)/((1 - sin φ2)(1 + sin φ1)))

    with h being the height of ξ in the triangle spanned by the three given points and
    φ1 and φ2 being the angles between h and the two sides extending from ξ, respectively.

    Note that the result is premultiplied by 4π!

    References:
    [1] S. Rjasanow. Vorkonditionierte iterative Auflösung von Randelementgleichungen
        für die Dirichlet-Aufgabe. Wissenschaftliche Schriftreihe der Technischen Uni-
        versität Karl-Marx-Stadt, 7/1990.

    @param ξ
        Observation point
    @param x1
        One node of the surface element
    @param x2
        x1's next neighbor in counterclockwise direction
    @param normal
        Unit normal vector of the surface element
    @return T
=#
function laplacepot{T}(::Type{InPlane}, ξ::Vector{T}, x1::Vector{T}, x2::Vector{T}, normal::Vector{T})
    # Construct triangle sides
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
    if max(0, h) < eps() ||            # ξ on the line through v1 and v2 (|φ1| = |φ2| = π/2)
        1 - abs(sinφ1) < eps() ||      # ξ on the line through v1 and v2 (φ1 = π/2)
        1 - abs(sinφ2) < eps() ||      # ξ on the line through v1 and v2 (φ2 = π/2)
        abs(sinφ1 - sinφ2) < eps()     # v1 and v2 on the line through ξ (φ1 = φ2 = 0)
        return zero(T)
    end

    # Since the observation point lies in the same plane as the surface element, we can decide
    # whether the result of this function is to be added or subtracted from the whole surface
    # triangle's Laplace potential by checking on which side of v the observation point lies.
    # This is equivalent to checking whether the normal of the triangle here and the one of
    # the surface element (which are obviously (anti)parallel) are oriented alike.
    #
    # TODO check why we have to multiply 0.5 here
    sign(u1, u2, normal) * .5 * h * log((1+sinφ2) * (1-sinφ1) / ((1-sinφ2) * (1+sinφ1)))
end

#=
    TODO
=#
function laplacepot{T}(::Type{InSpace}, ξ::Vector{T}, x1::Vector{T}, x2::Vector{T}, normal::Vector{T})
    zero(T)
end

#=
    TODO
=#
laplacepot{T, L <: ObservationPosition}(ξloc::Type{L}, ξ::Vector{T}, elem::Element{T}) = begin
    laplacepot(ξloc, ξ, elem.v1, elem.v2, elem.normal) +
    laplacepot(ξloc, ξ, elem.v2, elem.v3, elem.normal) +
    laplacepot(ξloc, ξ, elem.v3, elem.v1, elem.normal)
end

#=
    TODO
=#
laplacepot_dn(::Type{InPlane}) = 0

#=
    TODO
=#
function laplacepot_dn(::Type{InPlane})
    zero(T)
end

#=
    Generates a potential matrix according to the given function f. Use the function aliases
    with "SingleLayer" or "DoubleLayer" as defined below.

    Note that the result is premultiplied by 4π!

    References:
    [1] S. Rjasanow. Vorkonditionierte iterative Auflösung von Randelementgleichungen für
        die Dirichlet-Aufgabe. Wissenschaftliche Schriftreihe der Technischen Universität
        Karl-Marx-Stadt, 7/1990.

    @param dest
                Destination matrix
    @param elements
                List of all surface elements
    @param f
                Supported functions: laplacepot, laplacepot_dn
    @param zerodiag
                Specifies whether the diagonal elements should be zero
=#
function laplacecoll!_{T}(dest::DenseArray{T,2}, elements::Vector{Element{T}}, f::Function, zerodiag::Bool=false)
    numelem = length(elements)
    @inbounds for eidx in 1:numelem, oidx in 1:numelem
        #TODO check whether zerodiag is necessary
        if zerodiag && eidx == oidx
            dest[oidx, eidx] = zero(T)
            continue
        end

        ξ = elements[oidx].center
        elem = elements[eidx]

        # Check if ξ lies in the same plane as the surface element
        ξloc = abs(dist(ξ, elem)) < eps() ? InPlane : InSpace

        dest[oidx, eidx] = f(ξloc, ξ, elem)
    end
end
laplacecoll!{T}(::Type{SingleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}) = laplacecoll!_(dest, elements, laplacepot)
laplacecoll!{T}(::Type{DoubleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}) = laplacecoll!_(dest, elements, laplacepot_dn, true)

end # module
