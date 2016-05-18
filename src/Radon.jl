module Radon

using ..ProteinES
using ..ProteinES: ddot
using Distances: euclidean

export laplacecoll!, regularyukawacoll!

#=
    Compute the Laplace potential 1 / |x-ξ|.

    Note that the result is premultiplied by 4π!

    @param x
        Integration variable
    @param ξ
        Observation point
    @return T
=#
function laplacepot{T}(x::DenseArray{T,1}, ξ::Vector{T})
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= 1e-10 && return zero(T)

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
laplacepot{T}(x::DenseArray{T,1}, ξ::Vector{T}, ::Vector{T}, ::Option{T}) = laplacepot(x, ξ)

#=
    Compute the normal derivative of the Laplace potential:
    - 1 / |x-ξ|^2   * (x-ξ) ⋅ n / |x-ξ|

    Note that the result is premultiplied by 4π!

    @param x
        Integration variable
    @param ξ
        Observation point
    @param normal
        Normal unit vector at x
    @return T
=#
function ∂ₙlaplacepot{T}(x::DenseArray{T,1}, ξ::Vector{T}, normal::Vector{T})
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= 1e-10 && return zero(T)

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
∂ₙlaplacepot{T}(x::DenseArray{T,1}, ξ::Vector{T}, normal::Vector{T}, ::Option{T}) = ∂ₙlaplacepot(x, ξ, normal)

#=
    Compute the regular part of the yukawa potential, that is, Yukawa minus Laplace:
    e^[-√(εΣ/ε∞)/λ * |x-ξ|] / |x-ξ|  -  1 / |x-ξ|

    Note that the result is premultiplied by 4π!

    @param x
        Integration variable
    @param ξ
        Observation point
    @param opt
        Constants to be used
    @return T
=#
function regularyukawapot{T}(x::DenseArray{T,1}, ξ::Vector{T}, opt::Option{T}=defaultopt(T))
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= 1e-10 && return -opt.yukawa

    scalednorm = opt.yukawa * rnorm

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
regularyukawapot{T}(x::DenseArray{T,1}, ξ::Vector{T}, ::Vector{T}, opt::Option{T}=defaultopt(T)) = regularyukawapot(x, ξ, opt)

#=
    Compute the normal derivative of the regular part of the yukawa potential, that is,
    Yukawa minus Laplace:
    d/dn [e^[-√(εΣ/ε∞)/λ * |x-ξ|] / |x-ξ|  -  1 / |x-ξ|]
    = [1 - (1 + √(εΣ/ε∞)/λ * |x-ξ|)e^(√(εΣ/ε∞)/λ * |x-ξ|)] / |x-ξ|²   * (x-ξ)⋅n / |x-ξ|
    = [1 - (1 + c)e^(-c)] / |x-ξ|²   * (x-ξ)⋅n / |x-ξ|
    with c ≔ √(εΣ/ε∞)/λ * |x-ξ|

    Note that the result is premultiplied by 4π!

    @param x
        Integration variable
    @param ξ
        Observation point
    @param normal
        Normal unit vector at x
    @param opt
        Constants to be used
    @return T
=#
function ∂ₙregularyukawapot{T}(x::Vector{T}, ξ::Vector{T}, normal::Vector{T}, opt::Option{T}=defaultopt(T))
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    rnorm = euclidean(x, ξ)

    # limit for |x-ξ| → 0
    rnorm <= 1e-10 && return zero(T)

    cosovernorm2 = ddot(x, ξ, normal) / rnorm^3
    scalednorm = opt.yukawa * rnorm

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

#=
    Seven-point Radon cubature for a given function and a list of triangles and
    observation points. If `dest` is a vector, the function values f for each
    surface triangle have to be specified, since each element of the vector
    represents the dot product of the corresponding coefficient matrix row and
    the `fvals` vector.

    If you intend computing single/double layer potentials with this function,
    you might want to use the shorthand signatures `laplacecoll!` and
    `regularyukawacoll` instead.

    Note that the result is premultiplied by 4π.

    References:
    [1] V. I. Krilov. Priblizhennoe vichislenie integralov. Moskva, Nauka, 1967.
    [2] J. Radon. Zur mechanischen Kubatur. Monatsh. für Math. 52(4): 286-300, 1948.

    @param dest
        Destination matrix/vector
    @param elements
        Surface elements
    @param ξlist
        Observation points
    @param solution
        Fundamental solution. Supported functions: regularyukawapot, ∂ₙregularyukawapot,
        laplacepot, ∂ₙlaplacepot
    @param fvals
        Function values of the elements
    @param opt
        Constants to be used
=#
function radoncoll!{T}(dest::Union{DenseArray{T,1}, DenseArray{T,2}}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, solution::Function, fvals::DenseArray{T,1}=T[], opt::Option{T}=defaultopt(T))
    #=== MEMORY-CRITICAL CODE! ===#
    isvec  = isa(dest, DenseArray{T,1})
    isvec && @assert length(fvals) == length(elements)
    isvec && @assert length(dest) == length(ξlist)
    isvec || @assert size(dest) == (length(ξlist), length(elements))

    # pre-allocate memory for cubature points
    qpts = quadraturepoints(Triangle, T)
    cubpts = [zeros(T, 3) for _ in 1:qpts.num]

    @inbounds for (eidx, elem) in enumerate(elements)
        u = elem.v2 - elem.v1
        v = elem.v3 - elem.v1
        area = 2 * elem.area

        # compute cubature points
        # devectorized version of cubpts = [u * qpts.x[i] + v * qpts.y[i] + elem.v1 for i in 1:7]
        for i in 1:qpts.num, j in 1:3
            cubpts[i][j] = qpts.x[i] * u[j] + qpts.y[i] * v[j] + elem.v1[j]
        end

        for (oidx, obs) in enumerate(ξlist)
            value = zero(T)
            for i in 1:qpts.num
                value += solution(cubpts[i], obs, elem.normal, opt) * qpts.weight[i]
            end
            isvec ?
                dest[oidx] += value * area * fvals[eidx] :
                dest[oidx, eidx] = value * area
        end
    end
    nothing
end

laplacecoll!{T}(::Type{SingleLayer}, dest::DenseArray{T,1}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, fvals::DenseArray{T, 1}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, ξlist, laplacepot, fvals, opt)
laplacecoll!{T}(::Type{DoubleLayer}, dest::DenseArray{T,1}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, fvals::DenseArray{T, 1}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, ξlist, ∂ₙlaplacepot, fvals, opt)
laplacecoll!{T}(::Type{SingleLayer}, dest::DenseArray{T,2}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, ξlist, laplacepot, T[], opt)
laplacecoll!{T}(::Type{DoubleLayer}, dest::DenseArray{T,2}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, ξlist, ∂ₙlaplacepot, T[], opt)

regularyukawacoll!{T}(::Type{SingleLayer}, dest::DenseArray{T,1}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, fvals::DenseArray{T, 1}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, ξlist, regularyukawapot, fvals, opt)
regularyukawacoll!{T}(::Type{DoubleLayer}, dest::DenseArray{T,1}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, fvals::DenseArray{T, 1}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, ξlist, ∂ₙregularyukawapot, fvals, opt)
regularyukawacoll!{T}(::Type{SingleLayer}, dest::DenseArray{T,2}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, ξlist, regularyukawapot, T[], opt)
regularyukawacoll!{T}(::Type{DoubleLayer}, dest::DenseArray{T,2}, elements::Vector{Triangle{T}}, ξlist::Vector{Vector{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, ξlist, ∂ₙregularyukawapot, T[], opt)

end # module
