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
    rnorm = vecnorm(x-ξ)

    # limit for |x-ξ| →    0
    if rnorm <= 0
        return zero(T)
    end

    # guard against small rnorm
    if rnorm < .1
        # use alternating series to approximate
        # 1/c = Σ((-1)^i * (x-1)^i) for i = 0 to ∞
        term = one(T)
        tolerance = 1e-16
        tsum = zero(T)
        for i in 1:15
            if abs(term) < tolerance
                break
            end

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
function laplacepot_dn{T}(x::DenseArray{T,1}, ξ::Vector{T}, normal::Vector{T})
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    r = x - ξ
    rnorm = vecnorm(r)

    # limit for |x-ξ| → 0
    if rnorm <= 0
        return zero(T)
    end

    # guard against small rnorm
    if rnorm < .1
        # use alternating series to approximate
        # -1/c^3 = -1/2 * Σ((-1)^i * (x-1)^i * (i+1) * (i+2)) for i = 0 to ∞
        term = convert(T, 2)
        tolerance = 1e-16
        tsum = zero(T)
        for i in 1:15
            if abs(term) < tolerance
                break
            end

            tsum += term * (i+1) * (i+2)
            term *= -(rnorm -1)
        end
        return -.5 * tsum * (r ⋅ normal)
    end
    
    -1 / rnorm^3 * (r ⋅ normal)
end
laplacepot_dn{T}(x::DenseArray{T,1}, ξ::Vector{T}, normal::Vector{T}, ::Option{T}) = laplacepot_dn(x, ξ, normal)

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
    rnorm = vecnorm(x-ξ)

    # limit for |x-ξ| → 0
    if rnorm <= 0
        return -opt.yukawa
    end

    scalednorm = opt.yukawa * rnorm

    # guard against cancellation
    if scalednorm < .1
        # use alternating series to approximate
        # e^(-c) - 1 = Σ((-c)^i / i!) for i=1 to ∞
        term = -scalednorm
        tolerance = 1e-16 * abs(term)
        tsum = zero(T)     # DON'T EVER USE 0 HERE! Time: x2, Memory: x3
        for i in 1:15
            if abs(term) <= tolerance
                break
            end

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
function regularyukawapot_dn{T}(x::Vector{T}, ξ::Vector{T}, normal::Vector{T}, opt::Option{T}=defaultopt(T))
    #=== TIME- AND MEMORY-CRITICAL CODE! ===#
    r = x - ξ
    rnorm = vecnorm(r)

    # limit for |x-ξ| → 0
    if rnorm <= 0
        return zero(T)
    end

    cosovernorm2 = (r ⋅ normal) / rnorm^3
    scalednorm = opt.yukawa * rnorm

    # guard against cancellation
    if scalednorm < .1
        # use alternating series to approximate
        # 1 - (c+1)e^(-c) = Σ((-c)^i * (i-1) / i!) for i=2 to ∞
        term = scalednorm * scalednorm / 2
        tolerance = 1e-16 * abs(term)
        tsum = zero(T)  # DON'T EVER USE 0 HERE!
        for i in 2:16
            if abs(term #=* (i-1)=#) <= tolerance
                continue
            end

            tsum += term * (i-1)
            term *= -scalednorm / (i+1)
        end
        return tsum * cosovernorm2
    end

    # no danger of cancellation
    (1 - (1 + scalednorm) * exp(-scalednorm)) * cosovernorm2
end

#=
    Radon cubature with seven points to generate a potential matrix according to the given
    function f. For easy setup, use the function aliases laplacecoll! and
    regularyukawacoll! with "SingleLayer" or "DoubleLayer" instead.

    References:
    [1] V. I. Krilov. Priblizhennoe vichislenie integralov. Moskva, Nauka, 1967.
    [2] J. Radon. Zur mechanischen Kubatur. Monatsh. für Math. 52(4): 286-300, 1948.

    @param dest
                Destination matrix
    @param elements
                List of all surface elements
    @param f
                Supported functions: regularyukawapot, regularyukawapot_dn, laplacepot, 
                laplacepot_dn
    @param opt
                Constants to be used
=#
function radoncoll!{T}(dest::DenseArray{T,2}, elements::Vector{Element{T}}, f::Function, opt::Option{T}=defaultopt(T))
    #=== MEMORY-CRITICAL CODE! ===#
    numelem = length(elements)
    @assert size(dest) == (numelem, numelem)
    
    const r15 = √15
    const ξ = (1/3, (6+r15)/21, (9-2r15)/21, (6+r15)/21, (6-r15)/21, (9+2r15)/21, (6-r15)/21)
    const η = (1/3, (9-2r15)/21, (6+r15)/21, (6+r15)/21, (9+2r15)/21, (6-r15)/21, (6-r15)/21)
    const μ = (9/80, (155+r15)/2400, (155+r15)/2400, (155+r15)/2400, (155-r15)/2400, (155-r15)/2400, (155-r15)/2400)
    
    # pre-allocate memory for cubature points
    cubpts = [zeros(T, 3) for _ in 1:7]
    
    @inbounds for eidx in 1:numelem
        elem = elements[eidx]
        u = elem.v2 - elem.v1
        v = elem.v3 - elem.v1
        area = 2. * elem.area

        # compute cubature points
        # devectorized version of cubpts = [u * ξ[i] + v * η[i] + elem.v1 for i in 1:7]
        for i in 1:7, j in 1:3
            cubpts[i][j] = ξ[i] * u[j] + η[i] * v[j] + elem.v1[j]
        end

        for oidx in 1:numelem
            obs = elements[oidx]
            value = zero(T)
            for i in 1:7
                value += f(cubpts[i], obs.center, elem.normal, opt) * μ[i]
            end
            dest[oidx, eidx] = value * area
        end
    end
    nothing
end

#=
    Compute the Dirichlet trace of the single or double layer potential of Laplace.

    Please note that, in the latter case, the relation of K (Kʸ) to the full double layer
    potential W (Wʸ) is given by
    [(γ₀ξ^{int} W)f](ξ) = [-1 + σ(ξ)]f(ξ) + [Kf](ξ)

    This function uses a Radon cubature with seven points to generate the regular part of
    the Yukawa potential matrix.

    Note that the result is premultiplied by 4π!

    @param _
                SingleLayer or DoubleLayer
    @param dest
                Destination matrix
    @param elements
                List of elements in the system
    @param opt
                Constants to be used
=#
laplacecoll!{T}(::Type{SingleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, laplacepot, opt)
laplacecoll!{T}(::Type{DoubleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, laplacepot_dn, opt)

#=
    Compute the Dirichlet trace of the single layer potential or the essential part of the
    double layer potential of Yukawa minus Laplace.

    Please note that, in the latter case, the relation of K (Kʸ) to the full double layer
    potential W (Wʸ) is given by
    [(γ₀ξ^{int} W)f](ξ) = [-1 + σ(ξ)]f(ξ) + [Kf](ξ)

    This function uses a Radon cubature with seven points to generate the regular part of
    the Yukawa potential matrix.

    Note that the result is premultiplied by 4π!

    @param _
                SingleLayer or DoubleLayer
    @param dest
                Destination matrix
    @param elements
                List of elements in the system
    @param opt
                Constants to be used
=#
regularyukawacoll!{T}(::Type{SingleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, regularyukawapot, opt)
regularyukawacoll!{T}(::Type{DoubleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, regularyukawapot_dn, opt)

