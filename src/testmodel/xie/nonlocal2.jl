# =========================================================================================
"""
    struct NonlocalXieModel2{T}
        radius ::T                 # radius of the origin-centered sphere
        charges::Vector{Charge{T}} # point charges in the sphere
        params ::Option{T}         # system constants
        len    ::Int               # number of terms to be computed
        M₁     ::Matrix{T}         # coefficients C₁ for each charge
        M₂     ::Matrix{T}         # coefficients C₂ for each charge
        M₃     ::Matrix{T}         # coefficients C₃ for each charge
    end

Representation of the second nonlocal Poisson dielectric model described in
[[Xie16]](@ref Bibliography). This model comprises a full
[`XieSphere`](@ref NESSie.TestModel.XieSphere) and the coefficients ``C_{in}`` with
``i = 1, 2, 3`` (cf. Eqs. (39a-c)) for each point charge in the model, which are used in
the computation of the electrostatic potentials.

# Constructor
    NonlocalXieModel2(model::XieSphere{T}, len::Int)

The model is created solely from the given [`XieSphere`](@ref NESSie.TestModel.XieSphere) and
the number of terms to be used to approximate the original infinite sum (Eqs. 38a-b).
"""
const NonlocalXieModel2{T} = XieTestModel{T, :NonlocalXieModel2}

@inline function NonlocalXieModel2(model::XieSphere{T}, len::Int) where T
    NonlocalXieModel2{T}(model, len)
end

function _xie_coefficients(::Type{NonlocalXieModel2{T}}, model::XieSphere{T}, len::Int) where T
    a = model.radius
    λ = model.params.λ
    εΩ = model.params.εΩ
    εΣ = model.params.εΣ
    ε∞ = model.params.ε∞
    κ  = λ \ √(εΣ/ε∞)

    C₁ = Array{T}(undef, len, length(model.charges))  # Eq. (39a)
    C₂ = Array{T}(undef, len, length(model.charges))  # Eq. (39b)
    C₃ = Array{T}(undef, len, length(model.charges))  # Eq. (39c)

    iₐ = spherical_besseli(len, a/λ)
    kₐ = spherical_besselk(len, a/λ)
    kₖ = spherical_besselk(len, a * κ)

    # Eq. (57)
    u = collect(T,
        -n * εΩ + (εΣ - ε∞) * n * (n + 1) * 2a / T(π) / λ * iₐ(n) * kₐ(n)
        for n in 0:len-1
    )

    # Eq. (59)
    s = collect(T,
        (εΣ - ε∞)/ε∞ * (2n + 1) * λ / a * iₐ(n) * kₖ(n) +
        εΣ / ε∞ * iₐ(n+1) * kₖ(n) +
        κ * λ * iₐ(n) * kₖ(n+1)
        for n in 0:len-1
    )

    # Eq. (60)
    w = collect(T,
        (εΣ - ε∞)/ε∞ * n * λ / a * iₐ(n) * kₖ(n) +
        εΣ / ε∞ * iₐ(n+1) * kₖ(n) +
        κ * λ * iₐ(n) * kₖ(n+1)
        for n in 0:len-1
    )

    # Eq. (58)
    e = εΣ .* collect(T, 1:len) .* w .- u .* s

    @inbounds for (qi, q) in enumerate(model.charges)
        r = _norm(q.pos)
        for n in 0:len-1
            sₙ = s[n+1]
            wₙ = w[n+1]
            eₙ = e[n+1]
            C₁[n+1, qi] = (2n + 1) / T(4π) / eₙ * r^n * wₙ
            C₂[n+1, qi] = -(2n + 1) * (n + 1) / T(4π) / a^(n+2) / eₙ * λ * r^n * iₐ(n)
            C₃[n+1, qi] = r^n / T(4π) / a^(2n + 1) * ((2n + 1) * sₙ / eₙ - 1 / εΩ)
        end
    end
    (C₁, C₂, C₃)
end
