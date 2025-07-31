# =========================================================================================
"""
    struct LocalXieModel{T}
        radius ::T                 # radius of the origin-centered sphere
        charges::Vector{Charge{T}} # point charges in the sphere
        params ::Option{T}         # system constants
        len    ::Int               # number of terms to be computed
        M₁     ::Matrix{T}         # coefficients A₁ = C₁ for each charge
        M₂     ::Matrix{T}         # coefficients A₂ = C₂ for each charge
        M₃     ::Matrix{T}         # coefficients A₃ = C₃ for each charge
    end

Representation of the local Poisson dielectric model described in
[[Xie16]](@ref Bibliography). This model comprises a full
[`XieSphere`](@ref NESSie.TestModel.XieSphere) and the coefficients used for
the computation of the electrostatic potentials.

# Constructor
    LocalXieModel(model::XieSphere{T}, len::Int)

The model is created solely from the given [`XieSphere`](@ref NESSie.TestModel.XieSphere) and
the number of terms to be used to approximate the original infinite sum (Eqs. 38a-b).
"""
const LocalXieModel{T} = XieTestModel{T, :LocalXieModel}

@inline function LocalXieModel(model::XieSphere{T}, len::Int) where T
    LocalXieModel{T}(model, len)
end

function _xie_coefficients(::Type{LocalXieModel{T}}, model::XieSphere{T}, len::Int) where T
    # coefficients are based on the second nonlocal test model w/ ε∞ = εΣ
    a = model.radius
    λ = model.params.λ
    εΩ = model.params.εΩ
    εΣ = model.params.εΣ

    C₁ = Array{T}(undef, len, length(model.charges))  # Eq. (39a)
    C₂ = Array{T}(undef, len, length(model.charges))  # Eq. (39b)
    C₃ = Array{T}(undef, len, length(model.charges))  # Eq. (39c)

    iₐ = spherical_besseli(len, a/λ)
    kₐ = spherical_besselk(len, a/λ)

    # Eq. (57)
    u = collect(T, -n * εΩ for n in 0:len-1)

    # Eqs. (59) and (60)
    s = collect(T, iₐ(n+1) * kₐ(n) + iₐ(n) * kₐ(n+1) for n in 0:len-1)

    # Eq. (58)
    e = (εΣ .* collect(T, 1:len) .- u) .* s

    @inbounds for (qi, q) in enumerate(model.charges)
        r = _norm(q.pos)
        for n in 0:len-1
            sₙ = s[n+1]
            wₙ = s[n+1]
            eₙ = e[n+1]
            C₁[n+1, qi] = (2n + 1) / T(4π) / eₙ * r^n * wₙ
            C₂[n+1, qi] = -(2n + 1) * (n + 1) / T(4π) / a^(n+2) / eₙ * λ * r^n * iₐ(n)
            C₃[n+1, qi] = r^n / T(4π) / a^(2n + 1) * ((2n + 1) * sₙ / eₙ - 1 / εΩ)
        end
    end
    (C₁, C₂, C₃)
end
