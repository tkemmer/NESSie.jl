# =========================================================================================
"""
    struct NonlocalXieModel1{T}
        radius ::T                 # radius of the origin-centered sphere
        charges::Vector{Charge{T}} # point charges in the sphere
        params ::Option{T}         # system constants
        len    ::Int               # number of terms to be computed
        A₁     ::Array{T, 2}       # coefficients A₁ for each charge
        A₂     ::Array{T, 2}       # coefficients A₂ for each charge
        A₃     ::Array{T, 2}       # coefficients A₃ for each charge
    end

Representation of the first nonlocal Poisson dielectric model described in
[[Xie16]](@ref Bibliography). This model comprises a full
[`XieModel`](@ref NESSie.TestModel.XieModel) and the coefficients ``A_{in}`` with
``i = 1, 2, 3`` (cf. Eqs. (20a-c)) for each point charge in the model, which are used in
the computation of the electrostatic potentials.

!!! note
    This type does not provide a trivial constructor.

# Constructor
    NonlocalXieModel1(
        model::XieModel{T},
        len  ::Int
    )

The model is created solely from the given [`XieModel`](@ref NESSie.TestModel.XieModel) and
the number of terms to be used to approximate the original infinite sum (Eq. 18). The
coefficient vectors are computed automatically via
[`coefficients`](@ref NESSie.TestModel.coefficients).
"""
struct NonlocalXieModel1{T}
    """radius of the origin-centered sphere"""
    radius::T
    """point charges in the sphere"""
    charges::Vector{Charge{T}}
    """system constants"""
    params::Option{T}
    """number of terms to be computed"""
    len::Int
    """coefficients A₁ₙ for each charge"""
    A₁::Array{T, 2}
    """coefficients A₂ₙ for each charge"""
    A₂::Array{T, 2}
    """coefficients A₃ₙ for each charge"""
    A₃::Array{T, 2}

    NonlocalXieModel1{T}(
        model::XieModel{T},
        len  ::Int
    ) where T = begin
        (A₁, A₂, A₃) = coefficients(model, len)
        new(model.radius, model.charges, model.params, len, A₁, A₂, A₃)
    end
end

@inline NonlocalXieModel1(
    model::XieModel{T},
    len  ::Int
) where T = NonlocalXieModel1{T}(model, len)


# =========================================================================================
"""
    function coefficients(
        model::XieModel{T},
        len  ::Int
    )

Computes the coefficients ``A_{in}`` with ``i=1, 2, 3`` for the given
[`XieModel`](@ref NESSie.TestModel.XieModel) and the desired number of terms.

# Return type
`Tuple{
    Array{T, 2},
    Array{T, 2},
    Array{T, 2}
}`
"""
function coefficients(model::XieModel{T}, len::Int) where T
    a  = model.radius
    λ  = model.params.λ
    εΩ = model.params.εΩ
    εΣ = model.params.εΣ
    ε∞ = model.params.ε∞
    κ  = λ \ √(εΣ/ε∞)
    c  = λ * (εΣ - ε∞) / (a^2 * ε∞)  # use as -c for ..(ε∞ - εΣ)..   (cf. A₁)

    A₁ = Array{T}(undef, len, length(model.charges))  # Eq. (20a)
    A₂ = Array{T}(undef, len, length(model.charges))  # Eq. (20b)
    A₃ = Array{T}(undef, len, length(model.charges))  # Eq. (20c)

    iₐ = spherical_besseli(len, a/λ)  # Eq. (A1)
    kₐ = spherical_besselk(len, a*κ)  # Eq. (A2)

    # Eq. (21)
    d  = [
            (
                (n * (2n + 1) * εΩ * c) / a^n * iₐ(n) * kₐ(n) +
                (n * εΩ + (n + 1) * εΣ) / a^(n + 1) *
                (εΣ/ε∞ * iₐ(n+1) * kₐ(n) + κ * λ * iₐ(n) * kₐ(n+1))
            )
            for n in 0:len-1
         ]

    # Eq. (22), using Eqs. (A5)-(A7)
    w  = [
            (
                εΣ/ε∞ * kₐ(n) * (-λ * (n + 1) / a * iₐ(n) + iₐ(n-1)) +
                κ * λ * iₐ(n) * ((n + 1) / (κ * a) * kₐ(n) + kₐ(n-1))
            )
            for n in 0:len-1
         ]

    @inbounds for (qi, q) in enumerate(model.charges)
        r  = _norm(q.pos)
        iᵣ = spherical_besseli(len, max(1e-10, r/λ))
        for n in 0:len-1
            dₙ = d[n+1]
            wₙ = w[n+1]
            A₁[n+1, qi] = (2n + 1) / (4π * dₙ) *
                          (r^n / a^(n+1) * wₙ - c * n * iᵣ(n) * kₐ(n))
            A₂[n+1, qi] = λ * (2n + 1) / (4π * εΩ * dₙ) / a^(n+3) *
                          ((εΣ - εΩ) * (n + 1) * r^n / a^n * iₐ(n) -
                          (n * (εΩ + εΣ) + εΣ) * iᵣ(n))
            A₃[n+1, qi] = A₁[n+1, qi] / a^(2n + 1) +
                          A₂[n+1, qi] * (ε∞ - εΣ) / (ε∞ * a^n) * kₐ(n) -
                          r^n / (4π * εΩ * a^(2n+1))
        end
    end
    (A₁, A₂, A₃)
end
