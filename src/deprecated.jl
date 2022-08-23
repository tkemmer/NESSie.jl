# =========================================================================================
# deprecated since v1.3

@noinline function Radon.regularyukawacoll(
              ::Type{SingleLayer},
        ξ     ::AbstractVector{T},
        elem  ::Triangle{T},
        yukawa::T
    ) where T
    Base.depwarn(
        "This variant of the function is deprecated, painfully slow, and will be removed " *
        "soon! Please use the overloaded version with `elem::TriangleQuad{T}` parameter " *
        "instead!",
        :regularyukawacoll_single_triangle
    )
    Radon.radoncoll(ξ, quadraturepoints([elem])[1], yukawa, Radon.regularyukawapot)
end

@noinline function Radon.regularyukawacoll(
              ::Type{DoubleLayer},
        ξ     ::AbstractVector{T},
        elem  ::Triangle{T},
        yukawa::T
    ) where T
    Base.depwarn(
        "This variant of the function is deprecated, painfully slow, and will be removed " *
        "soon! Please use the overloaded version with `elem::TriangleQuad{T}` parameter " *
        "instead!",
        :regularyukawacoll_double_triangle
    )
    Radon.radoncoll(ξ, quadraturepoints([elem])[1], yukawa, Radon.∂ₙregularyukawapot)
end
