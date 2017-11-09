# =========================================================================================
"""
    mutable struct XieModel{T}
        radius ::T                                 # radius of the origin-centered sphere
        charges::Vector{Charge{T}}                 # point charges in the sphere
        params ::Option{T}         = defaultopt(T) # system constants
    end

System model of a dielectric sphere containing multiple point charges. On construction, the
given point charge model will be translated and rescaled to fit inside an origin-centered
sphere with the specified radius.
"""
mutable struct XieModel{T}
    "radius of the origin-centered sphere"
    radius::T
    "point charges in the sphere"
    charges::Vector{Charge{T}}
    "system constants"
    params::Option{T}

    XieModel{T}(
        radius ::T,
        charges::Vector{Charge{T}},
        params ::Option{T} = defaultopt(T)
    ) where T = new(radius, scalemodel(charges, radius), params)
end

XieModel(
    radius ::T,
    charges::Vector{Charge{T}},
    params ::Option{T} = defaultopt(T)
) where T = XieModel{T}(radius, charges, params)


# =========================================================================================
"""
    function scalemodel{T}(
        charges::Vector{Charge{T}},
        radius ::T
    )

Translates and rescales the given point charge model to fit inside an origin-centered
sphere with the specified radius. More specifically, the function will center the given
model at the origin and scaled in a way such that the outermost point charge will be
located at 80% `radius` distance from the origin.

# Return type
`Vector{Charge{T}}`
"""
function scalemodel(charges::Vector{Charge{T}}, radius::T) where T
    # center model
    cpos   = sum.([extrema(q.pos[i] for q in charges) for i in 1:3])/2
    newpos = [q.pos - cpos for q in charges]

    # compute and apply scaling factor
    sf = .8radius / maximum(vecnorm(pos) for pos in newpos)
    scale!(newpos, sf)

    [Charge{T}(pos, q.val) for (pos, q) in zip(newpos, charges)]
end
