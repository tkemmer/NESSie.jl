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

# Special contructors
```julia
XieModel{T}(
    radius ::T,
    charges::Vector{Charge{T}},
    params ::Option{T}          = defaultopt(T);
    # kwargs
    compat ::Bool               = false
)
```
`compat` enables the compatibility mode and scales the model exactly like the reference
implementation ([[Xie16]](@ref Bibliography)). Use this flag if you intend to compare
the results to the reference.
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
        params ::Option{T} = defaultopt(T);
        compat ::Bool = false
    ) where T = new(radius, scalemodel(charges, radius, compat=compat), params)
end

@inline XieModel(
    radius ::T,
    charges::Vector{Charge{T}},
    params ::Option{T} = defaultopt(T);
    compat ::Bool = false
) where T = XieModel{T}(radius, charges, params, compat=compat)


# =========================================================================================
"""
    function scalemodel{T}(
        charges::Vector{Charge{T}},
        radius ::T;
        # kwargs
        compat ::Bool              = false
    )

Translates and rescales the given point charge model to fit inside an origin-centered
sphere with the specified radius. More specifically, the function will center the given
model at the origin and scaled in a way such that the outermost point charge will be
located at 80% `radius` distance from the origin.

# Arguments
 * `compat` Enables compatibility mode and scales the model exactly like the reference
   implementation ([[Xie16]](@ref Bibliography)). Use this flag if you intend to compare
   the results to the reference.

# Return type
`Vector{Charge{T}}`
"""
function scalemodel(
        charges::Vector{Charge{T}},
        radius ::T;
        compat ::Bool               = false
    ) where T
    # center model
    cpos   = sum.(extrema(q.pos[i] for q in charges) for i in 1:3)/2
    newpos = [q.pos .- cpos for q in charges]

    # compute and apply scaling factor
    h = x -> compat ? ceil(x) : x  # in compat mode, scaling factor is rounded up
    sf = .8radius / h(maximum(_norm(pos) for pos in newpos))
    rmul!(newpos, sf)

    [Charge{T}(pos, q.val) for (pos, q) in zip(newpos, charges)]
end
