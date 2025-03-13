# =========================================================================================
"""
    mutable struct XieSphere{T}
        radius ::T                                 # radius of the origin-centered sphere
        charges::Vector{Charge{T}}                 # point charges in the sphere
        params ::Option{T}         = defaultopt(T) # system constants
    end

System model of a dielectric sphere containing multiple point charges. On construction, the
given point charge model will be translated and rescaled to fit inside an origin-centered
sphere with the specified radius.

# Special contructors
```julia
XieSphere(
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
@auto_hash_equals mutable struct XieSphere{T}
    "radius of the origin-centered sphere"
    radius::T
    "point charges in the sphere"
    charges::Vector{Charge{T}}
    "system constants"
    params::Option{T}

    XieSphere{T}(
        radius ::T,
        charges::Vector{Charge{T}},
        params ::Option{T} = defaultopt(T);
        compat ::Bool = false
    ) where T = new(radius, scalemodel(charges, radius, compat=compat), params)
end

@inline XieSphere(
    radius ::T,
    charges::Vector{Charge{T}},
    params ::Option{T} = defaultopt(T);
    compat ::Bool = false
) where T = XieSphere{T}(radius, charges, params, compat=compat)

@inline function Base.show(io::IO, ::MIME"text/plain", xie::XieSphere)
    show(io, xie)
end

@inline function Base.show(io::IO, xie::XieSphere)
    print(io,
        "$(typeof(xie))",
        "(charges = ", length(xie.charges),
        ", radius = $(xie.radius))"
    )
end


# =========================================================================================
"""
    struct XieTestModel{T, M} end

Common abstraction for all Xie sphere-based test models.
"""
@auto_hash_equals struct XieTestModel{T, M}
    """radius of the origin-centered sphere"""
    radius::T
    """point charges in the sphere"""
    charges::Vector{Charge{T}}
    """system constants"""
    params::Option{T}
    """number of terms to be computed"""
    len::Int
    """coefficients A₁ₙ or C₁ₙ for each charge"""
    M₁::Matrix{T}
    """coefficients A₂ₙ or C₂ₙ for each charge"""
    M₂::Matrix{T}
    """coefficients A₃ₙ or C₃ₙ for each charge"""
    M₃::Matrix{T}

    @inline function XieTestModel{T, M}(
        model::XieSphere{T},
        len::Int
    ) where {T, M}
        new(
            model.radius,
            model.charges,
            model.params,
            len,
            _xie_coefficients(XieTestModel{T, M}, model, len)...
        )
    end
end

@inline function Base.show(io::IO, ::MIME"text/plain", xie::XieTestModel)
    show(io, xie)
end

@inline function Base.show(io::IO, xie::XieTestModel)
    print(io,
        "$(typeof(xie))",
        "(charges = ", length(xie.charges),
        ", radius = $(xie.radius)",
        ", len = $(xie.len))"
    )
end


# =========================================================================================
"""
    function _xie_coefficients(
        ::XieTestModel{T}
        model::XieSphere{T},
        len  ::Int
    )

Depending on the concrete test model, computes the coefficients ``A_{in}`` or ``C_{in}``
with ``i=1, 2, 3`` for the given [`XieSphere`](@ref NESSie.TestModel.XieSphere) and the
desired number of terms.

# Return type
`NTuple{3, Matrix{T}}`
"""
function _xie_coefficients end


# =========================================================================================
"""
    function scalemodel(
        charges::Vector{Charge{T}},
        radius ::T;
        # kwargs
        compat ::Bool              = false
    )

Translates and rescales the given point charge model to fit inside an origin-centered
sphere with the specified radius. More specifically, the model is centered at the origin and,
if it contains more than one charge, scaled in a way such that the outermost point charge will
be located at 80% `radius` distance from the origin.

# Arguments
 * `compat` Enables compatibility mode and scales the model exactly like the reference
   implementation ([[Xie16]](@ref Bibliography)). Use this flag if you intend to compare
   the results to the reference.

# Return type
[`Vector{Charge{T}}`](@ref Charge)
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
    sf = h(maximum(_norm(pos) for pos in newpos))
    sf > 0 && rmul!(newpos, .8radius / sf)

    [Charge{T}(pos, q.val) for (pos, q) in zip(newpos, charges)]
end


# =========================================================================================
"""
    Model(xie::XieSphere)

Converts the given [Xie sphere](@ref XieSphere) into a triangle-based model, using
[Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl) for the mesh generation.

# Supported keyword arguments
 - `lc_min::Real = 0.12` corresponds to Gmsh's "Mesh.CharacteristicLengthMin"
 - `lc_max::Real = 0.13` corresponds to Gmsh's "Mesh.CharacteristicLengthMax"

# Return type
[`Model{T, Triangle{T}}`](@ref)
"""
@inline function NESSie.Model(xie::XieSphere{T}; lc_min::Real = 0.12, lc_max::Real = 0.13) where T
    model = _generate_sphere(T, zeros(T, 3), xie.radius; lc_min = lc_min, lc_max = lc_max)
    model.charges = xie.charges
    model.params  = xie.params
    model
end
