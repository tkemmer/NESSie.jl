module NESSieSciMLBaseExt

using NESSie.BEM:
    _ImplicitTriangleMatrix,
    _ImplicitTriangleQuadMatrix,
    LocalSystemMatrix,
    NonlocalSystemMatrix
using SciMLBase: AbstractNoTimeSolution

# Avoid method ambiguity (since SciMLBase v3.15.0)
# https://github.com/SciML/SciMLBase.jl/pull/1365
@static if hasmethod(Base.:*, Tuple{AbstractMatrix, AbstractNoTimeSolution})

    function Base.:*(
        A::_ImplicitTriangleMatrix{T},
        x::AbstractNoTimeSolution{T, 1}
    ) where T
        A * x.u
    end

    function Base.:*(
        A::_ImplicitTriangleQuadMatrix{T},
        x::AbstractNoTimeSolution{T, 1}
    ) where T
        A * x.u
    end

    function Base.:*(
        A::LocalSystemMatrix{T},
        x::AbstractNoTimeSolution{T, 1}
    ) where T
        A * x.u
    end

    function Base.:*(
        A::NonlocalSystemMatrix{T},
        x::AbstractNoTimeSolution{T, 1}
    ) where T
        A * x.u
    end
end

end # module
