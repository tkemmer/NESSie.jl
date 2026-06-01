module NESSieSciMLBaseExt

using NESSie.BEM:
    _ImplicitTriangleMatrix,
    _ImplicitTriangleQuadMatrix,
    LocalSystemMatrix,
    NonlocalSystemMatrix
using SciMLBase: AbstractNoTimeSolution

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

end # module
