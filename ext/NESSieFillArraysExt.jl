module NESSieFillArraysExt

using NESSie.BEM:
    _ImplicitTriangleMatrix,
    _ImplicitTriangleQuadMatrix,
    LocalSystemMatrix,
    NonlocalSystemMatrix
using FillArrays: AbstractZeros

function Base.:*(
    A::_ImplicitTriangleMatrix{T},
    x::AbstractZeros{T, 1}
) where T
    zeros(T, size(A, 1))
end

function Base.:*(
    A::_ImplicitTriangleQuadMatrix{T},
    x::AbstractZeros{T, 1}
) where T
    zeros(T, size(A, 1))
end

function Base.:*(
    A::LocalSystemMatrix{T},
    x::AbstractZeros{T, 1}
) where T
    zeros(T, size(A, 1))
end

function Base.:*(
    A::NonlocalSystemMatrix{T},
    x::AbstractZeros{T, 1}
) where T
    zeros(T, size(A, 1))
end

end # module
