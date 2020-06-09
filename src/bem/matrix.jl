struct Kfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T} end
const  KMat{T} = InteractionMatrix{T, Vector{T}, Triangle{T}, Kfun{T}}

function (::Kfun{T})(ξ::Vector{T}, elem::Triangle{T}) where T
    Rjasanow.laplacecoll(DoubleLayer, ξ, elem)
end


struct Kyfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T}
    yuk::T
end
const  KyMat{T} = InteractionMatrix{T, Vector{T}, Triangle{T}, Kyfun{T}}

function (A::Kyfun{T})(ξ::Vector{T}, elem::Triangle{T}) where T
    Radon.regularyukawacoll(DoubleLayer, ξ, elem, A.yuk)
end


struct Vfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T} end
const  VMat{T} = InteractionMatrix{T, Vector{T}, Triangle{T}, Vfun{T}}

function (::Vfun{T})(ξ::Vector{T}, elem::Triangle{T}) where T
    Rjasanow.laplacecoll(SingleLayer, ξ, elem)
end


struct Vyfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T}
    yuk::T
end
const  VyMat{T} = InteractionMatrix{T, Vector{T}, Triangle{T}, Vyfun{T}}

function (A::Vyfun{T})(ξ::Vector{T}, elem::Triangle{T}) where T
    Radon.regularyukawacoll(SingleLayer, ξ, elem, A.yuk)
end


function _solve_linear_system(A::AbstractArray{T, 2}, b::AbstractArray{T, 1}) where T
    gmres(A, b,
        verbose=true,
        restart=min(200, size(A, 2)),
        Pl=DiagonalPreconditioner(A)
    )
end

function Base.:*(
    A  ::InteractionMatrix{T, Vector{T}, Triangle{T}},
    x  ::AbstractArray{T, 1}
) where T
    dst = zeros(T, size(A, 1))
    Threads.@threads for i in 1:size(A, 1)
        s = zero(T)
        for j in 1:size(A, 2)
            s += A[i, j] * x[j]
        end
        dst[i] = s
    end
    dst
end

function Base.:*(
    A  ::InteractionMatrix{T, Vector{T}, Triangle{T}},
    B  ::AbstractArray{T, 2}
) where T
    m, n, p = (size(A)..., size(B, 2))
    dst = zeros(T, m, p)
    Threads.@threads for i in 1:m
        s = zeros(T, p)
        for j in 1:n
            s += A[i, j] .* B[j, :]
        end
        dst[i, :] .= s
    end
    dst
end

function LinearAlgebra.mul!(
    dst::AbstractArray{T, 1},
    A  ::InteractionMatrix{T, Vector{T}, Triangle{T}},
    B  ::AbstractArray{T}
) where T
    dst .= A * B
end
