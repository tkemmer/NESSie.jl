# =========================================================================================
"""
    struct Kfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T} end

Interaction function for an implicit representation of potential matrix `K`.
"""
struct Kfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T} end

@inline function (f::Kfun{T})(ξ::Vector{T}, elem::Triangle{T}; kwargs...) where T
    Rjasanow.laplacecoll(DoubleLayer, ξ, elem; kwargs...)
end


# =========================================================================================
"""
    struct Kyfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T}
        yuk::T   # exponent of the Yukawa operator's fundamental solution
    end

Interaction function for an implicit representation of potential matrix `Kʸ`.
"""
struct Kyfun{T} <: InteractionFunction{Vector{T}, TriangleQuad{T}, T}
    """Exponent of the Yukawa operator's fundamental solution"""
    yuk::T
end

@inline function (f::Kyfun{T})(ξ::Vector{T}, elem::TriangleQuad{T}) where T
    Radon.regularyukawacoll(DoubleLayer, ξ, elem, f.yuk)
end


# =========================================================================================
"""
    struct Vfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T} end

Interaction function for an implicit representation of potential matrix `V`.
"""
struct Vfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T} end

@inline function (f::Vfun{T})(ξ::Vector{T}, elem::Triangle{T}; kwargs...) where T
    Rjasanow.laplacecoll(SingleLayer, ξ, elem; kwargs...)
end


# =========================================================================================
"""
    struct Vyfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T}
        yuk::T   # exponent of the Yukawa operator's fundamental solution
    end

Interaction function for an implicit representation of potential matrix `Vʸ`.
"""
struct Vyfun{T} <: InteractionFunction{Vector{T}, TriangleQuad{T}, T}
    """Exponent of the Yukawa operator's fundamental solution"""
    yuk::T
end

@inline function (f::Vyfun{T})(ξ::Vector{T}, elem::TriangleQuad{T}) where T
    Radon.regularyukawacoll(SingleLayer, ξ, elem, f.yuk)
end


# =========================================================================================
"""
    _get_laplace_matrices(
        Ξ       ::Vector{Vector{T}},
        elements::Vector{Triangle{T}}
    )

Constructs and returns implicit representations for the single- and double-layer Laplace
potential matrices `V` and `K` for the given observation points and surface elements.

# Return type
[`Tuple{InteractionMatrix{T}, InteractionMatrix{T}}`]
(https://tkemmer.github.io/ImplicitArrays.jl/stable/#ImplicitArrays.InteractionMatrix)
"""
@inline function _get_laplace_matrices(
    Ξ       ::Vector{Vector{T}},
    elements::Vector{Triangle{T}}
) where T
    V = InteractionMatrix(Ξ, elements, Vfun{T}())
    K = InteractionMatrix(Ξ, elements, Kfun{T}())
    (V, K)
end


# =========================================================================================
"""
    _get_yukawa_matrices(
        Ξ       ::Vector{Vector{T}},
        elements::Vector{Triangle{T}},
        yuk     ::T
    )

Constructs and returns implicit representations for the single- and double-layer Yukawa
potential matrices `V^Y` and `K^Y` for the given observation points and surface elements.

# Return type
[`Tuple{InteractionMatrix{T}, InteractionMatrix{T}}`]
(https://tkemmer.github.io/ImplicitArrays.jl/stable/#ImplicitArrays.InteractionMatrix)
"""
@inline function _get_yukawa_matrices(
    Ξ::Vector{Vector{T}},
    elements::Vector{Triangle{T}},
    yuk::T
) where T
    pelm = TriangleQuad.(elements)
    Vy   = InteractionMatrix(Ξ, pelm, Vyfun{T}(yuk))
    Ky   = InteractionMatrix(Ξ, pelm, Kyfun{T}(yuk))
    (Vy, Ky)
end


# =========================================================================================
"""
    _solve_linear_system(
        A::AbstractArray{T, 2},
        b::AbstractArray{T, 1}
    )

Solves the linear system `Ax = b` from `A` and `b` using a Jacobi-preconditioned GMRES solver.

# Return type
`Array{T, 1}`

# Supported keyword arguments
Everything from [`IterativeSolvers.gmres`]
(https://iterativesolvers.julialinearalgebra.org/stable/linear_systems/gmres/) except for `log`.
"""
@inline function _solve_linear_system(A::AbstractArray{T, 2}, b::AbstractArray{T, 1}; kwargs...) where T
    gmres(A, b;
        verbose=true,
        restart=min(200, size(A, 2)),
        Pl=DiagonalPreconditioner(A),
        kwargs...,
        log = false # true would change the return type of gmres
    )
end

const _ImplicitBEMMatrix{T} = Union{
    InteractionMatrix{T, Vector{T}, Triangle{T}},
    InteractionMatrix{T, Vector{T}, TriangleQuad{T}}
}

function _mul!(
    dst::AbstractVector{T},
    A::_ImplicitBEMMatrix{T},
    x::AbstractVector{T};
    kwargs...
) where T
    for i in 1:size(A, 1)
        s = zero(T)
        for j in 1:size(A, 2)
            s += getindex(A, i, j; kwargs...) * x[j]
        end
        dst[i] = s
    end
    dst
end

function _mul!(
    dst::AbstractMatrix{T},
    A::_ImplicitBEMMatrix{T},
    B::AbstractMatrix{T};
    kwargs...
) where T
    m, n, p = (size(A)..., size(B, 2))
    for i in 1:m
        for k in 1:p
            s = zero(T)
            for j in 1:n
                s += getindex(A, i, j; kwargs...) * B[j, k]
            end
            dst[i, k] = s
        end
    end
    dst
end

function Base.:*(
    A::InteractionMatrix{T, Vector{T}, Triangle{T}},
    x::AbstractVector{T}
) where T
    dst = zeros(T, size(A, 1))
    tasks = map(index_chunks(eachrow(A); n = Threads.nthreads())) do idx
        Threads.@spawn _mul!(view(dst, idx), InteractionMatrix(A, idx, :), x; dat = Vector{T}(undef, 12))
    end
    wait.(tasks)
    dst
end

function Base.:*(
    A::InteractionMatrix{T, Vector{T}, TriangleQuad{T}},
    x::AbstractVector{T}
) where T
    dst = zeros(T, size(A, 1))
    tasks = map(index_chunks(eachrow(A); n = Threads.nthreads())) do idx
        Threads.@spawn _mul!(view(dst, idx), InteractionMatrix(A, idx, :), x)
    end
    wait.(tasks)
    dst
end

function Base.:*(
    A::InteractionMatrix{T, Vector{T}, Triangle{T}},
    B::SubArray{T, 2}
) where T
    dst = zeros(T, size(A, 1), size(B, 2))
    tasks = map(index_chunks(eachrow(A); n = Threads.nthreads())) do idx
        Threads.@spawn _mul!(view(dst, idx, :), InteractionMatrix(A, idx, :), B; dat = Vector{T}(undef, 12))
    end
    wait.(tasks)
    dst
end

function Base.:*(
    A::InteractionMatrix{T, Vector{T}, TriangleQuad{T}},
    B::SubArray{T, 2}
) where T
    dst = zeros(T, size(A, 1), size(B, 2))
    tasks = map(index_chunks(eachrow(A); n = Threads.nthreads())) do idx
        Threads.@spawn _mul!(view(dst, idx, :), InteractionMatrix(A, idx, :), B)
    end
    wait.(tasks)
    dst
end

@inline function LinearAlgebra.mul!(
    dst::AbstractArray{T, 1},
    A  ::_ImplicitBEMMatrix{T},
    B  ::AbstractArray{T}
) where T
    dst .= A * B
end
