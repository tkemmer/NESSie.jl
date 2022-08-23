# =========================================================================================
"""
    struct Kfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T} 
        dat::Vector{T}   # pre-allocated vector for internal use
    end

Interaction function for an implicit representation of potential matrix `K`.

# Special constructors
```julia
Kfun{T}()
```
Automatically initializes the internal data vector. Replaces the default constructor.
"""
struct Kfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T} 
    """Pre-allocated vector for internal use"""
    dat::Vector{T}

    Kfun{T}() where T = new(Vector{T}(undef, 12 * Threads.nthreads()))
end

@inline function (f::Kfun{T})(ξ::Vector{T}, elem::Triangle{T}) where T
    t = Threads.threadid()
    Rjasanow.laplacecoll(DoubleLayer, ξ, elem; dat=view(f.dat, (t-1)*12+1:t*12))
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
    struct Vfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T}
        dat::Vector{T}   # pre-allocated vector for internal use
    end

Interaction function for an implicit representation of potential matrix `V`.

# Special constructors
```julia
Vfun{T}()
```
Automatically initializes the internal data vector. Replaces the default constructor.
"""
struct Vfun{T} <: InteractionFunction{Vector{T}, Triangle{T}, T} 
    """Pre-allocated vector for internal use"""
    dat::Vector{T}

    Vfun{T}() where T = new(Vector{T}(undef, 12 * Threads.nthreads()))
end

@inline function (f::Vfun{T})(ξ::Vector{T}, elem::Triangle{T}) where T
    t = Threads.threadid()
    Rjasanow.laplacecoll(SingleLayer, ξ, elem; dat=view(f.dat, (t-1)*12+1:t*12))
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
    _get_laplace_matrices{T}(
        Ξ       ::Vector{Vector{T}},
        elements::Vector{Triangle{T}}
    )

Constructs and returns implicit representations for the single- and double-layer Laplace
potential matrices `V` and `K` for the given observation points and surface elements.

# Return type
`Tuple{InteractionMatrix{T}, InteractionMatrix{T}}`
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
    _get_yukawa_matrices{T}(
        Ξ       ::Vector{Vector{T}},
        elements::Vector{Triangle{T}},
        yuk     ::T
    )

Constructs and returns implicit representations for the single- and double-layer Yukawa
potential matrices `V^Y` and `K^Y` for the given observation points and surface elements.

# Return type
`Tuple{InteractionMatrix{T}, InteractionMatrix{T}}`
"""
@inline function _get_yukawa_matrices(
    Ξ::Vector{Vector{T}},
    elements::Vector{Triangle{T}},
    yuk::T
) where T
    pelm = quadraturepoints(elements)
    Vy   = InteractionMatrix(Ξ, pelm, Vyfun{T}(yuk))
    Ky   = InteractionMatrix(Ξ, pelm, Kyfun{T}(yuk))
    (Vy, Ky)
end


# =========================================================================================
"""
    _solve_linear_system{T}(
        A::AbstractArray{T, 2},
        b::AbstractArray{T, 1}
    )

Solves the linear system `Ax = b` from `A` and `b` using a Jacobi-preconditioned GMRES solver.

# Return type
`Array{T, 1}`
"""
@inline function _solve_linear_system(A::AbstractArray{T, 2}, b::AbstractArray{T, 1}) where T
    gmres(A, b,
        verbose=true,
        restart=min(200, size(A, 2)),
        Pl=DiagonalPreconditioner(A)
    )
end

function Base.:*(
    A  ::InteractionMatrix{T},
    x  ::AbstractVector{T}
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
    A  ::InteractionMatrix{T},
    B  ::AbstractMatrix{T}
) where T
    m, n, p = (size(A)..., size(B, 2))
    dst = zeros(T, m, p)
    Threads.@threads for i in 1:m
        for k in 1:p
            s = zero(T)
            for j in 1:n
                s += A[i, j] * B[j, k]
            end
            dst[i, k] = s
        end
    end
    dst
end

@inline function LinearAlgebra.mul!(
    dst::AbstractArray{T, 1},
    A  ::InteractionMatrix{T, Vector{T}, Triangle{T}},
    B  ::AbstractArray{T}
) where T
    dst .= A * B
end
