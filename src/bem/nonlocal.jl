# =========================================================================================
"""
    struct NonlocalBEMResult{T, E} <: BEMResult{T, E}
        model::Model{T, E}
        u    ::SubArray{T,1}   # [γ₀int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
        q    ::SubArray{T,1}   # [γ₁int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
        w    ::SubArray{T,1}   # [γ₀ext(Ψ)](ξ)     ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
        umol ::Vector{T}       # [γ₀int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
        qmol ::Vector{T}       # [γ₁int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0
    end

Result data of the nonlocal solving process to be used for potential computation and
post-processing, with `Ξ` being the list of observation points, that is, the set of
triangle centroids.
"""
struct NonlocalBEMResult{T, E} <: BEMResult{T, E}
    model::Model{T, E}
    u::SubArray{T,1}
    q::SubArray{T,1}
    w::SubArray{T,1}
    umol::Vector{T}
    qmol::Vector{T}
end


# =========================================================================================
# Documented in bem/local.jl
function solve(
         ::Type{NonlocalES},
    model::Model{T, Triangle{T}}
) where T
    # convenient access
    elements = model.elements
    εΩ       = model.params.εΩ
    εΣ       = model.params.εΣ
    ε∞       = model.params.ε∞
    yuk      = yukawa(model.params)

    # create system matrix
    numelem = length(elements)
    m = zeros(T, 3 * numelem, 3 * numelem)

    # convenient access to 9 blocks of the system matrix
    m11 = view(m,          1:numelem,           1:numelem )
    m12 = view(m,          1:numelem,   1+numelem:2numelem)
    m13 = view(m,          1:numelem,  1+2numelem:3numelem)
    m21 = view(m,  1+numelem:2numelem,          1:numelem )
    m22 = view(m,  1+numelem:2numelem,  1+numelem:2numelem)
#   m23 = view(m,  1+numelem:2numelem, 1+2numelem:3numelem)
#   m31 = view(m, 1+2numelem:3numelem,          1:numelem )
    m32 = view(m, 1+2numelem:3numelem,  1+numelem:2numelem)
    m33 = view(m, 1+2numelem:3numelem, 1+2numelem:3numelem)

    # initialize the system matrix;
    # since all other components of the system matrix will be premultiplied by 4π,
    # do the same for σ here
    pluseye!(m11, T(4π * σ))
    pluseye!(m21, T(4π * σ))
    pluseye!(m33, T(4π * σ))

    # compute molecular potential for the point charges;
    # molecular potentials are initially premultiplied by 4π⋅ε0⋅εΩ
    umol = εΩ .\   φmol(model)
    qmol = εΩ .\ ∂ₙφmol(model)

    # create right hand side
    rhs = zeros(T, 3 * numelem)

    # convenient access to the first block of rhs
    β = view(rhs, 1:numelem)

    # initialize rhs;
    # again, we apply a prefactor of 4π to σ to match the other components of the vector
    copyto!(β, umol)
    rmul!(β, -T(4π * σ))

    # create list of observation points
    Ξ = [e.center for e in elements]

    #=
        generate and apply Kʸ-K
    =#
    buffer = Array{T}(undef, numelem, numelem)
    Radon.regularyukawacoll!(DoubleLayer, buffer, elements, Ξ, yuk)

    # β += (1-εΩ/εΣ)(Kʸ-K)umol
    gemv!(1-εΩ/εΣ, buffer, umol, β)

    # m11 -= Kʸ-K
    axpy!(-1, buffer, m11)

    # m13 += ε∞/εΣ * (Kʸ-K)
    axpy!(ε∞/εΣ, buffer, m13)

    #=
        generate and apply Vʸ-V
    =#
    Radon.regularyukawacoll!(SingleLayer, buffer, elements, Ξ, yuk)

    # β += (εΩ/εΣ - εΩ/ε∞)(Vʸ-V)qmol
    gemv!(εΩ * (1/εΣ - 1/ε∞), buffer, qmol, β)

    # m12 += (εΩ/ε∞ - εΩ/εΣ)(Vʸ-V)
    axpy!(εΩ * (1/ε∞ - 1/εΣ), buffer, m12)

    #=
        generate and apply K
    =#
    Rjasanow.laplacecoll!(DoubleLayer, buffer, elements, Ξ)

    # β += K⋅umol
    gemv!(one(T), buffer, umol, β)

    # m11 -= K
    axpy!(-1, buffer, m11)

    # m21 += K
    axpy!(1, buffer, m21)

    # m33 -= K
    axpy!(-1, buffer, m33)

    #=
        generate and apply V
    =#
    Rjasanow.laplacecoll!(SingleLayer, buffer, elements, Ξ)

    # β -= εΩ/ε∞ * V * qmol
    gemv!(-εΩ/ε∞, buffer, qmol, β)

    # m12 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m12)

    # m22 -= V
    axpy!(-1, buffer, m22)

    # m32 += εΩ/ε∞ * V
    axpy!(εΩ/ε∞, buffer, m32)

    # solve system
    cauchy = m \ rhs

    NonlocalBEMResult(
        model,
        view(cauchy, 1:          numelem),
        view(cauchy, 1+numelem: 2numelem),
        view(cauchy, 1+2numelem:3numelem),
        umol,
        qmol
    )
end


# =========================================================================================
"""
    struct NonlocalSystemMatrix{T} <: AbstractArray{T, 2}
        V     ::InteractionMatrix{T}   # single-layer Laplace
        K     ::InteractionMatrix{T}   # double-layer Laplace
        Vy    ::InteractionMatrix{T}   # single-layer Yukawa
        Ky    ::InteractionMatrix{T}   # double-layer Yukawa
        params::Option{T}              # system constants
    end

Implicit representation of the nonlocal BEM system matrix.

# Special constructors
```julia

```
"""
struct NonlocalSystemMatrix{T} <: AbstractArray{T, 2}
    """Single-layer Laplace potentials"""
    V::InteractionMatrix{T, Vector{T}, Triangle{T}, Vfun{T}}
    """Double-layer Laplace potentials"""
    K::InteractionMatrix{T, Vector{T}, Triangle{T}, Kfun{T}}
    """Single-layer Yukawa potentials"""
    Vy::InteractionMatrix{T, Vector{T}, TriangleQuad{T}, Vyfun{T}}
    """Double-layer Yukawa potentials"""
    Ky::InteractionMatrix{T, Vector{T}, TriangleQuad{T}, Kyfun{T}}
    """System constants"""
    params::Option{T}

    function NonlocalSystemMatrix{T}(
        Ξ       ::Vector{Vector{T}},
        elements::Vector{Triangle{T}},
        params  ::Option{T}
    ) where T
        V, K = _get_laplace_matrices(Ξ, elements)
        Vy, Ky = _get_yukawa_matrices(Ξ, elements, yukawa(params))
        new(V, K, Vy, Ky, params)
    end
end

Base.size(A::NonlocalSystemMatrix{T}) where T = 3 .* size(A.K)

function LinearAlgebra.diag(A::NonlocalSystemMatrix{T}, k::Int = 0) where T
    k != 0 && error("diag not defined for k != 0 on ", typeof(A))
    σ  = T(2π) .* ones(T, size(A.K, 1))
    [σ .- diag(A.Ky); diag(A.V); σ]
end

function Base.:*(
    A::NonlocalSystemMatrix{T},
    x::AbstractArray{T, 1}
) where T
    εΩ  = A.params.εΩ
    εΣ  = A.params.εΣ
    ε∞  = A.params.ε∞
    numelem = size(A.K, 2)

    x1 = view(x, 1:numelem)
    x2 = view(x, numelem+1:2numelem)
    x3 = view(x, 2numelem+1:3numelem)

    Kx  = A.K * [x1 x3]
    Vx2 = A.V * x2
    σx1 = T(2π) .* x1

    [
        A.Ky * ((ε∞/εΣ) .* x3 .- x1) .- Kx[:,1] .+ (εΩ/ε∞-εΩ/εΣ) .* 
            (A.Vy * x2) .+ ((εΩ/ε∞) .* Vx2) .+ σx1;
        Kx[:,1] .- Vx2 .+ σx1;
        (εΩ/ε∞) .* Vx2 .- Kx[:,2] .+ (T(2π) * x3)
    ]
end

@inline function LinearAlgebra.mul!(
    Y::AbstractArray{T, 1},
    A::NonlocalSystemMatrix{T},
    v::AbstractArray{T, 1}
) where T
    Y .= A * v
end


# =========================================================================================
# Documented in bem/local.jl
function solve_implicit(
         ::Type{NonlocalES},
    model::Model{T, Triangle{T}}
) where T
    # observation points
    Ξ = [e.center for e in model.elements]

    # shortcuts
    εΩ  = model.params.εΩ
    εΣ  = model.params.εΣ
    ε∞  = model.params.ε∞
    numelem = length(model.elements)

    # compute molecular potential for the point charges;
    # molecular potentials are initially premultiplied by 4π⋅ε0⋅εΩ
    umol = εΩ .\   φmol(model, tolerance=_etol(T))
    qmol = εΩ .\ ∂ₙφmol(model)

    # create nonlocal system
    A = NonlocalSystemMatrix{T}(Ξ, model.elements, model.params)
    b = A.K * umol .+ (1 - εΩ/εΣ) .* (A.Ky * umol) .- T(2π) .* umol .- (εΩ/ε∞) .* 
        (A.V * qmol) .+ (εΩ/εΣ - εΩ/ε∞) .* (A.Vy * qmol)

    cauchy = _solve_linear_system(A, b)

    NonlocalBEMResult(
        model,
        view(cauchy, 1:          numelem),
        view(cauchy, 1+numelem: 2numelem),
        view(cauchy, 1+2numelem:3numelem),
        umol,
        qmol
    )
end
