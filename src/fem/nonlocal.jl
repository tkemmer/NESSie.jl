#=
    Returns the gradients of the linear basis functions for the
    given tetrahedron.

    Note that the results are premultiplied by the determinant!

    @param elem
        A tetrahedron
    @return (determinant, [part. derivatives])
=#
function basisfunctions{T}(elem::Tetrahedron{T})
    # map onto reference tetrahedron
    # x = Aα + b ⇔ α = A^(-1)⋅(x-b)
    a  = [elem.v2 - elem.v1 elem.v3 - elem.v1 elem.v4 - elem.v1]
    aᵢ = [(a[:,2] × a[:,3])'; (a[:,3] × a[:,1])'; (a[:,1] × a[:,2])'] # /det(a) # TODO devectorize
    bᵢ = aᵢ * elem.v1
    d = det(a)
    #d = a[:,1] ⋅ (a[:,2] × a[:,3]) # TODO check! (this is different to det(a))

    if(d < 0)
        println(d)
    end

    @assert d != 0

    # compute base functions
    #f = Function[x::Vector{T} -> aᵢ[i,:] ⋅ x - bᵢ[i] for i in 1:3]
    #insert!(f, 1, x::Vector{T} -> d - f[1](x) - f[2](x) - f[3](x))

    # compute partial derivatives
    ∇f = Vector{T}[aᵢ[i,:] for i in 1:3]
    insert!(∇f, 1, -∇f[1] - ∇f[2] - ∇f[3])

    (d, ∇f)
end

#=
    TODO
=#
function reverseprojection{T}(elem::Tetrahedron{T})
    α -> [elem.v2 - elem.v1 elem.v3 - elem.v1 elem.v4 - elem.v1] * α + elem.v1 # TODO devectorize
end

#=
    Returns the (local) stiffness matrix for the gradients of the
    linear basis functions of a single tetrahedron. The matrix does
    only represent the nodes present in the given element (all other
    entries would be zero anyway).

    ∫𝛁fᵢ⋅𝛁fⱼdx (for all nodes i, j of the element)
    = 𝛁fᵢ⋅𝛁fⱼ∫dx (since our base functions are linear)
    = 𝛁fᵢ⋅𝛁fⱼ⋅|det(A)|/6 (since our elements are tetrahedra)

    Note that the stiffness matrix is premultiplied by the absolute
    value of the determinant!

    @param ∇f
        Basis function gradients of the tetrahedron nodes v1 to v4
    @return Symmetric{T, Array{T,2}}
=#
# TODO devectorize, deduplicate (v is symmetric!)
localstiffnessmatrix_k{T}(𝛁f::Vector{T}) = Symmetric(6 \ [𝛁f[row]⋅𝛁f[col] for row in 1:length(𝛁f), col in 1:length(𝛁f)], :U)

#=
    TODO
=#
function stiffnessmatrix{T}(elements::Vector{Tetrahedron{T}}, revidx::Dict{UInt,UInt}, f::Function, domain::Symbol=:all)
    m = spzeros(T, length(revidx), length(revidx))
    for elem in elements
        # only process elements of the given domain (if domain != :all)
        domain != :all && elem.domain != domain && continue

        mlocal = f(elem)
        idx    = [revidx[object_id(e)] for e in (elem.v1, elem.v2, elem.v3, elem.v4)] # TODO pre-allocate

        # TODO m is symmetric!
        @inbounds for row in 1:4, col in 1:4
            m[idx[row], idx[col]] += mlocal[row,col]
        end
    end
    m
end
stiffnessmatrix_k{T}(elements::Vector{Tetrahedron{T}}, revidx::Dict{UInt,UInt}, domain::Symbol=:all) = stiffnessmatrix(elements, revidx, elem -> localstiffnessmatrix_k(elem), domain) # TODO
stiffnessmatrix_v{T}(elements::Vector{Tetrahedron{T}}, revidx::Dict{UInt,UInt}, domain::Symbol=:all) = stiffnessmatrix(elements, revidx, elem -> [2 1 1 1; 1 2 1 1; 1 1 2 1; 1 1 1 2], domain) # TODO

function localstiffnessmatrix_k{T}(elem::Tetrahedron{T})
    d, ∇f = basisfunctions(elem)
    #println(abs(d) \ Symmetric(6 \ [∇f[row]⋅∇f[col] for row in 1:length(∇f), col in 1:length(∇f)], :U))

    len = length(∇f)
    res = Array{T}(len, len)
    @inbounds for row in 1:len
        for col in row:len
            res[row, col] = ∇f[row]⋅∇f[col]
        end
    end
    scale!(res, 1 / (abs(d) * 6)) # FIXME pure evil!

    #=
    @inbounds for row in 1:len
        for col in row:len
            if(abs(res[row, col]) < 1e-10)
                res[row, col] = zero(T)
            end
        end
    end
    =#
    if(sum(res) > 100)
        println(res)
    end


    Symmetric(res, :U)
end

#=
    TODO
=#
function quadrature_ρ{T}(elem::Tetrahedron{T}, rproj::Function, charges::Vector{Charge{T}})
    qpts = quadraturepoints(Tetrahedron, T)
    # basis functions evaluated at cubature point
    const bt = Vector{T}[[.25, 1/6, .5, 1/6, 1/6], qpts.x, qpts.y, qpts.z]

    ρlocal = zeros(T, 4)
    q      = φmol(Vector{T}[rproj([qpts.x[i], qpts.y[i], qpts.z[i]]) for i in 1:qpts.num], charges)
    for node in 1:4, cpt in 1:qpts.num
        ρlocal[node] += qpts.weight[cpt] * bt[node][cpt] * q[cpt]
    end
    ρlocal
end

#=
    TODO
=#
function rhs_ρ{T}(elements::Vector{Tetrahedron{T}}, charges::Vector{Charge{T}}, revidx::Dict{UInt,UInt}, domain::Symbol=:all)
    ρ = zeros(T, length(revidx))

    for elem in elements
        # only process elements of the given domain (if domain != :all)
        domain != :all && elem.domain != domain && continue

        rproj  = reverseprojection(elem)
        ρlocal = quadrature_ρ(elem, rproj, charges)
        idx    = [revidx[object_id(e)] for e in (elem.v1, elem.v2, elem.v3, elem.v4)]

        @inbounds for row in 1:4
            ρ[idx[row]] += ρlocal[row]
        end
    end
    ρ
end

#=
    TODO
    merge with quadrature_ρ
=#
function quadrature_γ{T}(elem::Tetrahedron{T}, rproj::Function, ∇f::Vector{Vector{T}}, charges::Vector{Charge{T}})
    qpts   = quadraturepoints(Tetrahedron, T)
    γlocal = zeros(T, 4)
    q      = ∇φmol(Vector{T}[rproj([qpts.x[i], qpts.y[i], qpts.z[i]]) for i in 1:qpts.num], charges)
    for node in 1:4, cpt in 1:qpts.num
        γlocal[node] += qpts.weight[cpt] * (∇f[node] ⋅ q[cpt])
    end
    γlocal
end

#=
    TODO merge with rhs_ρ
=#
function rhs_γ{T}(elements::Vector{Tetrahedron{T}}, charges::Vector{Charge{T}}, revidx::Dict{UInt,UInt}, domain::Symbol=:all)
    γ = zeros(T, length(revidx))

    for elem in elements
        # only process elements of the given domain (if domain != :all)
        domain != :all && elem.domain != domain && continue

        rproj  = reverseprojection(elem)
        d, ∇f  = basisfunctions(elem)
        γlocal = d \ quadrature_γ(elem, rproj, ∇f, charges) # FIXME pure evil
        idx    = [revidx[object_id(e)] for e in (elem.v1, elem.v2, elem.v3, elem.v4)]

        @inbounds for row in 1:4
            γ[idx[row]] += γlocal[row]
        end
    end
    γ
end

#=
    TODO
=#
function espotential{T}(nodes::Vector{Vector{T}}, elements::Vector{Tetrahedron{T}}, charges::Vector{Charge{T}}, opt::Option{T}=defaultopt(T))
    # TODO reuse memory
    revidx = reverseindex(nodes)
    numnodes = length(nodes)
    v  = stiffnessmatrix_v(elements, revidx) / 120 # TODO
    kΩ = stiffnessmatrix_k(elements, revidx, :Ω)
    kΣ = stiffnessmatrix_k(elements, revidx, :Σ)

    @assert size(kΣ) == (numnodes, numnodes)


    ρ = rhs_ρ(elements, charges, revidx)

    #=
        1. Compute u₀
    =#
    println("Step 1: Solving for u0...")

    # stiffness matrix M₀
    m0 = spzeros(T, numnodes, numnodes)

    # m0 = λ² * K + V
    copy!(m0, kΩ)
    axpy!(1, kΣ, m0)
    scale!(m0, opt.λ^2)
    axpy!(1, v, m0)

    u0 = m0 \ ρ
    @assert isa(u0, Vector{T})
    @assert length(u0) == numnodes

    #=
        2. Assemble system matrix M
    =#
    println("Step 2: Assembling system matrix M...")
    m = spzeros(T, 2 * numnodes, 2 * numnodes)

    # convenient access to the matrix blocks
    m11 = view(m, 1:numnodes,                    1: numnodes)
    m12 = view(m, 1:numnodes,           1+numnodes:2numnodes)
    m21 = view(m, 1+numnodes:2numnodes,          1: numnodes)
    m22 = view(m, 1+numnodes:2numnodes, 1+numnodes:2numnodes)

    # m11 = ε∞ * KΣ + εΩ * KΩ
    axpy!(opt.ε∞, kΣ, m11)
    axpy!(opt.εΩ, kΩ, m11)

    # m12 = (εΣ - ε∞)KΣ
    axpy!(opt.εΣ - opt.ε∞, kΣ, m12)

    # m21 = -V
    axpy!(-1, v, m21)

    # m22 = λ² * K + V = M₀
    copy!(m22, m0)

    #=
        3. Compute (Ψ,u₁)
    =#
    println("Step 3: Solving for Ψ...")
    rhs = zeros(T, 2 * numnodes)
    β   = view(rhs, 1:numnodes)
    γΣ  = rhs_γ(elements, charges, revidx, :Σ) # TODO

    # β = (ε∞ - εΣ) * (KΣ * u₀) + (εΩ - ε∞) * γΣ
    #gemv!(opt.ε∞ - opt.εΣ, vΣ, u0, β) # BLAS does not support sparse matrices :\
    axpy!(opt.ε∞ - opt.εΣ, kΣ * u0, β) # TODO in-place solution
    axpy!(opt.εΩ - opt.ε∞, γΣ, β)

    Ψu1 = m \ rhs

    #=
        4. Compute (Φ*,u₂)
        TODO check if this can be ignored as long as we don't have ions in the solvent
    =#
    println("Step 4: Solving for Φ...")
    Φu2 = m \ zeros(T, 2 * numnodes)

    println("Step 5: Aggregating results...")
    ec / (4π * ε0 * opt.εΩ) * (view(Ψu1, 1:numnodes) + φmol(nodes, charges))
end
