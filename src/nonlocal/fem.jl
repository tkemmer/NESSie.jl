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
    aᵢ = [(a[:,2] × a[:,3])'; (a[:,3] × a[:,1])'; (a[:,1] × a[:,2])'] # /det(a)
    bᵢ = aᵢ * elem.v1
    d = det(a)
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
localstiffness{T}(𝛁f::Vector{T}) = Symmetric(6 \ [𝛁f[row]⋅𝛁f[col] for row in 1:length(𝛁f), col in 1:length(𝛁f)], :U)

#=
    TODO
=#
function globalstiffness{T}(nodes::Vector{Vector{T}}, elements::Vector{Tetrahedron{T}}, opt::Option{T}=defaultopt(T))
    revidx = reverseindex(nodes)
    cm = spzeros(T, length(nodes), length(nodes))
    c  = spzeros(T, length(nodes), length(nodes))

    for elem in elements
        v = [revidx[object_id(v)] for v in (elem.v1, elem.v2, elem.v3, elem.v4)]
        @inbounds for row in 1:4, col in 1:4
            # TODO
            cm[v[row], v[col]] += row == col ? 1/60 : 1/120
        end
    end

    scale!(cm, opt.εΣ / (opt.λ * opt.λ * opt.ε∞))

    for elem in elements
        d, ∇f = basisfunctions(elem)
        k = abs(d) \ localstiffness(∇f)
        v = [revidx[object_id(v)] for v in (elem.v1, elem.v2, elem.v3, elem.v4)]
        @inbounds for row in 1:4, col in 1:4
            cm[v[row], v[col]] += k[row,col]
            c[v[row], v[col]]  += k[row,col]
        end
    end

    (cm, c)
end

#=
    Computes and returns the discretized charge density distribution
    at the node positions.

    ρ(r) = 1/(4π*ε0) * Σ qᵢ/|rᵢ-r|

    Note that the result is premultiplied by 4π*ε0.
=#
function chargedensity{T}(nodes::Vector{Vector{T}}, charges::Vector{Charge{T}})
    [sum([charge.val / euclidean(charge.pos, node) for charge in charges]) for node in nodes]
end

#=
    TODO
=#
function espotential{T}(nodes::Vector{Vector{T}}, elements::Vector{Tetrahedron{T}}, charges::Vector{Charge{T}}, opt::Option{T}=defaultopt(T))
    cm, c = globalstiffness(nodes, elements, opt)
    p = chargedensity(nodes, charges)
    #scale!(p, opt.λ^2 / ε0 / opt.ε∞)

    # Φ = W
    Φ = cm \ p

    # Φ *= (εΣ - ε∞)
    scale!(Φ, opt.εΣ - opt.ε∞)

    # Φ += ε∞ * V
    axpy!(opt.ε∞, c \ p, Φ)

    # Φ *= 1/(λ² * εΣ)
    #scale!(Φ, 1 / (opt.λ^2 * opt.εΣ))
    scale!(Φ, ε0^2 * 4π * 1e9 / opt.λ^4 / 1.6) #TODO check
    Φ
end
