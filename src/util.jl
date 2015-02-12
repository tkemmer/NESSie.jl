#=
	Computes the element properties, that is, centroid, normal, distance to origin, and
	area.

	@param elem
				Element of interest
=#
function props!{T}(elem::Element{T})
	# reject degenerate triangles
	@assert !isdegenerate(elem) "Degenerate triangle $(elem)"

	# compute centroid
	elem.center = 3 \ (elem.v1 + elem.v2 + elem.v3)

	# compute normal
	elem.normal = (elem.v2 - elem.v1) × (elem.v3 - elem.v1)
	vnorm = vecnorm(elem.normal)
	elem.normal /= vnorm

	# compute distance to origin
	elem.distorig = - elem.normal ⋅	elem.v1	

	# compute area
	elem.area = 2 \ vnorm
	nothing
end

#=
	Initializes the given matrix m with αI, with I being an identity
	matrix with the same dimensions as m.

	@param m
				Corresponding matrix
	@param α
				Coefficient of the identity matrix
=#
function eye!{T}(m::Union(DenseArray{T,2}, SubArray{T,2}), α::Number=one(T))
	fill!(m, zero(T))
	α = convert(T, α)
	@inbounds for i in 1:min(size(m)...)
		m[i, i] = α
	end
	nothing
end

#=
	Tests whether the triangle with the given nodes is degenerate.

	@param v1
				First node of the triangle
	@param v2
				Second node of the triangle
	@param v3
				Third node of the triangle
	@return bool
=#
function isdegenerate{T <: FloatingPoint}(v1::Vector{T}, v2::Vector{T}, v3::Vector{T})
	@assert length(v1) == length(v2) == length(v3) == 3
	u1 = v2 - v1
	u2 = v3 - v1
	cosine = u1 ⋅ u2 / vecnorm(u1) / vecnorm(u2)
	v1 == v2 || v1 == v3 || v2 == v3 || 1 - abs(cosine) <= eps(T)
end
isdegenerate{T}(elem::Element{T}) = isdegenerate(elem.v1, elem.v2, elem.v3)

# Convenience aliases
gemv!{T}(α::T, m::Union(DenseArray{T,2}, SubArray{T,2}), v::Vector{T}, dest::Union(DenseArray{T,1}, SubArray{T,1})) = gemv!(α, m, v, one(T), dest)
gemv!{T}(α::T, m::Union(DenseArray{T,2}, SubArray{T,2}), v::Vector{T}, β::T, dest::Union(DenseArray{T,1}, SubArray{T,1})) = gemv!('N', α, m, v, β, dest)

