module NonLocalBEM

import ArrayViews: view
import Base.LinAlg.BLAS: gemv!, axpy!

export
	# types.jl
	Element,
	Charge,
	Option,
	SingleLayer,
	DoubleLayer,

	# hmo.jl
	readhmo,
	readhmo_nodes,
	readhmo_elements,
	readhmo_charges,

	# this file
	eye!,
	isdegenerate,
	compute_props!,
	compute_singularpot,
	compute_regularyukawapot,
	compute_regularyukawapot_dn,
	compute_regularyukawacoll!,
	compute_cauchy

include("types.jl")
include("hmo.jl")

const defaultopt = Option(2., 78., 1.8, 20.)

#=
	Computes the element properties, that is, centroid, normal, distance to origin, and
	area.

	@param elem
				Element of interest
=#
function compute_props!{T}(elem::Element{T})
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
	Computes the molecular potential (and the normal derivative) of the given system of
	point charges in a structureless medium.

	Note that the results are premultiplied by 4π!

	@param elements
				List of elements in the system
	@param charges
				List of charges in the system
	@param opt
				Constants to be used, including the dielectric constant of the solute
	@return (Vector{T}, Vector{T})
=#
function compute_singularpot{T}(elements::Vector{Element{T}}, charges::Vector{Charge{T}}, opt::Option{T}=defaultopt)
	umol = T[]; qmol = T[]
	for elem in elements
		push!(umol, 0)
		push!(qmol, 0)
		for charge in charges
			r = elem.center - charge.pos
			rnorm = vecnorm(r)

			umol[end] += charge.val / rnorm
			qmol[end] -= charge.val * (r ⋅ elem.normal) / rnorm^3
		end
	end
	(opt.εΩ \ umol, opt.εΩ \ qmol)
end

#=
	Compute the regular part of the yukawa potential, that is, Yukawa minus Laplace:
	e^[-√(εΣ/ε∞)/λ * |x-ξ|] / |x-ξ|  -  1 / |x-ξ|

	Note that the result is premultiplied by 4π!

	@param x
				Integration variable
	@param ξ
				Observation point
	@param opt
				Constants to be used
	@return T
=#
function compute_regularyukawapot{T}(x::DenseArray{T,1}, ξ::Vector{T}, opt::Option{T}=defaultopt)
	#=== TIME- AND MEMORY-CRITICAL CODE! ===#
	rnorm = vecnorm(x-ξ)

	# limit for |x-ξ| → 0
	if rnorm <= 0
		return -opt.yukawa
	end

	scalednorm = opt.yukawa * rnorm

	# guard against cancellation
	if scalednorm < .1
		# use alternating series to approximate
		# e^(-c) - 1 = Σ((-c)^i / i!) for i=1 to ∞
		term = -scalednorm
		tolerance = 1e-16 * abs(term)
		tsum = zero(T) 	# DON'T EVER USE 0 HERE! Time: x2, Memory: x3
		for i in 1:15
			if abs(term) <= tolerance
				continue
			end

			tsum += term
			term *= -scalednorm / (i+1)
		end
		return tsum / rnorm
	end

	# no danger of cancellation
	(exp(-scalednorm) - 1) / rnorm
end
compute_regularyukawapot{T}(x::DenseArray{T,1}, ξ::Vector{T}, ::Vector{T}, opt::Option{T}=defaultopt) = compute_regularyukawapot(x, ξ, opt)

#=
	Compute the normal derivative of the regular part of the yukawa potential, that is, 
	Yukawa minus Laplace:
	d/dn [e^[-√(εΣ/ε∞)/λ * |x-ξ|] / |x-ξ|  -  1 / |x-ξ|]
	= [1 - (1 + √(εΣ/ε∞)/λ * |x-ξ|)e^(√(εΣ/ε∞)/λ * |x-ξ|)] / |x-ξ|²   * (x-ξ)⋅n / |x-ξ|
	= [1 - (1 + c)e^(-c)] / |x-ξ|²   * (x-ξ)⋅n / |x-ξ|
	with c ≔ √(εΣ/ε∞)/λ * |x-ξ|

	Note that the result is premultiplied by 4π!

	@param x
				Integration variable
	@param ξ
				Observation point
	@param normal
				Unit normal vector at x
	@param opt
				Constants to be used
	@return T
=#
function compute_regularyukawapot_dn{T}(x::Vector{T}, ξ::Vector{T}, normal::Vector{T}, opt::Option{T}=defaultopt)
	#=== TIME- AND MEMORY-CRITICAL CODE! ===#
	r = x - ξ
	rnorm = vecnorm(r)

	# limit for |x-ξ| → 0
	if rnorm <= 0
		return zero(T)
	end

	cosovernorm2 = ((r / rnorm) ⋅ normal) / rnorm / rnorm
	scalednorm = opt.yukawa * rnorm

	# guard against cancellation
	if scalednorm < .1
		# use alternating series to approximate
		# 1 - (c+1)e^(-c) = Σ((-c)^i * (i-1) / i!) for i=2 to ∞
		term = scalednorm * scalednorm / 2
		tolerance = 1e-16 * abs(term)
		tsum = zero(T)  # DON'T EVER USE 0 HERE!
		for i in 2:16
			if abs(term * (i-1)) <= tolerance
				continue
			end

			tsum += term * (i-1)
			term *= -scalednorm / i
		end
		return tsum * cosovernorm2
	end

	# no danger of cancellation
	(1 - (1 + scalednorm) * exp(-scalednorm)) * cosovernorm2
end

#=
	Compute the Dirichlet trace of the single layer potential or the essential part of the
	double layer potential of Yukawa minus Laplace, depending on the function given. Use
	the function aliases with "SingleLayer" or "DoubleLayer" as defined	below.

	Please note that, in the latter case, the relation of K (Kʸ) to the full double layer
	potential W (Wʸ) is given by
	[(γ₀ξ^{int} W)f](ξ) = [-1 + σ(ξ)]f(ξ) + [Kf](ξ)

	This function uses a Radon cubature with seven points to generate the regular part of
	the Yukawa potential matrix. References:

	[1] V. I. Krilov. Priblizhennoe vichislenie integralov. Moskva, Nauka, 1967.
	[2] J. Radon. Zur mechanischen Kubatur. Monatsh. für Math. 52(4): 286-300, 1948.

	@param dest
				Destination matrix
	@param elements
				List of all surface elements
	@param f
				compute_regularyukawapot or compute_regularyukawapot_dn for single or
				double layer computation, respectively
	@param opt
				Constants to be used
=#
function compute_regularyukawacoll!_{T}(dest::DenseArray{T,2}, elements::Vector{Element{T}}, f::Function, opt::Option{T}=defaultopt)
	#=== MEMORY-CRITICAL CODE! ===#
	numelem = length(elements)
	@assert size(dest) == (numelem, numelem)
	
	const r15 = √15
	const ξ = (1/3, (6+r15)/21, (9-2r15)/21, (6+r15)/21, (6-r15)/21, (9+2r15)/21, (6-r15)/21)
	const η = (1/3, (9-2r15)/21, (6+r15)/21, (6+r15)/21, (9+2r15)/21, (6-r15)/21, (6-r15)/21)
	const μ = (9/80, (155+r15)/2400, (155+r15)/2400, (155+r15)/2400, (155-r15)/2400, (155-r15)/2400, (155-r15)/2400)
	
	# pre-allocate memory for cubature points
	cubpts = [zeros(T, 3) for _ in 1:7]
	# pre-compute μ * centroid for all elements and μ
	center = [[elem.center * μ[i] for i in 1:7] for elem in elements]

	@inbounds for eidx in 1:numelem
		elem = elements[eidx]
		u = elem.v2 - elem.v1
		v = elem.v3 - elem.v1
		area = 2. * elem.area
		
		# compute cubature points
		# devectorized version of cubpts = [u * ξ[i] + v * η[i] + elem.v1 for i in 1:7]
		for i in 1:7, j in 1:3
			cubpts[i][j] = ξ[i] * u[j] + η[i] * v[j] + elem.v1[j]
		end

		for oidx in 1:numelem
			obs = elements[oidx]
			value = zero(T)
			for i in 1:7
				#value += f(cubpts[i], obs.center * μ[i], elem.normal, opt)
				value += f(cubpts[i], center[eidx][i], elem.normal, opt) 
			end
			dest[oidx, eidx] = value * area
		end
	end
	nothing
end
compute_regularyukawacoll!{T}(::Type{SingleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}; opt::Option{T}=defaultopt) = compute_regularyukawacoll!_(dest, elements, compute_regularyukawapot, opt)
compute_regularyukawacoll!{T}(::Type{DoubleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}; opt::Option{T}=defaultopt) = compute_regularyukawacoll!_(dest, elements, compute_regularyukawapot_dn, opt)

#=
	TODO
=#
function compute_cauchy{T}(elements::Vector{Element{T}}, charges::Vector{Charge{T}}, opt::Option{T}=defaultopt)
	# convient access to constants
	const εΩ = opt.εΩ
	const εΣ = opt.εΣ
	const ε∞ = opt.ε∞

	# create system matrix
	numelem = length(elements)
	m = zeros(T, 3 * numelem, 3 * numelem)

	# convenient access to 9 blocks of the system matrix
	m11 = view(m,          1:numelem,           1:numelem )
	m21 = view(m,          1:numelem,   1+numelem:2numelem)
	m31 = view(m,          1:numelem,  1+2numelem:3numelem)
	m12 = view(m,  1+numelem:2numelem,          1:numelem )
	m22 = view(m,  1+numelem:2numelem,  1+numelem:2numelem)
	m32 = view(m,  1+numelem:2numelem, 1+2numelem:3numelem)
	m13 = view(m, 1+2numelem:3numelem,          1:numelem )
	m23 = view(m, 1+2numelem:3numelem,  1+numelem:2numelem)
	m33 = view(m, 1+2numelem:3numelem, 1+2numelem:3numelem)

	# initialize the system matrix
	eye!(m11, 2π)
	eye!(m21, 2π)
	eye!(m33, 2π)

	# compute molecular potential for the point charges
	umol, qmol = compute_singularpot(elements, charges, opt)

	# create right hand side
	rhs = zeros(T, 3 * numelem)
	
	# convenient access to the first block of rhs
	β = view(rhs, 1:numelem)
	
	# initialize rhs
	copy!(β, umol)
	scale!(β, -2π)


	#=
		generate and apply Kʸ-K
	=#
	buffer = Array(T, numelem, numelem)
	compute_regularyukawacoll!(DoubleLayer, buffer, elements)

	# β += (1-εΩ/εΣ)(Kʸ-K)umol
	gemv!(1-εΩ/εΣ, buffer, umol, β)

	# m11 -= Kʸ-K
	axpy!(-1., buffer, m11)

	# m13 += ε∞/εΣ * (Kʸ-K)
	axpy!(ε∞/εΣ, buffer, m13)
	

	#=
		generate and apply Vʸ-V
	=#
	compute_regularyukawacoll!(SingleLayer, buffer, elements)

	# β += (εΩ/εΣ - εΩ/ε∞)(Vʸ-V)qmol
	gemv!(εΩ * (1/εΣ - 1/ε∞), buffer, qmol, β)

	# m12 += (εΩ/ε∞ - εΩ/εΣ)(Vʸ-V)
	axpy!(εΩ * (1/ε∞ - 1/εΣ), buffer, m12)

	
	#=
		generate and apply K
	=#

	# TODO
	#…

	#println([m11 m12 m13; m21 m22 m23; m31 m32 m33]) 
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
function eye!{T}(m::DenseArray{T,2}, α::Number=one(T))
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
gemv!{T}(α::T, m::DenseArray{T,2}, v::Vector{T}, dest::DenseArray{T,1}) = gemv!(α, m, v, one(T), dest)
gemv!{T}(α::T, m::DenseArray{T,2}, v::Vector{T}, β::T, dest::DenseArray{T,1}) = gemv!('N', α, m, v, β, dest)

end # module

