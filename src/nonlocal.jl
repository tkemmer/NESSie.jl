module NonLocalBEM

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

	# rjasanow.jl
	rjasanowcoll!,
	rjasanowsinglepot,
	rjasanowdoublepot,

	# this file
	defaultopt,
	eye!,
	isdegenerate,
	props!,
	singularpot,
	laplacepot,
	laplacepot_dn,
	regularyukawapot,
	regularyukawapot_dn,
	radoncoll!,
	laplacecoll!,
	regularyukawacoll!,
	cauchy

include("types.jl")
include("hmo.jl")
include("rjasanow.jl")

# Default options
const defaultopt64 = Option(2., 78., 1.8, 20.)
const defaultopt32 = Option(2f0, 78f0, 1.8f0, 20f0)
defaultopt(::Type{Float64}) = defaultopt64
defaultopt(::Type{Float32}) = defaultopt32

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
function singularpot{T}(elements::Vector{Element{T}}, charges::Vector{Charge{T}}, opt::Option{T}=defaultopt(T))
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
	Compute the Laplace potential 1 / |x-ξ|.

	Note that the result is premultiplied by 4π!

	@param x
				Integration variable
	@param ξ
				Observation point
	@return T
=#
function laplacepot{T}(x::DenseArray{T,1}, ξ::Vector{T})
	#=== TIME- AND MEMORY-CRITICAL CODE! ===#
	rnorm = vecnorm(x-ξ)

	# limit for |x-ξ| →	0
	if rnorm <= 0
		return zero(T)
	end

	# guard against small rnorm
	if rnorm < .1
		# use alternating series to approximate
		# 1/c = Σ((-1)^i * (x-1)^i) for i = 0 to ∞
		term = one(T)
		tolerance = 1e-16
		tsum = zero(T)
		for i in 1:15
			if abs(term) < tolerance
				break
			end

			tsum += term
			term *= -(rnorm - 1)
		end
		return tsum
	end

	1 / rnorm
end
laplacepot{T}(x::DenseArray{T,1}, ξ::Vector{T}, ::Vector{T}, ::Option{T}) = laplacepot(x, ξ)

#=
	Compute the normal derivative of the Laplace potential:
	- 1 / |x-ξ|^2   * (x-ξ) ⋅ n / |x-ξ|

	Note that the result is premultiplied by 4π!

	@param x
				Integration variable
	@param ξ
				Observation point
	@param normal
				Normal unit vector at x
	@return T
=#
function laplacepot_dn{T}(x::DenseArray{T,1}, ξ::Vector{T}, normal::Vector{T})
	#=== TIME- AND MEMORY-CRITICAL CODE! ===#
	r = x - ξ
	rnorm = vecnorm(r)

	# limit for |x-ξ| → 0
	if rnorm <= 0
		return zero(T)
	end

	# guard against small rnorm
	if rnorm < .1
		# use alternating series to approximate
		# -1/c^3 = -1/2 * Σ((-1)^i * (x-1)^i * (i+1) * (i+2)) for i = 0 to ∞
		term = convert(T, 2)
		tolerance = 1e-16
		tsum = zero(T)
		for i in 1:15
			if abs(term) < tolerance
				break
			end

			tsum += term * (i+1) * (i+2)
			term *= -(rnorm -1)
		end
		return -.5 * tsum * (r ⋅ normal)
	end
	
	-1 / rnorm^3 * (r ⋅ normal)
end
laplacepot_dn{T}(x::DenseArray{T,1}, ξ::Vector{T}, normal::Vector{T}, ::Option{T}) = laplacepot_dn(x, ξ, normal)

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
function regularyukawapot{T}(x::DenseArray{T,1}, ξ::Vector{T}, opt::Option{T}=defaultopt(T))
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
				break
			end

			tsum += term
			term *= -scalednorm / (i+1)
		end
		return tsum / rnorm
	end

	# no danger of cancellation
	(exp(-scalednorm) - 1) / rnorm
end
regularyukawapot{T}(x::DenseArray{T,1}, ξ::Vector{T}, ::Vector{T}, opt::Option{T}=defaultopt(T)) = regularyukawapot(x, ξ, opt)

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
				Normal unit vector at x
	@param opt
				Constants to be used
	@return T
=#
function regularyukawapot_dn{T}(x::Vector{T}, ξ::Vector{T}, normal::Vector{T}, opt::Option{T}=defaultopt(T))
	#=== TIME- AND MEMORY-CRITICAL CODE! ===#
	r = x - ξ
	rnorm = vecnorm(r)

	# limit for |x-ξ| → 0
	if rnorm <= 0
		return zero(T)
	end

	cosovernorm2 = (r ⋅ normal) / rnorm^3
	scalednorm = opt.yukawa * rnorm

	# guard against cancellation
	if scalednorm < .1
		# use alternating series to approximate
		# 1 - (c+1)e^(-c) = Σ((-c)^i * (i-1) / i!) for i=2 to ∞
		term = scalednorm * scalednorm / 2
		tolerance = 1e-16 * abs(term)
		tsum = zero(T)  # DON'T EVER USE 0 HERE!
		for i in 2:16
			if abs(term #=* (i-1)=#) <= tolerance
				continue
			end

			tsum += term * (i-1)
			term *= -scalednorm / (i+1)
		end
		return tsum * cosovernorm2
	end

	# no danger of cancellation
	(1 - (1 + scalednorm) * exp(-scalednorm)) * cosovernorm2
end

#=
	Radon cubature with seven points to generate a potential matrix according to the given
	function f. For easy setup, use the function aliases laplacecoll! and
	regularyukawacoll! with "SingleLayer" or "DoubleLayer" instead.

	References:
	[1] V. I. Krilov. Priblizhennoe vichislenie integralov. Moskva, Nauka, 1967.
	[2] J. Radon. Zur mechanischen Kubatur. Monatsh. für Math. 52(4): 286-300, 1948.

	@param dest
				Destination matrix
	@param elements
				List of all surface elements
	@param f
				Supported functions: regularyukawapot, regularyukawapot_dn, laplacepot, 
				laplacepot_dn
	@param opt
				Constants to be used
=#
function radoncoll!{T}(dest::DenseArray{T,2}, elements::Vector{Element{T}}, f::Function, opt::Option{T}=defaultopt(T))
	#=== MEMORY-CRITICAL CODE! ===#
	numelem = length(elements)
	@assert size(dest) == (numelem, numelem)
	
	const r15 = √15
	const ξ = (1/3, (6+r15)/21, (9-2r15)/21, (6+r15)/21, (6-r15)/21, (9+2r15)/21, (6-r15)/21)
	const η = (1/3, (9-2r15)/21, (6+r15)/21, (6+r15)/21, (9+2r15)/21, (6-r15)/21, (6-r15)/21)
	const μ = (9/80, (155+r15)/2400, (155+r15)/2400, (155+r15)/2400, (155-r15)/2400, (155-r15)/2400, (155-r15)/2400)
	
	# pre-allocate memory for cubature points
	cubpts = [zeros(T, 3) for _ in 1:7]
	
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
				value += f(cubpts[i], obs.center, elem.normal, opt) * μ[i]
			end
			dest[oidx, eidx] = value * area
		end
	end
	nothing
end

#=
	Compute the Dirichlet trace of the single or double layer potential of Laplace.

	Please note that, in the latter case, the relation of K (Kʸ) to the full double layer
	potential W (Wʸ) is given by
	[(γ₀ξ^{int} W)f](ξ) = [-1 + σ(ξ)]f(ξ) + [Kf](ξ)

	This function uses a Radon cubature with seven points to generate the regular part of
	the Yukawa potential matrix.

	Note that the result is premultiplied by 4π!

	@param _
				SingleLayer or DoubleLayer
	@param dest
				Destination matrix
	@param elements
				List of elements in the system
	@param opt
				Constants to be used
=#
laplacecoll!{T}(::Type{SingleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, laplacepot, opt)
laplacecoll!{T}(::Type{DoubleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, laplacepot_dn, opt)

#=
	Compute the Dirichlet trace of the single layer potential or the essential part of the
	double layer potential of Yukawa minus Laplace.

	Please note that, in the latter case, the relation of K (Kʸ) to the full double layer
	potential W (Wʸ) is given by
	[(γ₀ξ^{int} W)f](ξ) = [-1 + σ(ξ)]f(ξ) + [Kf](ξ)

	This function uses a Radon cubature with seven points to generate the regular part of
	the Yukawa potential matrix.

	Note that the result is premultiplied by 4π!

	@param _
				SingleLayer or DoubleLayer
	@param dest
				Destination matrix
	@param elements
				List of elements in the system
	@param opt
				Constants to be used
=#
regularyukawacoll!{T}(::Type{SingleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, regularyukawapot, opt)
regularyukawacoll!{T}(::Type{DoubleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}, opt::Option{T}=defaultopt(T)) = radoncoll!(dest, elements, regularyukawapot_dn, opt)

#=
	Computes the full cauchy data on the surface of the biomolecule.

	Note that the result is premultiplied by 4π!

	@param elements
				List of surface elements
	@param charges
				List of charges in the biomolecule
	@param opt
				Constants to be used
	@return Vector{T}
=#
function cauchy{T}(elements::Vector{Element{T}}, charges::Vector{Charge{T}}, opt::Option{T}=defaultopt(T))
	# convient access to constants
	const εΩ = opt.εΩ
	const εΣ = opt.εΣ
	const ε∞ = opt.ε∞

	# create system matrix
	numelem = length(elements)
	m = zeros(T, 3 * numelem, 3 * numelem)

	# convenient access to 9 blocks of the system matrix
	m11 = sub(m,          1:numelem,           1:numelem )
	m12 = sub(m,          1:numelem,   1+numelem:2numelem)
	m13 = sub(m,          1:numelem,  1+2numelem:3numelem)
	m21 = sub(m,  1+numelem:2numelem,          1:numelem )
	m22 = sub(m,  1+numelem:2numelem,  1+numelem:2numelem)
	m23 = sub(m,  1+numelem:2numelem, 1+2numelem:3numelem)
	m31 = sub(m, 1+2numelem:3numelem,          1:numelem )
	m32 = sub(m, 1+2numelem:3numelem,  1+numelem:2numelem)
	m33 = sub(m, 1+2numelem:3numelem, 1+2numelem:3numelem)

	# initialize the system matrix
	eye!(m11, 2π)
	eye!(m21, 2π)
	eye!(m33, 2π)

	# compute molecular potential for the point charges
	umol, qmol = singularpot(elements, charges, opt)

	# create right hand side
	rhs = zeros(T, 3 * numelem)
	
	# convenient access to the first block of rhs
	β = sub(rhs, 1:numelem)
	
	# initialize rhs
	copy!(β, umol)
	scale!(β, -2π)

	#=
		generate and apply Kʸ-K
	=#
	buffer = Array(T, numelem, numelem)
	regularyukawacoll!(DoubleLayer, buffer, elements)

	# β += (1-εΩ/εΣ)(Kʸ-K)umol
	gemv!(1-εΩ/εΣ, buffer, umol, β)

	# m11 -= Kʸ-K
	axpy!(-1., buffer, m11)

	# m13 += ε∞/εΣ * (Kʸ-K)
	axpy!(ε∞/εΣ, buffer, m13)

	#=
		generate and apply Vʸ-V
	=#
	regularyukawacoll!(SingleLayer, buffer, elements)

	# β += (εΩ/εΣ - εΩ/ε∞)(Vʸ-V)qmol
	gemv!(εΩ * (1/εΣ - 1/ε∞), buffer, qmol, β)

	# m12 += (εΩ/ε∞ - εΩ/εΣ)(Vʸ-V)
	axpy!(εΩ * (1/ε∞ - 1/εΣ), buffer, m12)
	
	#=
		generate and apply K
	=#
	rjasanowcoll!(DoubleLayer, buffer, elements)

	# β += K
	gemv!(1., buffer, umol, β)

	# m11 -= K
	axpy!(-1., buffer, m11)

	# m21 += K
	axpy!(1., buffer, m21)

	# m33 -= K
	axpy!(-1., buffer, m33)

	#=
		generate and apply V
	=#
	rjasanowcoll!(SingleLayer, buffer, elements)

	# β -= εΩ/ε∞ * V * qmol
	gemv!(-εΩ/ε∞, buffer, qmol, β)

	# m12 += εΩ/ε∞ * V
	axpy!(εΩ/ε∞, buffer, m12)

	# m22 -= V
	axpy!(-1., buffer, m22)

	# m32 += εΩ/ε∞ * V
	axpy!(εΩ/ε∞, buffer, m32)

	# solve system
	m\rhs
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

end # module

