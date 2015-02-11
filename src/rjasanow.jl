import Base: log, asin, vecnorm

#=
	Rjasanow's implementation to generate a potential matrix according to the given
	function f. Use the function aliases with "SingleLayer" or "DoubleLayer" as defined
	below.

	Note that the result is premultiplied by 4π!

	@param dest
				Destination matrix
	@param elements
				List of all surface elements
	@param f
				Supported functions: compute_rjasanowsinglepot, compute_rjasanowdoublepot
	@param zerodiag
				Specifies whether the diagonal elements should be zero
=#
function compute_rjasanowcoll!_{T}(dest::DenseArray{T,2}, elements::Vector{Element{T}}, f::Function, zerodiag::Bool=false)
	numelem = length(elements)
	gridsize = zero(T)
	for elem in elements, node in (elem.v1, elem.v2, elem.v3)
		gridsize = max(gridsize, vecnorm(node)^2)
	end
	gridsize = √(gridsize)

	@inbounds for eidx in 1:numelem, oidx in 1:numelem
		dest[oidx, eidx] = zerodiag && eidx == oidx ? zero(T) : f(elements[eidx], elements[oidx].center, gridsize)
	end
end
compute_rjasanowcoll!{T}(::Type{SingleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}) = compute_rjasanowcoll!_(dest, elements, compute_rjasanowsinglepot)
compute_rjasanowcoll!{T}(::Type{DoubleLayer}, dest::DenseArray{T,2}, elements::Vector{Element{T}}) = compute_rjasanowcoll!_(dest, elements, compute_rjasanowdoublepot, true)

#=
	Rjasanow's implementation to compute the Laplace potential.

	Note that the result is premultiplied by 4π!

	@param elem
				Integration element
	@param obs
				Observation point
	@param gridsize
				Size of the grid
	@return T
=#
function compute_rjasanowsinglepot{T}(elem::Element{T}, obs::Vector{T}, gridsize::T)
	t = elem.normal .* elem.normal
	kom = t[2] > t[1] ? (t[3] > t[2] ? 3 : 2) : (t[3] > t[1] ? 3 : 1)
	ivord = elem.normal[kom] > 0 ? 1 : -1
	d = obs ⋅ elem.normal + elem.distorig
	
	if abs(d) < 1e-15 * gridsize
		return rjasanow_intebn(obs, elem.v1, elem.v2, kom, ivord) +
		       rjasanow_intebn(obs, elem.v2, elem.v3, kom, ivord) +
			   rjasanow_intebn(obs, elem.v3, elem.v1, kom, ivord)
	end
	
	ym = obs - d*elem.normal
	d = abs(d)
	rjasanow_intrau(ym, elem.v1, elem.v2, kom, ivord, gridsize, d) +
	rjasanow_intrau(ym, elem.v2, elem.v3, kom, ivord, gridsize, d) +
	rjasanow_intrau(ym, elem.v3, elem.v1, kom, ivord, gridsize, d)
end

#=
	Rjasanow's implementation to compute the normal derivative of the Laplace potential.

	Note that the result is premultiplied by 4π!

	@param elem
				Integration element
	@param obs
				Observation point
	@param gridsize
				Size of the grid
	@return T
=#
function compute_rjasanowdoublepot{T}(elem::Element{T}, obs::Vector{T}, gridsize::T)
	t = elem.normal .* elem.normal
	kom = t[2] > t[1] ? (t[3] > t[2] ? 3 : 2) : (t[3] > t[1] ? 3 : 1)
	ivord = elem.normal[kom] > 0 ? 1 : -1
	d = obs ⋅ elem.normal + elem.distorig
	
	if abs(d) < 1e-15 * gridsize
		return zero(T)
	end

	ym = obs - d*elem.normal
	rjasanow_intrau_neu(ym, elem.v1, elem.v2, kom, ivord, gridsize, d) +
	rjasanow_intrau_neu(ym, elem.v2, elem.v3, kom, ivord, gridsize, d) +
	rjasanow_intrau_neu(ym, elem.v3, elem.v1, kom, ivord, gridsize, d)
end

#=
	???
=#
function rjasanow_sihdrk{T}(x::Vector{T}, y::Vector{T}, z::Vector{T}, kom::Int)
	u1 = y - x
	d1 = vecnorm(u1)^2
	u2 = z - x
	d2 = vecnorm(u2)^2
	u  = z - y
	d3 = vecnorm(u)^2
	@assert d3 != 0
	
	if d1 == 0 || d2 == 0
		return (zero(T), zero(T), zero(T), 0)::(T,T,T,Int)
	end

	sif1 = max(-one(T), min(u1 ⋅ u / √(d1*d3), one(T)))
	sif2 = max(-one(T), min(u2 ⋅ u / √(d2*d3), one(T)))
	h = √(d2*(one(T)-sif2^2))
	s = kom == 1 ? u1[2]*u2[3] - u1[3]*u2[2] : 
	    kom == 2 ? u1[3]*u2[1] - u1[1]*u2[3] :
		           u1[1]*u2[2] - u1[2]*u2[1]
	ivora = s > 0 ?  1 :
	        s < 0 ? -1 :
			         0
	(sif1, sif2, h, ivora)::(T,T,T,Int)
end

#=
	???
=#
function rjasanow_intebn{T}(x::Vector{T}, y::Vector{T}, z::Vector{T}, kom::Int, ivord::Int)
	sif1, sif2, h, ivora = rjasanow_sihdrk(x, y, z, kom)
	if h == 0 || sif1 == sif2 || abs(sif1) == 1 || abs(sif2) == 1
		return zero(T)
	end
	(ivord == ivora ? .5 : -.5) * h * log((1+sif2) * (1-sif1) / ((1-sif2) * (1+sif1)))
end

#=
	???
=#
function rjasanow_intrau{T}(x::Vector{T}, y::Vector{T}, z::Vector{T}, kom::Int, ivord::Int, gridsize::T, d::T)
	sif1, sif2, h, ivora = rjasanow_sihdrk(x, y, z, kom)
	if abs(h) < 1e-15 * gridsize || abs(sif1-sif2) < 1e-15 || abs(abs(sif1)-1) < 1e-15 || abs(abs(sif2)-1) < 1e-15
		return zero(T)
	end
	q = 1/√(d^2 + h^2)
	qk = h * q
	q = d * q
	qq = q^2
	w1 = √(1 - qq*sif1^2)
	w2 = √(1 - qq*sif2^2)
	rr = .5 * h * log((w2 + qk*sif2) * (w1 - qk*sif1) / ((w2 - qk*sif2) * (w1 + qk*sif1))) + 
	     d * (asin(q * sif2) - asin(q * sif1) - asin(sif2) + asin(sif1))
	ivora == ivord ? rr : -rr
end

#=
	???
=#
function rjasanow_intrau_neu{T}(x::Vector{T}, y::Vector{T}, z::Vector{T}, kom::Int, ivord::Int, gridsize::T, d::T)
	sif1, sif2, h, ivora = rjasanow_sihdrk(x, y, z, kom)
	if abs(h) < 1e-15 * gridsize || abs(sif1-sif2) < 1e-15 || abs(abs(sif1)-1) < 1e-15 || abs(abs(sif2)-1) < 1e-15
		return zero(T)
	end
	akappa = abs(d) / √(d^2 + h^2)
	sint = asin(sif2) - asin(sif1) - asin(akappa * sif2) + asin(akappa * sif1)
	if d < 0
		sint = -sint
	end
	ivord == ivora ? sint : -sint
end

