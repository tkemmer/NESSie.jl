include("nonlocal.jl")

using NonLocalBEM
import Base: min, max

fh = open(ARGS[1])
const nodes, elements, charges = readhmo(fh, dtype=Float64)
close(fh)

@printf("%5d nodes\n", length(nodes))
@printf("%5d elements\n", length(elements))
@printf("%5d charges\n", length(charges))

map(props!, elements)
#@time cauchy(elements, charges)

function foo{T}(elem::Vector{Element{T}})
	#elem = elem[1:10]
	numelem = length(elem)
	buffer1 = Array(T, numelem, numelem)
	buffer2 = Array(T, numelem, numelem)

 	ptype = SingleLayer
	println("=== Radon ===")
	@time laplacecoll!(ptype, buffer1, elem)
	printm(buffer1)
	@printf "Min: %d\t Max: %f\tMean: %f\n\n" min(buffer1) max(buffer1) (sum(buffer1)/length(buffer1))
	
	println("=== Rjasanow ===")
	@time rjasanowcoll!(ptype, buffer2, elem)
	printm(buffer2)
	@printf "Min: %d\t Max: %f\tMean: %f\n\n" min(buffer1) max(buffer1) (sum(buffer1)/length(buffer1))

	b = buffer2-buffer1
	println("=== Δ ===")
	printm(b)

	b = b .* b
	mse = length(b) \ sum(b)
	println("\nFull result matrix")
	@printf "MSE:  %.8f\n" mse
	@printf "RMSE: %.8f\n\n" √mse

	mse = numelem \ sum(vec(b)[1:numelem+1:length(b)])
	println("Diagonal only")
	@printf "MSE:  %.8f\n" mse
	@printf "RMSE: %.8f\n\n" √mse
	
	mse = (numelem * (numelem-1)) \ sum(deleteat!(vec(b), 1:numelem+1:length(b)))
	println("Result matrix w/o diagonal")
	@printf "MSE:  %.8f\n" mse
	@printf "RMSE: %.8f\n\n" √mse

	nothing
end

function printm{T}(m::DenseArray{T,2}, mx::Int=10)
	for i in 1:mx
		for j in 1:mx
			@printf "%+3.8f\t" m[i,j]
		end
		println()
	end
end

function min{T}(m::DenseArray{T,2})
	mn = Inf
	for el in m
		mn = min(mn, el)
	end
	mn
end

function max{T}(m::DenseArray{T,2})
	mx = -Inf
	for el in m
		mx = max(mx, el)
	end
	mx
end

#for i in 1:10
#	println(elements[i].distorig)
#end
println()
foo(elements)
