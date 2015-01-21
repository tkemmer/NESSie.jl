include("nonlocal.jl")

using NonLocalBEM

fh = open(ARGS[1])
const nodes, elements, charges = readhmo(fh, dtype=Float64)
close(fh)

@printf("%5d nodes\n", length(nodes))
@printf("%5d elements\n", length(elements))
@printf("%5d charges\n", length(charges))

map(compute_props!, elements)
@time compute_cauchy(elements, charges)


