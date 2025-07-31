bem_suite = SUITE["BEM"]["Post-processing"]

for T in (Float32, Float64)
    model = Format.readhmo(nessie_data_path("benchmark/2lzx-5k.hmo"), T)
    bem_local = BEM.solve(LocalES, model; method = :blas)
    bem_nonlocal = BEM.solve(NonlocalES, model; method = :blas)

    bem_suite["rfenergy"]["local ($T)"] = @benchmarkable rfenergy($bem_local)
    bem_suite["rfenergy"]["nonlocal ($T)"] = @benchmarkable rfenergy($bem_nonlocal)

    gridsize = 20
    x = range(-15one(T), 30one(T); length=gridsize)
    y = range(-15one(T), 30one(T); length=gridsize)
    z = range(-15one(T), 30one(T); length=gridsize)
    Ξ = [[xi, yi, zi] for xi in x for yi in y for zi in z]

    bem_suite["rfpotential"]["local ($T)"] = @benchmarkable rfpotential($Ξ, $bem_local)
    bem_suite["rfpotential"]["nonlocal ($T)"] = @benchmarkable rfpotential($Ξ, $bem_nonlocal)
    bem_suite["espotential"]["local ($T)"] = @benchmarkable espotential($Ξ, $bem_local)
    bem_suite["espotential"]["nonlocal ($T)"] = @benchmarkable espotential($Ξ, $bem_nonlocal)
end
