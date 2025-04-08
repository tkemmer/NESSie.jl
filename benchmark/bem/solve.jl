bem_suite = SUITE["BEM"]["Solvers"]

for T in (Float32, Float64), method in (:gmres, :blas)
    model = Model(TestModel.bornion("Ca", T))

    bem_suite["$method"]["local ($T)"] = @benchmarkable BEM.solve(LocalES, $model; method = $method)
    bem_suite["$method"]["nonlocal ($T)"] = @benchmarkable BEM.solve(NonlocalES, $model; method = $method)
end
