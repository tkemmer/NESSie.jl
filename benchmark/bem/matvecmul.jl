bem_suite = SUITE["BEM"]["Matrix-vector mult"]

for T in (Float32, Float64)
    model   = Format.readhmo(nessie_data_path("benchmark/2lzx-5k.hmo"), T)
    numelem = length(model.elements)
    yuk     = NESSie.yukawa(model.params)
    Ξ       = [e.center for e in model.elements]
    x       = ones(T, numelem)
    x3      = ones(T, 3numelem)
    V, K    = BEM._get_laplace_matrices(Ξ, model.elements)
    Vy, Ky  = BEM._get_yukawa_matrices(Ξ, model.elements, yuk)
    loc     = BEM.LocalSystemMatrix(K, model.params)
    nloc    = BEM.NonlocalSystemMatrix(Ξ, model.elements, model.params)

    bem_suite["V ($T)"] = @benchmarkable $V * $x
    bem_suite["K ($T)"] = @benchmarkable $K * $x
    bem_suite["Vy ($T)"] = @benchmarkable $Vy * $x
    bem_suite["Ky ($T)"] = @benchmarkable $Ky * $x
    bem_suite["local ($T))"] = @benchmarkable $loc * $x
    bem_suite["nonlocal ($T)"] = @benchmarkable $nloc * $x3
end
