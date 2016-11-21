function φΩ{T}(nodes::Vector{Vector{T}}, elements::Vector{Triangle{T}}, charges::Vector{Charge{T}}, LaplaceMod::Module=Rjasanow, opt::Option{T}=defaultopt(T))
    c = cauchy(elements, charges, LaplaceMod, opt)
    numelem = length(elements)
    φ = zeros(T, length(nodes))

    # convenient access
    γ0intφstar = c[1:numelem]
    γ1intφstar = c[1+numelem:2numelem]

    LaplaceMod.laplacecoll!(DoubleLayer, φ, elements, nodes, γ0intφstar)
    scale!(φ, -1)
    LaplaceMod.laplacecoll!(SingleLayer, φ, elements, nodes, γ1intφstar)
    scale!(φ, 2)
    T(1.69e-9 / 4π / ε0) * (4π \ φ + opt.εΩ \ φmol(nodes, charges))
end
