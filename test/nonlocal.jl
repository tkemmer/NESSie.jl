# check singularpot
for T in (Float64, Float32)
    # empty lists
    @test singularpot(Triangle{T}[], Charge{T}[]) == (T[], T[])
    # empty charge list
    elements = [Triangle(map(T, [0, 0, 0]), map(T, [0, 0, 3]), map(T, [0, 3, 0])),
                Triangle(map(T, [3, 0, 0]), map(T, [0, 4, 0]), map(T, [0, 0, 5]))]
    map(props!, elements)
    (umol, qmol) = singularpot(elements, Charge{T}[])
    @test isa(umol, Vector{T}) && isa(qmol, Vector{T})
    @test umol == qmol == [0, 0]
    # single charge (trivial options)
    opt = Option(one(T), zeros(T, 3)...)
    charges = [Charge(T, 0, 0, 0, √2)]
    (umol, qmol) = singularpot(elements, charges, opt)
    @test isa(umol, Vector{T}) && isa(qmol, Vector{T})
    @test_approx_eq umol [1, .6]
    @test_approx_eq qmol [0, -162/25/√769]
    # single charge (non-trivial options)
    opt = Option(convert(T, 2), zeros(T, 3)...)
    (umol, qmol) = singularpot(elements, charges, opt)
    @test isa(umol, Vector{T}) && isa(qmol, Vector{T})
    @test_approx_eq umol [.5, .3]
    @test_approx_eq qmol [0, -162/50/√769]
    # multiple charges (non-trivial options)
    push!(charges, Charge(T, 1, 1, 1, -√5))
    (umol, qmol) = singularpot(elements, charges, opt)
    @test isa(umol, Vector{T}) && isa(qmol, Vector{T})
    @test_approx_eq umol [2 \ (1-√5), -1.2]
    @test_approx_eq qmol [2 \ √5, 1593/50/√769]
end

# TODO cauchy
