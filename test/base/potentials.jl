using NESSie: φmol, ∂ₙφmol

context("φmol and ∂ₙφmol") do
    # TODO check list variants of the functions
    for T in testtypes
        # empty lists
        model = Model{T, Triangle{T}}()
        @fact φmol(model) --> T[]
        # empty charge list
        elements = [Triangle(T[0, 0, 0], T[0, 0, 3], T[0, 3, 0]),
                    Triangle(T[3, 0, 0], T[0, 4, 0], T[0, 0, 5])]
        model = Model(Vector{T}[], elements)
        umol = φmol(model)
        @fact typeof(umol) --> Vector{T}
        @fact umol --> zeros(T, 2)
        qmol = ∂ₙφmol(model)
        @fact typeof(qmol) --> Vector{T}
        @fact qmol --> zeros(T, 2)
        # single charge
        charges = [Charge(T[0, 0, 0, √2]...)]
        model = Model(Vector{T}[], elements, charges)
        umol = φmol(model)
        @fact typeof(umol) --> Vector{T}
        @fact umol --> roughly(T[1, .6])
        qmol = ∂ₙφmol(model)
        @fact typeof(qmol) --> Vector{T}
        @fact qmol --> roughly(T[0, -162/25/√769])
        # multiple charges
        push!(model.charges, Charge(T[1, 1, 1, -√5]...))
        umol = φmol(model)
        @fact typeof(umol) --> Vector{T}
        @fact umol --> roughly(T[1 - √5, -2.4])
        qmol = ∂ₙφmol(model)
        @fact typeof(qmol) --> Vector{T}
        @fact qmol --> roughly(T[√5, 1593/25/√769])
        # singular case
        umol = φmol(zeros(T, 3), [Charge(T[0, 0, 0, √2]...)])
        @fact typeof(umol) --> T
        @fact umol --> roughly(T(√2 / 1e-10))
        umol = φmol(zeros(T, 3), [Charge(T[0, 0, 0, √2]...)], tolerance=T(1e-12))
        @fact typeof(umol) --> T
        @fact umol --> roughly(T(√2 / 1e-12))
    end
end

@pending ∇φmol --> :nothing
