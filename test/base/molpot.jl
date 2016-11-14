using ProteinES: φmol, ∂ₙφmol

context("φmol and ∂ₙφmol") do
    for T in testtypes
        # empty lists
        @fact φmol(Vector{T}[], Charge{T}[]) --> T[]
        # empty charge list
        elements = [Triangle(T[0, 0, 0], T[0, 0, 3], T[0, 3, 0]),
                    Triangle(T[3, 0, 0], T[0, 4, 0], T[0, 0, 5])]
        map(props!, elements)
        umol = φmol(elements, Charge{T}[])
        @fact typeof(umol) --> Vector{T}
        @fact umol --> zeros(T, 2)
        qmol = ∂ₙφmol(elements, Charge{T}[])
        @fact typeof(qmol) --> Vector{T}
        @fact qmol --> zeros(T, 2)
        # single charge
        charges = [Charge(T, 0, 0, 0, √2)]
        umol = φmol(elements, charges)
        @fact typeof(umol) --> Vector{T}
        @fact umol --> roughly(T[1, .6])
        qmol = ∂ₙφmol(elements, charges)
        @fact typeof(qmol) --> Vector{T}
        @fact qmol --> roughly(T[0, -162/25/√769])
        # multiple charges
        push!(charges, Charge(T, 1, 1, 1, -√5))
        umol = φmol(elements, charges)
        @fact typeof(umol) --> Vector{T}
        @fact umol --> roughly(T[1 - √5, -2.4])
        qmol = ∂ₙφmol(elements, charges)
        @fact typeof(qmol) --> Vector{T}
        @fact qmol --> roughly(T[√5, 1593/25/√769])
    end
end

@pending ∇φmol --> :nothing
