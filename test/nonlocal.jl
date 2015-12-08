context("singularpot") do
    for T in testtypes
        # empty lists
        @fact singularpot(Triangle{T}[], Charge{T}[]) --> (T[], T[])
        # empty charge list
        elements = [Triangle(map(T, [0, 0, 0]), map(T, [0, 0, 3]), map(T, [0, 3, 0])),
                    Triangle(map(T, [3, 0, 0]), map(T, [0, 4, 0]), map(T, [0, 0, 5]))]
        map(props!, elements)
        (umol, qmol) = singularpot(elements, Charge{T}[])
        @fact isa(umol, Vector{T}) --> true
        @fact isa(qmol, Vector{T}) --> true
        @fact umol --> [0, 0]
        @fact qmol --> [0, 0]
        # single charge (trivial options)
        opt = Option(one(T), zeros(T, 3)...)
        charges = [Charge(T, 0, 0, 0, √2)]
        (umol, qmol) = singularpot(elements, charges, opt)
        @fact isa(umol, Vector{T}) --> true
        @fact isa(qmol, Vector{T}) --> true
        @fact umol --> roughly(map(T, [1, .6]))
        @fact qmol --> roughly(map(T, [0, -162/25/√769]))
        # single charge (non-trivial options)
        opt = Option(T(2), zeros(T, 3)...)
        (umol, qmol) = singularpot(elements, charges, opt)
        @fact isa(umol, Vector{T}) --> true
        @fact isa(qmol, Vector{T}) --> true
        @fact umol --> roughly(map(T, [.5, .3]))
        @fact qmol --> roughly(map(T, [0, -162/50/√769]))
        # multiple charges (non-trivial options)
        push!(charges, Charge(T, 1, 1, 1, -√5))
        (umol, qmol) = singularpot(elements, charges, opt)
        @fact isa(umol, Vector{T}) --> true
        @fact isa(qmol, Vector{T}) --> true
        @fact umol --> roughly(map(T, [2 \ (1-√5), -1.2]))
        @fact qmol --> roughly(map(T, [2 \ √5, 1593/50/√769]))
    end
end

@pending cauchy --> :nothing
