using ProteinES.Nonlocal
using ProteinES.Nonlocal: basisfunctions, localstiffness

context("basisfunctions") do
    for T in testtypes
        # reference tetrahedron
        elem = Tetrahedron(T[0, 0, 0], T[1, 0, 0], T[0, 1, 0], T[0, 0, 1])
        d, ∇f = basisfunctions(elem)
        @fact typeof(d) --> T
        @fact d --> 1
        @fact typeof(∇f) --> Vector{Vector{T}}
        @fact ∇f[1] --> d * T[-1, -1, -1]
        @fact ∇f[2] --> d * T[1, 0, 0]
        @fact ∇f[3] --> d * T[0, 1, 0]
        @fact ∇f[4] --> d * T[0, 0, 1]

        # tetrahedron to be mapped onto the reference
        elem = Tetrahedron(T[1, 1, 1], T[3, 1, 1], T[1, 3, 1], T[1, 1, 3])
        d, ∇f = basisfunctions(elem)
        @fact typeof(d) --> T
        @fact d --> 8
        @fact typeof(∇f) --> Vector{Vector{T}}
        @fact ∇f[1] --> d * T[-.5, -.5, -.5]
        @fact ∇f[2] --> d * T[.5, 0, 0]
        @fact ∇f[3] --> d * T[0, .5, 0]
        @fact ∇f[4] --> d * T[0, 0, .5]

        # tetrahedron crossing all quadrants
        elem = Tetrahedron(T[1, -1, -1], T[0, 1, -1], T[-1, -1, -1], T[0, 0, 1])
        d, ∇f = basisfunctions(elem)
        @fact typeof(d) --> T
        @fact d --> 8
        @fact typeof(∇f) --> Vector{Vector{T}}
        @fact ∇f[1] --> d * T[.5, -.25, -.125]
        @fact ∇f[2] --> d * T[0, .5, -.25]
        @fact ∇f[3] --> d * T[-.5, -.25, -.125]
        @fact ∇f[4] --> d * T[0, 0, .5]

        # TODO negative determinant
    end
end

context("localstiffness") do
    for T in testtypes
        # reference tetrahedron
        k = localstiffness(Vector{T}[T[-1, -1, -1], T[1, 0, 0], T[0, 1, 0], T[0, 0, 1]])
        @fact typeof(k) --> Symmetric{T, Array{T, 2}}
        @fact k --> Symmetric(6 \ T[3 -1 -1 -1; 0 1 0 0; 0 0 1 0; 0 0 0 1], :U)

        # tetrahedron to be mapped onto the reference
        k = localstiffness(Vector{T}[T[-4, -4, -4], T[4, 0, 0], T[0, 4, 0], T[0, 0, 4]])
        @fact typeof(k) --> Symmetric{T, Array{T, 2}}
        @fact k --> Symmetric(6 \ T[48 -16 -16 -16; 0 16 0 0; 0 0 16 0; 0 0 0 16], :U)
    end
end

@pending globalstiffness --> :nothing
@pending chargedensity --> :nothing
@pending espotential --> :nothing
