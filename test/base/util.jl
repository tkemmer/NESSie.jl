import ProteinES: eye!, isdegenerate, seek, reverseindex, unpack, vertexnormals, cos, cathetus, sign, distance

context("eye!") do
    for T in (Int, testtypes...)
        m = zeros(T, 3, 3)
        eye!(m)
        @fact typeof(m) --> Array{T, 2}
        @fact size(m) --> (3, 3)
        @fact m --> [1 0 0; 0 1 0; 0 0 1]
        eye!(m, 2)
        @fact typeof(m) --> Array{T, 2}
        @fact size(m) --> (3, 3)
        @fact m --> [2 0 0; 0 2 0; 0 0 2]
        m = zeros(T, 2, 3)
        eye!(m)
        @fact typeof(m) --> Array{T, 2}
        @fact size(m) --> (2, 3)
        @fact m --> [1 0 0; 0 1 0]
        eye!(m, 2)
        @fact typeof(m) --> Array{T, 2}
        @fact size(m) --> (2, 3)
        @fact m --> [2 0 0; 0 2 0]
        m = zeros(T, 3, 2)
        eye!(m)
        @fact typeof(m) --> Array{T, 2}
        @fact size(m) --> (3, 2)
        @fact m --> [1 0; 0 1; 0 0]
        eye!(m, 2)
        @fact typeof(m) --> Array{T, 2}
        @fact size(m) --> (3, 2)
        @fact m --> [2 0; 0 2; 0 0]
    end
end

context("props! and isdegenerate") do
    @fact_throws AssertionError isdegenerate([1., 1.], [1., 1., 1.], [1., 1., 1.])
    @fact_throws AssertionError isdegenerate([1., 1., 1.], [1., 1.], [1., 1., 1.])
    @fact_throws AssertionError isdegenerate([1., 1., 1.], [1., 1., 1.], [1., 1.])
    for T in testtypes
        # degenerate triangles
        elem = Triangle(T[0, 0, 0], T[0, 0, 0], T[1, 0, 0])
        @fact_throws AssertionError props!(elem)
        elem = Triangle(T[0, 0, 0], T[0, 0, 1], T[0, 0, 0])
        @fact_throws AssertionError props!(elem)
        elem = Triangle(T[1, 0, 0], T[0, 0, 0], T[0, 0, 0])
        @fact_throws AssertionError props!(elem)
        elem = Triangle(T[0, 1, 0], T[0, 2, 0], T[0, 3, 0])
        @fact_throws AssertionError props!(elem)
        # simple 2D triangle
        elem = Triangle(T[0, 0, 0], T[0, 0, 3], T[0, 3, 0])
        props!(elem)
        @fact typeof(elem.center) --> Vector{T}
        @fact length(elem.center) --> 3
        @fact typeof(elem.normal) --> Vector{T}
        @fact length(elem.normal) --> 3
        @fact typeof(elem.distorig) --> T
        @fact typeof(elem.area) --> T
        @fact elem.center --> roughly([0, 1, 1])
        @fact elem.normal --> roughly([-1, 0, 0])
        @fact elem.distorig --> roughly(0)
        @fact elem.area --> roughly(4.5)
        # simple 3D triangle
        elem = Triangle(T[3, 0, 0], T[0, 4, 0], T[0, 0, 5])
        props!(elem)
        @fact typeof(elem.center) --> Vector{T}
        @fact length(elem.center) --> 3
        @fact typeof(elem.normal) --> Vector{T}
        @fact length(elem.normal) --> 3
        @fact typeof(elem.distorig) --> T
        @fact typeof(elem.area) --> T
        @fact elem.center --> roughly(map(T, ([1, 4/3, 5/3])))
        @fact elem.normal --> roughly(map(T, (√769 \ [20, 15, 12])))
        @fact elem.distorig --> roughly(map(T, 60/√769))
        @fact elem.area --> roughly(T(√769/2))
    end
end

context("seek") do
    fname, fh = mktemp()
    write(fh, "Lorem\nipsum\ndolor\nsit\namet.\n")
    seekstart(fh)
    seek(fh, "Foo")
    @fact eof(fh) --> true
    seekstart(fh)
    seek(fh, "Foo", false)
    @fact eof(fh) --> true
    seekstart(fh)
    seek(fh, "dolor")
    @fact readline(fh) --> "sit\n"
    seekstart(fh)
    seek(fh, "dolor", false)
    @fact readline(fh) --> "dolor\n"
    close(fh)
    rm(fname)
end

context("reverseindex") do
    for T in testtypes
        v1 = T[1, 2, 3]
        v2 = T[4, 5, 6]
        v3 = T[7, 8, 9]
        d = reverseindex(Vector{T}[])
        @fact typeof(d) --> Dict{UInt, UInt}
        @fact d --> Dict()
        d = reverseindex(Vector{T}[v1, v2, v3])
        @fact typeof(d) --> Dict{UInt, UInt}
        @fact length(d) --> 3
        @fact d[object_id(v1)] --> 1
        @fact d[object_id(v2)] --> 2
        @fact d[object_id(v3)] --> 3
        d = reverseindex(Vector{T}[v1, v1, v2])
        @fact typeof(d) --> Dict{UInt, UInt}
        @fact length(d) --> 2
        @fact d[object_id(v1)] --> 2
        @fact d[object_id(v2)] --> 3
        @fact_throws KeyError d[object_id(v3)]
    end
end

context("unpack") do
    for T in testtypes
        d = unpack(Vector{T}[T[1], T[2], T[3]], 0)
        @fact typeof(d) --> Vector{T}
        @fact d --> []
        d = unpack(Vector{T}[T[1], T[2], T[3]], 1)
        @fact typeof(d) --> Vector{T}
        @fact d --> [1, 2, 3]
        d = unpack(Vector{T}[T[1, 2], T[3, 4]], 1)
        @fact typeof(d) --> Vector{T}
        @fact d --> [1, 3]
        d = unpack(Vector{T}[T[1, 2], T[3, 4]], 2)
        @fact typeof(d) --> Vector{T}
        @fact d --> [1, 2, 3, 4]
        d = unpack(Vector{T}[T[1, 2, 3], T[4, 5, 6]])
        @fact typeof(d) --> Vector{T}
        @fact d --> [1, 2, 3, 4, 5, 6]
        d = unpack(Vector{T}[T[1, 2, 3, 4, 5, 6]], 6)
        @fact typeof(d) --> Vector{T}
        @fact d --> [1, 2, 3, 4, 5, 6]
        @fact_throws BoundsError unpack(Vector{T}[T[1], T[2], T[3]], 2)
    end
end

context("vertexnormals") do
    for T in testtypes
        nodes    = Vector{T}[T[0, 0, 0], T[0, 0, 3], T[0, 3, 0], T[1, -3, 3]]
        elements = [Triangle(nodes[1], nodes[2], nodes[3]),
                    Triangle(nodes[1], nodes[4], nodes[2])]
        map(props!, elements)
        d = vertexnormals(Vector{T}[], Triangle{T}[])
        @fact typeof(d) --> Vector{Vector{T}}
        @fact d --> []
        d = vertexnormals(Vector{T}[nodes[1], nodes[2], nodes[3]], Triangle{T}[elements[1]])
        @fact typeof(d) --> Vector{Vector{T}}
        @fact length(d) --> 3
        @fact d[1] --> roughly([-1, 0, 0])
        @fact d[2] --> roughly([-1, 0, 0])
        @fact d[3] --> roughly([-1, 0, 0])
        d = vertexnormals(Vector{T}[nodes[1], nodes[2], nodes[4]], Triangle{T}[elements[2]])
        @fact typeof(d) --> Vector{Vector{T}}
        @fact length(d) --> 3
        @fact d[1] --> roughly(map(T, √90 \ [-9, -3, 0]))
        @fact d[2] --> roughly(map(T, √90 \ [-9, -3, 0]))
        @fact d[3] --> roughly(map(T, √90 \ [-9, -3, 0]))
        d = vertexnormals(nodes, elements)
        @fact typeof(d) --> Vector{Vector{T}}
        @fact length(d) --> 4
        @fact d[1] --> roughly(map(T, √360 \ [-9 - √90, -3, 0]))
        @fact d[2] --> roughly(map(T, √360 \ [-9 - √90, -3, 0]))
        @fact d[3] --> roughly([-1, 0, 0])
        @fact d[4] --> roughly(map(T, √90 \ [-9, -3, 0]))
        d = vertexnormals(nodes, elements, true)
        @fact typeof(d) --> Vector{Vector{T}}
        @fact length(d) --> 4
        @fact d[1] --> roughly(map(T, -√360 \ [-9 - √90, -3, 0]))
        @fact d[2] --> roughly(map(T, -√360 \ [-9 - √90, -3, 0]))
        @fact d[3] --> roughly([1, 0, 0])
        @fact d[4] --> roughly(map(T, -√90 \ [-9, -3, 0]))
    end
end

@pending cos --> :nothing
@pending cathetus --> :nothing
@pending sign --> :nothing
@pending distance --> :nothing
