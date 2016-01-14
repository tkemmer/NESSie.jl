import NonLocalES: eye!, isdegenerate, seek, reverseindex, unpack, vertexnormals, cos, cathetus, sign, distance
import JSON: parse

context("eye!") do
    for T in (Int, testtypes...)
        m = zeros(T, 3, 3)
        eye!(m)
        @fact isa(m, Array{T, 2}) --> true
        @fact size(m) --> (3, 3)
        @fact m --> [1 0 0; 0 1 0; 0 0 1]
        eye!(m, 2)
        @fact isa(m, Array{T, 2}) --> true
        @fact size(m) --> (3, 3)
        @fact m --> [2 0 0; 0 2 0; 0 0 2]
        m = zeros(T, 2, 3)
        eye!(m)
        @fact isa(m, Array{T, 2}) --> true
        @fact size(m) --> (2, 3)
        @fact m --> [1 0 0; 0 1 0]
        eye!(m, 2)
        @fact isa(m, Array{T, 2}) --> true
        @fact size(m) --> (2, 3)
        @fact m --> [2 0 0; 0 2 0]
        m = zeros(T, 3, 2)
        eye!(m)
        @fact isa(m, Array{T, 2}) --> true
        @fact size(m) --> (3, 2)
        @fact m --> [1 0; 0 1; 0 0]
        eye!(m, 2)
        @fact isa(m, Array{T, 2}) --> true
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
        @fact isa(elem.center, Vector{T}) --> true
        @fact length(elem.center) --> 3
        @fact isa(elem.normal, Vector{T}) --> true
        @fact length(elem.normal) --> 3
        @fact isa(elem.distorig, T) --> true
        @fact isa(elem.area, T) --> true
        @fact elem.center --> roughly([0, 1, 1])
        @fact elem.normal --> roughly([-1, 0, 0])
        @fact elem.distorig --> roughly(0)
        @fact elem.area --> roughly(4.5)
        # simple 3D triangle
        elem = Triangle(T[3, 0, 0], T[0, 4, 0], T[0, 0, 5])
        props!(elem)
        @fact isa(elem.center, Vector{T}) --> true
        @fact length(elem.center) --> 3
        @fact isa(elem.normal, Vector{T}) --> true
        @fact length(elem.normal) --> 3
        @fact isa(elem.distorig, T) --> true
        @fact isa(elem.area, T) --> true
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
        @fact isa(d, Dict{UInt, UInt}) --> true
        @fact d --> Dict()
        d = reverseindex(Vector{T}[v1, v2, v3])
        @fact isa(d, Dict{UInt, UInt}) --> true
        @fact length(d) --> 3
        @fact d[object_id(v1)] --> 1
        @fact d[object_id(v2)] --> 2
        @fact d[object_id(v3)] --> 3
        d = reverseindex(Vector{T}[v1, v1, v2])
        @fact isa(d, Dict{UInt, UInt}) --> true
        @fact length(d) --> 2
        @fact d[object_id(v1)] --> 2
        @fact d[object_id(v2)] --> 3
        @fact_throws KeyError d[object_id(v3)]
    end
end

context("unpack") do
    for T in testtypes
        d = unpack(Vector{T}[T[1], T[2], T[3]], 0)
        @fact isa(d, Vector{T}) --> true
        @fact d --> []
        d = unpack(Vector{T}[T[1], T[2], T[3]], 1)
        @fact isa(d, Vector{T}) --> true
        @fact d --> [1, 2, 3]
        d = unpack(Vector{T}[T[1, 2], T[3, 4]], 1)
        @fact isa(d, Vector{T}) --> true
        @fact d --> [1, 3]
        d = unpack(Vector{T}[T[1, 2], T[3, 4]], 2)
        @fact isa(d, Vector{T}) --> true
        @fact d --> [1, 2, 3, 4]
        d = unpack(Vector{T}[T[1, 2, 3], T[4, 5, 6]])
        @fact isa(d, Vector{T}) --> true
        @fact d --> [1, 2, 3, 4, 5, 6]
        d = unpack(Vector{T}[T[1, 2, 3, 4, 5, 6]], 6)
        @fact isa(d, Vector{T}) --> true
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
        @fact isa(d, Vector{Vector{T}}) --> true
        @fact d --> []
        d = vertexnormals(Vector{T}[nodes[1], nodes[2], nodes[3]], Triangle{T}[elements[1]])
        @fact isa(d, Vector{Vector{T}}) --> true
        @fact length(d) --> 3
        @fact d[1] --> roughly([-1, 0, 0])
        @fact d[2] --> roughly([-1, 0, 0])
        @fact d[3] --> roughly([-1, 0, 0])
        d = vertexnormals(Vector{T}[nodes[1], nodes[2], nodes[4]], Triangle{T}[elements[2]])
        @fact isa(d, Vector{Vector{T}}) --> true
        @fact length(d) --> 3
        @fact d[1] --> roughly(map(T, √90 \ [-9, -3, 0]))
        @fact d[2] --> roughly(map(T, √90 \ [-9, -3, 0]))
        @fact d[3] --> roughly(map(T, √90 \ [-9, -3, 0]))
        d = vertexnormals(nodes, elements)
        @fact isa(d, Vector{Vector{T}}) --> true
        @fact length(d) --> 4
        @fact d[1] --> roughly(map(T, √360 \ [-9 - √90, -3, 0]))
        @fact d[2] --> roughly(map(T, √360 \ [-9 - √90, -3, 0]))
        @fact d[3] --> roughly([-1, 0, 0])
        @fact d[4] --> roughly(map(T, √90 \ [-9, -3, 0]))
        d = vertexnormals(nodes, elements, true)
        @fact isa(d, Vector{Vector{T}}) --> true
        @fact length(d) --> 4
        @fact d[1] --> roughly(map(T, -√360 \ [-9 - √90, -3, 0]))
        @fact d[2] --> roughly(map(T, -√360 \ [-9 - √90, -3, 0]))
        @fact d[3] --> roughly([1, 0, 0])
        @fact d[4] --> roughly(map(T, -√90 \ [-9, -3, 0]))
    end
end

context("xml3d_mesh") do
    for T in testtypes
        # empty system
        js = parse(xml3d_mesh(Vector{T}[], Triangle{T}[]))
        @fact js["format"] --> "xml3d-json"
        @fact haskey(js, "version") --> true
        @fact js["data"]["index"]["type"] --> "int"
        @fact length(js["data"]["index"]["seq"]) --> 1
        @fact js["data"]["index"]["seq"][1]["value"] --> []
        @fact js["data"]["position"]["type"] --> "float3"
        @fact length(js["data"]["position"]["seq"]) --> 1
        @fact js["data"]["position"]["seq"][1]["value"] --> []
        @fact js["data"]["normal"]["type"] --> "float3"
        @fact length(js["data"]["normal"]["seq"]) --> 1
        @fact js["data"]["normal"]["seq"][1]["value"] --> []
        # small system
        nodes    = Vector{T}[T[0, 0, 0], T[0, 0, 3], T[0, 3, 0], T[1, -3, 3]]
        elements = [Triangle(nodes[1], nodes[2], nodes[3]),
                    Triangle(nodes[1], nodes[4], nodes[2])]
        map(props!, elements)
        js = parse(xml3d_mesh(nodes, elements))
        @fact js["format"] --> "xml3d-json"
        @fact haskey(js, "version") --> true
        @fact js["data"]["index"]["type"] --> "int"
        @fact length(js["data"]["index"]["seq"]) --> 1
        @fact js["data"]["index"]["seq"][1]["value"] --> [0, 1, 2, 0, 3, 1]
        @fact js["data"]["position"]["type"] --> "float3"
        @fact length(js["data"]["position"]["seq"]) --> 1
        @fact js["data"]["position"]["seq"][1]["value"] --> [0, 0, 0, 0, 0, 3, 0, 3, 0, 1, -3, 3]
        @fact js["data"]["normal"]["type"] --> "float3"
        @fact length(js["data"]["normal"]["seq"]) --> 1
        @fact map(T, js["data"]["normal"]["seq"][1]["value"]) --> roughly(map(T,
            [√360 \ [-9 - √90, -3, 0]; √360 \ [-9 - √90, -3, 0]; [-1, 0, 0]; √90 \ [-9, -3, 0]]
        ))
        # small system, inverted normals
        js = parse(xml3d_mesh(nodes, elements, true))
        @fact js["format"] --> "xml3d-json"
        @fact haskey(js, "version") --> true
        @fact js["data"]["index"]["type"] --> "int"
        @fact length(js["data"]["index"]["seq"]) --> 1
        @fact js["data"]["index"]["seq"][1]["value"] --> [0, 1, 2, 0, 3, 1]
        @fact js["data"]["position"]["type"] --> "float3"
        @fact length(js["data"]["position"]["seq"]) --> 1
        @fact js["data"]["position"]["seq"][1]["value"] --> [0, 0, 0, 0, 0, 3, 0, 3, 0, 1, -3, 3]
        @fact js["data"]["normal"]["type"] --> "float3"
        @fact length(js["data"]["normal"]["seq"]) --> 1
        @fact map(T, js["data"]["normal"]["seq"][1]["value"]) --> roughly(map(T,
            [-√360 \ [-9 - √90, -3, 0]; -√360 \ [-9 - √90, -3, 0]; [1, 0, 0]; -√90 \ [-9, -3, 0]]
        ))
    end
end

@pending cos --> :nothing
@pending cathetus --> :nothing
@pending sign --> :nothing
@pending distance --> :nothing
