using ProteinES: eye!, pluseye!, isdegenerate, seek, reverseindex, unpack, vertexnormals, cos, cathetus, sign, distance

context("eye! and pluseye!") do
    for T in (Int, testtypes...)
        m = -ones(T, 3, 3)
        pluseye!(m)
        @fact typeof(m) --> Array{T, 2}
        @fact size(m) --> (3, 3)
        @fact m --> [0 -1 -1; -1 0 -1; -1 -1 0]
        pluseye!(m, 2)
        @fact typeof(m) --> Array{T, 2}
        @fact size(m) --> (3, 3)
        @fact m --> [2 -1 -1; -1 2 -1; -1 -1 2]
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
    @fact_throws AssertionError isdegenerate(Triangle([1., 1.], [1., 1., 1.], [1., 1., 1.]))
    @fact_throws AssertionError isdegenerate(Triangle([1., 1., 1.], [1., 1.], [1., 1., 1.]))
    @fact_throws AssertionError isdegenerate(Triangle([1., 1., 1.], [1., 1., 1.], [1., 1.]))
    for T in testtypes
        # degenerate triangles
        @fact_throws AssertionError Triangle(T[0, 0, 0], T[0, 0, 0], T[1, 0, 0])
        @fact_throws AssertionError Triangle(T[0, 0, 0], T[0, 0, 1], T[0, 0, 0])
        @fact_throws AssertionError Triangle(T[1, 0, 0], T[0, 0, 0], T[0, 0, 0])
        @fact_throws AssertionError Triangle(T[0, 1, 0], T[0, 2, 0], T[0, 3, 0])
        @fact_throws AssertionError Triangle(T[0, 1, 0], T[1e-11, 1, 0], T[0, 3, 0])
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
    @fact readline(fh) --> "sit"
    seekstart(fh)
    seek(fh, "dolor", false)
    @fact readline(fh) --> "dolor"
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
        d = unpack(Vector{T}[])
        @fact typeof(d) --> Vector{T}
        @fact d --> []
        d = unpack(Vector{T}[T[1], T[2], T[3]])
        @fact typeof(d) --> Vector{T}
        @fact d --> [1, 2, 3]
        d = unpack(Vector{T}[T[1, 2], T[3, 4]])
        @fact typeof(d) --> Vector{T}
        @fact d --> [1, 2, 3, 4]
        d = unpack(Vector{T}[T[1, 2, 3], T[4, 5, 6]])
        @fact typeof(d) --> Vector{T}
        @fact d --> [1, 2, 3, 4, 5, 6]
        d = unpack(Vector{T}[T[1, 2, 3, 4, 5, 6]])
        @fact typeof(d) --> Vector{T}
        @fact d --> [1, 2, 3, 4, 5, 6]
    end
end

context("vertexnormals") do
    for T in testtypes
        nodes    = Vector{T}[T[0, 0, 0], T[0, 0, 3], T[0, 3, 0], T[1, -3, 3]]
        elements = [Triangle(nodes[1], nodes[2], nodes[3]),
                    Triangle(nodes[1], nodes[4], nodes[2])]
        map(props!, elements)
        model = SurfaceModel{T}()
        d = vertexnormals(model)
        @fact typeof(d) --> Vector{Vector{T}}
        @fact d --> []
        model.nodes = Vector{T}[nodes[1], nodes[2], nodes[3]]
        model.elements = Triangle{T}[elements[1]]
        d = vertexnormals(model)
        @fact typeof(d) --> Vector{Vector{T}}
        @fact length(d) --> 3
        @fact d[1] --> roughly([-1, 0, 0])
        @fact d[2] --> roughly([-1, 0, 0])
        @fact d[3] --> roughly([-1, 0, 0])
        model.nodes = Vector{T}[nodes[1], nodes[2], nodes[4]]
        model.elements = Triangle{T}[elements[2]]
        d = vertexnormals(model)
        @fact typeof(d) --> Vector{Vector{T}}
        @fact length(d) --> 3
        @fact d[1] --> roughly(map(T, √90 \ [-9, -3, 0]))
        @fact d[2] --> roughly(map(T, √90 \ [-9, -3, 0]))
        @fact d[3] --> roughly(map(T, √90 \ [-9, -3, 0]))
        model.nodes = nodes
        model.elements = elements
        d = vertexnormals(model)
        @fact typeof(d) --> Vector{Vector{T}}
        @fact length(d) --> 4
        @fact d[1] --> roughly(map(T, √360 \ [-9 - √90, -3, 0]))
        @fact d[2] --> roughly(map(T, √360 \ [-9 - √90, -3, 0]))
        @fact d[3] --> roughly([-1, 0, 0])
        @fact d[4] --> roughly(map(T, √90 \ [-9, -3, 0]))
    end
end

context("cos") do
    for T in testtypes
        u = T[1, 0, 0]
        @fact typeof(cos(u, T[1, 0, 0])) --> T
        @fact cos(u, T[1, 0, 0]) --> roughly(one(T))
        @fact cos(u, T[0, 1, 0]) --> roughly(zero(T))
        @fact cos(u, T[0, 1, 0]) --> roughly(zero(T))
        @fact cos(u, T[-1, 0, 0]) --> roughly(-one(T))
        @fact cos(u, T[0, -1, 0]) --> roughly(zero(T))
        @fact cos(u, T[1, 1, 0]) --> roughly(T(cos(π/4)))
        @fact cos(u, T[-1, 1, 0]) --> roughly(T(cos(3π/4)))
        @fact cos(u, T[1, 0, 1]) --> roughly(T(cos(π/4)))
    end
end

context("meshunion") do
    for T in testtypes
        # empty system
        empty = VolumeModel{T}()
        res = meshunion(empty, empty)
        @fact typeof(res.nodes) --> Vector{Vector{T}}
        @fact typeof(res.elements) --> Vector{Tetrahedron{T}}
        @fact length(res.nodes) --> 0
        @fact length(res.elements) --> 0

        # small system
        nodesΩ = Vector{T}[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
        modelΩ = VolumeModel{T}(
            nodesΩ,
            [Tetrahedron(nodesΩ...), Tetrahedron(nodesΩ...)]
        )
        nodesΣ = Vector{T}[[0, 0, 0], [-1, 0, 0], [0, -1, 0], [0, 0, -1]]
        modelΣ = VolumeModel{T}(
            nodesΣ,
            [Tetrahedron(nodesΣ...)]
        )
        res = meshunion(modelΩ, modelΣ)
        oids = Set([object_id(e) for e in res.nodes])
        @fact typeof(res.nodes) --> Vector{Vector{T}}
        @fact typeof(res.elements) --> Vector{Tetrahedron{T}}
        @fact length(res.nodes) --> 7
        @fact length(res.elements) --> 3
        for node in nodesΩ ∪ nodesΣ
            @fact node ∈ res.nodes --> true
        end
        @fact modelΩ.elements[1] ∈ res.elements --> true
        @fact modelΩ.elements[2] ∈ res.elements --> true
        @fact modelΣ.elements[1] ∈ res.elements --> true
        @fact modelΣ.elements[1].v1 --> exactly(nodesΩ[1])
        for elem in res.elements
            for v in (elem.v1, elem.v2, elem.v3, elem.v4)
                v ∈ oids
            end
        end
    end
end

context("obspoints_line") do
    for T in testtypes
        Ξ = [ξ for ξ in obspoints_line(T[1, 1, 1], T[1, 1, 1], 3)]
        @fact typeof(Ξ) --> Vector{Vector{T}}
        @fact length(Ξ) --> 3
        for ξ in Ξ
            @fact ξ --> T[1, 1, 1]
        end
        Ξ = [ξ for ξ in obspoints_line(T[0, 0, 0], T[2, 2, 2], 3)]
        @fact typeof(Ξ) --> Vector{Vector{T}}
        @fact length(Ξ) --> 3
        @fact Ξ[1] --> T[0, 0, 0]
        @fact Ξ[2] --> T[1, 1, 1]
        @fact Ξ[3] --> T[2, 2, 2]
        Ξ = [ξ for ξ in obspoints_line(T[2, 2, 2], T[0, 0, 0], 3)]
        @fact typeof(Ξ) --> Vector{Vector{T}}
        @fact length(Ξ) --> 3
        @fact Ξ[3] --> T[0, 0, 0]
        @fact Ξ[2] --> T[1, 1, 1]
        @fact Ξ[1] --> T[2, 2, 2]
    end
end

context("obspoints_plane") do
    for T in testtypes
        Ξ = [collect(Ξ) for Ξ in obspoints_plane(T[1, 1, 1], T[1, 1, 1], T[1, 1, 1], 3, 3)]
        @fact typeof(Ξ) --> Vector{Vector{Vector{T}}}
        @fact length(Ξ) --> 3
        for row in Ξ
            @fact row --> [T[1, 1, 1], T[1, 1, 1], T[1, 1, 1]]
        end
        Ξ = [collect(Ξ) for Ξ in obspoints_plane(T[1, 0, 0], T[0, 0, -1], T[0, 1, 0], 3, 3)]
        @fact typeof(Ξ) --> Vector{Vector{Vector{T}}}
        @fact length(Ξ) --> 3
        @fact Ξ[1] --> [T[0, 0, -1], T[0, 0.5, -0.5], T[0, 1, 0]]
        @fact Ξ[2] --> [T[0.5, 0, -0.5], T[0.5, 0.5, 0], T[0.5, 1, 0.5]]
        @fact Ξ[3] --> [T[1, 0, 0], T[1, 0.5, 0.5], T[1, 1, 1]]
        Ξ = [collect(Ξ) for Ξ in obspoints_plane(T[1, 0, 0], T[1, 1, 1], T[0, 1, 0], 3, 3)]
        @fact typeof(Ξ) --> Vector{Vector{Vector{T}}}
        @fact length(Ξ) --> 3
        @fact Ξ[1] --> [T[1, 1, 1], T[0.5, 1, 0.5], T[0, 1, 0]]
        @fact Ξ[2] --> [T[1, 0.5, 0.5], T[0.5, 0.5, 0], T[0, 0.5, -0.5]]
        @fact Ξ[3] --> [T[1, 0, 0], T[0.5, 0, -0.5], T[0, 0, -1]]
    end
end

@pending cathetus --> :nothing
@pending sign --> :nothing
@pending distance --> :nothing
@pending ddot --> :nothing
