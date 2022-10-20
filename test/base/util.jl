using NESSie: _cos, _dot, _norm, _seek, eye!, pluseye!, isdegenerate, reverseindex, unpack,
    vertexnormals, distance
using LinearAlgebra: ⋅, norm

@testset "_cos" begin
    for T in testtypes
        u = T[1, 0, 0]
        @test typeof(_cos(u, T[1, 0, 0])) == T
        @test _cos(u, T[1, 0, 0]) ≈ one(T)
        @test _cos(u, T[0, 1, 0]) ≈ zero(T)
        @test _cos(u, T[-1, 0, 0]) ≈ -one(T)
        @test _cos(u, T[0, -1, 0]) ≈ zero(T)
        @test _cos(u, T[1, 1, 0]) ≈ T(cos(π/4))
        @test _cos(u, T[-1, 1, 0]) ≈ T(cos(3π/4))
        @test _cos(u, T[1, 0, 1]) ≈ T(cos(π/4))
    end
end

@testset "_dot" begin
    for T in testtypes
        x1 = T[1, 2, 3]
        x2 = T[4, 6, 8]
        x3 = T[-1, -1, -1]

        @test typeof(_dot(x1, x2)) == T
        @test _dot(x1, x2) ≈ x1 ⋅ x2
        @test _dot(x2, x1) ≈ x2 ⋅ x1
        @test _dot(x1, x3) ≈ x1 ⋅ x3
        @test _dot(x3, x1) ≈ x3 ⋅ x1
        @test _dot(x1, zeros(T, 3)) ≈ zero(T)
        @test _dot(zeros(T, 3), x1) ≈ zero(T)
    end
end

@testset "_norm" begin
    for T in testtypes
        x1 = T[1, 2, 3]
        x2 = T[4, 6, 8]
        x3 = T[-1, -1, -1]

        @test typeof(_norm(x1)) == T
        @test _norm(x1) ≈ norm(x1)
        @test _norm(x2) ≈ norm(x2)
        @test _norm(x3) ≈ norm(x3)
        @test _norm(zeros(T, 3)) ≈ zero(T)
    end
end

@testset "_seek" begin
    fname, fh = mktemp()
    write(fh, "Lorem\nipsum\ndolor\nsit\namet.\n")
    seekstart(fh)
    _seek(fh, "Foo")
    @test eof(fh)
    seekstart(fh)
    _seek(fh, "Foo", false)
    @test eof(fh)
    seekstart(fh)
    _seek(fh, "dolor")
    @test readline(fh) == "sit"
    seekstart(fh)
    _seek(fh, "dolor", false)
    @test readline(fh) == "dolor"
    close(fh)
    rm(fname)
end

@testset "distance" begin
    for T in testtypes
        elem = Triangle(T[0, 0, 0], T[0, 1, 0], T[0, 0, 1])

        @test typeof(distance(elem.v1, elem)) == T
        @test distance(elem.v1, elem) ≈ zero(T)
        @test distance(elem.v2, elem) ≈ zero(T)
        @test distance(elem.v3, elem) ≈ zero(T)
        @test distance(T[1, 0, 0], elem) ≈ one(T)
        @test distance(T[1, 1, 0], elem) ≈ one(T)
        @test distance(T[1, 0, 1], elem) ≈ one(T)
        @test distance(T[-1, 0, 0], elem) ≈ -1 * one(T)
        @test distance(T[-2, 1, 0], elem) ≈ -2 * one(T)
        @test distance(T[-3, 0, 1], elem) ≈ -3 * one(T)
    end
end

@testset "eye! and pluseye!" begin
    for T in (Int, testtypes...)
        m = -ones(T, 3, 3)
        pluseye!(m)
        @test typeof(m) == Array{T, 2}
        @test size(m) == (3, 3)
        @test m == [0 -1 -1; -1 0 -1; -1 -1 0]
        pluseye!(m, 2)
        @test typeof(m) == Array{T, 2}
        @test size(m) == (3, 3)
        @test m == [2 -1 -1; -1 2 -1; -1 -1 2]
        eye!(m)
        @test typeof(m) == Array{T, 2}
        @test size(m) == (3, 3)
        @test m == [1 0 0; 0 1 0; 0 0 1]
        eye!(m, 2)
        @test typeof(m) == Array{T, 2}
        @test size(m) == (3, 3)
        @test m == [2 0 0; 0 2 0; 0 0 2]
        m = zeros(T, 2, 3)
        eye!(m)
        @test typeof(m) == Array{T, 2}
        @test size(m) == (2, 3)
        @test m == [1 0 0; 0 1 0]
        eye!(m, 2)
        @test typeof(m) == Array{T, 2}
        @test size(m) == (2, 3)
        @test m == [2 0 0; 0 2 0]
        m = zeros(T, 3, 2)
        eye!(m)
        @test typeof(m) == Array{T, 2}
        @test size(m) == (3, 2)
        @test m == [1 0; 0 1; 0 0]
        eye!(m, 2)
        @test typeof(m) == Array{T, 2}
        @test size(m) == (3, 2)
        @test m == [2 0; 0 2; 0 0]
    end
end

@testset "props and isdegenerate" begin
    @test_throws AssertionError isdegenerate(Triangle([1., 1.], [1., 1., 1.], [1., 1., 1.]))
    @test_throws AssertionError isdegenerate(Triangle([1., 1., 1.], [1., 1.], [1., 1., 1.]))
    @test_throws AssertionError isdegenerate(Triangle([1., 1., 1.], [1., 1., 1.], [1., 1.]))
    for T in testtypes
        # degenerate triangles
        @test_throws AssertionError Triangle(T[0, 0, 0], T[0, 0, 0], T[1, 0, 0])
        @test_throws AssertionError Triangle(T[0, 0, 0], T[0, 0, 1], T[0, 0, 0])
        @test_throws AssertionError Triangle(T[1, 0, 0], T[0, 0, 0], T[0, 0, 0])
        @test_throws AssertionError Triangle(T[0, 1, 0], T[0, 2, 0], T[0, 3, 0])
        @test_throws AssertionError Triangle(T[0, 1, 0], T[1e-11, 1, 0], T[0, 3, 0])
        # simple 2D triangle
        elem = Triangle(T[0, 0, 0], T[0, 0, 3], T[0, 3, 0])
        @test typeof(elem.center) == Vector{T}
        @test length(elem.center) == 3
        @test typeof(elem.normal) == Vector{T}
        @test length(elem.normal) == 3
        @test typeof(elem.distorig) == T
        @test typeof(elem.area) == T
        @test elem.center ≈ [0, 1, 1]
        @test elem.normal ≈ [-1, 0, 0]
        @test elem.distorig ≈ 0
        @test elem.area ≈ 4.5
        # simple 3D triangle
        elem = Triangle(T[3, 0, 0], T[0, 4, 0], T[0, 0, 5])
        @test typeof(elem.center) == Vector{T}
        @test length(elem.center) == 3
        @test typeof(elem.normal) == Vector{T}
        @test length(elem.normal) == 3
        @test typeof(elem.distorig) == T
        @test typeof(elem.area) == T
        @test elem.center ≈ map(T, ([1, 4/3, 5/3]))
        @test elem.normal ≈ map(T, (√769 \ [20, 15, 12]))
        @test elem.distorig ≈ map(T, 60/√769)
        @test elem.area ≈ T(√769/2)
    end
end

@testset "reverseindex" begin
    for T in testtypes
        v1 = T[1, 2, 3]
        v2 = T[4, 5, 6]
        v3 = T[7, 8, 9]
        d = reverseindex(Vector{T}[])
        @test typeof(d) == Dict{UInt, UInt}
        @test d == Dict()
        d = reverseindex(Vector{T}[v1, v2, v3])
        @test typeof(d) == Dict{UInt, UInt}
        @test length(d) == 3
        @test d[objectid(v1)] == 1
        @test d[objectid(v2)] == 2
        @test d[objectid(v3)] == 3
        d = reverseindex(Vector{T}[v1, v1, v2])
        @test typeof(d) == Dict{UInt, UInt}
        @test length(d) == 2
        @test d[objectid(v1)] == 2
        @test d[objectid(v2)] == 3
        @test_throws KeyError d[objectid(v3)]
    end
end

@testset "unpack" begin
    for T in testtypes
        d = unpack(Vector{T}[])
        @test typeof(d) == Vector{T}
        @test d == []
        d = unpack(Vector{T}[T[1], T[2], T[3]])
        @test typeof(d) == Vector{T}
        @test d == [1, 2, 3]
        d = unpack(Vector{T}[T[1, 2], T[3, 4]])
        @test typeof(d) == Vector{T}
        @test d == [1, 2, 3, 4]
        d = unpack(Vector{T}[T[1, 2, 3], T[4, 5, 6]])
        @test typeof(d) == Vector{T}
        @test d == [1, 2, 3, 4, 5, 6]
        d = unpack(Vector{T}[T[1, 2, 3, 4, 5, 6]])
        @test typeof(d) == Vector{T}
        @test d == [1, 2, 3, 4, 5, 6]
    end
end

@testset "vertexnormals" begin
    for T in testtypes
        nodes    = Vector{T}[T[0, 0, 0], T[0, 0, 3], T[0, 3, 0], T[1, -3, 3]]
        elements = [Triangle(nodes[1], nodes[2], nodes[3]),
                    Triangle(nodes[1], nodes[4], nodes[2])]
        d = vertexnormals(Model{T, Triangle{T}}())
        @test typeof(d) == Vector{Vector{T}}
        @test d == []
        d = vertexnormals(Model(
            Vector{T}[nodes[1], nodes[2], nodes[3]],
            Triangle{T}[elements[1]]
        ))
        @test typeof(d) == Vector{Vector{T}}
        @test length(d) == 3
        @test d[1] ≈ [-1, 0, 0]
        @test d[2] ≈ [-1, 0, 0]
        @test d[3] ≈ [-1, 0, 0]
        d = vertexnormals(Model(
            Vector{T}[nodes[1], nodes[2], nodes[4]],
            Triangle{T}[elements[2]]
        ))
        @test typeof(d) == Vector{Vector{T}}
        @test length(d) == 3
        @test d[1] ≈ map(T, √90 \ [-9, -3, 0])
        @test d[2] ≈ map(T, √90 \ [-9, -3, 0])
        @test d[3] ≈ map(T, √90 \ [-9, -3, 0])
        d = vertexnormals(Model(nodes, elements))
        @test typeof(d) == Vector{Vector{T}}
        @test length(d) == 4
        @test d[1] ≈ map(T, √360 \ [-9 - √90, -3, 0])
        @test d[2] ≈ map(T, √360 \ [-9 - √90, -3, 0])
        @test d[3] ≈ [-1, 0, 0]
        @test d[4] ≈ map(T, √90 \ [-9, -3, 0])
    end
end

@testset "meshunion" begin
    for T in testtypes
        # empty system
        empty = Model{T, Tetrahedron{T}}()
        res = meshunion(empty, empty)
        @test typeof(res.nodes) == Vector{Vector{T}}
        @test typeof(res.elements) == Vector{Tetrahedron{T}}
        @test length(res.nodes) == 0
        @test length(res.elements) == 0

        # small system
        nodesΩ = Vector{T}[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
        modelΩ = Model(
            nodesΩ,
            [Tetrahedron(nodesΩ...), Tetrahedron(nodesΩ...)]
        )
        nodesΣ = Vector{T}[[0, 0, 0], [-1, 0, 0], [0, -1, 0], [0, 0, -1]]
        modelΣ = Model(
            nodesΣ,
            [Tetrahedron(nodesΣ...)]
        )
        res = meshunion(modelΩ, modelΣ)
        oids = Set([objectid(e) for e in res.nodes])
        @test typeof(res) == Model{T, Tetrahedron{T}}
        @test length(res.nodes) == 7
        @test length(res.elements) == 3
        for node in nodesΩ ∪ nodesΣ
            @test node ∈ res.nodes
        end
        @test res.elements[1] == modelΩ.elements[1]
        @test res.elements[2] == modelΩ.elements[2]
        @test res.elements[3].v1 === nodesΩ[1]
        @test res.elements[3].v2 === nodesΣ[2]
        @test res.elements[3].v3 === nodesΣ[3]
        @test res.elements[3].v4 === nodesΣ[4]
        for elem in res.elements
            for v in (elem.v1, elem.v2, elem.v3, elem.v4)
                v ∈ oids
            end
        end
    end
end

@testset "obspoints_line" begin
    for T in testtypes
        Ξ = [ξ for ξ in obspoints_line(T[1, 1, 1], T[1, 1, 1], 3)]
        @test typeof(Ξ) == Vector{Vector{T}}
        @test length(Ξ) == 3
        for ξ in Ξ
            @test ξ == T[1, 1, 1]
        end
        Ξ = [ξ for ξ in obspoints_line(T[0, 0, 0], T[2, 2, 2], 3)]
        @test typeof(Ξ) == Vector{Vector{T}}
        @test length(Ξ) == 3
        @test Ξ[1] == T[0, 0, 0]
        @test Ξ[2] == T[1, 1, 1]
        @test Ξ[3] == T[2, 2, 2]
        Ξ = [ξ for ξ in obspoints_line(T[2, 2, 2], T[0, 0, 0], 3)]
        @test typeof(Ξ) == Vector{Vector{T}}
        @test length(Ξ) == 3
        @test Ξ[3] == T[0, 0, 0]
        @test Ξ[2] == T[1, 1, 1]
        @test Ξ[1] == T[2, 2, 2]
    end
end

@testset "obspoints_plane" begin
    for T in testtypes
        Ξ = [collect(Ξ) for Ξ in obspoints_plane(T[1, 1, 1], T[1, 1, 1], T[1, 1, 1], 3, 3)]
        @test typeof(Ξ) == Vector{Vector{Vector{T}}}
        @test length(Ξ) == 3
        for row in Ξ
            @test row == [T[1, 1, 1], T[1, 1, 1], T[1, 1, 1]]
        end
        Ξ = [collect(Ξ) for Ξ in obspoints_plane(T[1, 0, 0], T[0, 0, -1], T[0, 1, 0], 3, 3)]
        @test typeof(Ξ) == Vector{Vector{Vector{T}}}
        @test length(Ξ) == 3
        @test Ξ[1] == [T[0, 0, -1], T[0, 0.5, -0.5], T[0, 1, 0]]
        @test Ξ[2] == [T[0.5, 0, -0.5], T[0.5, 0.5, 0], T[0.5, 1, 0.5]]
        @test Ξ[3] == [T[1, 0, 0], T[1, 0.5, 0.5], T[1, 1, 1]]
        Ξ = [collect(Ξ) for Ξ in obspoints_plane(T[1, 0, 0], T[1, 1, 1], T[0, 1, 0], 3, 3)]
        @test typeof(Ξ) == Vector{Vector{Vector{T}}}
        @test length(Ξ) == 3
        @test Ξ[1] == [T[1, 1, 1], T[0.5, 1, 0.5], T[0, 1, 0]]
        @test Ξ[2] == [T[1, 0.5, 0.5], T[0.5, 0.5, 0], T[0, 0.5, -0.5]]
        @test Ξ[3] == [T[1, 0, 0], T[0.5, 0, -0.5], T[0, 0, -1]]
    end
end

@test_skip _sign
@test_skip cathetus
@test_skip ddot
