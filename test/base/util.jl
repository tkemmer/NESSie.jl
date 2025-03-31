@testitem "Utilities" begin
    include("../testsetup.jl")

    import GeometryBasics
    using LinearAlgebra: ⋅, norm

    @testset "_cos" begin
        using NESSie: _cos

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
        using NESSie: _dot

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
        using NESSie: _norm

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
        using NESSie: _seek

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
        using NESSie: distance

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

    @testset "_pluseye!" begin
        using NESSie: _pluseye!

        for T in (Int, testtypes...)
            m = -ones(T, 3, 3)
            @test _pluseye!(m) === m
            @test typeof(m) == Array{T, 2}
            @test size(m) == (3, 3)
            @test m == [0 -1 -1; -1 0 -1; -1 -1 0]
            @test _pluseye!(m, T(2)) === m
            @test typeof(m) == Array{T, 2}
            @test size(m) == (3, 3)
            @test m == [2 -1 -1; -1 2 -1; -1 -1 2]
        end
    end

    @testset "props and isdegenerate" begin
        using NESSie: props, isdegenerate

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

    @testset "_reverseindex" begin
        using NESSie: _reverseindex

        for T in testtypes
            v1 = T[1, 2, 3]
            v2 = T[4, 5, 6]
            v3 = T[7, 8, 9]
            d = _reverseindex(Vector{T}[])
            @test typeof(d) == IdDict{Vector{T}, Int}
            @test d == IdDict()
            d = _reverseindex(Vector{T}[v1, v2, v3])
            @test typeof(d) == IdDict{Vector{T}, Int}
            @test length(d) == 3
            @test d[v1] == 1
            @test d[v2] == 2
            @test d[v3] == 3
            d = _reverseindex(Vector{T}[v1, v1, v2])
            @test typeof(d) == IdDict{Vector{T}, Int}
            @test length(d) == 2
            @test d[v1] == 2
            @test d[v2] == 3
            @test_throws KeyError d[v3]
        end
    end

    @testset "vertexnormals" begin
        using NESSie: vertexnormals

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

    @testset "guess_domain" begin
        for T in testtypes
            nodes    = Vector{T}[T[0, 0, 0], T[0, 0, 3], T[0, 3, 0], T[3, 0, 0]]
            elements = [Triangle(nodes[1], nodes[2], nodes[3]),
                        Triangle(nodes[1], nodes[4], nodes[2]),
                        Triangle(nodes[1], nodes[3], nodes[4]),
                        Triangle(nodes[2], nodes[4], nodes[3])]

            let model = Model(nodes, elements)
                for node in model.nodes
                    @test guess_domain(node, model) === :Γ
                end
                for elem in model.elements
                    @test guess_domain(elem.center, model) === :Γ
                end
                @test guess_domain(T[1.5, 0, 0], model) === :Γ
                @test guess_domain(T[0.00001, 0.00001, 0.00001], model; surface_margin = T(1e-4)) === :Γ
                @test guess_domain(T[0.00001, 0.00001, 0.00001], model) === :Ω
                @test guess_domain(T[0.001, 0.001, 0.001], model) === :Ω
                @test guess_domain(T[3//4, 3//4, 3//4], model) === :Ω
                @test guess_domain(T[-0.00001, -0.00001, -0.00001], model; surface_margin = T(1e-4)) === :Γ
                @test guess_domain(T[-0.00001, -0.00001, -0.00001], model) === :Σ
                @test guess_domain(T[-0.001, -0.001, -0.001], model) === :Σ
                @test guess_domain(T[-1, -1, -1], model) === :Σ
            end
        end
    end

    @testset "GeometryBasics.mesh" begin
        for T in testtypes
            let model = Model{T, NESSie.Triangle{T}}()
                gbm = GeometryBasics.mesh(model)
                @test gbm isa GeometryBasics.Mesh{3, T, GeometryBasics.TriangleFace{Int}}
                @test isempty(gbm.position)
                @test isempty(gbm.faces)
            end

            let model = Format.readoff(NESSie._data_path("born/na.off"), T)
                nodes = Set((v...,) for v in model.nodes)
                elems = Set((e.v1..., e.v2..., e.v3...) for e in model.elements)

                gbm = GeometryBasics.mesh(model)
                @test gbm isa GeometryBasics.Mesh{3, T, GeometryBasics.TriangleFace{Int}}
                @test length(gbm.position) == length(model.nodes)
                @test Set((p...,) for p in gbm.position) == nodes
                @test length(gbm.faces) == length(model.elements)
                @test Set(Tuple(c for p in gbm.position[[f...]] for c in p) for f in gbm.faces) == elems
            end
        end
    end

    @testset "NESSie.Model" begin
        for T in testtypes
            let model = Model{T, Triangle{T}}()
                model2 = Model(GeometryBasics.mesh(model))
                @test model2 isa Model{T, Triangle{T}}
                @test isempty(model2.nodes)
                @test isempty(model2.elements)
                @test isempty(model2.charges)
                @test model2.params == defaultopt(T)
            end

            let model = Format.readoff(NESSie._data_path("born/na.off"), T)
                model.charges = Format.readpqr(NESSie._data_path("born/na.pqr"), T)
                model.params  = Option{T}(1, 2, 3, 4)

                model2 = Model(GeometryBasics.mesh(model); charges = model.charges, params = model.params)
                @test model2 isa Model{T, Triangle{T}}
                @test model2.nodes == model.nodes
                @test model2.elements == model.elements
                @test model2.charges == model.charges
                @test model2.params == model.params
            end
        end
    end

    @testset "GeometryBasics.Rect" begin
        for T in testtypes
            let model = Model{T, NESSie.Triangle{T}}()
                box = GeometryBasics.Rect(model)
                @test box isa GeometryBasics.Rect3{T}
                @test box.origin == zeros(T, 3)
                @test box.widths == zeros(T, 3)

                box = GeometryBasics.Rect(model; padding = one(T))
                @test box isa GeometryBasics.Rect3{T}
                @test box.origin == -ones(T, 3)
                @test box.widths == 2 .* ones(T, 3)
            end

            let model = Format.readoff(NESSie._data_path("born/na.off"), T)
                box = GeometryBasics.Rect(model)
                @test box isa GeometryBasics.Rect3{T}
                @test box.origin ≈ fill(T(-1.0049999952), 3)
                @test box.widths ≈ fill(T(2 * 1.0049999952), 3)

                box = GeometryBasics.Rect(model; padding = one(T))
                @test box isa GeometryBasics.Rect3{T}
                @test box.origin ≈ fill(T(-2.0049999952), 3)
                @test box.widths ≈ fill(T(2 * 1.0049999952 + 2), 3)
            end
        end
    end

    @test_skip _data_path
    @test_skip _sign
    @test_skip cathetus
    @test_skip ddot
end
