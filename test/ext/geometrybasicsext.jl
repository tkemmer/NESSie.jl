@testitem "GeometryBasicsExt" begin
    include("../testsetup.jl")

    using GeometryBasics

    @testset "GeometryBasics.mesh" begin
        for T in testtypes
            let model = Model{T, NESSie.Triangle{T}}()
                gbm = GeometryBasics.mesh(model)
                @test gbm isa GeometryBasics.Mesh{3, T, TriangleFace{Int}}
                @test isempty(gbm.position)
                @test isempty(gbm.faces)
            end

            let model = Format.readoff(NESSie._data_path("born/na.off"), T)
                nodes = Set((v...,) for v in model.nodes)
                elems = Set((e.v1..., e.v2..., e.v3...) for e in model.elements)

                gbm = GeometryBasics.mesh(model)
                @test gbm isa GeometryBasics.Mesh{3, T, TriangleFace{Int}}
                @test length(gbm.position) == length(model.nodes)
                @test Set((p...,) for p in gbm.position) == nodes
                @test length(gbm.faces) == length(model.elements)
                @test Set(Tuple(c for p in gbm.position[[f...]] for c in p) for f in gbm.faces) == elems
            end
        end
    end

    @testset "NESSie.Model" begin
        for T in testtypes
            let model = Model{T, NESSie.Triangle{T}}()
                model2 = Model(GeometryBasics.mesh(model))
                @test model2 isa Model{T, NESSie.Triangle{T}}
                @test isempty(model2.nodes)
                @test isempty(model2.elements)
                @test isempty(model2.charges)
                @test model2.params == defaultopt(T)
            end

            let model = Format.readoff(NESSie._data_path("born/na.off"), T)
                model.charges = Format.readpqr(NESSie._data_path("born/na.pqr"), T)
                model.params  = Option{T}(1, 2, 3, 4)

                model2 = Model(GeometryBasics.mesh(model); charges = model.charges, params = model.params)
                @test model2 isa Model{T, NESSie.Triangle{T}}
                @test model2.nodes == model.nodes
                @test model2.elements == model.elements
                @test model2.charges == model.charges
                @test model2.params == model.params
            end
        end
    end
end
