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

            let model = Format.readoff("../../data/born/na.off", T)
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
end
