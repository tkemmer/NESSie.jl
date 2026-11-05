@testset "readobj" begin
    for T in testtypes
        # empty file
        let model = Format.readobj(nessie_data_path("test/empty"), T)
            @test model isa Model{T, Triangle{T}}
            @test length(model.nodes) == 0
            @test length(model.elements) == 0
        end

        # small file
        let model = Format.readobj(nessie_data_path("test/simple.obj"), T)
            @test model isa Model{T}
            nodes, elements = (model.nodes, model.elements)

            @test length(nodes) == 3
            @test length(elements) == 2

            @test nodes[1] == T[1, 0, 0]
            @test nodes[2] == T[0, 1, 0]
            @test nodes[3] == T[0, 0, 1]
            @test elements[1].v1 === nodes[1]
            @test elements[1].v2 === nodes[2]
            @test elements[1].v3 === nodes[3]
            @test elements[2].v1 === nodes[3]
            @test elements[2].v2 === nodes[1]
            @test elements[2].v3 === nodes[2]
        end
    end
end

@testset "writeobj" begin
    @test_skip writeobj
end
