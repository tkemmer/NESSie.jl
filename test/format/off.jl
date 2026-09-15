@testset "readoff" begin
    for T in testtypes
        # empty file
        let model = Format.readoff(nessie_data_path("test/empty"), T)
            @test model isa Model{T}
            @test length(model.nodes) == 0
            @test length(model.elements) == 0
        end

        # small file
        let model = Format.readoff(nessie_data_path("test/simple.off"), T)
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

            # input stream
            model2 = open(io -> Format.readoff(io, T), nessie_data_path("test/simple.off"))
            @test model2 == model
        end
    end
end

@testset "writeoff" begin
    for T in testtypes
        # empty model
        model = Model{T, Triangle{T}}()
        mktemp() do fname, fh
            Format.writeoff(fh, model)
            model2 = Format.readoff(fname, T)
            @test model == model2
        end

        # small model
        nodes = Vector{T}[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
        elements = Triangle{T}[
            Triangle(nodes[1], nodes[2], nodes[3]),
            Triangle(nodes[1], nodes[3], nodes[4])
        ]
        model = Model(nodes, elements)

        mktemp() do fname, fh
            Format.writeoff(fh, model)
            model2 = Format.readoff(fname, T)
            @test model == model2
        end

        mktemp() do fname, _
            Format.writeoff(fname, model)
            model2 = Format.readoff(fname, T)
            @test model == model2
        end
    end
end
