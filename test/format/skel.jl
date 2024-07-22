@testitem "Format: SKEL/.skel" begin
    include("../testsetup.jl")

    @testset "writeskel (surface model)" begin
        for T in testtypes
            nodes = Vector{T}[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
            model = Model(
                nodes,
                Triangle{T}[
                    Triangle(nodes[1], nodes[2], nodes[3]),
                    Triangle(nodes[1], nodes[2], nodes[4])
                ]
            )
            (fname, fh) = mktemp()
            try
                Format.writeskel(fh, model)
                seek(fh, 0)
                @test readline(fh) == "SKEL"
                @test split(readline(fh)) == ["4", "2"]
                for node in nodes
                    @test [parse(T, e) for e in split(readline(fh))] == node
                end
                @test split(readline(fh)) == ["4", "0", "1", "2", "0"]
                @test split(readline(fh)) == ["4", "0", "1", "3", "0"]
                @test strip(readline(fh)) == ""
                @test eof(fh)
            finally
                close(fh)
                rm(fname)
            end
        end
    end

    @testset "writeskel (volume model)" begin
        for T in testtypes
            nodes = Vector{T}[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [-2, 0, 0]]
            model = Model(
                nodes,
                Tetrahedron{T}[
                    Tetrahedron(nodes[1], nodes[2], nodes[3], nodes[4]),
                    Tetrahedron(nodes[1], nodes[2], nodes[3], nodes[5])
                ]
            )
            (fname, fh) = mktemp()
            try
                Format.writeskel(fh, model)
                seek(fh, 0)
                @test readline(fh) == "SKEL"
                @test split(readline(fh)) == ["5", "2"]
                for node in nodes
                    @test [parse(T, e) for e in split(readline(fh))] == node
                end
                @test split(readline(fh)) == ["7", "0", "1", "2", "3", "1", "3", "0"]
                @test split(readline(fh)) == ["7", "0", "1", "2", "4", "1", "4", "0"]
                @test strip(readline(fh)) == ""
                @test eof(fh)
            finally
                close(fh)
                rm(fname)
            end
        end
    end
end
