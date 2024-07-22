@testitem "Format: GAMer/.off" begin
    include("../testsetup.jl")

    testfiles = ((mktemp()..., (0,0)), # empty file
                (mktemp()..., (3,2))) # dummy file
    write(testfiles[2][2], """
    OFF
    3 2 5
    1.0000000000e+00     0.0000000000e+00     0.0000000000e+00
    0.0000000000e+00     1.0000000000e+00     0.0000000000e+00
    0.0000000000e+00     0.0000000000e+00     1.0000000000e+00
    3 0 1 2
    3 2 0 1
    """)

    @testset "readoff" begin
        try
            for (fname, fh, len) in testfiles, fn in ("single", "bulk"), T in testtypes
                seekstart(fh)
                nodes = (); elements = ()
                if fn == "single"
                    # check single methods
                    readline(fh)
                    readline(fh)
                    nodes = Format.readoff_nodes(fh, len[1], T)
                    elements = Format.readoff_elements(fh, len[2], nodes)
                elseif fn == "bulk"
                    # check bulk method
                    model = Format.readoff(fh, T)
                    nodes, elements = (model.nodes, model.elements)
                end
                # check return types
                @test typeof(nodes) == Vector{Vector{T}}
                @test typeof(elements) == Vector{Triangle{T}}

                # check lengths
                @test length(nodes) == len[1]
                @test length(elements) == len[2]

                # check content
                if fname == testfiles[2][1]
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
        finally
            [begin close(fh); rm(fn); end for (fn, fh, _) in testfiles]
        end
    end
end
