@testitem "Format: HyperMesh/.hmo" begin
    include("../testsetup.jl")

    testfiles = ((mktemp()..., (0,0,0)), # empty file
                (mktemp()..., (3,2,2))) # dummy file
    write(testfiles[2][2], """
    # stuffing
    stuffing

    BEG_NODL_DATA
        stuffing
        x 0.0 0.0 0.0
        x 1.0 2.0 3.0
        x -2.00001 1.337 42.0
    END_NODL_DATA
    stuffing

    BEG_ELEM_DATA
        stuffing
        x x x 1 2 3
        x x x 2 3 1
    END_ELEM_DATA
    stuffing

    BEG_CHARGE_DATA
        stuffing
        x 1.0 2.0 3.0 999.0
        x -2.00001 1.337 42.0 -999.9
    END_CHARGE_DATA

    """)

    @testset "readhmo" begin
        try
            for (fname, fh, len) in testfiles, fn in ("single", "bulk"), T in testtypes
                seekstart(fh)
                nodes = (); elements = (); charges = ()
                if fn == "single"
                    # check single methods (atstart=false)
                    nodes = Format.readhmo_nodes(fh, T)
                    elements = Format.readhmo_elements(fh, nodes)
                    charges = Format.readhmo_charges(fh, T)
                elseif fn == "bulk"
                    # check bulk method (implicit single methods with atstart=true)
                    model = Format.readhmo(fh, T)
                    nodes, elements, charges = (model.nodes, model.elements, model.charges)
                end
                # check return types
                @test typeof(nodes) == Vector{Vector{T}}
                @test typeof(elements) == Vector{Triangle{T}}
                @test typeof(charges) == Vector{Charge{T}}
                # check lengths
                @test length(nodes) == len[1]
                @test length(elements) == len[2]
                @test length(charges) == len[3]

                # check content
                if fname == testfiles[2][1]
                    @test nodes[1] == T[0, 0, 0]
                    @test nodes[2] == T[1, 2, 3]
                    @test nodes[3] == T[-2.00001, 1.337, 42]
                    @test elements[1].v1 === nodes[1]
                    @test elements[1].v2 === nodes[2]
                    @test elements[1].v3 === nodes[3]
                    @test elements[2].v1 === nodes[2]
                    @test elements[2].v2 === nodes[3]
                    @test elements[2].v3 === nodes[1]
                    @test charges[1].pos == T[1, 2, 3]
                    @test charges[1].val == T(999)
                    @test charges[2].pos == T[-2.00001, 1.337, 42]
                    @test charges[2].val == T(-999.9)
                end
            end
        finally
            [begin close(fh); rm(fn); end for (fn, fh, _) in testfiles]
        end
    end

    @testset "writehmo" begin
        @test_skip writehmo
    end
end
