const tempfile = mktemp()
const testtypes = (Float64, Float32)
const testfiles = ((tempfile[1], tempfile[2], (0,0,0)),   # empty file
                   ("hmo.hmo", open("hmo.hmo"), (2,2,2))) # dummy file

try
    for (fname, fh, len) in testfiles, fn in ("single", "bulk"), dtype in testtypes
        nodes = (); elements = (); charges = ()
        if fn == "single"
            # check single methods (atstart=false)
            nodes = readhmo_nodes(fh, false, dtype=dtype)
            elements = readhmo_elements(fh, nodes, false, dtype=dtype)
            charges = readhmo_charges(fh, false, dtype=dtype)
        elseif fn == "bulk"
            # check bulk method (implicit single methods with atstart=true)
            nodes, elements, charges = readhmo(fh, dtype=dtype)
        end
        # check return types
        @test isa(nodes, Vector{Vector{dtype}})
        @test isa(elements, Vector{Element{dtype}})
        @test isa(charges, Vector{Charge{dtype}})
        # check lengths
        @test length(nodes) == len[1]
        @test length(elements) == len[2]
        @test length(charges) == len[3]
        # check content
        if fname == "hmo.hmo"
            @test nodes[1] == dtype[1, 2, 3]
            @test nodes[2] == dtype[-2.00001, 1.337, 42]
            @test elements[1].v1 === elements[1].v3 === elements[2].v2 === nodes[1]
            @test elements[1].v2 === elements[2].v1 === elements[2].v3 === nodes[2]
            @test charges[1].pos == dtype[1, 2, 3] && charges[1].val == convert(dtype, 999)
            @test charges[2].pos == dtype[-2.00001, 1.337, 42] && charges[2].val == convert(dtype, -999.9)
        end
        seek(fh,0)
    end
finally
    [close(fh) for (_, fh, __) in testfiles]
    rm(tempfile[1])
end

