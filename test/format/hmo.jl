using ProteinES.Format
using ProteinES.Format: readhmo_nodes, readhmo_elements, readhmo_charges

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

context("readhmo") do
    try
        for (fname, fh, len) in testfiles, fn in ("single", "bulk"), T in testtypes
            seekstart(fh)
            nodes = (); elements = (); charges = ()
            if fn == "single"
                # check single methods (atstart=false)
                nodes = readhmo_nodes(fh, T)
                elements = readhmo_elements(fh, nodes)
                charges = readhmo_charges(fh, T)
            elseif fn == "bulk"
                # check bulk method (implicit single methods with atstart=true)
                model = readhmo(fh, T)
                nodes, elements, charges = (model.nodes, model.elements, model.charges)
            end
            # check return types
            @fact typeof(nodes) --> Vector{Vector{T}}
            @fact typeof(elements) --> Vector{Triangle{T}}
            @fact typeof(charges) --> Vector{Charge{T}}
            # check lengths
            @fact length(nodes) --> len[1]
            @fact length(elements) --> len[2]
            @fact length(charges) --> len[3]

            # check content
            if fname == testfiles[2][1]
                @fact nodes[1] --> T[0, 0, 0]
                @fact nodes[2] --> T[1, 2, 3]
                @fact nodes[3] --> T[-2.00001, 1.337, 42]
                @fact elements[1].v1 --> exactly(nodes[1])
                @fact elements[1].v2 --> exactly(nodes[2])
                @fact elements[1].v3 --> exactly(nodes[3])
                @fact elements[2].v1 --> exactly(nodes[2])
                @fact elements[2].v2 --> exactly(nodes[3])
                @fact elements[2].v3 --> exactly(nodes[1])
                @fact charges[1].pos --> T[1, 2, 3]
                @fact charges[1].val --> T(999)
                @fact charges[2].pos --> T[-2.00001, 1.337, 42]
                @fact charges[2].val --> T(-999.9)
            end
        end
    finally
        [begin close(fh); rm(fn); end for (fn, fh, _) in testfiles]
    end
end
