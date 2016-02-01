using ProteinES.IO
import ProteinES.IO: readoff_nodes, readoff_elements

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

context("readoff") do
    try
        for (fname, fh, len) in testfiles, fn in ("single", "bulk"), T in testtypes
            seekstart(fh)
            nodes = (); elements = ()
            if fn == "single"
                # check single methods
                readline(fh)
                readline(fh)
                nodes = readoff_nodes(fh, len[1], T)
                elements = readoff_elements(fh, len[2], nodes, T)
            elseif fn == "bulk"
                # check bulk method
                nodes, elements = readoff(fh, T)
            end
            # check return types
            @fact isa(nodes, Vector{Vector{T}}) --> true
            @fact isa(elements, Vector{Triangle{T}}) --> true

            # check lengths
            @fact length(nodes) --> len[1]
            @fact length(elements) --> len[2]

            # check content
            if fname == testfiles[2][1]
                @fact nodes[1] --> T[1, 0, 0]
                @fact nodes[2] --> T[0, 1, 0]
                @fact nodes[3] --> T[0, 0, 1]
                @fact elements[1].v1 --> exactly(nodes[1])
                @fact elements[1].v2 --> exactly(nodes[2])
                @fact elements[1].v3 --> exactly(nodes[3])
                @fact elements[2].v1 --> exactly(nodes[3])
                @fact elements[2].v2 --> exactly(nodes[1])
                @fact elements[2].v3 --> exactly(nodes[2])
            end
        end
    finally
        [begin close(fh); rm(fn); end for (fn, fh, _) in testfiles]
    end
end
