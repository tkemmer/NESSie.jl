using ProteinES.Format
using ProteinES.Format: readmsms_nodes, readmsms_elements

testfiles = ((mktemp()..., mktemp()..., (0,0)), # empty file
             (mktemp()..., mktemp()..., (3,2))) # dummy file

write(testfiles[2][2], """
# MSMS solvent excluded surface vertices
#vertex #sphere density probe_r
  x  x  x  x
  0.0 0.0 0.0 x x x x x x
  1.0 2.0 3.0 x x x x x x
  -2.00001 1.337 42.0 x x x x x x
""")

write(testfiles[2][4], """
# MSMS solvent excluded surface faces
#faces  #sphere density probe_r
  x  x  x  x
  1 2 3 x x
  2 3 1 x x
""")

context("readmsms") do
    try
        for (fnamev, fv, fnamef, ff, len) in testfiles, fn in ("single", "bulk"), T in testtypes
            seekstart(fv)
            seekstart(ff)
            nodes = (); elements = ()
            if fn == "single"
                # check single methods
                nodes = readmsms_nodes(fv, T)
                elements = readmsms_elements(ff, nodes)
            elseif fn == "bulk"
                # check bulk method (implicit single methods with atstart=true)
                model = readmsms(fv, ff, T)
                nodes, elements = (model.nodes, model.elements)
            end
            # check return types
            @fact typeof(nodes) --> Vector{Vector{T}}
            @fact typeof(elements) --> Vector{Triangle{T}}
            # check lengths
            @fact length(nodes) --> len[1]
            @fact length(elements) --> len[2]

            # check content
            if fnamev == testfiles[2][1]
                @fact nodes[1] --> T[0, 0, 0]
                @fact nodes[2] --> T[1, 2, 3]
                @fact nodes[3] --> T[-2.00001, 1.337, 42]
                @fact elements[1].v1 --> exactly(nodes[1])
                @fact elements[1].v2 --> exactly(nodes[2])
                @fact elements[1].v3 --> exactly(nodes[3])
                @fact elements[2].v1 --> exactly(nodes[2])
                @fact elements[2].v2 --> exactly(nodes[3])
                @fact elements[2].v3 --> exactly(nodes[1])
            end
        end
    finally
        [begin close(fv); rm(fnamev); close(ff); rm(fnamef) end for (fnamev, fv, fnamef, ff, _) in testfiles]
    end
end
