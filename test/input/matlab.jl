import NonLocalES: readmatlab_nodes, readmatlab_elements

testfiles = ((mktemp()..., (0,0)), # empty file
             (mktemp()..., (4,2))) # dummy file
write(testfiles[2][2], """
mcsf_begin=1;

      dim = 3;
    dimii = 3;
 vertices = 4;
simplices = 2;
vert=[
%------------------------------------------------------------------------------
%  Node-ID  Chrt        X-Coordinate        Y-coordinate        Z-coordinate
%---------  ----     ----------------     ----------------     ----------------
        0    2      0.0000000000e+00     0.0000000000e+00     0.0000000000e+00
        1    2      1.0000000000e+00     0.0000000000e+00     0.0000000000e+00
        2    2      0.0000000000e+00     1.0000000000e+00     0.0000000000e+00
        3    2      0.0000000000e+00     0.0000000000e+00     1.0000000000e+00
];
simp=[
%------------------------------------------------------------------------------------------
%  Simp-ID Grp    Mat          Face-Types                      Vertex-Numbers
%--------- ---    ---    ---------------------  -------------------------------------------
        0  0      0      0     0     0     0           0          1          2          3
        1  0      1      0     0     0     0           3          2          1          0
];
mcsf_end=1;""")

try
    for (fname, fh, len) in testfiles, fn in ("single", "bulk"), T in testtypes
        seekstart(fh)
        nodes = (); elements = ()
        if fn == "single"
            # check single methods
            nodes = readmatlab_nodes(fh, T)
            elements = readmatlab_elements(fh, nodes, T)
        elseif fn == "bulk"
            # check bulk method
            nodes, elements = readmatlab(fh, T)
        end
        # check return types
        @fact isa(nodes, Vector{Vector{T}}) --> true
        @fact isa(elements, Vector{Tetrahedron{T}}) --> true
        # check lengths
        @fact length(nodes) --> len[1]
        @fact length(elements) --> len[2]

        # check content
        if fname == testfiles[2][1]
            @fact nodes[1] --> T[0, 0, 0]
            @fact nodes[2] --> T[1, 0, 0]
            @fact nodes[3] --> T[0, 1, 0]
            @fact nodes[4] --> T[0, 0, 1]
            @fact elements[1].v1 --> exactly(nodes[1])
            @fact elements[1].v2 --> exactly(nodes[2])
            @fact elements[1].v3 --> exactly(nodes[3])
            @fact elements[1].v4 --> exactly(nodes[4])
            @fact elements[2].v1 --> exactly(nodes[4])
            @fact elements[2].v2 --> exactly(nodes[3])
            @fact elements[2].v3 --> exactly(nodes[2])
            @fact elements[2].v4 --> exactly(nodes[1])
        end
    end
finally
    [begin close(fh); rm(fn); end for (fn, fh, _) in testfiles]
end