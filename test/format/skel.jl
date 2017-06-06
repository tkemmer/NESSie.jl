using NESSie.Format

context("writeskel (surface model)") do
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
            writeskel(fh, model)
            seek(fh, 0)
            @fact readline(fh) --> "SKEL"
            @fact split(readline(fh)) --> ["4", "2"]
            for node in nodes
                @fact [parse(T, e) for e in split(readline(fh))] --> node
            end
            @fact split(readline(fh)) --> ["4", "0", "1", "2", "0"]
            @fact split(readline(fh)) --> ["4", "0", "1", "3", "0"]
            @fact strip(readline(fh)) --> ""
            @fact eof(fh) --> true
        finally
            close(fh)
            rm(fname)
        end
    end
end

context("writeskel (volume model)") do
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
            writeskel(fh, model)
            seek(fh, 0)
            @fact readline(fh) --> "SKEL"
            @fact split(readline(fh)) --> ["5", "2"]
            for node in nodes
                @fact [parse(T, e) for e in split(readline(fh))] --> node
            end
            @fact split(readline(fh)) --> ["7", "0", "1", "2", "3", "1", "3", "0"]
            @fact split(readline(fh)) --> ["7", "0", "1", "2", "4", "1", "4", "0"]
            @fact strip(readline(fh)) --> ""
            @fact eof(fh) --> true
        finally
            close(fh)
            rm(fname)
        end
    end
end
