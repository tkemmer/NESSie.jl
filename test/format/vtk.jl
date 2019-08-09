using NESSie.Format
using LightXML: parse_string, root, name, child_elements, attribute, content

@testset "writevtk (surface model)" begin
    for T in testtypes
        # empty system
        model = Model{T, Triangle{T}}()
        (fname, fh) = mktemp()
        try
            xroot = root(parse_string(readback(fh -> writevtk(fh, model))))
            @test name(xroot) == "VTKFile"
            @test attribute(xroot, "type") == "PolyData"
            xgrid = collect(child_elements(xroot))
            @test length(xgrid) == 1
            @test name(xgrid[1]) == "PolyData"
            xgrid = collect(child_elements(xgrid[1]))
            @test length(xgrid) == 1
            @test name(xgrid[1]) == "Piece"
            @test attribute(xgrid[1], "NumberOfPoints") == "0"
            @test attribute(xgrid[1], "NumberOfVerts") == "0"
            @test attribute(xgrid[1], "NumberOfLines") == "0"
            @test attribute(xgrid[1], "NumberOfStrips") == "0"
            @test attribute(xgrid[1], "NumberOfPolys") == "0"
            xpiece = collect(child_elements(xgrid[1]))
            @test length(xpiece) == 2
            @test name(xpiece[1]) == "Points"
            xdata = collect(child_elements(xpiece[1]))
            @test length(xdata) == 1
            @test name(xdata[1]) == "DataArray"
            @test attribute(xdata[1], "type") == "$T"
            @test attribute(xdata[1], "NumberOfComponents") == "3"
            @test attribute(xdata[1], "format") == "ascii"
            @test strip(content(xdata[1])) == ""
            @test name(xpiece[2]) == "Polys"
            xdata = collect(child_elements(xpiece[2]))
            @test length(xdata) == 2
            @test name(xdata[1]) == "DataArray"
            @test attribute(xdata[1], "type") == "Int32"
            @test attribute(xdata[1], "Name") == "connectivity"
            @test attribute(xdata[1], "format") == "ascii"
            @test strip(content(xdata[1])) == ""
            @test name(xdata[2]) == "DataArray"
            @test attribute(xdata[2], "type") == "Int32"
            @test attribute(xdata[2], "Name") == "offsets"
            @test attribute(xdata[2], "format") == "ascii"
            @test strip(content(xdata[2])) == ""
        finally
            close(fh)
            rm(fname)
        end

        # small system
        nodes = Vector{T}[[0, 0, 0], [0, 0, 3], [0, 3, 0], [1, -3, 3]]
        model = Model(
            nodes,
            Triangle{T}[
                Triangle(nodes[1], nodes[2], nodes[3]),
                Triangle(nodes[1], nodes[4], nodes[2])
            ]
        )
        (fname, fh) = mktemp()
        try
            xroot = root(parse_string(readback(fh -> writevtk(fh, model))))
            @test name(xroot) == "VTKFile"
            @test attribute(xroot, "type") == "PolyData"
            xgrid = collect(child_elements(xroot))
            @test length(xgrid) == 1
            @test name(xgrid[1]) == "PolyData"
            xgrid = collect(child_elements(xgrid[1]))
            @test length(xgrid) == 1
            @test name(xgrid[1]) == "Piece"
            @test attribute(xgrid[1], "NumberOfPoints") == "4"
            @test attribute(xgrid[1], "NumberOfVerts") == "0"
            @test attribute(xgrid[1], "NumberOfLines") == "0"
            @test attribute(xgrid[1], "NumberOfStrips") == "0"
            @test attribute(xgrid[1], "NumberOfPolys") == "2"
            xpiece = collect(child_elements(xgrid[1]))
            @test length(xpiece) == 2
            @test name(xpiece[1]) == "Points"
            xdata = collect(child_elements(xpiece[1]))
            @test length(xdata) == 1
            @test name(xdata[1]) == "DataArray"
            @test attribute(xdata[1], "type") == "$T"
            @test attribute(xdata[1], "NumberOfComponents") == "3"
            @test attribute(xdata[1], "format") == "ascii"
            @test [parse(T, e) for e in split(strip(content(xdata[1])))] == T[0, 0, 0, 0, 0, 3, 0, 3, 0, 1, -3, 3]
            @test name(xpiece[2]) == "Polys"
            xdata = collect(child_elements(xpiece[2]))
            @test length(xdata) == 2
            @test name(xdata[1]) == "DataArray"
            @test attribute(xdata[1], "type") == "Int32"
            @test attribute(xdata[1], "Name") == "connectivity"
            @test attribute(xdata[1], "format") == "ascii"
            @test [parse(Int, e) for e in split(strip(content(xdata[1])))] == [0, 1, 2, 0, 3, 1]
            @test name(xdata[2]) == "DataArray"
            @test attribute(xdata[2], "type") == "Int32"
            @test attribute(xdata[2], "Name") == "offsets"
            @test attribute(xdata[2], "format") == "ascii"
            @test [parse(Int, e) for e in split(strip(content(xdata[2])))] == [3, 6]
        finally
            close(fh)
            rm(fname)
        end
    end
end

@testset "writevtk (volume model)" begin
    for T in testtypes
        # empty system
        model = Model{T, Tetrahedron{T}}()
        (fname, fh) = mktemp()
        try
            xroot = root(parse_string(readback(fh -> writevtk(fh, model))))
            @test name(xroot) == "VTKFile"
            @test attribute(xroot, "type") == "UnstructuredGrid"
            xgrid = collect(child_elements(xroot))
            @test length(xgrid) == 1
            @test name(xgrid[1]) == "UnstructuredGrid"
            xgrid = collect(child_elements(xgrid[1]))
            @test length(xgrid) == 1
            @test name(xgrid[1]) == "Piece"
            @test attribute(xgrid[1], "NumberOfPoints") == "0"
            @test attribute(xgrid[1], "NumberOfCells") == "0"
            xpiece = collect(child_elements(xgrid[1]))
            @test length(xpiece) == 2
            @test name(xpiece[1]) == "Points"
            xdata = collect(child_elements(xpiece[1]))
            @test length(xdata) == 1
            @test name(xdata[1]) == "DataArray"
            @test attribute(xdata[1], "type") == "$T"
            @test attribute(xdata[1], "NumberOfComponents") == "3"
            @test attribute(xdata[1], "format") == "ascii"
            @test strip(content(xdata[1])) == ""
            @test name(xpiece[2]) == "Cells"
            xdata = collect(child_elements(xpiece[2]))
            @test length(xdata) == 3
            @test name(xdata[1]) == "DataArray"
            @test attribute(xdata[1], "type") == "Int32"
            @test attribute(xdata[1], "Name") == "connectivity"
            @test attribute(xdata[1], "format") == "ascii"
            @test strip(content(xdata[1])) == ""
            @test name(xdata[2]) == "DataArray"
            @test attribute(xdata[2], "type") == "Int32"
            @test attribute(xdata[2], "Name") == "offsets"
            @test attribute(xdata[2], "format") == "ascii"
            @test strip(content(xdata[2])) == ""
            @test name(xdata[3]) == "DataArray"
            @test attribute(xdata[3], "type") == "UInt8"
            @test attribute(xdata[3], "Name") == "types"
            @test attribute(xdata[3], "format") == "ascii"
            @test strip(content(xdata[3])) == ""
        finally
            close(fh)
            rm(fname)
        end

        # small system
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
            xroot = root(parse_string(readback(fh -> writevtk(fh, model))))
            @test name(xroot) == "VTKFile"
            @test attribute(xroot, "type") == "UnstructuredGrid"
            xgrid = collect(child_elements(xroot))
            @test length(xgrid) == 1
            @test name(xgrid[1]) == "UnstructuredGrid"
            xgrid = collect(child_elements(xgrid[1]))
            @test length(xgrid) == 1
            @test name(xgrid[1]) == "Piece"
            @test attribute(xgrid[1], "NumberOfPoints") == "5"
            @test attribute(xgrid[1], "NumberOfCells") == "2"
            xpiece = collect(child_elements(xgrid[1]))
            @test length(xpiece) == 2
            @test name(xpiece[1]) == "Points"
            xdata = collect(child_elements(xpiece[1]))
            @test length(xdata) == 1
            @test name(xdata[1]) == "DataArray"
            @test attribute(xdata[1], "type") == "$T"
            @test attribute(xdata[1], "NumberOfComponents") == "3"
            @test attribute(xdata[1], "format") == "ascii"
            @test [parse(T, e) for e in split(strip(content(xdata[1])))] == T[0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, -2, 0, 0]
            @test name(xpiece[2]) == "Cells"
            xdata = collect(child_elements(xpiece[2]))
            @test length(xdata) == 3
            @test name(xdata[1]) == "DataArray"
            @test attribute(xdata[1], "type") == "Int32"
            @test attribute(xdata[1], "Name") == "connectivity"
            @test attribute(xdata[1], "format") == "ascii"
            @test [parse(Int, e) for e in split(strip(content(xdata[1])))] == [0, 1, 2, 3, 0, 1, 2, 4]
            @test name(xdata[2]) == "DataArray"
            @test attribute(xdata[2], "type") == "Int32"
            @test attribute(xdata[2], "Name") == "offsets"
            @test attribute(xdata[2], "format") == "ascii"
            @test [parse(Int, e) for e in split(strip(content(xdata[2])))] == [4, 8]
            @test name(xdata[3]) == "DataArray"
            @test attribute(xdata[3], "type") == "UInt8"
            @test attribute(xdata[3], "Name") == "types"
            @test attribute(xdata[3], "format") == "ascii"
            @test [parse(Int, e) for e in split(strip(content(xdata[3])))] == [10, 10]
        finally
            close(fh)
            rm(fname)
        end
    end
end
