using NESSie.Format
using JSON3
using LightXML: parse_string, root, name, child_elements, attribute, content

@testset "writexml3d_json (nodes)" begin
    for T in testtypes
        # empty system
        js = JSON3.read(readback(fh -> writexml3d_json(fh, Vector{T}[])))
        @test js["format"] == "xml3d-json"
        @test haskey(js, "version") == true
        @test js["data"]["position"]["type"] == "float3"
        @test length(js["data"]["position"]["seq"]) == 1
        @test js["data"]["position"]["seq"][1]["value"] == []
        # small system
        nodes = Vector{T}[T[0, 0, 0], T[0, 0, 3], T[0, 3, 0], T[1, -3, 3]]
        js = JSON3.read(readback(fh -> writexml3d_json(fh, nodes)))
        @test js["format"] == "xml3d-json"
        @test haskey(js, "version") == true
        @test js["data"]["position"]["type"] == "float3"
        @test length(js["data"]["position"]["seq"]) == 1
        @test js["data"]["position"]["seq"][1]["value"] == [0, 0, 0, 0, 0, 3, 0, 3, 0, 1, -3, 3]
    end
end

@testset "writexml3d_json (surface model)" begin
    for T in testtypes
        # empty system
        js = JSON3.read(readback(fh -> writexml3d_json(fh, Model{T, Triangle{T}}())))
        @test js["format"] == "xml3d-json"
        @test haskey(js, "version") == true
        @test js["data"]["index"]["type"] == "int"
        @test length(js["data"]["index"]["seq"]) == 1
        @test js["data"]["index"]["seq"][1]["value"] == []
        @test js["data"]["position"]["type"] == "float3"
        @test length(js["data"]["position"]["seq"]) == 1
        @test js["data"]["position"]["seq"][1]["value"] == []
        @test js["data"]["normal"]["type"] == "float3"
        @test length(js["data"]["normal"]["seq"]) == 1
        @test js["data"]["normal"]["seq"][1]["value"] == []
        # small system
        nodes    = Vector{T}[T[0, 0, 0], T[0, 0, 3], T[0, 3, 0], T[1, -3, 3]]
        elements = [Triangle(nodes[1], nodes[2], nodes[3]),
                    Triangle(nodes[1], nodes[4], nodes[2])]
        js = JSON3.read(readback(fh -> writexml3d_json(fh, Model(nodes, elements))))
        @test js["format"] == "xml3d-json"
        @test haskey(js, "version") == true
        @test js["data"]["index"]["type"] == "int"
        @test length(js["data"]["index"]["seq"]) == 1
        @test js["data"]["index"]["seq"][1]["value"] == [0, 1, 2, 0, 3, 1]
        @test js["data"]["position"]["type"] == "float3"
        @test length(js["data"]["position"]["seq"]) == 1
        @test js["data"]["position"]["seq"][1]["value"] == [0, 0, 0, 0, 0, 3, 0, 3, 0, 1, -3, 3]
        @test js["data"]["normal"]["type"] == "float3"
        @test length(js["data"]["normal"]["seq"]) == 1
        @test map(T, js["data"]["normal"]["seq"][1]["value"]) ≈ map(T,
            [√360 \ [-9 - √90, -3, 0]; √360 \ [-9 - √90, -3, 0]; [-1, 0, 0]; √90 \ [-9, -3, 0]]
        )
    end
end

@testset "writexml3d_xml" begin
    for T in testtypes
        # empty system
        xroot = root(parse_string(readback(fh -> writexml3d_xml(fh, Vector{T}[]))))
        @test name(xroot) == "xml3d"
        xchildren = collect(child_elements(xroot))
        @test length(xchildren) == 1
        @test name(xchildren[1]) == "data"
        @test attribute(xchildren[1], "id") == "mesh"
        xchildren = collect(child_elements(xchildren[1]))
        @test length(xchildren) == 1
        @test name(xchildren[1]) == "float3"
        @test strip(content(xchildren[1])) == ""
        # small system
        nodes    = Vector{T}[T[0, 0, 0], T[0, 0, 3], T[0, 3, 0], T[1, -3, 3]]
        xroot = root(parse_string(readback(fh -> writexml3d_xml(fh, nodes))))
        @test name(xroot) == "xml3d"
        xchildren = collect(child_elements(xroot))
        @test length(xchildren) == 1
        @test name(xchildren[1]) == "data"
        @test attribute(xchildren[1], "id") == "mesh"
        xchildren = collect(child_elements(xchildren[1]))
        @test length(xchildren) == 1
        @test name(xchildren[1]) == "float3"
        @test T[parse(T, e) for e in split(strip(content(xchildren[1])))] == T[0, 0, 0, 0, 0, 3, 0, 3, 0, 1, -3, 3]
    end
end
