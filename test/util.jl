import JSON: parse

# check eye!
for dtype in (Int, Float64, Float32)
    m = zeros(dtype, 3, 3)
    eye!(m)
    @test isa(m, Array{dtype, 2})
    @test size(m) == (3, 3)
    @test m == [1 0 0; 0 1 0; 0 0 1]
    eye!(m, 2)
    @test isa(m, Array{dtype, 2})
    @test size(m) == (3, 3)
    @test m == [2 0 0; 0 2 0; 0 0 2]
    m = zeros(dtype, 2, 3)
    eye!(m)
    @test isa(m, Array{dtype, 2})
    @test size(m) == (2, 3)
    @test m == [1 0 0; 0 1 0]
    eye!(m, 2)
    @test isa(m, Array{dtype, 2})
    @test size(m) == (2, 3)
    @test m == [2 0 0; 0 2 0]
    m = zeros(dtype, 3, 2)
    eye!(m)
    @test isa(m, Array{dtype, 2})
    @test size(m) == (3, 2)
    @test m == [1 0; 0 1; 0 0]
    eye!(m, 2)
    @test isa(m, Array{dtype, 2})
    @test size(m) == (3, 2)
    @test m == [2 0; 0 2; 0 0]
end

# check props! and isdegenerate
@test_throws AssertionError isdegenerate([1., 1.], [1., 1., 1.], [1., 1., 1.])
@test_throws AssertionError isdegenerate([1., 1., 1.], [1., 1.], [1., 1., 1.])
@test_throws AssertionError isdegenerate([1., 1., 1.], [1., 1., 1.], [1., 1.])
for dtype in (Float64, Float32)
    # degenerate triangles
    elem = Element(dtype[0, 0, 0], dtype[0, 0, 0], dtype[1, 0, 0])
    @test_throws AssertionError props!(elem)
    elem = Element(dtype[0, 0, 0], dtype[0, 0, 1], dtype[0, 0, 0])
    @test_throws AssertionError props!(elem)
    elem = Element(dtype[1, 0, 0], dtype[0, 0, 0], dtype[0, 0, 0])
    @test_throws AssertionError props!(elem)
    elem = Element(dtype[0, 1, 0], dtype[0, 2, 0], dtype[0, 3, 0])
    @test_throws AssertionError props!(elem)
    # simple 2D triangle
    elem = Element(dtype[0, 0, 0], dtype[0, 0, 3], dtype[0, 3, 0])
    props!(elem)
    @test isa(elem.center, Vector{dtype})
    @test length(elem.center) == 3
    @test isa(elem.normal, Vector{dtype})
    @test length(elem.normal) == 3
    @test isa(elem.distorig, dtype)
    @test isa(elem.area, dtype)
    @test_approx_eq elem.center [0, 1, 1]
    @test_approx_eq elem.normal [-1, 0, 0]
    @test_approx_eq elem.distorig 0.
    @test_approx_eq elem.area 4.5
    # simple 3D triangle
    elem = Element(dtype[3, 0, 0], dtype[0, 4, 0], dtype[0, 0, 5])
    props!(elem)
    @test isa(elem.center, Vector{dtype})
    @test length(elem.center) == 3
    @test isa(elem.normal, Vector{dtype})
    @test length(elem.normal) == 3
    @test isa(elem.distorig, dtype)
    @test isa(elem.area, dtype)
    @test_approx_eq elem.center [1., 4/3, 5/3]
    @test_approx_eq elem.normal (√769 \ [20, 15, 12])
    @test_approx_eq elem.distorig 60/√769
    @test_approx_eq elem.area √769/2
end

# check indexmap
for dtype in (Float64, Float32)
    v1 = dtype[1, 2, 3]
    v2 = dtype[4, 5, 6]
    v3 = dtype[7, 8, 9]
    d = indexmap(Vector{dtype}[])
    @test isa(d, Dict{Ptr{dtype}, Int})
    @test d == Dict()
    d = indexmap(Vector{dtype}[v1, v2, v3])
    @test isa(d, Dict{Ptr{dtype}, Int})
    @test length(d) == 3
    @test d[pointer(v1)] == 1 && d[pointer(v2)] == 2 && d[pointer(v3)] == 3
    d = indexmap(Vector{dtype}[v1, v1, v2])
    @test isa(d, Dict{Ptr{dtype}, Int})
    @test length(d) == 2
    @test d[pointer(v1)] == 2 && d[pointer(v2)] == 3
    @test_throws KeyError d[pointer(v3)]
end

# check unpack
for dtype in (Float64, Float32)
    d = unpack(Vector{dtype}[dtype[1], dtype[2], dtype[3]], 0)
    @test isa(d, Vector{dtype})
    @test d == []
    d = unpack(Vector{dtype}[dtype[1], dtype[2], dtype[3]], 1)
    @test isa(d, Vector{dtype})
    @test d == [1, 2, 3]
    d = unpack(Vector{dtype}[dtype[1, 2], dtype[3, 4]], 1)
    @test isa(d, Vector{dtype})
    @test d == [1, 3]
    d = unpack(Vector{dtype}[dtype[1, 2], dtype[3, 4]], 2)
    @test isa(d, Vector{dtype})
    @test d == [1, 2, 3, 4]
    d = unpack(Vector{dtype}[dtype[1, 2, 3], dtype[4, 5, 6]])
    @test isa(d, Vector{dtype})
    @test d == [1, 2, 3, 4, 5, 6]
    d = unpack(Vector{dtype}[dtype[1, 2, 3, 4, 5, 6]], 6)
    @test isa(d, Vector{dtype})
    @test d == [1, 2, 3, 4, 5, 6]
    @test_throws BoundsError unpack(Vector{dtype}[dtype[1], dtype[2], dtype[3]], 2)
end

# check vertexnormals
for dtype in (Float64, Float32)
    nodes    = Vector{dtype}[dtype[0, 0, 0], dtype[0, 0, 3], dtype[0, 3, 0], dtype[1, -3, 3]]
    elements = [Element(nodes[1], nodes[2], nodes[3]),
                Element(nodes[1], nodes[4], nodes[2])]
    map(props!, elements)
    d = vertexnormals(Vector{dtype}[], Element{dtype}[])
    @test isa(d, Vector{Vector{dtype}})
    @test d == []
    d = vertexnormals(Vector{dtype}[nodes[1], nodes[2], nodes[3]], Element{dtype}[elements[1]])
    @test isa(d, Vector{Vector{dtype}})
    @test length(d) == 3
    @test_approx_eq d[1] [-1, 0, 0]
    @test_approx_eq d[2] [-1, 0, 0]
    @test_approx_eq d[3] [-1, 0, 0]
    d = vertexnormals(Vector{dtype}[nodes[1], nodes[2], nodes[4]], Element{dtype}[elements[2]])
    @test isa(d, Vector{Vector{dtype}})
    @test length(d) == 3
    @test_approx_eq d[1] √90 \ [-9, -3, 0]
    @test_approx_eq d[2] √90 \ [-9, -3, 0]
    @test_approx_eq d[3] √90 \ [-9, -3, 0]
    d = vertexnormals(nodes, elements)
    @test isa(d, Vector{Vector{dtype}})
    @test length(d) == 4
    @test_approx_eq d[1] √360 \ [-9 - √90, -3, 0]
    @test_approx_eq d[2] √360 \ [-9 - √90, -3, 0]
    @test_approx_eq d[3] [-1, 0, 0]
    @test_approx_eq d[4] √90 \ [-9, -3, 0]
    d = vertexnormals(nodes, elements, true)
    @test isa(d, Vector{Vector{dtype}})
    @test length(d) == 4
    @test_approx_eq d[1] -√360 \ [-9 - √90, -3, 0]
    @test_approx_eq d[2] -√360 \ [-9 - √90, -3, 0]
    @test_approx_eq d[3] [1, 0, 0]
    @test_approx_eq d[4] -√90 \ [-9, -3, 0]
end

# check xml3d_mesh
for dtype in (Float64, Float32)
    # empty system
    js = parse(xml3d_mesh(Vector{dtype}[], Element{dtype}[]))
    @test js["format"] == "xml3d-json"
    @test haskey(js, "version")
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
    nodes    = Vector{dtype}[dtype[0, 0, 0], dtype[0, 0, 3], dtype[0, 3, 0], dtype[1, -3, 3]]
    elements = [Element(nodes[1], nodes[2], nodes[3]),
                Element(nodes[1], nodes[4], nodes[2])]
    map(props!, elements)
    js = parse(xml3d_mesh(nodes, elements))
    @test js["format"] == "xml3d-json"
    @test haskey(js, "version")
    @test js["data"]["index"]["type"] == "int"
    @test length(js["data"]["index"]["seq"]) == 1
    @test js["data"]["index"]["seq"][1]["value"] == [0, 1, 2, 0, 3, 1]
    @test js["data"]["position"]["type"] == "float3"
    @test length(js["data"]["position"]["seq"]) == 1
    @test js["data"]["position"]["seq"][1]["value"] == [0, 0, 0, 0, 0, 3, 0, 3, 0, 1, -3, 3]
    @test js["data"]["normal"]["type"] == "float3"
    @test length(js["data"]["normal"]["seq"]) == 1
    @test_approx_eq(
        convert(Vector{dtype}, js["data"]["normal"]["seq"][1]["value"]),
        [√360 \ [-9 - √90, -3, 0]; √360 \ [-9 - √90, -3, 0]; [-1, 0, 0]; √90 \ [-9, -3, 0]]
    )
    # small system, inverted normals
    js = parse(xml3d_mesh(nodes, elements, true))
    @test js["format"] == "xml3d-json"
    @test haskey(js, "version")
    @test js["data"]["index"]["type"] == "int"
    @test length(js["data"]["index"]["seq"]) == 1
    @test js["data"]["index"]["seq"][1]["value"] == [0, 1, 2, 0, 3, 1]
    @test js["data"]["position"]["type"] == "float3"
    @test length(js["data"]["position"]["seq"]) == 1
    @test js["data"]["position"]["seq"][1]["value"] == [0, 0, 0, 0, 0, 3, 0, 3, 0, 1, -3, 3]
    @test js["data"]["normal"]["type"] == "float3"
    @test length(js["data"]["normal"]["seq"]) == 1
    @test_approx_eq(
        convert(Vector{dtype}, js["data"]["normal"]["seq"][1]["value"]),
        [-√360 \ [-9 - √90, -3, 0]; -√360 \ [-9 - √90, -3, 0]; [1, 0, 0]; -√90 \ [-9, -3, 0]]
    )
end
