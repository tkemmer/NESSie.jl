import JSON: parse

# check eye!
for T in (Int, Float64, Float32)
    m = zeros(T, 3, 3)
    eye!(m)
    @test isa(m, Array{T, 2})
    @test size(m) == (3, 3)
    @test m == [1 0 0; 0 1 0; 0 0 1]
    eye!(m, 2)
    @test isa(m, Array{T, 2})
    @test size(m) == (3, 3)
    @test m == [2 0 0; 0 2 0; 0 0 2]
    m = zeros(T, 2, 3)
    eye!(m)
    @test isa(m, Array{T, 2})
    @test size(m) == (2, 3)
    @test m == [1 0 0; 0 1 0]
    eye!(m, 2)
    @test isa(m, Array{T, 2})
    @test size(m) == (2, 3)
    @test m == [2 0 0; 0 2 0]
    m = zeros(T, 3, 2)
    eye!(m)
    @test isa(m, Array{T, 2})
    @test size(m) == (3, 2)
    @test m == [1 0; 0 1; 0 0]
    eye!(m, 2)
    @test isa(m, Array{T, 2})
    @test size(m) == (3, 2)
    @test m == [2 0; 0 2; 0 0]
end

# check props! and isdegenerate
@test_throws AssertionError isdegenerate([1., 1.], [1., 1., 1.], [1., 1., 1.])
@test_throws AssertionError isdegenerate([1., 1., 1.], [1., 1.], [1., 1., 1.])
@test_throws AssertionError isdegenerate([1., 1., 1.], [1., 1., 1.], [1., 1.])
for T in (Float64, Float32)
    # degenerate triangles
    elem = Triangle(T[0, 0, 0], T[0, 0, 0], T[1, 0, 0])
    @test_throws AssertionError props!(elem)
    elem = Triangle(T[0, 0, 0], T[0, 0, 1], T[0, 0, 0])
    @test_throws AssertionError props!(elem)
    elem = Triangle(T[1, 0, 0], T[0, 0, 0], T[0, 0, 0])
    @test_throws AssertionError props!(elem)
    elem = Triangle(T[0, 1, 0], T[0, 2, 0], T[0, 3, 0])
    @test_throws AssertionError props!(elem)
    # simple 2D triangle
    elem = Triangle(T[0, 0, 0], T[0, 0, 3], T[0, 3, 0])
    props!(elem)
    @test isa(elem.center, Vector{T})
    @test length(elem.center) == 3
    @test isa(elem.normal, Vector{T})
    @test length(elem.normal) == 3
    @test isa(elem.distorig, T)
    @test isa(elem.area, T)
    @test_approx_eq elem.center [0, 1, 1]
    @test_approx_eq elem.normal [-1, 0, 0]
    @test_approx_eq elem.distorig 0.
    @test_approx_eq elem.area 4.5
    # simple 3D triangle
    elem = Triangle(T[3, 0, 0], T[0, 4, 0], T[0, 0, 5])
    props!(elem)
    @test isa(elem.center, Vector{T})
    @test length(elem.center) == 3
    @test isa(elem.normal, Vector{T})
    @test length(elem.normal) == 3
    @test isa(elem.distorig, T)
    @test isa(elem.area, T)
    @test_approx_eq elem.center [1., 4/3, 5/3]
    @test_approx_eq elem.normal (√769 \ [20, 15, 12])
    @test_approx_eq elem.distorig 60/√769
    @test_approx_eq elem.area √769/2
end

# check indexmap
for T in (Float64, Float32)
    v1 = T[1, 2, 3]
    v2 = T[4, 5, 6]
    v3 = T[7, 8, 9]
    d = indexmap(Vector{T}[])
    @test isa(d, Dict{Ptr{T}, Int})
    @test d == Dict()
    d = indexmap(Vector{T}[v1, v2, v3])
    @test isa(d, Dict{Ptr{T}, Int})
    @test length(d) == 3
    @test d[pointer(v1)] == 1 && d[pointer(v2)] == 2 && d[pointer(v3)] == 3
    d = indexmap(Vector{T}[v1, v1, v2])
    @test isa(d, Dict{Ptr{T}, Int})
    @test length(d) == 2
    @test d[pointer(v1)] == 2 && d[pointer(v2)] == 3
    @test_throws KeyError d[pointer(v3)]
end

# check unpack
for T in (Float64, Float32)
    d = unpack(Vector{T}[T[1], T[2], T[3]], 0)
    @test isa(d, Vector{T})
    @test d == []
    d = unpack(Vector{T}[T[1], T[2], T[3]], 1)
    @test isa(d, Vector{T})
    @test d == [1, 2, 3]
    d = unpack(Vector{T}[T[1, 2], T[3, 4]], 1)
    @test isa(d, Vector{T})
    @test d == [1, 3]
    d = unpack(Vector{T}[T[1, 2], T[3, 4]], 2)
    @test isa(d, Vector{T})
    @test d == [1, 2, 3, 4]
    d = unpack(Vector{T}[T[1, 2, 3], T[4, 5, 6]])
    @test isa(d, Vector{T})
    @test d == [1, 2, 3, 4, 5, 6]
    d = unpack(Vector{T}[T[1, 2, 3, 4, 5, 6]], 6)
    @test isa(d, Vector{T})
    @test d == [1, 2, 3, 4, 5, 6]
    @test_throws BoundsError unpack(Vector{T}[T[1], T[2], T[3]], 2)
end

# check vertexnormals
for T in (Float64, Float32)
    nodes    = Vector{T}[T[0, 0, 0], T[0, 0, 3], T[0, 3, 0], T[1, -3, 3]]
    elements = [Triangle(nodes[1], nodes[2], nodes[3]),
                Triangle(nodes[1], nodes[4], nodes[2])]
    map(props!, elements)
    d = vertexnormals(Vector{T}[], Triangle{T}[])
    @test isa(d, Vector{Vector{T}})
    @test d == []
    d = vertexnormals(Vector{T}[nodes[1], nodes[2], nodes[3]], Triangle{T}[elements[1]])
    @test isa(d, Vector{Vector{T}})
    @test length(d) == 3
    @test_approx_eq d[1] [-1, 0, 0]
    @test_approx_eq d[2] [-1, 0, 0]
    @test_approx_eq d[3] [-1, 0, 0]
    d = vertexnormals(Vector{T}[nodes[1], nodes[2], nodes[4]], Triangle{T}[elements[2]])
    @test isa(d, Vector{Vector{T}})
    @test length(d) == 3
    @test_approx_eq d[1] √90 \ [-9, -3, 0]
    @test_approx_eq d[2] √90 \ [-9, -3, 0]
    @test_approx_eq d[3] √90 \ [-9, -3, 0]
    d = vertexnormals(nodes, elements)
    @test isa(d, Vector{Vector{T}})
    @test length(d) == 4
    @test_approx_eq d[1] √360 \ [-9 - √90, -3, 0]
    @test_approx_eq d[2] √360 \ [-9 - √90, -3, 0]
    @test_approx_eq d[3] [-1, 0, 0]
    @test_approx_eq d[4] √90 \ [-9, -3, 0]
    d = vertexnormals(nodes, elements, true)
    @test isa(d, Vector{Vector{T}})
    @test length(d) == 4
    @test_approx_eq d[1] -√360 \ [-9 - √90, -3, 0]
    @test_approx_eq d[2] -√360 \ [-9 - √90, -3, 0]
    @test_approx_eq d[3] [1, 0, 0]
    @test_approx_eq d[4] -√90 \ [-9, -3, 0]
end

# check xml3d_mesh
for T in (Float64, Float32)
    # empty system
    js = parse(xml3d_mesh(Vector{T}[], Triangle{T}[]))
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
    nodes    = Vector{T}[T[0, 0, 0], T[0, 0, 3], T[0, 3, 0], T[1, -3, 3]]
    elements = [Triangle(nodes[1], nodes[2], nodes[3]),
                Triangle(nodes[1], nodes[4], nodes[2])]
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
        convert(Vector{T}, js["data"]["normal"]["seq"][1]["value"]),
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
        convert(Vector{T}, js["data"]["normal"]["seq"][1]["value"]),
        [-√360 \ [-9 - √90, -3, 0]; -√360 \ [-9 - √90, -3, 0]; [1, 0, 0]; -√90 \ [-9, -3, 0]]
    )
end

# TODO check cos
# TODO check cathetus
# TODO check sign
# TODO check distance
