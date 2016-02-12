using ProteinES.IO
import JSON: parse

context("xml3dmesh") do
    for T in testtypes
        # empty system
        js = parse(xml3dmesh(Vector{T}[], Triangle{T}[]))
        @fact js["format"] --> "xml3d-json"
        @fact haskey(js, "version") --> true
        @fact js["data"]["index"]["type"] --> "int"
        @fact length(js["data"]["index"]["seq"]) --> 1
        @fact js["data"]["index"]["seq"][1]["value"] --> []
        @fact js["data"]["position"]["type"] --> "float3"
        @fact length(js["data"]["position"]["seq"]) --> 1
        @fact js["data"]["position"]["seq"][1]["value"] --> []
        @fact js["data"]["normal"]["type"] --> "float3"
        @fact length(js["data"]["normal"]["seq"]) --> 1
        @fact js["data"]["normal"]["seq"][1]["value"] --> []
        # small system
        nodes    = Vector{T}[T[0, 0, 0], T[0, 0, 3], T[0, 3, 0], T[1, -3, 3]]
        elements = [Triangle(nodes[1], nodes[2], nodes[3]),
                    Triangle(nodes[1], nodes[4], nodes[2])]
        map(props!, elements)
        js = parse(xml3dmesh(nodes, elements))
        @fact js["format"] --> "xml3d-json"
        @fact haskey(js, "version") --> true
        @fact js["data"]["index"]["type"] --> "int"
        @fact length(js["data"]["index"]["seq"]) --> 1
        @fact js["data"]["index"]["seq"][1]["value"] --> [0, 1, 2, 0, 3, 1]
        @fact js["data"]["position"]["type"] --> "float3"
        @fact length(js["data"]["position"]["seq"]) --> 1
        @fact js["data"]["position"]["seq"][1]["value"] --> [0, 0, 0, 0, 0, 3, 0, 3, 0, 1, -3, 3]
        @fact js["data"]["normal"]["type"] --> "float3"
        @fact length(js["data"]["normal"]["seq"]) --> 1
        @fact map(T, js["data"]["normal"]["seq"][1]["value"]) --> roughly(map(T,
            [√360 \ [-9 - √90, -3, 0]; √360 \ [-9 - √90, -3, 0]; [-1, 0, 0]; √90 \ [-9, -3, 0]]
        ))
        # small system, inverted normals
        js = parse(xml3dmesh(nodes, elements, true))
        @fact js["format"] --> "xml3d-json"
        @fact haskey(js, "version") --> true
        @fact js["data"]["index"]["type"] --> "int"
        @fact length(js["data"]["index"]["seq"]) --> 1
        @fact js["data"]["index"]["seq"][1]["value"] --> [0, 1, 2, 0, 3, 1]
        @fact js["data"]["position"]["type"] --> "float3"
        @fact length(js["data"]["position"]["seq"]) --> 1
        @fact js["data"]["position"]["seq"][1]["value"] --> [0, 0, 0, 0, 0, 3, 0, 3, 0, 1, -3, 3]
        @fact js["data"]["normal"]["type"] --> "float3"
        @fact length(js["data"]["normal"]["seq"]) --> 1
        @fact map(T, js["data"]["normal"]["seq"][1]["value"]) --> roughly(map(T,
            [-√360 \ [-9 - √90, -3, 0]; -√360 \ [-9 - √90, -3, 0]; [1, 0, 0]; -√90 \ [-9, -3, 0]]
        ))
    end
end