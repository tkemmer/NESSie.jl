#=
    Returns a mesh representation of the given system in a XML3D-specific JSON format.

    @param nodes
        List of nodes
    @param elements
        List of surface elements
    @return ASCIIString
=#
function xml3dmesh{T}(nodes::Vector{Vector{T}}, elements::Vector{Triangle{T}}, invertnormals::Bool=false)
    revidx = reverseindex(nodes)
    json(Dict(
        "format" => "xml3d-json",
        "version" => "0.4.0",
        "data" => Dict(
            "index" => Dict(
                "type" => "int",
                "seq" => [Dict{ASCIIString, Vector{Int}}("value" => [revidx[object_id(n)]-1 for n in unpack([Vector{T}[o.v1, o.v2, o.v3] for o in elements])])]
            ),
            "position" => Dict(
                "type" => "float3",
                "seq" => [Dict{ASCIIString, Vector{Float64}}("value" => unpack(nodes))]
            ),
            "normal" => Dict(
                "type" => "float3",
                "seq" => [Dict{ASCIIString, Vector{Float64}}("value" => unpack(vertexnormals(nodes, elements, invertnormals)))]
            )
        )
    ))
end
