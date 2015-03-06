#=
    Returns a mesh representation of the given system in a XML3D-specific JSON format.

    @param nodes
        List of nodes
    @param elements
        List of surface elements
    @return ASCIIString
=#
function xml3d_mesh{T}(nodes::Vector{Vector{T}}, elements::Vector{Element{T}}, invertnormals::Bool=false)
    idx = indexmap(nodes)
    json(Dict(
        "format" => "xml3d-json",
        "version" => "0.4.0",
        "data" => Dict(
            "index" => Dict(
                "type" => "int",
                "seq" => [Dict{ASCIIString, Vector{Int}}("value" => [idx[pointer(n)]-1 for n in unpack([Vector{T}[o.v1, o.v2, o.v3] for o in elements])])]
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
