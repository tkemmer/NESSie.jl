#=
    Returns a point cloud representation of the given system in a XML3D-specific JSON format.

    @param nodes
        List of nodes
    @return String
=#
function xml3djson{T}(nodes::Vector{Vector{T}})
    json(Dict(
        "format" => "xml3d-json",
        "version" => "0.4.0",
        "data" => Dict(
            "position" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}("value" => unpack(nodes))]
            )
        )
    ))
end

#=
    Returns a mesh representation of the given system in a XML3D-specific JSON format.

    @param nodes
        List of nodes
    @param elements
        List of surface elements
    @return String
=#
function xml3djson{T}(nodes::Vector{Vector{T}}, elements::Vector{Triangle{T}}, invertnormals::Bool=false)
    revidx = reverseindex(nodes)
    json(Dict(
        "format" => "xml3d-json",
        "version" => "0.4.0",
        "data" => Dict(
            "index" => Dict(
                "type" => "int",
                "seq" => [Dict{String, Vector{Int}}("value" => [revidx[object_id(n)]-1 for n in unpack([Vector{T}[o.v1, o.v2, o.v3] for o in elements])])]
            ),
            "position" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}("value" => unpack(nodes))]
            ),
            "normal" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}("value" => unpack(vertexnormals(nodes, elements, invertnormals)))]
            )
        )
    ))
end

#=
    Returns a point cloud representation of the given system in a XML3D-specific XML format.

    @param nodes
        List of nodes
    @return String
=#
function xml3dmesh{T}(nodes::Vector{Vector{T}})
    xdoc = XMLDocument()
    xroot = create_root(xdoc, "xml3d")
    set_attribute(xroot, "xmlns", "http://www.xml3d.org/2009/xml3d")
    xmesh = new_child(xroot, "data")
    set_attribute(xmesh, "id", "mesh")
    xpos = new_child(xmesh, "float3")
    set_attribute(xpos, "name", "position")
    add_text(xpos, join(unpack(nodes), " "))
    string(xdoc)
end
