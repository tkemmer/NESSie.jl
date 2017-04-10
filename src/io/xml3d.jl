#=
    Returns a point cloud representation of the given system in a XML3D-specific JSON format.

    @param fname/stream
        Path or handle to (writable) JSON file
    @param nodes
        List of nodes
    @return String
=#
function writexml3d_json{T}(stream::IOStream, nodes::Vector{Vector{T}})
    println(stream, json(Dict(
        "format" => "xml3d-json",
        "version" => "0.4.0",
        "data" => Dict(
            "position" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}("value" => unpack(nodes))]
            )
        )
    )))
end
writexml3d_json{T}(fname::String, nodes::Vector{Vector{T}}) = open(fh -> writexml3d_json(fh, nodes), fname, "w")

#=
    Returns a mesh representation of the given system in a XML3D-specific JSON format.

    @param fname/stream
        Path or handle to (writable) JSON file
    @param model
        Surface model
    @return String
=#
function writexml3d_json{T}(stream::IOStream, model::SurfaceModel{T}, invertnormals::Bool=false)
    revidx = reverseindex(model.nodes)
    println(stream, json(Dict(
        "format" => "xml3d-json",
        "version" => "0.4.0",
        "data" => Dict(
            "index" => Dict(
                "type" => "int",
                "seq" => [Dict{String, Vector{Int}}("value" => [revidx[object_id(n)]-1 for n in unpack([Vector{T}[o.v1, o.v2, o.v3] for o in model.elements])])]
            ),
            "position" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}("value" => unpack(model.nodes))]
            ),
            "normal" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}("value" => unpack(vertexnormals(model.nodes, model.elements, invertnormals)))]
            )
        )
    )))
end
writexml3d_json{T}(fname::String, model::SurfaceModel{T}, invertnormals::Bool=false) =
    open(fh -> writexml3d_json(fh, model, invertnormals), fname, "w")

#=
    Returns a point cloud representation of the given system in a XML3D-specific XML format.

    @param fname/stream
        Path or handle to (writable) JSON file
    @param nodes
        List of nodes
    @return String
=#
function writexml3d_xml{T}(stream::IOStream, nodes::Vector{Vector{T}})
    xdoc = XMLDocument()
    xroot = create_root(xdoc, "xml3d")
    set_attribute(xroot, "xmlns", "http://www.xml3d.org/2009/xml3d")
    xmesh = new_child(xroot, "data")
    set_attribute(xmesh, "id", "mesh")
    xpos = new_child(xmesh, "float3")
    set_attribute(xpos, "name", "position")
    add_text(xpos, join(unpack(nodes), " "))
    println(stream, string(xdoc))
end
writexml3d_xml{T}(fname::String, nodes::Vector{Vector{T}}) = open(fh -> writexml3d_xml(fh, nodes), fname, "w")
