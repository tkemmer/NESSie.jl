# =========================================================================================
"""
    writexml3d_json(
        stream::IOStream,
        model ::Union{Vector{Vector{T}}, Model{T, Triangle{T}}}
    )

Creates a XML3D-specific JSON file either from a given collection of nodes (representing the
latter as point cloud) or from a given surface model.

# Specification
<https://github.com/xml3d/xml3d.js/wiki/External-resources>

# Return type
`Nothing`

# Alias

    writexml3d_json(
        fname::String,
        nodes::Union{Vector{Vector{T}}), Model{T, Triangle{T}}}
    )

Creates the JSON file by name rather than `IOStream` object.
"""
function writexml3d_json(
        stream::IOStream,
        nodes ::Vector{Vector{T}}
    ) where T
    JSON3.write(stream, Dict(
        "format" => "xml3d-json",
        "version" => "0.4.0",
        "data" => Dict(
            "position" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}("value" => unpack(nodes))]
            )
        )
    ))
    nothing
end

function writexml3d_json(
        stream::IOStream,
        model ::Model{T, Triangle{T}}
    ) where T
    revidx = reverseindex(model.nodes)
    JSON3.write(stream, Dict(
        "format" => "xml3d-json",
        "version" => "0.4.0",
        "data" => Dict(
            "index" => Dict(
                "type" => "int",
                "seq" => [Dict{String, Vector{Int}}(
                            "value" => [revidx[objectid(n)]-1 for n
                                in unpack([Vector{T}[o.v1, o.v2, o.v3] for o
                                in model.elements])]
                          )]
            ),
            "position" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}("value" => unpack(model.nodes))]
            ),
            "normal" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}(
                            "value" => unpack(vertexnormals(model))
                          )]
            )
        )
    ))
    nothing
end

@inline function writexml3d_json(
        fname::String,
        model::Union{Vector{Vector{T}},Model{T, Triangle{T}}}
    ) where T
    open(fh -> writexml3d_json(fh, model), fname, "w")
end


# =========================================================================================
"""
    writexml3d_xml(
        stream::IOStream,
        nodes ::Vector{Vector{T}}
    )

Creates a XML3D-specific XML file from a given collection of nodes, representing the latter
as point cloud.

# Specification
<https://github.com/xml3d/xml3d.js/wiki/External-resources>

# Return type
`Nothing`

# Alias

    writexml3d_xml(
        fname::String,
        nodes::Vector{Vector{T}}
    )

Creates the XML file by name rather than `IOStream` object.
"""
function writexml3d_xml(
        stream::IOStream,
        nodes ::Vector{Vector{T}}
    ) where T
    xdoc = XMLDocument()
    xroot = create_root(xdoc, "xml3d")
    set_attribute(xroot, "xmlns", "http://www.xml3d.org/2009/xml3d")
    xmesh = new_child(xroot, "data")
    set_attribute(xmesh, "id", "mesh")
    xpos = new_child(xmesh, "float3")
    set_attribute(xpos, "name", "position")
    add_text(xpos, join(unpack(nodes), " "))
    println(stream, string(xdoc))
    nothing
end

@inline function writexml3d_xml(
        fname::String,
        nodes::Vector{Vector{T}}
    ) where T
    open(fh -> writexml3d_xml(fh, nodes), fname, "w")
end
