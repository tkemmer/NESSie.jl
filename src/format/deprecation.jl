# deprecated in v1.6 (to be removed in v2.0)
export writexml3d_json

@inline function writexml3d_json(
        stream::IOStream,
        model::Union{Vector{Vector{T}}, Model{T, Triangle{T}}}
    ) where T
    Base.depwarn(
        "`writexml3d_json` is deprecated and will be removed in the next major version (v2.0). " *
        "Please update your code to use alternative functionality before upgrading.",
        :writexml3d_json
    )
    # Call original implementation
    _writexml3d_json(stream, model)
end

@inline function writexml3d_json(
        fname::AbstractString,
        model::Union{Vector{Vector{T}}, Model{T, Triangle{T}}}
    ) where T
    Base.depwarn(
        "`writexml3d_json` is deprecated and will be removed in the next major version (v2.0). " *
        "Please update your code to use alternative functionality before upgrading.",
        :writexml3d_json
    )
    # Call original implementation
    _writexml3d_json(fname, model)
end

# deprecated in v1.6 (to be removed in v2.0)
export writexml3d_xml

@inline function writexml3d_xml(
        stream::IOStream,
        nodes::Vector{Vector{T}}
    ) where T
    Base.depwarn(
        "`writexml3d_xml` is deprecated and will be removed in the next major version (v2.0). " *
        "Please update your code to use alternative functionality before upgrading.",
        :writexml3d_xml
    )
    # Call original implementation
    _writexml3d_xml(stream, nodes)
end

@inline function writexml3d_xml(
        fname::AbstractString,
        nodes::Vector{Vector{T}}
    ) where T
    Base.depwarn(
        "`writexml3d_xml` is deprecated and will be removed in the next major version (v2.0). " *
        "Please update your code to use alternative functionality before upgrading.",
        :writexml3d_xml
    )
    # Call original implementation
    _writexml3d_xml(fname, nodes)
end

function _writexml3d_json(
        stream::IOStream,
        nodes ::Vector{Vector{T}}
    ) where T
    JSON.json(stream, Dict(
        "format" => "xml3d-json",
        "version" => "0.4.0",
        "data" => Dict(
            "position" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}("value" => collect(T, Iterators.flatten(nodes)))]
            )
        )
    ))
    nothing
end

function _writexml3d_json(
        stream::IOStream,
        model ::Model{T, Triangle{T}}
    ) where T
    revidx = _reverseindex(model.nodes)
    JSON.json(stream, Dict(
        "format" => "xml3d-json",
        "version" => "0.4.0",
        "data" => Dict(
            "index" => Dict(
                "type" => "int",
                "seq" => [Dict{String, Vector{Int}}(
                            "value" => [revidx[n]-1 for n
                                in Iterators.flatten(Vector{T}[o.v1, o.v2, o.v3] for o
                                in model.elements)]
                          )]
            ),
            "position" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}(
                            "value" => collect(T, Iterators.flatten(model.nodes))
                          )]
            ),
            "normal" => Dict(
                "type" => "float3",
                "seq" => [Dict{String, Vector{Float64}}(
                            "value" => collect(T, Iterators.flatten(vertexnormals(model)))
                          )]
            )
        )
    ))
    nothing
end

@inline function _writexml3d_json(
        fname::AbstractString,
        model::Union{Vector{Vector{T}},Model{T, Triangle{T}}}
    ) where T
    open(fh -> writexml3d_json(fh, model), fname, "w")
end

function _writexml3d_xml(
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
    add_text(xpos, join(Iterators.flatten(nodes), " "))
    println(stream, string(xdoc))
    nothing
end

@inline function _writexml3d_xml(
        fname::AbstractString,
        nodes::Vector{Vector{T}}
    ) where T
    open(fh -> writexml3d_xml(fh, nodes), fname, "w")
end
