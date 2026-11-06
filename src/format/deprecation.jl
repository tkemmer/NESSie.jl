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

# deprecated in v1.6 (to be removed in v2.0)
export writeskel

function writeskel(
        stream::IOStream,
        model ::Model{T}
    ) where T
    Base.depwarn(
        "`writeskel` is deprecated and will be removed in the next major version (v2.0). " *
        "Please update your code to use alternative functionality before upgrading.",
        :writeskel
    )
    # Call original implementation
    _writeskel(stream, model)
end

function _writeskel(
        stream::IOStream,
        model ::Model{T, Triangle{T}}
    ) where T
    println(stream, "SKEL")
    println(stream, "$(length(model.nodes))\t$(length(model.elements))")
    for node in model.nodes
        println(stream, "$(node[1])\t$(node[2])\t$(node[3])")
    end
    revidx = _reverseindex(model.nodes)
    for elem in model.elements
        (v1, v2, v3) = (revidx[elem.v1] - 1,
                        revidx[elem.v2] - 1,
                        revidx[elem.v3] - 1)
        println(stream, "4\t$v1\t$v2\t$v3\t$v1")
    end
    nothing
end

function _writeskel(
        stream::IOStream,
        model ::Model{T, Tetrahedron{T}}
    ) where T
    println(stream, "SKEL")
    println(stream, "$(length(model.nodes))\t$(length(model.elements))")
    for node in model.nodes
        println(stream, "$(node[1])\t$(node[2])\t$(node[3])")
    end
    revidx = _reverseindex(model.nodes)
    for elem in model.elements
        (v1, v2, v3, v4) = (revidx[elem.v1] - 1,
                            revidx[elem.v2] - 1,
                            revidx[elem.v3] - 1,
                            revidx[elem.v4] - 1)
        println(stream, "7\t$v1\t$v2\t$v3\t$v4\t$v2\t$v4\t$v1")
    end
    nothing
end

@inline function writeskel(
        fname::AbstractString,
        model::M
    ) where {T, M <: Model{T}}
    Base.depwarn(
        "`writeskel` is deprecated and will be removed in the next major version (v2.0). " *
        "Please update your code to use alternative functionality before upgrading.",
        :writeskel
    )
    open(fh -> writeskel(fh, model), fname, "w")
end

# deprecated in v1.6 (to be removed in v2.0)
export readmcsf

function readmcsf(
        stream::IOStream,
              ::Type{T}=Float64;
        domain::Symbol=:none
    ) where T <: AbstractFloat
    Base.depwarn(
        "`readmcsf` is deprecated and will be removed in the next major version (v2.0). " *
        "Please update your code to use alternative functionality before upgrading.",
        :readmcsf
    )
    nodes = readmcsf_nodes(stream, T)
    Model(nodes, readmcsf_elements(stream, nodes, domain=domain))
end

@inline function readmcsf(
        fname ::AbstractString,
              ::Type{T}=Float64;
        domain::Symbol=:none
    ) where T
    Base.depwarn(
        "`readmcsf` is deprecated and will be removed in the next major version (v2.0). " *
        "Please update your code to use alternative functionality before upgrading.",
        :readmcsf
    )
    open(fh -> readmcsf(fh, T, domain=domain), fname)
end

@inline function readmcsf(
        fnameΩ::AbstractString,
        fnameΣ::AbstractString,
              ::Type{T}=Float64
    ) where T
    Base.depwarn(
        "`readmcsf` is deprecated and will be removed in the next major version (v2.0). " *
        "Please update your code to use alternative functionality before upgrading.",
        :readmcsf
    )
    meshunion(readmcsf(fnameΩ, T, domain=:Ω), readmcsf(fnameΣ, T, domain=:Σ))
end

function readmcsf_nodes(
        stream::IOStream,
              ::Type{T}=Float64
    ) where T <: AbstractFloat
    nodes = Vector{T}[]
    _seek(stream, "vert=[")
    for line in eachline(stream)
        startswith(line, "%") && continue        # skip comments
        startswith(line, "];") && break          # all nodes read
        push!(nodes, [parse(T, e) for e in split(line)[3:5]])
    end
    nodes
end

function readmcsf_elements(
        stream::IOStream,
        nodes ::Vector{Vector{T}};
        domain::Symbol=:none
    ) where T <: AbstractFloat
    elements = Tetrahedron{T}[]
    _seek(stream, "simp=[")
    for line in eachline(stream)
        startswith(line, "%") && continue        # skip comments
        startswith(line, "];") && break          # all simplices read
        push!(
            elements,
            Tetrahedron([nodes[parse(Int, e)+1] for e in split(line)[8:end]]..., domain)
        )
    end
    elements
end

# deprecated in v1.6 (to be removed in v2.0)
function writevtk(
        stream::IOStream,
        model ::Model{T, Tetrahedron{T}}
    ) where T
    Base.depwarn(
        "`writevtk` for volume models is deprecated and will be removed in the next major " *
        "version (v2.0). Please update your code to use alternative functionality before upgrading.",
        :writevtk
    )
    xdoc = XMLDocument()
    revidx = _reverseindex(model.nodes)

    # VTKFile
    xroot = create_root(xdoc, "VTKFile")
    set_attribute(xroot, "type", "UnstructuredGrid")

    # UnstructuredGrid/Piece
    xpiece = new_child(new_child(xroot, "UnstructuredGrid"), "Piece")
    set_attribute(xpiece, "NumberOfPoints", length(model.nodes))
    set_attribute(xpiece, "NumberOfCells",  length(model.elements))

    # Points
    xpoints = new_child(new_child(xpiece, "Points"), "DataArray")
    set_attribute(xpoints, "type", "$T")
    set_attribute(xpoints, "NumberOfComponents", "3")
    set_attribute(xpoints, "format", "ascii")
    add_text(xpoints, join(Iterators.flatten(model.nodes), " "))

    # Cells
    xcells = new_child(xpiece, "Cells")

    # Cells/connectivity
    xconn = new_child(xcells, "DataArray")
    set_attribute(xconn, "type", "Int32")
    set_attribute(xconn, "Name", "connectivity")
    set_attribute(xconn, "format", "ascii")
    add_text(xconn, join([revidx[n]-1 for n in
        Iterators.flatten(Vector{T}[o.v1, o.v2, o.v3, o.v4] for o in model.elements)], " "))

    # Cells/offsets
    xoffs = new_child(xcells, "DataArray")
    set_attribute(xoffs, "type", "Int32")
    set_attribute(xoffs, "Name", "offsets")
    set_attribute(xoffs, "format", "ascii")
    add_text(xoffs, join([4i for i in 1:length(model.elements)], " "))

    # Cells/types
    xtypes = new_child(xcells, "DataArray")
    set_attribute(xtypes, "type", "UInt8")
    set_attribute(xtypes, "Name", "types")
    set_attribute(xtypes, "format", "ascii")
    add_text(xtypes, "10 "^length(model.elements))

    # create .vtu file
    println(stream, string(xdoc))
    nothing
end
