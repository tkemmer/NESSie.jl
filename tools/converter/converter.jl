include("../../src/ProteinES.jl")

using ProteinES
using ProteinES.IO

abstract MeshType
type Nodes <: MeshType end
type Surface <: MeshType end
type Volume <: MeshType end

const formats = Dict{String, String}(
    "hmo"  =>  "hmo",
    "m"    =>  "mcsf",
    ""     =>  "msms",
    "off"  =>  "off"
)

const fin = Dict{String, Tuple{DataType, Function, String}}(
    "hmo"   =>  (Surface, readhmo, "HyperMesh/.hmo"),
    "mcsf"  =>  (Volume, readmcsf, "GAMer/.m"),
    "msms"  =>  (Surface, readmsms, "MSMS/.{face,vert}"),
    "off"   =>  (Surface, readoff, "GAMer/.off")
)

const fout = Dict{String, Tuple{DataType, Function, String}}(
    "nodes.json"    =>  (Nodes, xml3djson, "XML3D/.json"),
    "nodes.xml"     =>  (Nodes, xml3dmesh, "XML3D/.xml"),
    "surface.json"  =>  (Surface, xml3djson, "XML3D/.json")
)

writeOutput{T}(::Type{Nodes}, f::Function, nodes::Vector{Vector{T}}, ::Union{Vector{Triangle{T}}, Vector{Tetrahedron{T}}}) = println(f(nodes))

writeOutput{T}(::Type{Surface}, f::Function, nodes::Vector{Vector{T}}, elements::Vector{Triangle{T}}) = begin
    map(props!, elements)
    println(f(nodes, elements))
end

writeOutput{T}(::Type{Surface}, f::Function, nodes::Vector{Vector{T}}, elements::Vector{Tetrahedron{T}})  = println("\033[1m\033[31mERROR: Volume-to-surface mesh conversion is not supported\033[0m")

iformat = ""; ifname = ""; oformat = ""

if length(ARGS) == 2
    ifname = ARGS[1]
    oformat = ARGS[2]
    spl = split(split(ifname, "/")[end], ".")
    iformat = length(spl) == 1 ? formats[""] : formats[spl[end]]
elseif length(ARGS) == 3
    iformat = ARGS[1]
    ifname = ARGS[2]
    oformat = ARGS[3]
else
    println("\n\033[1mProteinES.jl file format converter\033[0m")
    println("==================================\n")
    println("\033[1mUsage:\033[0m\n")
    println("\tconverter.jl [<format in>] <file in> <format out>\n")
    println("If not specified, the input format will be guessed by file extension.\nMulti-file input, as provided by some file formats (e.g., MSMS), has\nto be specified without file extension. For instance,")
    println("\n\tconverter.jl msms mymesh surface.json\n")
    println("will read the MSMS-generated surface mesh from both files mymesh.face\nand mymesh.vert and convert it into a XML3D-compatible JSON format.")
    println("\n\033[1mAvailable input formats:\033[0m\n")
    for (format, (mesh, _, desc)) in fin
        println("\t\033[1m$format\033[0m\t\t$mesh\t\t$desc")
    end
    println("\n\033[1mAvailable output formats:\033[0m\n")
    for (format, (mesh, _, desc)) in fout
        println("\t\033[1m$format\033[0m\t$mesh\t\t$desc")
    end
    println()
    exit(1)
end

writeOutput(fout[oformat][1:2]..., fin[iformat][2](ifname)...)
