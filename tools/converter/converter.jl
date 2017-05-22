include("../../src/ProteinES.jl")

using ProteinES
using ProteinES.IO

abstract type MeshType end
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
    "hmo"   =>  (Surface, readhmo,  "HyperMesh/.hmo"),
    "mcsf"  =>  (Volume,  readmcsf, "GAMer/.m"),
    "msms"  =>  (Surface, readmsms, "MSMS/.{face,vert}"),
    "off"   =>  (Surface, readoff,  "GAMer/.off")
)

const fout = Dict{String, Tuple{DataType, Function, String}}(
    "nodes.json"    =>  (Nodes,   writexml3d_json, "XML3D/.json"),
    "nodes.xml"     =>  (Nodes,   writexml3d_xml,  "XML3D/.xml"),
    "surface.json"  =>  (Surface, writexml3d_json, "XML3D/.json"),
    "surface.skel"  =>  (Surface, writeskel,       "SKEL/.skel"),
    "surface.vtp"   =>  (Surface, writevtk,        "VTK/.vtp"),
    "volume.skel"   =>  (Volume,  writeskel,       "SKEL/.skel"),
    "volume.vtu"    =>  (Volume,  writevtk,        "VTK/.vtu")
)

iformat = ""; ifname = ""; oformat = ""; ofname = ""

if length(ARGS) == 3
    ifname  = ARGS[1]
    oformat = ARGS[2]
    ofname  = ARGS[3]
    spl = split(split(ifname, "/")[end], ".")
    iformat = length(spl) == 1 ? formats[""] : formats[spl[end]]
elseif length(ARGS) == 4
    iformat = ARGS[1]
    ifname  = ARGS[2]
    oformat = ARGS[3]
    ofname  = ARGS[4]
else
    println("\n\e[1mProteinES.jl file format converter\e[0m")
    println("==================================\n")
    println("\e[1mUsage:\e[0m\n")
    println("\tconverter.jl [<format in>] <file in> <format out> <file out>\n")
    println("If not specified, the input format will be guessed by file extension.")
    println("Multi-file input, as provided by some file formats (e.g., MSMS), has")
    println("to be specified without file extension. For instance,")
    println("\n\tconverter.jl msms mymesh surface.json mymesh.json\n")
    println("will read the MSMS-generated surface mesh from both files mymesh.face")
    println("and mymesh.vert and convert it into a XML3D-compatible JSON format.")
    println("\n\e[1mAvailable input formats:\e[0m\n")
    for (format, (mesh, _, desc)) in fin
        println("\t\e[1m$format\e[0m\t\t$mesh\t\t$desc")
    end
    println("\n\e[1mAvailable output formats:\e[0m\n")
    for (format, (mesh, _, desc)) in fout
        println("\t\e[1m$format\e[0m\t$mesh\t\t$desc")
    end
    println()
    exit(1)
end

mtype, f = fout[oformat][1:2]
model    = fin[iformat][2](ifname)
f(ofname, mtype === Nodes ? model.nodes : model)
