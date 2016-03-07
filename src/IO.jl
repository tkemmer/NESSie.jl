module IO

using ..ProteinES
using ..ProteinES: reverseindex, unpack, vertexnormals
using JSON: json
using LightXML: XMLDocument, create_root, new_child, add_text, set_attribute

export readhmo, readmatlab, readoff, readpqr, xml3djson, xml3dmesh

include("io/hmo.jl")
include("io/matlab.jl")
include("io/off.jl")
include("io/pqr.jl")
include("io/xml3d.jl")

end # module
