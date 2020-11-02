var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#NESSie.jl-1",
    "page": "Home",
    "title": "NESSie.jl",
    "category": "section",
    "text": "Nonlocal Electrostatics in Structured Solvents"
},

{
    "location": "index.html#Manual-1",
    "page": "Home",
    "title": "Manual",
    "category": "section",
    "text": "Introduction"
},

{
    "location": "man/introduction.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "man/introduction.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": "Electrostatic interactions are a major contributor to protein-protein and protein-ligand interactions. In contrast to other molecular interaction components, they can be significant over medium to long distances and are thus crucial for molecular visibility. One major challenge in this context is the treatment of the solvent the molecules are immersed in, e.g., water in a biological context. Strong simplifications of the structure of such polarizable and highly structured solvents are commonplace to achieve the required computational efficiency, but invariably lead to inaccuracies."
},

{
    "location": "man/references.html#",
    "page": "Bibliography",
    "title": "Bibliography",
    "category": "page",
    "text": ""
},

{
    "location": "man/references.html#Bibliography-1",
    "page": "Bibliography",
    "title": "Bibliography",
    "category": "section",
    "text": "[Åqv90] J. Åqvist, Ion-water interaction potentials derived from free energy pertubation simulations. J. Phys. Chem. 94: 8021, 1990.\n[Kea86] P. Keast, Moderate degree tetrahedral quadrature formulas. CMAME 55: 339-348, 1986.\n[Rad48] J. Radon, Zur mechanischen Kubatur (in German). Monatsh. für Math. 52(4): 286-300, 1948.\n[Rja90] S. Rjasanow, Vorkonditionierte iterative Auflösung von Randelementgleichungen für die Dirichlet-Aufgabe (in German). Wissenschaftliche Schriftreihe der Technischen Universität Karl-Marx-Stadt, 7/1990.\n[Ste03] O. Steinbach, Numerische Näherungsverfahren für elliptische Randwertprobleme - Finite Elemente und Randelemente (in German). Advances in Numerical Matheamtics. Teubner Verlag/GWV Fachverlage GmbH, Wiesbaden, 2003.\n[Xie16] D. Xie, H. W. Volkmer, and J. Ying, Analytical solutions of nonlocal Poisson dielectric models with multiple point charges inside a dielectric sphere. Physical Review E 93(4): 043304, 2016."
},

{
    "location": "formats/input.html#",
    "page": "Input formats",
    "title": "Input formats",
    "category": "page",
    "text": ""
},

{
    "location": "formats/input.html#Input-formats-1",
    "page": "Input formats",
    "title": "Input formats",
    "category": "section",
    "text": "    CurrentModule = NESSie.FormatCurrently supported input file formats with different models:File type Surface model Volume model Charges included\nHMO ✓  ✓\nMcsf  ✓ \nMSMS ✓  \nOFF ✓  \nPQR   ✓ (charges only)"
},

{
    "location": "formats/input.html#NESSie.Format.readhmo",
    "page": "Input formats",
    "title": "NESSie.Format.readhmo",
    "category": "Function",
    "text": "readhmo{T <: AbstractFloat}(\n    stream::IOStream,\n          ::Type{T}=Float64\n)\n\nReads a complete surface model from the given HMO file.\n\nReturn type\n\nModel{T, Triangle{T}}\n\nAlias\n\nreadhmo{T}(fname::String, ::Type{T}=Float64)\n\nReads the model using a file name rather than a IOStream object.\n\n\n\n"
},

{
    "location": "formats/input.html#HMO-1",
    "page": "Input formats",
    "title": "HMO",
    "category": "section",
    "text": "    readhmo"
},

{
    "location": "formats/input.html#NESSie.Format.readmcsf",
    "page": "Input formats",
    "title": "NESSie.Format.readmcsf",
    "category": "Function",
    "text": "readmcsf{T <: AbstractFloat}(\n    stream::IOStream,\n          ::Type{T}=Float64;\n    # kwargs\n    domain::Symbol=:none\n)\n\nReads a volume model from the given GAMer-generated mcsf file.\n\nnote: Note\nThis file type does not support charge models! Hence, the charge list of the returning Model object is empty and has to be set separately.\n\nArguments\n\ndomain Element domain (solute :Ω, solvent :Σ, or :none)\n\nReturn type\n\nModel{T, Tetrahedron{T}}\n\nAliases\n\nreadmcsf{T}(\n    fname ::String,\n          ::Type{T}=Float64;\n    # kwargs\n    domain::Symbol=:none\n)\n\nReads the model using a file name rather than a IOStream object.\n\nreadmcsf{T}(\n    fnameΩ::String,\n    fnameΣ::String,\n          ::Type{T}=Float64\n)\n\nReads the contents of two separate files (given by name), sets the element domains to :Ω or :Σ, respectively, and returns a single (merged) Model object.\n\n\n\n"
},

{
    "location": "formats/input.html#Mcsf-1",
    "page": "Input formats",
    "title": "Mcsf",
    "category": "section",
    "text": "    readmcsf"
},

{
    "location": "formats/input.html#NESSie.Format.readmsms",
    "page": "Input formats",
    "title": "NESSie.Format.readmsms",
    "category": "Function",
    "text": "readmsms{T <: AbstractFloat}(\n    vertstream::IOStream,\n    facestream::IOStream,\n              ::Type{T}=Float64\n)\n\nReads a surface model from the given MSMS-generated .face and .vert files.\n\nnote: Note\nThis file type does not support charge models! Hence, the charge list of the returning Model object is empty and has to be set separately.\n\nReturn type\n\nModel{T, Triangle{T}}\n\nAlias\n\nreadmsms{T}(fname::String, ::Type{T}=Float64)\n\nReads the model using a common file name prefix (fname.{vert,face}) for both files rather than IOStream objects.\n\n\n\n"
},

{
    "location": "formats/input.html#MSMS-1",
    "page": "Input formats",
    "title": "MSMS",
    "category": "section",
    "text": "    readmsms"
},

{
    "location": "formats/input.html#NESSie.Format.readoff",
    "page": "Input formats",
    "title": "NESSie.Format.readoff",
    "category": "Function",
    "text": "readoff{T <: AbstractFloat}(\n    stream::IOStream,\n          ::Type{T}=Float64\n)\n\nReads a surface model from the given OFF file.\n\nnote: Note\nThis file type does not support charge models! Hence, the charge list of the returning Model object is empty and has to be set separately.\n\nSpecification\n\nhttp://www.geomview.org/docs/html/OFF.html\n\nReturn type\n\nModel{T}\n\nAlias\n\nreadoff{T}(fname::String, ::Type{T}=Float64)\n\nReads the model using a file name rather than a IOStream object.\n\n\n\n"
},

{
    "location": "formats/input.html#OFF-1",
    "page": "Input formats",
    "title": "OFF",
    "category": "section",
    "text": "    readoff"
},

{
    "location": "formats/input.html#NESSie.Format.readpqr",
    "page": "Input formats",
    "title": "NESSie.Format.readpqr",
    "category": "Function",
    "text": "readpqr{T <: AbstractFloat}(\n    stream::IOStream,\n          ::Type{T}=Float64\n)\n\nReads a charge model from the given PQR file.\n\nReturn type\n\nVector{Charge{T}}\n\nAlias\n\nreadpqr{T}(fname::String, ::Type{T}=Float64)\n\nReads the charge model using a file name rather than a IOStream object.\n\n\n\n"
},

{
    "location": "formats/input.html#PQR-1",
    "page": "Input formats",
    "title": "PQR",
    "category": "section",
    "text": "    readpqr"
},

{
    "location": "formats/output.html#",
    "page": "Output formats",
    "title": "Output formats",
    "category": "page",
    "text": ""
},

{
    "location": "formats/output.html#Output-formats-1",
    "page": "Output formats",
    "title": "Output formats",
    "category": "section",
    "text": "    CurrentModule = NESSie.FormatCurrently supported output file formats with different models:File type Point cloud Surface model Volume model\nSKEL  ✓ ✓\nVTK  ✓ ✓\nXML3D/JSON ✓ ✓ \nXML3D/XML ✓  "
},

{
    "location": "formats/output.html#NESSie.Format.writeskel",
    "page": "Output formats",
    "title": "NESSie.Format.writeskel",
    "category": "Function",
    "text": "writeskel{T, M <: Model{T}}(\n    stream::IOStream,\n    model ::M\n)\n\nCreates a SKEL file from a given surface or volume model, representing the model as a collection of points and polylines.\n\nSpecification\n\nhttp://www.geomview.org/docs/html/SKEL.html\n\nReturn type\n\nVoid\n\nAlias\n\nwriteskel{T, M <: Model{T}}(\n    fname::String,\n    model::M\n)\n\nCreates the SKEL file by name rather than IOStream object.\n\n\n\n"
},

{
    "location": "formats/output.html#SKEL-1",
    "page": "Output formats",
    "title": "SKEL",
    "category": "section",
    "text": "    writeskel"
},

{
    "location": "formats/output.html#NESSie.Format.writevtk",
    "page": "Output formats",
    "title": "NESSie.Format.writevtk",
    "category": "Function",
    "text": "function writevtk{T, M <: Model{T}}(\n    stream::IOStream,\n    model ::M\n)\n\nCreates a VTK-compatible output file from a given surface or volume model. The exact file type is determined by the given model:\n\nModel type Resulting file type\nSurface model VTK PolyData\nVolume model VTK UnstructuredGrid\n\nSpecification\n\nhttp://www.vtk.org/VTK/img/file-formats.pdf\n\nReturn type\n\nVoid\n\nAlias\n\nwritevtk{T, M <: Model{T}}(\n    fname::String,\n    model::M\n)\n\nCreates the VTK file by name rather than IOStream object.\n\n\n\n"
},

{
    "location": "formats/output.html#VTK-1",
    "page": "Output formats",
    "title": "VTK",
    "category": "section",
    "text": "    writevtk"
},

{
    "location": "formats/output.html#NESSie.Format.writexml3d_json",
    "page": "Output formats",
    "title": "NESSie.Format.writexml3d_json",
    "category": "Function",
    "text": "writexml3d_json{T}(\n    stream::IOStream,\n    model ::Union{Vector{Vector{T}}, Model{T, Triangle{T}}}\n)\n\nCreates a XML3D-specific JSON file either from a given collection of nodes (representing the latter as point cloud) or from a given surface model.\n\nSpecification\n\nhttps://github.com/xml3d/xml3d.js/wiki/External-resources\n\nReturn type\n\nVoid\n\nAlias\n\nwritexml3d_json{T}(\n    fname::String,\n    nodes::Union{Vector{Vector{T}}), Model{T, Triangle{T}}}\n)\n\nCreates the JSON file by name rather than IOStream object.\n\n\n\n"
},

{
    "location": "formats/output.html#NESSie.Format.writexml3d_xml",
    "page": "Output formats",
    "title": "NESSie.Format.writexml3d_xml",
    "category": "Function",
    "text": "writexml3d_xml{T}(\n    stream::IOStream,\n    nodes ::Vector{Vector{T}}\n)\n\nCreates a XML3D-specific XML file from a given collection of nodes, representing the latter as point cloud.\n\nSpecification\n\nhttps://github.com/xml3d/xml3d.js/wiki/External-resources\n\nReturn type\n\nVoid\n\nAlias\n\nwritexml3d_xml{T}(\n    fname::String,\n    nodes::Vector{Vector{T}}\n)\n\nCreates the XML file by name rather than IOStream object.\n\n\n\n"
},

{
    "location": "formats/output.html#XML3D-1",
    "page": "Output formats",
    "title": "XML3D",
    "category": "section",
    "text": "    writexml3d_json\n    writexml3d_xml"
},

{
    "location": "lib/constants.html#",
    "page": "Constants",
    "title": "Constants",
    "category": "page",
    "text": ""
},

{
    "location": "lib/constants.html#Constants-1",
    "page": "Constants",
    "title": "Constants",
    "category": "section",
    "text": ""
},

{
    "location": "lib/constants.html#NESSie.ε0",
    "page": "Constants",
    "title": "NESSie.ε0",
    "category": "Constant",
    "text": "Vacuum permittivity\n\nUnit\n\nfracFm\n\n\n\n"
},

{
    "location": "lib/constants.html#Global-constants-1",
    "page": "Constants",
    "title": "Global constants",
    "category": "section",
    "text": "    ε0"
},

{
    "location": "lib/constants.html#NESSie.Option",
    "page": "Constants",
    "title": "NESSie.Option",
    "category": "Type",
    "text": "struct Option{T <: AbstractFloat}\n    εΩ::T       # dielectric constant of the solute\n    εΣ::T       # dielectric constant of the solvent\n    ε∞::T       # large-scale (bulk) solvent response\n    λ ::T       # correlation length scale [λ] = Å\nend\n\nSystem parameters\n\n\n\n"
},

{
    "location": "lib/constants.html#NESSie.defaultopt",
    "page": "Constants",
    "title": "NESSie.defaultopt",
    "category": "Function",
    "text": "defaultopt(T::Type{Float64} = Float64)\ndefaultopt(T::Type{Float32})\n\nDefault system parameters\n\nReturn type\n\nOption{T}\n\nDefault values\n\n_ = 2\n_ = 78\n_ = 18\n = 20 \n\n\n\n"
},

{
    "location": "lib/constants.html#System-constants-1",
    "page": "Constants",
    "title": "System constants",
    "category": "section",
    "text": "    Option\n    defaultopt"
},

{
    "location": "lib/electrostatics.html#",
    "page": "Electrostatics",
    "title": "Electrostatics",
    "category": "page",
    "text": ""
},

{
    "location": "lib/electrostatics.html#Electrostatics-1",
    "page": "Electrostatics",
    "title": "Electrostatics",
    "category": "section",
    "text": ""
},

{
    "location": "lib/electrostatics.html#NESSie.PotentialType",
    "page": "Electrostatics",
    "title": "NESSie.PotentialType",
    "category": "Type",
    "text": "abstract type PotentialType end\nstruct SingleLayer <: PotentialType end\nstruct DoubleLayer <: PotentialType end\n\nEnum-like representation of single and double layer potentials\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#NESSie.SingleLayer",
    "page": "Electrostatics",
    "title": "NESSie.SingleLayer",
    "category": "Type",
    "text": "abstract type PotentialType end\nstruct SingleLayer <: PotentialType end\nstruct DoubleLayer <: PotentialType end\n\nEnum-like representation of single and double layer potentials\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#NESSie.DoubleLayer",
    "page": "Electrostatics",
    "title": "NESSie.DoubleLayer",
    "category": "Type",
    "text": "abstract type PotentialType end\nstruct SingleLayer <: PotentialType end\nstruct DoubleLayer <: PotentialType end\n\nEnum-like representation of single and double layer potentials\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#Potential-types-1",
    "page": "Electrostatics",
    "title": "Potential types",
    "category": "section",
    "text": "    PotentialType\n    SingleLayer\n    DoubleLayer"
},

{
    "location": "lib/electrostatics.html#NESSie.LocalityType",
    "page": "Electrostatics",
    "title": "NESSie.LocalityType",
    "category": "Type",
    "text": "abstract type LocalityType end\nstruct NonlocalES <: LocalityType end\nstruct LocalES    <: LocalityType end\n\nEnum-like representation of locality assumption:\n\nLocal electrostatics: Complete independence of solvent molecules\nNonlocal electrostatics: Allow solvent molecule correlation effects (with area-of-effect radius λ)\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#NESSie.LocalES",
    "page": "Electrostatics",
    "title": "NESSie.LocalES",
    "category": "Type",
    "text": "abstract type LocalityType end\nstruct NonlocalES <: LocalityType end\nstruct LocalES    <: LocalityType end\n\nEnum-like representation of locality assumption:\n\nLocal electrostatics: Complete independence of solvent molecules\nNonlocal electrostatics: Allow solvent molecule correlation effects (with area-of-effect radius λ)\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#NESSie.NonlocalES",
    "page": "Electrostatics",
    "title": "NESSie.NonlocalES",
    "category": "Type",
    "text": "abstract type LocalityType end\nstruct NonlocalES <: LocalityType end\nstruct LocalES    <: LocalityType end\n\nEnum-like representation of locality assumption:\n\nLocal electrostatics: Complete independence of solvent molecules\nNonlocal electrostatics: Allow solvent molecule correlation effects (with area-of-effect radius λ)\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#Locality-assumption-1",
    "page": "Electrostatics",
    "title": "Locality assumption",
    "category": "section",
    "text": "    LocalityType\n    LocalES\n    NonlocalES"
},

{
    "location": "lib/electrostatics.html#Potentials-1",
    "page": "Electrostatics",
    "title": "Potentials",
    "category": "section",
    "text": ""
},

{
    "location": "lib/electrostatics.html#NESSie.φmol",
    "page": "Electrostatics",
    "title": "NESSie.φmol",
    "category": "Function",
    "text": "φmol{T}(\n    ξ        ::Vector{T},\n    charges  ::Vector{Charge{T}};\n    # kwargs\n    tolerance::T                 = T(1e-10)\n)\n\nComputes and returns the molecular potential of the given system of point charges in a structureless medium for the given observation point ξ:\n\n_mol() = frac14 _0 _ sum_i fracqr-\n\nIf r- is smaller than the given tolerance, the value is replaced by tolerance for the affected charge.\n\nnote: Note\nThe return value is premultiplied by 4    _\n\nReturn type\n\nT\n\nAliases\n\nφmol{T}(\n    Ξ        ::Vector{Vector{T}},\n    charges  ::Vector{Charge{T}};\n    # kwargs\n    tolerance::T                 = T(1e-10)\n)\n\nComputes the molecular potentials for a list of observation points.\n\nφmol{T}(\n    model    ::Model{T, Triangle{T}};\n    # kwargs\n    tolerance::T                     = T(1e-10)\n)\n\nComputes the molecular potentials for the given surface model, using each triangle center as observation point.\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#NESSie.∂ₙφmol",
    "page": "Electrostatics",
    "title": "NESSie.∂ₙφmol",
    "category": "Function",
    "text": "∂ₙφmol{T}(ξ::Triangle{T}, charges::Vector{Charge{T}})\n\nComputes and returns the normal derivative of the given system's molecular potential in a structureless medium, using the given triangle's center as observation point and the triangle's normal as reference normal.\n\n_mol() = -frac14 _0 _ sum_i fracqr- (r-)  n\n\nnote: Note\nThe return value is premultiplied by 4    _\n\nReturn type\n\nT\n\nAliases\n\n∂ₙφmol{T}(model::Model{T, Triangle{T}})\n\nComputes the normal derivatives of the molecular potentials for the given surface model, using each triangle center and normal as observation point.\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#NESSie.∇φmol",
    "page": "Electrostatics",
    "title": "NESSie.∇φmol",
    "category": "Function",
    "text": "∇φmol{T}(ξ::Vector{T}, charges::Vector{Charge{T}})\n\nComputes and returns the gradient of the given system's molecular potential in a structureless medium for the given observation point ξ.\n\n_mol() = -frac14 _0 _ sum_i fracqr- (r-)\n\nnote: Note\nThe return value is premultiplied by 4    _\n\nReturn type\n\nVector{T}\n\nAliases\n\n∇φmol{T}(Ξ::Vector{Vector{T}})\n\nComputes the molecular potential gradients for a list of observation points.\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#Molecular-potentials-1",
    "page": "Electrostatics",
    "title": "Molecular potentials",
    "category": "section",
    "text": "    φmol\n    ∂ₙφmol\n    NESSie.∇φmol"
},

{
    "location": "lib/electrostatics.html#NESSie.BEM.φΩ",
    "page": "Electrostatics",
    "title": "NESSie.BEM.φΩ",
    "category": "Function",
    "text": "φΩ{T, B <: BEMResult{T}}(\n    Ξ         ::Vector{Vector{T}},\n    bem       ::B\n)\n\nComputes the local or nonlocal interior electrostatic potential _ for the given set of observation points Ξ.\n\nwarning: Warning\nThis function does not verify whether all points in Ξ are located in !\n\nUnit\n\nV = fracCF\n\nReturn type\n\nVector{T}\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#NESSie.TestModel.φΩ",
    "page": "Electrostatics",
    "title": "NESSie.TestModel.φΩ",
    "category": "Function",
    "text": "φΩ{T, L <: LocalityType}(\n       ::Type{L},\n    ξ  ::Vector{T},\n    ion::BornIon{T},\n    opt::Option{T} = defaultopt(T)\n)\n\nComputes the interior local or nonlocal electrostatic potential _ for the given observation point .\n\nUnit\n\nV = fracCF\n\nReturn type\n\nT\n\nwarning: Warning\nThis function does not verify whether ξ is located inside of the sphere!\n\n\n\nfunction φΩ{T}(\n    ξ    ::Vector{T},\n    model::NonlocalXieModel1{T}\n)\n\nComputes the interior nonlocal electrostatic potential _ for the given observation point .\n\nUnit\n\nV = fracCF\n\nReturn type\n\nT\n\nwarning: Warning\nThis function does not verify whether ξ is located inside of the sphere!\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#Interior-potentials-1",
    "page": "Electrostatics",
    "title": "Interior potentials",
    "category": "section",
    "text": "    NESSie.BEM.φΩ\n    NESSie.TestModel.φΩ"
},

{
    "location": "lib/electrostatics.html#NESSie.BEM.φΣ",
    "page": "Electrostatics",
    "title": "NESSie.BEM.φΣ",
    "category": "Function",
    "text": "φΣ{T, B <: BEMResult{T}}(\n    Ξ         ::Vector{Vector{T}},\n    bem       ::B\n)\n\nComputes the local or nonlocal exterior electrostatic potential _ for the given set of observation points Ξ.\n\nwarning: Warning\nThis function does not verify whether all points in Ξ are located in !\n\nUnit\n\nV = fracCF\n\nReturn type\n\nVector{T}\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#NESSie.TestModel.φΣ",
    "page": "Electrostatics",
    "title": "NESSie.TestModel.φΣ",
    "category": "Function",
    "text": "φΣ{T, L <: LocalityType}(\n       ::Type{L},\n    ξ  ::Vector{T},\n    ion::BornIon{T},\n    opt::Option{T} = defaultopt(T)\n)\n\nComputes the exterior local or nonlocal electrostatic potential _ for the given observation point .\n\nUnit\n\nV = fracCF\n\nReturn type\n\nT\n\nwarning: Warning\nThis function does not verify whether ξ is located outside of the sphere!\n\n\n\nfunction φΣ{T}(\n    ξ    ::Vector{T},\n    model::NonlocalXieModel1{T}\n)\n\nComputes the exterior nonlocal electrostatic potential _ for the given observation point .\n\nUnit\n\nV = fracCF\n\nReturn type\n\nT\n\nwarning: Warning\nThis function does not verify whether ξ is located outside of the sphere!\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#Exterior-potentials-1",
    "page": "Electrostatics",
    "title": "Exterior potentials",
    "category": "section",
    "text": "    NESSie.BEM.φΣ\n    NESSie.TestModel.φΣ"
},

{
    "location": "lib/electrostatics.html#NESSie.BEM.rfenergy",
    "page": "Electrostatics",
    "title": "NESSie.BEM.rfenergy",
    "category": "Function",
    "text": "rfenergy{T, B <: BEMResult{T}}(\n    bem       ::B\n)\n\nComputes the local or nonlocal reaction field energy W* as\n\nW^* = ^*  quad d\n\nwhere ^* is the reaction field and  is the corresponding charge distribution.\n\nUnit\n\nfrackJmol\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "lib/electrostatics.html#Energies-1",
    "page": "Electrostatics",
    "title": "Energies",
    "category": "section",
    "text": "    NESSie.BEM.rfenergy"
},

{
    "location": "lib/models.html#",
    "page": "Models",
    "title": "Models",
    "category": "page",
    "text": ""
},

{
    "location": "lib/models.html#Models-1",
    "page": "Models",
    "title": "Models",
    "category": "section",
    "text": ""
},

{
    "location": "lib/models.html#NESSie.Element",
    "page": "Models",
    "title": "NESSie.Element",
    "category": "Type",
    "text": "abstract type Element{T <: AbstractFloat} end\nabstract type SurfaceElement{T} <: Element{T} end\nabstract type VolumeElement{T}  <: Element{T} end\n\nAbstract base types for all elements.\n\n\n\n"
},

{
    "location": "lib/models.html#NESSie.Triangle",
    "page": "Models",
    "title": "NESSie.Triangle",
    "category": "Type",
    "text": "struct Triangle{T} <: SurfaceElement{T}\n    v2      ::Vector{T}   # position of the second node\n    v1      ::Vector{T}   # position of the first node\n    v3      ::Vector{T}   # position of the third node\n    center  ::Vector{T}   # centroid of the triangle\n    normal  ::Vector{T}   # normal vector of the triangle\n    area    ::T           # area of the triangle\n    distorig::T           # distance to the origin\nend\n\nRepresentation of a single surface triangle.\n\nSpecial constructors\n\nTriangle{T}(\n    v1  ::Vector{T},\n    v2  ::Vector{T},\n    v3  ::Vector{T}\n)\n\nMost commonly used constructor variant for creating a triangle by only specifying its nodes. The remaining member variables will automatically be computed via props.\n\n\n\n"
},

{
    "location": "lib/models.html#NESSie.Tetrahedron",
    "page": "Models",
    "title": "NESSie.Tetrahedron",
    "category": "Type",
    "text": "struct Tetrahedron{T} <: VolumeElement{T}\n    v1::Vector{T}       # position of the first node\n    v2::Vector{T}       # position of the second node\n    v3::Vector{T}       # position of the third node\n    v4::Vector{T}       # position of the fourth node\n    domain::Symbol      # element domain (solvent :Σ, solute :Ω, or :none)\nend\n\nRepresentation of a single tetrahedron.\n\nSpecial constructors\n\nTetrahedron{T}(\n    v1::Vector{T},\n    v2::Vector{T},\n    v3::Vector{T},\n    v4::Vector{T}\n)\n\nSets domain to :none.\n\n\n\n"
},

{
    "location": "lib/models.html#Elements-1",
    "page": "Models",
    "title": "Elements",
    "category": "section",
    "text": "    Element\n    Triangle\n    Tetrahedron"
},

{
    "location": "lib/models.html#NESSie.Charge",
    "page": "Models",
    "title": "NESSie.Charge",
    "category": "Type",
    "text": "struct Charge{T <: AbstractFloat}\n    pos::Vector{T}  # position of the charge\n    val::T          # charge value\nend\n\nRepresentation of a single point charge.\n\nSpecial constructors\n\nCharge{T}(\n    posx::T,\n    posy::T,\n    posz::T,\n    val ::T\n)\n\nConstructor variant with flat argument list for pos.\n\n\n\n"
},

{
    "location": "lib/models.html#Charge-models-1",
    "page": "Models",
    "title": "Charge models",
    "category": "section",
    "text": "    Charge"
},

{
    "location": "lib/models.html#NESSie.Model",
    "page": "Models",
    "title": "NESSie.Model",
    "category": "Type",
    "text": "mutable struct Model{T, E <: Element{T}}\n    nodes   ::Vector{Vector{T}  = Vector{T}[]    # mesh nodes\n    elements::Vector{E}         = E[]            # mesh elements\n    charges ::Vector{Charge{T}} = Charge{T}[]    # point charges in the molecule\n    params  ::Option{T}         = defaultopt(T)  # system constants\nend\n\nSystem model representing a biomelecule in solvation, including a collection of point charges in the molecule and a set of system constants. The system can either be represented as a surface model (e.g., a collection of molecule surface triangles) or as a volume model (e.g., a collection of tetrahedra for the molecule and its surrounding space).\n\n\n\n"
},

{
    "location": "lib/models.html#NESSie.TestModel.BornIon",
    "page": "Models",
    "title": "NESSie.TestModel.BornIon",
    "category": "Type",
    "text": "struct BornIon{T <: AbstractFloat}\n    charge::Charge{T}  # point charge at the sphere's center\n    radius::T          # sphere radius in Å\nend\n\nSingle Born ion, that is, a monoatomic ion represented as a spherically symmetric domain with a single point charge located at its center (_ = 1).\n\nSpecial constructors\n\nBornIon{T}(charge::T, radius::T)\n\nCenters the sphere at (0 0 0)^T.\n\n\n\n"
},

{
    "location": "lib/models.html#NESSie.TestModel.XieModel",
    "page": "Models",
    "title": "NESSie.TestModel.XieModel",
    "category": "Type",
    "text": "mutable struct XieModel{T}\n    radius ::T                                 # radius of the origin-centered sphere\n    charges::Vector{Charge{T}}                 # point charges in the sphere\n    params ::Option{T}         = defaultopt(T) # system constants\nend\n\nSystem model of a dielectric sphere containing multiple point charges. On construction, the given point charge model will be translated and rescaled to fit inside an origin-centered sphere with the specified radius.\n\nSpecial contructors\n\nXieModel{T}(\n    radius ::T,\n    charges::Vector{Charge{T}},\n    params ::Option{T}          = defaultopt(T);\n    # kwargs\n    compat ::Bool               = false\n)\n\ncompat enables the compatibility mode and scales the model exactly like the reference implementation ([Xie16]). Use this flag if you intend to compare the results to the reference.\n\n\n\n"
},

{
    "location": "lib/models.html#NESSie.TestModel.NonlocalXieModel1",
    "page": "Models",
    "title": "NESSie.TestModel.NonlocalXieModel1",
    "category": "Type",
    "text": "struct NonlocalXieModel1{T}\n    radius ::T                 # radius of the origin-centered sphere\n    charges::Vector{Charge{T}} # point charges in the sphere\n    params ::Option{T}         # system constants\n    len    ::Int               # number of terms to be computed\n    A₁     ::Array{T, 2}       # coefficients A₁ for each charge\n    A₂     ::Array{T, 2}       # coefficients A₂ for each charge\n    A₃     ::Array{T, 2}       # coefficients A₃ for each charge\nend\n\nRepresentation of the first nonlocal Poisson dielectric model described in [Xie16]. This model comprises a full XieModel and the coefficients A_in with i = 1 2 3 (cf. Eqs. (20a-c)) for each point charge in the model, which are used in the computation of the electrostatic potentials.\n\nnote: Note\nThis type does not provide a trivial constructor.\n\nConstructor\n\nNonlocalXieModel1{T}(\n    model::XieModel{T},\n    len  ::Int\n)\n\nThe model is created solely from the given XieModel and the number of terms to be used to approximate the original infinite sum (Eq. 18). The coefficient vectors are computed automatically via coefficients.\n\n\n\n"
},

{
    "location": "lib/models.html#System-models-1",
    "page": "Models",
    "title": "System models",
    "category": "section",
    "text": "    Model\n    NESSie.TestModel.BornIon\n    NESSie.TestModel.XieModel\n    NESSie.TestModel.NonlocalXieModel1"
},

{
    "location": "lib/quadrature.html#",
    "page": "Quadrature",
    "title": "Quadrature",
    "category": "page",
    "text": ""
},

{
    "location": "lib/quadrature.html#Quadrature-1",
    "page": "Quadrature",
    "title": "Quadrature",
    "category": "section",
    "text": ""
},

{
    "location": "lib/quadrature.html#NESSie.QuadraturePoints",
    "page": "Quadrature",
    "title": "NESSie.QuadraturePoints",
    "category": "Type",
    "text": "abstract type QuadraturePoints{T <: AbstractFloat} end\n\nAbstract base type for all sets of quadrature points\n\n\n\n"
},

{
    "location": "lib/quadrature.html#NESSie.QuadPts2D",
    "page": "Quadrature",
    "title": "NESSie.QuadPts2D",
    "category": "Type",
    "text": "struct QuadPts2D{T} <: QuadraturePoints{T}\n    num   ::Int         # number of points\n    x     ::Vector{T}   # x values\n    y     ::Vector{T}   # y values\n    weight::Vector{T}   # weights\nend\n\nQuadrature points and weights for triangles\n\n\n\n"
},

{
    "location": "lib/quadrature.html#NESSie.QuadPts3D",
    "page": "Quadrature",
    "title": "NESSie.QuadPts3D",
    "category": "Type",
    "text": "struct QuadPts3D{T} <: QuadraturePoints{T}\n    num   ::Int         # number of points\n    x     ::Vector{T}   # x values\n    y     ::Vector{T}   # y values\n    z     ::Vector{T}   # z values\n    weight::Vector{T}   # weights\nend\n\nQuadrature points and weights for tetrahedra\n\n\n\n"
},

{
    "location": "lib/quadrature.html#Quadrature-points-1",
    "page": "Quadrature",
    "title": "Quadrature points",
    "category": "section",
    "text": "    QuadraturePoints\n    QuadPts2D\n    QuadPts3D"
},

{
    "location": "lib/quadrature.html#NESSie.quadraturepoints",
    "page": "Quadrature",
    "title": "NESSie.quadraturepoints",
    "category": "Function",
    "text": "quadraturepoints(::Type{Tetrahedron{Float64}})\nquadraturepoints(::Type{Tetrahedron{Float32}})\nquadraturepoints(::Type{Triangle{Float64}})\nquadraturepoints(::Type{Triangle{Float32}})\n\nGenerator function for quadrature points:\n\nTriangles: 7 points per element [Rad48]\nTetrahedra: 5 points per element [Kea86]\n\nReturn type\n\nQuadPts2D or QuadPts3D\n\n\n\n"
},

{
    "location": "lib/quadrature.html#Generators-1",
    "page": "Quadrature",
    "title": "Generators",
    "category": "section",
    "text": "    quadraturepoints"
},

{
    "location": "lib/quadrature.html#NESSie.Rjasanow.laplacecoll!",
    "page": "Quadrature",
    "title": "NESSie.Rjasanow.laplacecoll!",
    "category": "Function",
    "text": "laplacecoll!{T, P <: PotentialType}(\n    ptype   ::Type{P},\n    dest    ::DenseArray{T,1},\n    elements::Vector{Triangle{T}},\n    Ξ       ::Vector{Vector{T}},\n    fvals   ::Union{DenseArray{T,1},SubArray{T,1}}\n)\n\nlaplacecoll!{T, P <: PotentialType}(\n    ptype   ::Type{P},\n    dest    ::DenseArray{T, 2},\n    elements::Vector{Triangle{T}},\n    Ξ       ::Vector{Vector{T}}\n)\n\nAnalytical solution for the single or double layer Laplace potential for a given list of triangles and observation points Ξ [Rja90].\n\nThe first version of this function uses a vector as destination dest, where each element represents the dot product of the corresponding coefficient matrix row and the fvals vector.\n\nnote: Note\nThe result is premultiplied by 4π.\n\nReturn type\n\nVoid\n\n\n\n"
},

{
    "location": "lib/quadrature.html#Laplace-potential-1",
    "page": "Quadrature",
    "title": "Laplace potential",
    "category": "section",
    "text": "    NESSie.Rjasanow.laplacecoll!"
},

{
    "location": "lib/quadrature.html#NESSie.Radon.regularyukawacoll!",
    "page": "Quadrature",
    "title": "NESSie.Radon.regularyukawacoll!",
    "category": "Function",
    "text": "regularyukawacoll!{T, P <: PotentialType}(\n            ::Type{P},\n    dest    ::DenseArray{T,1},\n    elements::Vector{Triangle{T}},\n    Ξ       ::Vector{Vector{T}},\n    yukawa  ::T,\n    fvals   ::Union{DenseArray{T,1},SubArray{T,1}}\n)\n\nregularyukawacoll!{T, P <: PotentialType}(\n            ::Type{P},\n    dest    ::DenseArray{T,2},\n    elements::Vector{Triangle{T}},\n    Ξ       ::Vector{Vector{T}},\n    yukawa  ::T\n)\n\nComputes the regular part of the single or double layer Yukawa potential (that is, Yukawa minus Laplace) using a seven-point Radon cubature [Rad48] for a given list of triangles and observation points Ξ.\n\nThe first version of this function uses a vector as destination dest, where each element represents the dot product of the corresponding coefficient matrix row and the fvals vector.\n\nnote: Note\nThe result is premultiplied by 4π.\n\nArguments\n\nyukawa Exponent of the Yukawa operator's fundamental solution\n\nReturn type\n\nVoid\n\n\n\n"
},

{
    "location": "lib/quadrature.html#Yukawa-potential-1",
    "page": "Quadrature",
    "title": "Yukawa potential",
    "category": "section",
    "text": "    NESSie.Radon.regularyukawacoll!"
},

{
    "location": "lib/solvers.html#",
    "page": "Solvers",
    "title": "Solvers",
    "category": "page",
    "text": ""
},

{
    "location": "lib/solvers.html#Solvers-1",
    "page": "Solvers",
    "title": "Solvers",
    "category": "section",
    "text": ""
},

{
    "location": "lib/solvers.html#NESSie.BEM.BEMResult",
    "page": "Solvers",
    "title": "NESSie.BEM.BEMResult",
    "category": "Type",
    "text": "abstract type BEMResult{T, E <: SurfaceElement{T}} end\n\nAbstract base type for all BEM solver results\n\n\n\n"
},

{
    "location": "lib/solvers.html#NESSie.BEM.LocalBEMResult",
    "page": "Solvers",
    "title": "NESSie.BEM.LocalBEMResult",
    "category": "Type",
    "text": "struct LocalBEMResult{T, E} <: BEMResult{T, E}\n    model::Model{T, E}\n    u    ::Vector{T}   # [γ₀int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0\n    q    ::Vector{T}   # [γ₁int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0\n    umol ::Vector{T}   # [γ₀int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0\n    qmol ::Vector{T}   # [γ₁int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0\nend\n\nResult data of the local solving process to be used for potential computation and post-processing, with Ξ being the list of observation points, that is, the set of triangle centroids.\n\n\n\n"
},

{
    "location": "lib/solvers.html#NESSie.BEM.NonlocalBEMResult",
    "page": "Solvers",
    "title": "NESSie.BEM.NonlocalBEMResult",
    "category": "Type",
    "text": "struct NonlocalBEMResult{T, E} <: BEMResult{T, E}\n    model::Model{T, E}\n    u    ::SubArray{T,1}   # [γ₀int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0\n    q    ::SubArray{T,1}   # [γ₁int(φ*)](ξ)    ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0\n    w    ::SubArray{T,1}   # [γ₀ext(Ψ)](ξ)     ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0\n    umol ::Vector{T}       # [γ₀int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0\n    qmol ::Vector{T}       # [γ₁int(φ*mol)](ξ) ∀ ξ ∈ Ξ; premultiplied by 4π⋅ε0\nend\n\nResult data of the nonlocal solving process to be used for potential computation and post-processing, with Ξ being the list of observation points, that is, the set of triangle centroids.\n\n\n\n"
},

{
    "location": "lib/solvers.html#NESSie.BEM.solve",
    "page": "Solvers",
    "title": "NESSie.BEM.solve",
    "category": "Function",
    "text": "solve{T, L <: LocalityType}(\n              ::L,\n    model     ::Model{T, Triangle{T}}\n)\n\nComputes the full local or nonlocal cauchy data on the surface of the biomolecule.\n\nReturn type\n\nLocalBEMResult{T, Triangle{T}} or NonlocalBEMResult{T, Triangle{T}}\n\n\n\n"
},

{
    "location": "lib/solvers.html#BEM-solver-1",
    "page": "Solvers",
    "title": "BEM solver",
    "category": "section",
    "text": "    CurrentModule = NESSie.BEM    BEMResult\n    LocalBEMResult\n    NonlocalBEMResult\n    solve"
},

{
    "location": "lib/util.html#",
    "page": "Utility functions",
    "title": "Utility functions",
    "category": "page",
    "text": ""
},

{
    "location": "lib/util.html#NESSie.TestModel.bornion",
    "page": "Utility functions",
    "title": "NESSie.TestModel.bornion",
    "category": "Function",
    "text": "bornion(name::String, ::Type{Float64} = Float64)\nbornion(name::String, ::Type{Float32})\n\nGenerator function for built-in Born ions:\n\nName Charge Radius [Åqv90]\nLi +1 0.645\nNa +1 1.005\nK +1 1.365\nRb +1 1.505\nCs +1 1.715\nMg +2 0.615\nCa +2 1.015\nSr +2 1.195\nBa +2 1.385\n\nReturn type\n\nBornIon\n\n\n\n"
},

{
    "location": "lib/util.html#NESSie.meshunion",
    "page": "Utility functions",
    "title": "NESSie.meshunion",
    "category": "Function",
    "text": "meshunion{T}(\n    model1::Model{T, Tetrahedron{T}},\n    model2::Model{T, Tetrahedron{T}}\n)\n\nMerges two volume models, e.g., the models of a protein and the solvent. Duplicate nodes (e.g., the nodes on the protein surface) are merged into a single node, duplicate elements and charges (if any) are retained as well as the system constants of the first model.\n\nReturn type\n\nModel{T, Tetrahedron{T}}\n\nnote: Note\nThis function assumes that there are no duplicates within either of the node lists!\n\n\n\n"
},

{
    "location": "lib/util.html#NESSie.obspoints_line",
    "page": "Utility functions",
    "title": "NESSie.obspoints_line",
    "category": "Function",
    "text": "obspoints_line{T}(\n    u::Vector{T},\n    v::Vector{T},\n    n::Int\n)\n\nGenerates n evenly distributed observation points along the line segment from u to v.\n\nReturn type\n\nFunction\n\nExample\n\nfor ξ in obspoints_line([0, 0, 0], [1, 1, 1], 10)\n    ...\nend\n\n\n\n"
},

{
    "location": "lib/util.html#NESSie.obspoints_plane",
    "page": "Utility functions",
    "title": "NESSie.obspoints_plane",
    "category": "Function",
    "text": "obspoints_plane{T}(\n    a  ::Vector{T},\n    b  ::Vector{T},\n    c  ::Vector{T},\n    nba::Int,\n    nbc::Int\n)\n\nGenerates nba ⋅ nbc evenly distributed observation points on the parallelogram with the sides overlineBA and overlineBC and a, b, c being the location vectors to the points A, B, and C, respectively.\n\nArguments\n\nnba   Number of observation points along overlineBA\nnbc   Number of observation points along overlineBC\n\nReturn type\n\nFunction\n\nExample\n\nfor Ξ in obspoints_plane(...)\n    for ξ in Ξ\n        ...\n    end\nend\n\n\n\n"
},

{
    "location": "lib/util.html#Utility-functions-1",
    "page": "Utility functions",
    "title": "Utility functions",
    "category": "section",
    "text": "    NESSie.TestModel.bornion\n    meshunion\n    obspoints_line\n    obspoints_plane"
},

{
    "location": "intern/constants.html#",
    "page": "Constants",
    "title": "Constants",
    "category": "page",
    "text": ""
},

{
    "location": "intern/constants.html#NESSie.σ",
    "page": "Constants",
    "title": "NESSie.σ",
    "category": "Constant",
    "text": "Geometric quantity ()\n\n() = lim_0 frac14 _r  r-= d = frac12\n\nfor almost all    (cf. [Ste03]).\n\n\n\n"
},

{
    "location": "intern/constants.html#NESSie.ec",
    "page": "Constants",
    "title": "NESSie.ec",
    "category": "Constant",
    "text": "10^10 times the elementary charge (for   m conversion)\n\nUnit\n\nC\n\n\n\n"
},

{
    "location": "intern/constants.html#NESSie.potprefactor",
    "page": "Constants",
    "title": "NESSie.potprefactor",
    "category": "Function",
    "text": "potprefactor(T::Type{Float64} = Float64)\npotprefactor(T::Type{Float32})\n\nCommon prefactor for all potentials _ and _:\n\nfrac1602  10^-1910^-10  4     1145  4\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "intern/constants.html#NESSie.yukawa",
    "page": "Constants",
    "title": "NESSie.yukawa",
    "category": "Function",
    "text": "yukawa{T}(opt::Option{T})\n\nExponent 1 for the fundamental solution of the yukawa operator\n\n = sqrtfrac__\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "intern/constants.html#int-constants-1",
    "page": "Constants",
    "title": "Constants",
    "category": "section",
    "text": "    CurrentModule = NESSie    σ\n    ec\n    potprefactor\n    yukawa"
},

{
    "location": "intern/input.html#",
    "page": "Input formats",
    "title": "Input formats",
    "category": "page",
    "text": ""
},

{
    "location": "intern/input.html#NESSie.Format.readhmo_nodes",
    "page": "Input formats",
    "title": "NESSie.Format.readhmo_nodes",
    "category": "Function",
    "text": "readhmo_nodes{T <: AbstractFloat}(\n    stream::IOStream,\n          ::Type{T}=Float64\n)\n\nReads all nodes from the given HMO file.\n\nReturn type\n\nVector{Vector{T}}\n\n\n\n"
},

{
    "location": "intern/input.html#NESSie.Format.readhmo_elements",
    "page": "Input formats",
    "title": "NESSie.Format.readhmo_elements",
    "category": "Function",
    "text": "readhmo_elements{T <: AbstractFloat}(\n    stream::IOStream,\n    nodes ::Vector{Vector{T}}\n)\n\nReads all elements from the given HMO file.\n\nReturn type\n\nVector{Triangle{T}}\n\n\n\n"
},

{
    "location": "intern/input.html#NESSie.Format.readhmo_charges",
    "page": "Input formats",
    "title": "NESSie.Format.readhmo_charges",
    "category": "Function",
    "text": "readhmo_charges{T <: AbstractFloat}(\n    stream::IOStream,\n          ::Type{T}=Float64\n)\n\nReads all charges from the given HMO file.\n\nReturn type\n\nVector{Charge{T}}\n\n\n\n"
},

{
    "location": "intern/input.html#NESSie.Format.readmcsf_nodes",
    "page": "Input formats",
    "title": "NESSie.Format.readmcsf_nodes",
    "category": "Function",
    "text": "readmcsf_nodes{T <: AbstractFloat}(\n    stream::IOStream,\n          ::Type{T}=Float64\n)\n\nReads all nodes from the given GAMer-generated mcsf file.\n\nReturn type\n\nVector{Vector{T}}\n\n\n\n"
},

{
    "location": "intern/input.html#NESSie.Format.readmcsf_elements",
    "page": "Input formats",
    "title": "NESSie.Format.readmcsf_elements",
    "category": "Function",
    "text": "readmcsf_elements{T <: AbstractFloat}(\n    stream::IOStream,\n    nodes ::Vector{Vector{T}};\n    # kwargs\n    domain::Symbol=:none\n)\n\nReads all elements from the given GAMer-generated mcsf file.\n\nReturn type\n\nVector{Tetrahedron{T}}\n\n\n\n"
},

{
    "location": "intern/input.html#NESSie.Format.readmsms_nodes",
    "page": "Input formats",
    "title": "NESSie.Format.readmsms_nodes",
    "category": "Function",
    "text": "readmsms_nodes{T <: AbstractFloat}(\n    stream::IOStream,\n          ::Type{T}=Float64\n)\n\nReads all nodes from the given MSMS-generated .vert file.\n\nReturn type\n\nVector{Vector{T}}\n\n\n\n"
},

{
    "location": "intern/input.html#NESSie.Format.readmsms_elements",
    "page": "Input formats",
    "title": "NESSie.Format.readmsms_elements",
    "category": "Function",
    "text": "readmsms_elements{T <: AbstractFloat}(\n    stream::IOStream,\n    nodes ::Vector{Vector{T}}\n)\n\nReads all elements from the given MSMS-generated .face file.\n\nReturn type\n\nVector{Triangle{T}}\n\n\n\n"
},

{
    "location": "intern/input.html#NESSie.Format.readoff_nodes",
    "page": "Input formats",
    "title": "NESSie.Format.readoff_nodes",
    "category": "Function",
    "text": "readoff_nodes{T <: AbstractFloat}(\n    stream::IOStream,\n    n     ::Int,\n          ::Type{T}=Float64\n)\n\nReads the first n nodes from the given OFF file.\n\nReturn type\n\nVector{Vector{T}}\n\n\n\n"
},

{
    "location": "intern/input.html#NESSie.Format.readoff_elements",
    "page": "Input formats",
    "title": "NESSie.Format.readoff_elements",
    "category": "Function",
    "text": "readoff_elements{T <: AbstractFloat}(\n    stream::IOStream,\n    n     ::Int,\n    nodes ::Vector{Vector{T}},\n          ::Type{T}=Float64\n)\n\nReads the first n elements from the given OFF file.\n\nReturn type\n\nVector{Triangle{T}}\n\n\n\n"
},

{
    "location": "intern/input.html#int-input-1",
    "page": "Input formats",
    "title": "Input formats",
    "category": "section",
    "text": "    CurrentModule = NESSie.Format    readhmo_nodes\n    readhmo_elements\n    readhmo_charges\n    readmcsf_nodes\n    readmcsf_elements\n    readmsms_nodes\n    readmsms_elements\n    readoff_nodes\n    readoff_elements"
},

{
    "location": "intern/quadrature.html#",
    "page": "Quadrature",
    "title": "Quadrature",
    "category": "page",
    "text": ""
},

{
    "location": "intern/quadrature.html#int-quadrature-1",
    "page": "Quadrature",
    "title": "Quadrature",
    "category": "section",
    "text": ""
},

{
    "location": "intern/quadrature.html#NESSie.Radon.radoncoll!",
    "page": "Quadrature",
    "title": "NESSie.Radon.radoncoll!",
    "category": "Function",
    "text": "radoncoll!{T}(\n        dest    ::DenseArray{T,1},\n        elements::Vector{Triangle{T}},\n        Ξ       ::Vector{Vector{T}},\n        solution::Function,\n        yukawa  ::T,\n        fvals   ::DenseArray{T,1}\n)\n\nradoncoll!{T}(\n        dest    ::DenseArray{T,2},\n        elements::Vector{Triangle{T}},\n        Ξ       ::Vector{Vector{T}},\n        solution::Function,\n        yukawa  ::T\n)\n\nSeven-point Radon cubature [Rad48] for a given function and a list of triangles and observation points Ξ. If dest is a vector, the function values f for each surface triangle have to be specified, since each element of the vector represents the dot product of the corresponding coefficient matrix row and the fvals vector.\n\nIf you intend computing single/double layer potentials with this function, you might want to use the shorthand signature regularyukawacoll! instead.\n\nnote: Note\nThe result is premultiplied by 4π.\n\nArguments\n\nsolution Fundamental solution; supported functions: regularyukawapot, ∂ₙregularyukawapot\nyukawa Exponent of the Yukawa operator's fundamental solution\n\nReturn type\n\nVoid\n\n\n\n"
},

{
    "location": "intern/quadrature.html#NESSie.Radon.regularyukawapot",
    "page": "Quadrature",
    "title": "NESSie.Radon.regularyukawapot",
    "category": "Function",
    "text": "regularyukawapot{T}(\n    x     ::DenseArray{T,1},\n    ξ     ::Vector{T},\n    yukawa::T,\n          ::Vector{T}=T[]\n)\n\nComputes the regular part of the Yukawa potential, that is, Yukawa minus Laplace:\n\nmathcalG^Y-mathcalG^L = frac14frace^-fracx -  - 1x-\n\nnote: Note\nThe result is premultiplied by 4π.\n\nArguments\n\nyukawa Exponent of the Yukawa operator's fundamental solution\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "intern/quadrature.html#NESSie.Radon.∂ₙregularyukawapot",
    "page": "Quadrature",
    "title": "NESSie.Radon.∂ₙregularyukawapot",
    "category": "Function",
    "text": "∂ₙregularyukawapot{T}(\n    x     ::Vector{T},\n    ξ     ::Vector{T},\n    yukawa::T,\n    normal::Vector{T}\n)\n\nComputes the normal derivative of the regular part of the Yukawa potential, that is, Yukawa minus Laplace:\n\nfracn frac14 frace^-fracx -  - 1x-\n= frac14 frac1 - (1 + ^-1 x - )e^-fracx - x-\nfrac(x - )  nx - \n\nnote: Note\nThe result is premultiplied by 4π.\n\nArguments\n\nyukawa Exponent of the Yukawa operator's fundamental solution\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "intern/quadrature.html#NESSie.Radon.setcubpts!",
    "page": "Quadrature",
    "title": "NESSie.Radon.setcubpts!",
    "category": "Function",
    "text": "setcubpts!{T}(\n    dest::Vector{Vector{T}},\n    qpts::QuadPts2D{T},\n    elem::Triangle{T}\n)\n\nPrepare cubature points for one surface element.\n\nReturn type\n\nVoid\n\n\n\n"
},

{
    "location": "intern/quadrature.html#Radon-cubature-1",
    "page": "Quadrature",
    "title": "Radon cubature",
    "category": "section",
    "text": "    CurrentModule = NESSie.Radon    radoncoll!\n    regularyukawapot\n    ∂ₙregularyukawapot\n    setcubpts!"
},

{
    "location": "intern/quadrature.html#NESSie.Rjasanow.ObservationPosition",
    "page": "Quadrature",
    "title": "NESSie.Rjasanow.ObservationPosition",
    "category": "Type",
    "text": "abstract type ObservationPosition end\nstruct InPlane <: ObservationPosition end\nstryct InSpace <: ObservationPosition end\n\nEnum-like representation of the obseration point's position relative to the corresponding surface element\n\n\n\n"
},

{
    "location": "intern/quadrature.html#NESSie.Rjasanow.InPlane",
    "page": "Quadrature",
    "title": "NESSie.Rjasanow.InPlane",
    "category": "Type",
    "text": "abstract type ObservationPosition end\nstruct InPlane <: ObservationPosition end\nstryct InSpace <: ObservationPosition end\n\nEnum-like representation of the obseration point's position relative to the corresponding surface element\n\n\n\n"
},

{
    "location": "intern/quadrature.html#NESSie.Rjasanow.InSpace",
    "page": "Quadrature",
    "title": "NESSie.Rjasanow.InSpace",
    "category": "Type",
    "text": "abstract type ObservationPosition end\nstruct InPlane <: ObservationPosition end\nstryct InSpace <: ObservationPosition end\n\nEnum-like representation of the obseration point's position relative to the corresponding surface element\n\n\n\n"
},

{
    "location": "intern/quadrature.html#NESSie.Rjasanow.laplacepot",
    "page": "Quadrature",
    "title": "NESSie.Rjasanow.laplacepot",
    "category": "Function",
    "text": "function laplacepot{T, P <: PotentialType}(\n    ptype::Type{P},\n    ξ    ::Vector{T},\n    elem ::Triangle{T},\n    dist ::T\n)\n\nComputes the single or double layer Laplace potential of the given triangle for the given observation point ξ. The latter needs to be projected onto the surface element plane [Rja90].\n\nnote: Note\nThe result is premultiplied by 4π.\n\nArguments\n\ndist Distance from the original ξ to the surface element plane\n\nReturn type\n\nT\n\n\n\nlaplacepot{T, P <: PotentialType}(\n    ptype::Type{P},\n    ξ::Vector{T},\n    x1::Vector{T},\n    x2::Vector{T},\n    normal::Vector{T},\n    dist::T\n)\n\nComputes the single or double layer Laplace potential of the triangle defined by the observation point ξ and two nodes x1 and x2 of the surface element. x2 is required to be x1's next neighbor in counterclockwise direction. Also, ξ needs to be projected onto the surface element plane  [Rja90].\n\nnote: Note\nThe result is premultiplied by 4π.\n\nArguments\n\nnormal Unit normal vector of the surface element\ndist Distance from the original ξ to the surface element plane\n\nReturn type\n\nT\n\n\n\nlaplacepot{T, P <: PotentialType, O <: ObservationPosition}(\n         ::Type{P},\n         ::Type{O},\n    sinφ1::T,\n    sinφ2::T,\n    h    ::T,\n    d    ::T\n)\n\nComputes the Laplace potential (or its normal derivative) of the triangle with the given height h at the observation point ξ (projected onto the surface element plane) and the sines of the angles  and  between h and the triangle sides extending from ξ [Rja90].\n\nnote: Note\nThe result is premultiplied by 4π.\n\nArguments\n\nd Distance from the original ξ to the surface element plane\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "intern/quadrature.html#NESSie.Rjasanow.logterm",
    "page": "Quadrature",
    "title": "NESSie.Rjasanow.logterm",
    "category": "Function",
    "text": "logterm{T}(χ2::T, sinφ::T)\n\nUtility function to compute\n\nfrac\nsqrt1 -  sin() + sqrt1 -  sin()\nsqrt1 -  sin() - sqrt1 -  sin()\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "intern/quadrature.html#NESSie.Rjasanow.projectξ",
    "page": "Quadrature",
    "title": "NESSie.Rjasanow.projectξ",
    "category": "Function",
    "text": "projectξ{T}(ξ::Vector{T}, elem::Triangle{T})\n\nProjects ξ onto the surface element plane. Returns the projection and the original distance to the plane.\n\nReturn type\n\nTuple{Vector{T}, T}\n\n\n\n"
},

{
    "location": "intern/quadrature.html#Rjasanow-1",
    "page": "Quadrature",
    "title": "Rjasanow",
    "category": "section",
    "text": "    CurrentModule = NESSie.Rjasanow    ObservationPosition\n    InPlane\n    InSpace\n    laplacepot\n    logterm\n    projectξ"
},

{
    "location": "intern/util.html#",
    "page": "Utility functions",
    "title": "Utility functions",
    "category": "page",
    "text": ""
},

{
    "location": "intern/util.html#NESSie.cathetus",
    "page": "Utility functions",
    "title": "NESSie.cathetus",
    "category": "Function",
    "text": "cathetus{T}(hyp::T, cosθ::T)\n\nComputes the cathetus c of a triangle given the hypotenuse h and the cosine of the exterior angle  between the hypotenuse and the other cathetus c.\n\nc = h  cos() \nh = c + c \n c = sqrth  (1 - cos())\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.TestModel.coefficients",
    "page": "Utility functions",
    "title": "NESSie.TestModel.coefficients",
    "category": "Function",
    "text": "function coefficients{T}(\n    model::XieModel{T},\n    len  ::Int\n)\n\nComputes the coefficients A_in with i=1 2 3 for the given XieModel and the desired number of terms.\n\nReturn type\n\nTuple{     Array{T, 2},     Array{T, 2},     Array{T, 2} }\n\n\n\n"
},

{
    "location": "intern/util.html#Base.cos",
    "page": "Utility functions",
    "title": "Base.cos",
    "category": "Function",
    "text": "cos{T}(\n    u    ::Vector{T},\n    v    ::Vector{T},\n    unorm::T=vecnorm(u),\n    vnorm::T=vecnorm(v)\n)\n\nComputes the cosine of the angle between the given vectors u and v with lengths unorm and vnorm, respectively.\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.ddot",
    "page": "Utility functions",
    "title": "NESSie.ddot",
    "category": "Function",
    "text": "ddot{T}(\n    u::Vector{T},\n    v::Vector{T},\n    n::Vector{T}\n)\n\nDevectorized computation of (u-v)⋅n.\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.distance",
    "page": "Utility functions",
    "title": "NESSie.distance",
    "category": "Function",
    "text": "distance{T}(\n    q   ::Vector{T},\n    elem::Triangle{T}\n)\n\nCalculates the (positive or negative) distance from the given point q to the plane the given triangle elem is located in.\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.eye!",
    "page": "Utility functions",
    "title": "NESSie.eye!",
    "category": "Function",
    "text": "eye!{T}(\n    m::Union{DenseArray{T,2}, SubArray{T,2}},\n    α::Number=one(T)\n)\n\nInitializes the given matrix m with αI, with I being an identity matrix with the same dimensions as m.\n\nReturn type\n\nVoid\n\nExample\n\njulia> m = 2 * ones(2, 2)\n2×2 Array{Float64,2}:\n 2.0  2.0\n 2.0  2.0\n\njulia> eye!(m); m\n2×2 Array{Float64,2}:\n 1.0  0.0\n 0.0  1.0\n\njulia> eye!(m, 2); m\n2×2 Array{Float64,2}:\n 2.0  0.0\n 0.0  2.0\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.isdegenerate",
    "page": "Utility functions",
    "title": "NESSie.isdegenerate",
    "category": "Function",
    "text": "isdegenerate{T}(elem::Triangle{T})\n\nTests whether the given triangle is degenerate.\n\nReturn type\n\nBool\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.TestModel.legendre",
    "page": "Utility functions",
    "title": "NESSie.TestModel.legendre",
    "category": "Function",
    "text": "legendre{T <: AbstractFloat}(\n    maxn::Int,\n    x   ::T\n)\n\nPrecomputes the Legendre polynomials P(x) for n = 0 1  textttmaxn-1 and returns a generic function to access the values.\n\nReturn type\n\n(generic function)\n\n(n::Int) -> Pₙ(x)   # return type: T\n\nExample\n\njulia> p = legendre(3, 0.5);\n\njulia> [p(n) for n in 0:2]    # [P₀(0.5), P₁(0.5), P₂(0.5)]\n3-element Array{Float64,1}:\n  1.0\n  0.5\n -0.125\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.pluseye!",
    "page": "Utility functions",
    "title": "NESSie.pluseye!",
    "category": "Function",
    "text": "pluseye!{T}(\n    m::Union{DenseArray{T,2}, SubArray{T,2}},\n    α::Number=one(T)\n)\n\nAdds α to all diagonal elements of matrix m.\n\nReturn type\n\nVoid\n\nExample\n\njulia> m = 2 * ones(2, 2)\n2×2 Array{Float64,2}:\n 2.0  2.0\n 2.0  2.0\n\njulia> pluseye!(m); m\n2×2 Array{Float64,2}:\n 3.0  2.0\n 2.0  3.0\n\njulia> pluseye!(m, 2); m\n2×2 Array{Float64,2}:\n 5.0  2.0\n 2.0  5.0\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.props",
    "page": "Utility functions",
    "title": "NESSie.props",
    "category": "Function",
    "text": "props{T}(\n    elem::Triangle{T}\n)\n\nComputes the given triangle's properties, that is, centroid, normal, distance to origin, and area. Returns the completely initialized Triangle as a copy.\n\nwarning: Warning\nThe given triangle remains unchanged!\n\nReturn type\n\nTriangle\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.reverseindex",
    "page": "Utility functions",
    "title": "NESSie.reverseindex",
    "category": "Function",
    "text": "reverseindex{T}(v::Vector{T})\n\nCreates a reverse index for the given vector v, that is, a dictionary linking the object IDs of the vector elements to the corresponding position in the vector.\n\nReturn type\n\nDict{UInt, UInt}\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.TestModel.scalemodel",
    "page": "Utility functions",
    "title": "NESSie.TestModel.scalemodel",
    "category": "Function",
    "text": "function scalemodel{T}(\n    charges::Vector{Charge{T}},\n    radius ::T;\n    # kwargs\n    compat ::Bool              = false\n)\n\nTranslates and rescales the given point charge model to fit inside an origin-centered sphere with the specified radius. More specifically, the function will center the given model at the origin and scaled in a way such that the outermost point charge will be located at 80% radius distance from the origin.\n\nArguments\n\ncompat Enables compatibility mode and scales the model exactly like the reference implementation ([Xie16]). Use this flag if you intend to compare the results to the reference.\n\nReturn type\n\nVector{Charge{T}}\n\n\n\n"
},

{
    "location": "intern/util.html#Base.seek",
    "page": "Utility functions",
    "title": "Base.seek",
    "category": "Function",
    "text": "seek(\n    fh         ::IOStream,\n    prefix     ::String,\n    skiptheline::Bool = true\n)\n\nFast-forwards an IOStream to the next line starting with the given prefix. In case there is no such line, the stream handle will be set to EOF.\n\nArguments\n\nskiptheline    If true, said line will also be skipped\n\nReturn type\n\nVoid\n\n\n\n"
},

{
    "location": "intern/util.html#Base.sign",
    "page": "Utility functions",
    "title": "Base.sign",
    "category": "Function",
    "text": "sign{T}(\n    u::Vector{T},\n    v::Vector{T},\n    n::Vector{T}\n)\n\nDetermines whether the normal vector of the plane specified by the vectors u and v has the same orientation as the given normal vector n. Returns 1 if both normals have the same orientation, 0 if at least one of the vectors is zero, and -1 otherwise.\n\nReturn type\n\nT\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.TestModel.spherical_besseli",
    "page": "Utility functions",
    "title": "NESSie.TestModel.spherical_besseli",
    "category": "Function",
    "text": "spherical_besseli{T <: AbstractFloat}(\n    maxn::Int,\n    r   ::T\n)\n\nPrecomputes the modified spherical Bessel function of the first kind i(r) for n=-1 0  textttmaxn and returns a generic function to access the values. i(r) is defined as\n\ni(r) = sqrtfrac2r I_n+05(r)\n\nwhere I_ is the modified Bessel function of the first kind [Xie16].\n\nReturn type\n\n(generic function)\n\n(n::Int) -> iₙ(r)   # return type: T\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.TestModel.spherical_besselk",
    "page": "Utility functions",
    "title": "NESSie.TestModel.spherical_besselk",
    "category": "Function",
    "text": "spherical_besselk{T <: AbstractFloat}(\n    maxn::Int,\n    r   ::T\n)\n\nPrecomputes the modified spherical Bessel function of the second kind k(r) for n=-1 0  textttmaxn and returns a generic function to access the values. k(r) is defined as\n\nk(r) = sqrtfrac2r K_n+05(r)\n\nwhere K_ is the modified Bessel function of the second kind [Xie16].\n\nReturn type\n\n(generic function)\n\n(n::Int) -> kₙ(r)   # return type: T\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.unpack",
    "page": "Utility functions",
    "title": "NESSie.unpack",
    "category": "Function",
    "text": "unpack{T}(data::Vector{Vector{T}})\n\nUnpacks the given vector of vectors into a single vector.\n\nReturn type\n\nVector{T}\n\nExample\n\njulia> unpack([[1, 2], [3]])\n3-element Array{Int64,1}:\n 1\n 2\n 3\n\n\n\n"
},

{
    "location": "intern/util.html#NESSie.vertexnormals",
    "page": "Utility functions",
    "title": "NESSie.vertexnormals",
    "category": "Function",
    "text": "vertexnormals{T}(model::Model{T, Triangle{T}})\n\nReturns a vector containing the normal vectors of the given model's triangles.\n\nReturn type\n\nVector{Vector{T}}\n\n\n\n"
},

{
    "location": "intern/util.html#int-util-1",
    "page": "Utility functions",
    "title": "Utility functions",
    "category": "section",
    "text": "    CurrentModule = NESSie\n    DocTestSetup = quote\n        using NESSie: eye!, pluseye!, unpack\n        using NESSie.TestModel:legendre\n    end    cathetus\n    NESSie.TestModel.coefficients\n    cos\n    ddot\n    distance\n    eye!\n    isdegenerate\n    NESSie.TestModel.legendre\n    pluseye!\n    props\n    reverseindex\n    NESSie.TestModel.scalemodel\n    seek\n    sign\n    NESSie.TestModel.spherical_besseli\n    NESSie.TestModel.spherical_besselk\n    unpack\n    vertexnormals    DocTestSetup = nothing"
},

]}
