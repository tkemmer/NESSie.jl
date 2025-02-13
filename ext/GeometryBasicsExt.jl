module GeometryBasicsExt

import GeometryBasics

using NESSie

"""
    GeometryBasics.mesh(model::Model{T, Triangle{T}})

Converts the given model into a GeometryBasics.jl-compatible mesh, e.g., for visualization
through [Makie.jl](https://docs.makie.org/stable/reference/plots/mesh).

# Return type
[`GeometryBasics.Mesh`]
(https://juliageometry.github.io/GeometryBasics.jl/stable/meshes/#GeometryBasics.Mesh-meshes)
"""
function GeometryBasics.mesh(model::Model{T, Triangle{T}}) where T
    ridx = NESSie._reverseindex(model.nodes)

    points = GeometryBasics.Point{3, T}.(model.nodes)
    faces = [
        GeometryBasics.TriangleFace(ridx[elem.v1], ridx[elem.v2], ridx[elem.v3])
        for elem in model.elements
    ]

    GeometryBasics.Mesh(points, faces)
end

"""
    function NESSie.Model(
        mesh::GeometryBasics.Mesh{3, T, <: GeometryBasics.NgonFace{3}};
        charges::Vector{Charge{T}} = Charge{T}[],
        params::Option{T} = defaultopt(T)
    )

Converts the given GeometryBasics.jl mesh into a triangle-based model.

# Return type
[`Model{T, Triangle{T}}`](@ref)
"""
function NESSie.Model(
    mesh::GeometryBasics.Mesh{3, T, <: GeometryBasics.NgonFace{3}};
    charges::Vector{Charge{T}} = Charge{T}[],
    params::Option{T} = defaultopt(T)
) where T
    model = Model{T, Triangle{T}}()
    model.nodes = mesh.position
    model.elements = collect(
        Triangle{T},
        Triangle(model.nodes[f[1]], model.nodes[f[2]], model.nodes[f[3]]) for f in mesh.faces
    )
    model.charges = charges
    model.params = params
    model
end

end # module
