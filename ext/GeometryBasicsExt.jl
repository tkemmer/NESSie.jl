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

end # module
