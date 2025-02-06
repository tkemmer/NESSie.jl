module GeometryBasicsExt

import GeometryBasics

using NESSie

function GeometryBasics.mesh(model::Model{T, Triangle{T}}) where T
    ridx = IdDict(n => i for (i, n) in enumerate(model.nodes))

    points = GeometryBasics.Point{3, T}.(model.nodes)
    faces = [
        GeometryBasics.TriangleFace(ridx[elem.v1], ridx[elem.v2], ridx[elem.v3])
        for elem in model.elements
    ]

    GeometryBasics.Mesh(points, faces)
end

end # module
