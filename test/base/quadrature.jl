context("quadraturepoints") do
    for T in testtypes
        @fact typeof(quadraturepoints(Triangle{T})) --> QuadPts2D{T}
        @fact typeof(quadraturepoints(Tetrahedron{T})) --> QuadPts3D{T}
    end
end
