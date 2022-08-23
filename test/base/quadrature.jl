@testset "quadraturepoints" begin
    for T in testtypes
        # simple generators
        @test typeof(quadraturepoints(Triangle{T})) == QuadPts2D{T}
        @test typeof(quadraturepoints(Tetrahedron{T})) == QuadPts3D{T}

        # on-triangle generators
        elements = Triangle{T}[]
        ret = quadraturepoints(elements)
        @test typeof(ret) == Vector{TriangleQuad{T}}
        @test isempty(ret)

        push!(elements, Triangle(T[0, 0, 0], T[0, 1, 0], T[0, 0, 1]))
        qpts = quadraturepoints(Triangle{T})
        ret  = quadraturepoints(elements)
        @test typeof(ret) == Vector{TriangleQuad{T}}
        @test length(ret) == 1
        @test ret[1].elem === elements[1]
        @test ret[1].weights === qpts.weight
        @test size(ret[1].qpts) == (3, qpts.num)
        for i in 1:qpts.num
            @test ret[1].qpts[:, i] ≈ T[0, qpts.x[i], qpts.y[i]]
        end

        push!(elements, Triangle(T[1, 1, 1], T[1, 1, 3], T[1, 3, 1]))
        ret  = quadraturepoints(elements)
        @test typeof(ret) == Vector{TriangleQuad{T}}
        @test length(ret) == 2
        @test ret[2].elem === elements[2]
        @test ret[2].weights === qpts.weight
        @test size(ret[2].qpts) == (3, qpts.num)
        for i in 1:qpts.num
            @test ret[2].qpts[:, i] ≈ T[1, 1 + 2 * qpts.y[i], 1 + 2 * qpts.x[i]]
        end
    end
end
