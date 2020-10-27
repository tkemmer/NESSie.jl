using NESSie.Format

@testset "readstl" begin
    @test_skip readstl
end

@testset "writestl" begin
    for T in testtypes
        # empty system
        model = Model{T, Triangle{T}}()
        (fname, fh) = mktemp()
        try
            writestl(fh, model)
            seek(fh, 0)
            @test [read(fh, UInt8) for _ in 1:80] ≈ zeros(Float32, 80)
            @test read(fh, UInt32) == 0
            @test_throws EOFError read(fh, UInt8)
        finally
            close(fh)
            rm(fname)
        end

        # small system
        nodes = Vector{T}[[0, 0, 0], [0, 0, 3], [0, 3, 0], [1, -3, 3]]
        model = Model(
            nodes,
            Triangle{T}[
                Triangle(nodes[1], nodes[2], nodes[3]),
                Triangle(nodes[1], nodes[4], nodes[2])
            ]
        )
        (fname, fh) = mktemp()
        try
            writestl(fh, model)
            seek(fh, 0)
            @test [read(fh, UInt8) for _ in 1:80] ≈ zeros(Float32, 80)
            @test read(fh, UInt32) == 2
            for elem in model.elements
                @test [read(fh, Float32) for _ in 1:3] ≈ elem.normal
                @test [read(fh, Float32) for _ in 1:3] ≈ elem.v1
                @test [read(fh, Float32) for _ in 1:3] ≈ elem.v2
                @test [read(fh, Float32) for _ in 1:3] ≈ elem.v3
                @test read(fh, UInt16) == zero(UInt16)
            end
            @test_throws EOFError read(fh, UInt8)
        finally
            close(fh)
            rm(fname)
        end
    end
end
