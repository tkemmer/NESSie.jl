@testitem "Format: Wavefront/.obj" begin
    include("../testsetup.jl")

    @testset "writeobj" begin
        @test_skip writeobj
    end
end
