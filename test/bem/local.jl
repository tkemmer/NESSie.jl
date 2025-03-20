@testitem "Local BEM" begin
    include("../testsetup.jl")

    @test_skip solve
    @test_skip espotential
    @test_skip rfenergy
end
