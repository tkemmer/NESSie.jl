@testitem "Nonlocal BEM" begin
    include("../testsetup.jl")

    @test_skip solve
    @test_skip φΩ
    @test_skip φΣ
    @test_skip φΓ
    @test_skip rfenergy
end
