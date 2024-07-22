@testitem "Constants" begin
    include("../testsetup.jl")

    @testset "potprefactor" begin
        using NESSie: potprefactor

        for T in testtypes
            @test typeof(potprefactor(T)) == T
        end
    end

    @testset "defaultopt" begin
        for T in testtypes
            @test typeof(defaultopt(T)) == Option{T}
        end
    end
end
