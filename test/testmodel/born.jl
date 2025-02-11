@testitem "Test model: Born" begin
    include("../testsetup.jl")

    using NESSie.TestModel

    @testset "defaultopt" begin
        for T in testtypes
            @test typeof(defaultopt(BornIon{T})) == Option{T}
        end
    end

    @testset "bornion" begin
        for T in testtypes
            @test typeof(bornion("Na", T)) == BornIon{T}
        end
    end

    @test_skip φΩ(LocalES)
    @test_skip φΣ(LocalES)
    @test_skip φΩ(NonlocalES)
    @test_skip φΣ(NonlocalES)
end
