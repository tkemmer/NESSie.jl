using NESSie: potprefactor

@testset "potprefactor" begin
    for T in testtypes
        @test typeof(potprefactor(T)) == T
    end
end

@testset "defaultopt" begin
    for T in testtypes
        @test typeof(defaultopt(T)) == Option{T}
    end
end
