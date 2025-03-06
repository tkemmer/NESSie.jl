@testitem "Test model: Born" begin
    include("../testsetup.jl")

    using NESSie.TestModel

    const ion_names = [
        "Li", "li", "Na", "na", "K", "k", "Rb", "rb", "Cs", "cs",
        "Mg", "mg", "Ca", "ca", "Sr", "sr", "Ba", "ba"
    ]

    @testset "defaultopt" begin
        for T in testtypes
            @test typeof(defaultopt(BornIon{T})) == Option{T}
        end
    end

    @testset "bornion" begin
        for T in testtypes, name in ion_names
            @test typeof(bornion(name, T)) == BornIon{T}
        end
    end

    @testset "Model" begin
        for T in testtypes
            let model = Model(bornion(first(ion_names), T))
                @test model isa Model{T, Triangle{T}}
                @test !isempty(model.nodes)
                @test !isempty(model.elements)
                @test !isempty(model.charges)
                @test model.params == defaultopt(BornIon{T})
            end
        end
    end

    @test_skip φΩ(LocalES)
    @test_skip φΣ(LocalES)
    @test_skip φΩ(NonlocalES)
    @test_skip φΣ(NonlocalES)
end
