using ProteinES: potprefactor

context("potprefactor") do
    for T in testtypes
        @fact typeof(potprefactor(T)) --> T
    end
end

context("defaultopt") do
    for T in testtypes
        @fact typeof(defaultopt(T)) --> Option{T}
    end
end
