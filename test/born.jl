using ProteinES.Born

context("bornion") do
    for T in testtypes
        @fact typeof(bornion(T, "Na")) --> BornIon{T}
    end
end

@pending φΩ(LocalES) --> :nothing
@pending φΣ(LocalES) --> :nothing
@pending φΩ(NonlocalES) --> :nothing
@pending φΣ(NonlocalES) --> :nothing
