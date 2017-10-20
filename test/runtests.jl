push!(LOAD_PATH,"../src/")

using NESSie
using FactCheck

# All tests. Will be used if no command line arguments are given.
const tests = [
    "base/constants.jl",
    "base/model.jl",
    "base/potentials.jl",
    "base/quadrature.jl",
    "base/util.jl",
    "bem/local.jl",
    "bem/nonlocal.jl",
    "format/hmo.jl",
    "format/mcsf.jl",
    "format/msms.jl",
    "format/off.jl",
    "format/pqr.jl",
    "format/skel.jl",
    "format/vtk.jl",
    "format/xml3d.jl",
    "testmodel/born.jl",
    "radon.jl",
    "rjasanow.jl"
]
const testtypes = (Float64, Float32)

function printusage()
    println("\n\e[1mNESSie.jl test suite\e[0m")
    println("====================\n")
    println("\e[1mUsage:\e[0m\n")
    println("\truntests.jl [--help] [<test>...]\n")
    println("If called without arguments, this program will run all available tests.")
    println("\n\e[1mOptions:\e[0m\n")
    println("\t\e[1m--help\e[0m\tPrint this message and exit")
    println("\n\e[1mAvailable tests:\e[0m\n")
    for test in tests
        println("\t$(test[1:end-3])")
    end
    println()
end

#=
    Validate args and return test files
=#
function checkargs(args)
    ret = String[]
    err = String[]
    for arg in ARGS
        if arg == "--help"
            printusage()
            exit(0)
        end
        if arg[1] == '-'
            push!(err, "\e[1;31mInvalid option \e[39m$arg\e[0m")
            continue
        end
        if arg in tests || "$(arg).jl" in tests
            push!(ret, endswith(arg, ".jl") ? arg : "$(arg).jl")
        else
            push!(err, "\e[1;31mInvalid test file \e[39m$arg\e[0m")
        end
    end
    if length(err) > 0
        println("\e[1;31mERROR:\e[0m\n")
        for e in err
            println(" * $e")
        end
        print("\n\e[1;31mRun\e[39m runtests.jl --help")
        println(" \e[31mfor a list of available options and tests.\e[0m")
        exit(1)
    end
    ret
end

# Run tests
include("testutils.jl")
for t in (length(ARGS) > 0 ? checkargs(ARGS) : tests)
    facts(t) do
        include(t)
    end
    println()
end
FactCheck.exitstatus()
