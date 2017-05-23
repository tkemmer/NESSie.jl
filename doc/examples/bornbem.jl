push!(LOAD_PATH,"../../src/")

#=
    bornbem.jl

    Comparison of local and nonlocal exterior potentials of Born ions.
=#

using ProteinES
using ProteinES.IO
using ProteinES.BEM
using ProteinES.Born
using PyPlot: figure, plot, show, legend, xlabel, ylabel, title, axvline

function plotexterior(name::String, maxdist::Float64=10., resolution::Int=100)
    name          = lowercase(name)
    born          = bornion(Float64, name)
    model         = readoff("../../data/born/$name.off")
    model.charges = readpqr("../../data/born/$name.pqr")
    model.params  = Option(1., 78., 1., 23.)
    Ξ             = collect(obspoints_line(
                        [0., 0., born.radius],
                        [0., 0., maxdist],
                        resolution)
                    )
    figure()
    title("$name ion (radius=$(born.radius) Å)")
    axvline(born.radius, color=(.5, .5, .5), linestyle="--", linewidth=3)

    xvals = [ξ[3] for ξ in Ξ]
    print("Local Born... ")
    plot(
        xvals,
        [Born.φΣ(LocalES, ξ, born, model.params) for ξ in Ξ],
        label="Born (local)"
    )

    print("done!\nNonlocal Born... ")
    plot(
        xvals,
        [Born.φΣ(NonlocalES, ξ, born, model.params) for ξ in Ξ],
        label="Born (nonlocal)"
    )

    print("done!\nLocal BEM... ")
    plot(
        xvals,
        BEM.φΣ(Ξ, solve(LocalES, model)),
        label="BEM (local)"
    )

    print("done!\nNonlocal BEM... ")
    plot(
        xvals,
        BEM.φΣ(Ξ, solve(NonlocalES, model)),
        label="BEM (nonlocal)"
    )
    println("done!")

    legend()
    xlabel("Distance from point charge [Å]")
    ylabel("Exterior potential [V]")
    show()
end

if length(ARGS) < 1
    println("\n\e[1mUsage\e[0m: julia bornbem.jl ION_NAME")
    println("\nwhere ION_NAME is a valid Born ion name (see documentation).")
    exit(1)
end

plotexterior(ARGS[1])
