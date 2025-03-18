#=
    bornbem.jl

    Comparison of local and nonlocal potentials of Born ions.
=#

using NESSie
using NESSie.BEM, NESSie.TestModel
using Plots

function plotpot(name::String, maxdist::Float64=10., resolution::Int=100)
    name  = titlecase(name)
    born  = bornion(name)
    model = Model(born)
    Ξ     = LinRange([0., 0., born.radius / 2], [0., 0., maxdist], resolution)

    plot_x = getindex.(Ξ, 3)
    p = vline([born.radius]; color = :black, linewidth = 1, linestyle = :solid, label = nothing)

    print("Local Born... ")
    plot!(p, plot_x, espotential(LocalES, Ξ, born); label="Born (local)",
        title  = "$name ion (radius=$(born.radius) Å)",
        xlabel = "Distance from point charge in Å",
        ylabel = "Electrostatic potential in V",
        yscale = :log10
    )

    print("done!\nNonlocal Born... ")
    plot!(p, plot_x, espotential(NonlocalES, Ξ, born);
        label = "Born (nonlocal)"
    )

    print("done!\nLocal BEM... ")
    plot!(p, plot_x, espotential(Ξ, solve(LocalES, model; method = :blas));
        label = "BEM (local)",
        marker = :circle,
        markersize = 2,
        linetype = :scatter
    )

    print("done!\nNonlocal BEM... ")
    plot!(p, plot_x, espotential(Ξ, solve(NonlocalES, model; method = :blas));
        label = "BEM (nonlocal)",
        marker = :plus,
        markersize = 2,
        markercolor = :black,
        linetype = :scatter
    )

    println("done!")
    display(p)
end

if length(ARGS) < 1
    println("\n\e[1mUsage\e[0m: julia -i bornbem.jl ION_NAME")
    println("\nwhere ION_NAME is a valid Born ion name (see documentation).")
    exit(1)
end

plotpot(ARGS[1])
