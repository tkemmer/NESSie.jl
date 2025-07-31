using NESSie
using NESSie.BEM, NESSie.TestModel
using Plots

#=
    xiebem.jl

    Comparison of nonlocal electrostatic and reaction field potentials for a Poisson test model
    generated from a given PQR file. The potentials are computed using the NESSie.TestModel and
    NESSie.BEM modules.
=#
if length(ARGS) > 1 || startswith(get(ARGS, 1, ""), "-")
    println("\n\e[1mUsage\e[0m: julia -i xiebem.jl [PQR_FILE]")
    exit(1)
end

const εΩ = 2.               # dielectric constant (protein)
const εΣ = 80.              # dielectric constant (solvent)
const ε∞ = 1.8              # dielectric bulk response (solvent)
const λ  = 15.              # correlation length (Å)
const radius = 1.           # test model sphere radius (Å)
const numiter = 20          # number of iterations (test model)
const lval = 2.             # grid size [-l, l]² in xz-plane
const grid = 40             # grid resolution (grid x grid)
const fpqr = "xie/2LZX.pqr" # PQR file (ignored when given as argument)

@inline function plotpotential(
    x   ::Vector{T},
    y   ::Vector{T},
    z   ::Vector{T};
    name::String = "",
    reac::Bool = false
) where T
    surface(x, y, z;
        cmap = :viridis,
        title = name,
        xlabel = "x in Å",
        ylabel = "y in Å",
        zlabel = reac ? "Reaction field potential in V" : "Electrostatic potential in V",
        colorbar = false
    )
end

# generate observation points
Ξ = collect(Iterators.flatten(obspoints_plane(
    [-lval,  lval, zero(lval)],
    [-lval, -lval, zero(lval)],
    [ lval, -lval, zero(lval)],
    grid, grid
)))
plot_x = getindex.(Ξ, 1)
plot_y = getindex.(Ξ, 2)

# Nonlocal Xie model 2
xie = TestModel.XieModel(
    radius,
    Format.readpqr(get(ARGS, 1, nessie_data_path(fpqr))),
    Option(εΩ, εΣ, ε∞, λ),
    compat=true
)
nlxie2 = NonlocalXieModel2(xie, numiter)
p1 = plotpotential(plot_x, plot_y, espotential(Ξ, nlxie2); name = "Nonlocal Poisson Test Model")
p2 = plotpotential(plot_x, plot_y, rfpotential(Ξ, nlxie2); reac = true)

# BEM
surf = Model(xie)

bem = solve(NonlocalES, surf; method = :blas)
p3 = plotpotential(plot_x, plot_y, espotential(Ξ, bem); name = "Nonlocal BEM")
p4 = plotpotential(plot_x, plot_y, rfpotential(Ξ, bem); reac = true)

# synchronize z-axes
zl13 = extrema([zlims(p1)..., zlims(p3)...])
zl24 = extrema([zlims(p2)..., zlims(p4)...])
zlims!(p1, zl13); zlims!(p3, zl13)
zlims!(p2, zl24); zlims!(p4, zl24)

display(plot(p1, p3, p2, p4; layout = (2, 2), size = (1000, 1000)))
