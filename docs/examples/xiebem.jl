push!(LOAD_PATH,"../../src/")

using NESSie
using NESSie.BEM
using NESSie.Format: readoff, readpqr
using NESSie.TestModel

using LinearAlgebra: norm
using PyPlot: cm_get_cmap, figure, plot_trisurf, show, title, xlabel, ylabel, zlabel, zlim

#=
    xiebem.jl

    Comparison of nonlocal reaction field potentials (or electrostatic potentials) for a
    Poisson test model generated from a given PQR file. The potentials are computed using
    the NESSie.TestModel and NESSie.BEM modules.
=#
if length(ARGS) > 1 || startswith(get(ARGS, 1, ""), "-")
    println("\n\e[1mUsage\e[0m: julia xiebem.jl [PQR_FILE]")
    exit(1)
end


const εΩ = 2.                           # dielectric constant (protein)
const εΣ = 80.                          # dielectric constant (solvent)
const ε∞ = 1.8                          # dielectric bulk response (solvent)
const λ  = 15.                          # correlation length (Å)
const radius = 1.                       # test model sphere radius (Å)
const numiter = 20                      # number of iterations (test model)
const reac = true                       # generate reaction field only
const lval = 2.                         # grid size [-l, l]² in xz-plane
const grid = 40                         # grid resolution (grid x grid)
const fpqr = "../../data/xie/2LZX.pqr"  # PQR file (ignored when given as argument)

function plotpotential(
        x   ::Vector{T},
        y   ::Vector{T},
        z   ::Vector{T},
        name::String = ""
    ) where T
    figure()
    plot_trisurf(x, y, z, cmap=cm_get_cmap("viridis"))
    length(name) > 0 && (title(name))
    xlabel("x component [\$\\AA\$]")
    ylabel("z component [\$\\AA\$]")
    zlabel(reac ? "Reaction field potential [V]" : "Electrostatic potential [V]")
    nothing
end

# generate observation points
Ξ = [ξ for Ξ in obspoints_plane(
        [-lval, zero(lval),  lval],
        [-lval, zero(lval), -lval],
        [ lval, zero(lval), -lval],
        grid, grid
    ) for ξ in Ξ]

# x- and y-axes of the plots
x = [ξ[1] for ξ in Ξ]
y = [ξ[3] for ξ in Ξ]

#=
    Nonlocal Xie model 1
=#
nlxie1 = NonlocalXieModel1(
            TestModel.XieModel(
                radius,
                readpqr(get(ARGS, 1, fpqr)),
                Option(εΩ, εΣ, ε∞, λ),
                compat=true),
            numiter
        )
molpot  = φmol(Ξ, nlxie1.charges)/4π/εΩ/ε0*NESSie.ec # molecular potential

z = [
        norm(ξ) < nlxie1.radius ?
        TestModel.φΩ(ξ, nlxie1) :
        TestModel.φΣ(ξ, nlxie1) for ξ in Ξ
    ]
reac && (z -= molpot)

zmin = (any(x -> x < 0, z) ? 1.2 : .8) * minimum(z)
zmax = (any(x -> x > 0, z) ? 1.2 : .8) * maximum(z)
plotpotential(x, y, z, "Nonlocal Poisson Test Model")

#=
    BEM
=#
surf = readoff("../../data/xie/unitsphere.off")
surf.charges = nlxie1.charges
surf.params  = nlxie1.params
bem = solve(NonlocalES, surf)

zΣ = BEM.φΣ(Ξ, bem)
zΩ = BEM.φΩ(Ξ, bem)
z  = [norm(ξ) < nlxie1.radius ? zΩ[i] : zΣ[i] for (i, ξ) in enumerate(Ξ)]
reac && (z -= molpot)

zmin = min(zmin, (any(x -> x < 0, z) ? 1.2 : .8) * minimum(z))
zmax = max(zmax, (any(x -> x > 0, z) ? 1.2 : .8) * maximum(z))
zlim((zmin, zmax))
plotpotential(x, y, z, "Nonlocal BEM")
zlim((zmin, zmax))

show()
