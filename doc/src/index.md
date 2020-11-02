# NESSie.jl
Nonlocal Electrostatics in Structured Solvents


## Motivation
Electrostatic interactions are a major contributor to protein-protein and protein-ligand interactions. In contrast to other molecular interaction components, they can be significant over medium to long distances and are thus crucial for molecular visibility. One major challenge in this context is the treatment of the solvent the molecules are immersed in, e.g., water in a biological context. Strong simplifications of the structure of such polarizable and highly structured solvents are commonplace to achieve the required computational efficiency, but invariably lead to inaccuracies.


## Usage example
The following Julia code shows how to compute and print the nonlocal reaction field energy
of a single Na+ ion (modelled as a sperical vacuum-filled system) in water:

```julia
using NESSie
using NESSie.BEM
using NESSie.Format: readoff, readpqr

# I. Create model
model           = readoff("data/born/na.off")
model.charges   = readpqr("data/born/na.pqr")
model.params.εΩ = 1   # dielectric constant for vacuum model
model.params.εΣ = 78  # dielectric constant for water

# II. Apply nonlocal solver
bem = solve(NonlocalES, model)

# III. Apply postprocessor
val = rfenergy(bem)
println("Reaction field energy: $val kJ/mol")
```
More examples are available in the `doc/examples/` directory.


## Citing
If you use `NESSie.jl` in your research, please cite the following publication:
> Kemmer, T, Rjasanow, S., Hildebrandt, A (2018). NESSie. jl - Efficient and Intuitive
> Finite Element and Boundary Element Methods for Nonlocal Protein Electrostatics in the
> Julia Language. Journal of Computational Science 28, 193-203.

With BibTeX, you can use the following entry:
```
@article{nessie-2018,
    author = {Kemmer, Thomas and Rjasanow, Sergej and Hildebrandt, Andreas},
    title = {{NESSie.jl -- Efficient and Intuitive Finite Element and Boundary Element Methods for Nonlocal Protein Electrostatics in the Julia Language}},
    year = {2018},
    journal = {Journal of Computational Science},
    volume = {28},
    pages = {193-203},
    doi = {10.1016/j.jocs.2018.08.008}
}
```


## [References](@id Bibliography)

 * **[Åqv90]**
   J. Åqvist, Ion-water interaction potentials derived from free energy pertubation
   simulations. J. Phys. Chem. 94: 8021, 1990.

 * **[Kea86]**
   P. Keast, Moderate degree tetrahedral quadrature formulas. CMAME 55: 339-348, 1986.

 * **[Rad48]**
   J. Radon, Zur mechanischen Kubatur (in German). Monatsh. für Math. 52(4): 286-300, 1948.

 * **[Rja90]**
   S. Rjasanow, Vorkonditionierte iterative Auflösung von Randelementgleichungen für die
   Dirichlet-Aufgabe (in German). Wissenschaftliche Schriftreihe der Technischen
   Universität Karl-Marx-Stadt, 7/1990.

 * **[Ste03]**
   O. Steinbach, Numerische Näherungsverfahren für elliptische Randwertprobleme - Finite
   Elemente und Randelemente (in German). Advances in Numerical Matheamtics. Teubner
   Verlag/GWV Fachverlage GmbH, Wiesbaden, 2003.

 * **[Xie16]**
   D. Xie, H. W. Volkmer, and J. Ying, Analytical solutions of nonlocal Poisson dielectric
   models with multiple point charges inside a dielectric sphere. Physical Review E 93(4):
   043304, 2016.
