# Nonlocal Electrostatics in Structured Solvents
[![](https://img.shields.io/github/workflow/status/tkemmer/NESSie.jl/CI?style=for-the-badge)](https://github.com/tkemmer/NESSie.jl/actions/workflows/CI.yml)
[![](https://img.shields.io/github/license/tkemmer/NESSie.jl?style=for-the-badge)](https://github.com/tkemmer/NESSie.jl/blob/master/LICENSE)
[![](https://img.shields.io/badge/docs-stable-blue.svg?style=for-the-badge)](https://tkemmer.github.io/NESSie.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg?style=for-the-badge)](https://tkemmer.github.io/NESSie.jl/dev)

`NESSie.jl` is a software package for the popular [Julia language](https://julialang.org),
providing efficient numerical solvers for local and nonlocal protein electrostatics
computations in structured solvents.


## Installation
This version of `NESSie.jl` requires Julia 1.0 or above. In the Julia shell, switch to the
`Pkg` shell by pressing `]` and enter the following command:

```sh
pkg> add https://github.com/tkemmer/NESSie.jl
```

Analogously, the `NESSie.jl` can be removed at any time using the following command in the
`Pkg` shell:
```sh
pkg> rm NESSie
```


## Usage

### Modules
All functionality is organized into separate modules:
 * **`NESSie          `** Models and utility functions
 * **`NESSie.BEM      `** Local and nonlocal BEM solvers
 * **`NESSie.Format   `** Input and output file formats
 * **`NESSie.TestModel`** Test models with analytical solutions


### Example
The following Julia code shows how to compute and print the nonlocal reaction field energy
of a single Na+ ion (modeled as a spherically-symmetric, vacuum-filled system) in water:

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


## Documentation
A detailed documentation of `NESSie.jl` is available
[online](https://tkemmer.github.io/NESSie.jl/dev/).

You can also build the same documentation manually using the following command in the
`doc/` directory:
```sh
shell> julia make.jl
```
Please note that this additionally requires the
[`Documenter.jl`](http://github.com/JuliaDocs/Documenter.jl) package to be installed. This can be achieved by using the following command in the `Pkg` shell:
```sh
pkg> add Documenter
```
After the build process has finished successfully, the documentation can be
found in the `doc/build/` directory.


## Testing
`NESSie.jl` provides tests for most of its functions. You can run the test suite with the
following command in the `Pkg` shell:
```sh
pkg> test NESSie
```


## Tools
`NESSie.jl` ships with the following extra tools:

 * **Converter.jl**, a simple Julia script to convert surface and volume meshes into
 different formats
 * **Mesher**, a wrapper application for the [GAMer](http://www.fetk.org/codes/gamer/)
 library to create surface and volume meshes from [PDB](https://www.rcsb.org/pdb) files

Please refer to the respective tool's README file for more details.


## Citing
If you use `NESSie.jl` in your research, please cite the following publication:
> Kemmer, T, Rjasanow, S., Hildebrandt, A (2018). NESSie. jl - Efficient and Intuitive
> Finite Element and Boundary Element Methods for Nonlocal Protein Electrostatics in the
> Julia Language. Journal of Computational Science 28, 193-203.

A pre-formatted citation for BibTeX can be found in [CITATION.bib](https://github.com/tkemmer/NESSie.jl/blob/master/CITATION.bib).
