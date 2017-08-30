# Nonlocal Electrostatics in Structured Solvents

`NESSie.jl` is a software package for the popular [Julia language](https://julialang.org), providing efficient numerical solvers for local and nonlocal protein electrostatics computations in structured solvents.


## Requirements
In addition to the Julia language itself (version 0.6), `NESSie.jl` depends on the following Julia packages:

 * [Distances.jl](http://github.com/JuliaStats/Distances.jl)
 * [JSON.jl](http://github.com/JuliaIO/JSON.jl)
 * [LightXML.jl](http://github.com/JuliaIO/LightXML.jl)
 * [Documenter.jl](http://github.com/JuliaDocs/Documenter.jl) (optional; required for documentation)
 * [FactCheck.jl](http://github.com/JuliaArchive/FactCheck.jl) (optional; required for tests)

All of these packages can be directly installed from Julia's interactive command line (invoked by typing `julia` into your terminal):
```julia
Pkg.add("Distances")
Pkg.add("JSON")
Pkg.add("LightXML")
Pkg.add("Documenter") # optional
Pkg.add("FactCheck")  # optional
```
The same commands can be used to check whether the respective package is already installed.


## Installation
`NESSie.jl` is currently under development and not yet included in the official [Julia package repositories](https://pkg.julialang.org/). However, you can download a copy of the package using [this link](https://github.com/tkemmer/NESSie.jl/archive/master.zip) or `git` on your terminal:

```sh
git clone https://github.com/tkemmer/nessie.jl
```

In order to use `NESSie.jl`, just tell Julia where to find it. This can be done by adding the following line either to Julia's startup file (`~/.juliarc.jl`) or to the top of your script:
```julia
push!(LOAD_PATH, "/path/to/nessie.jl/src")
```

Alternatively, Julia can be started with an additional argument pointing to the `src/` directory of your local `NESSie.jl` copy:
```sh
julia -L /path/to/nessie.jl/src
```


## Usage

### Modules
All functionality is organized into separate modules:
 * **`NESSie         `** Models and utility functions
 * **`NESSie.BEM     `** Local and nonlocal BEM solvers
 * **`NESSie.FEM     `** Nonlocal FEM solver (experimental; only available in `fem` branch)
 * **`NESSie.Format  `** Input and output file formats


### Example
The following Julia code shows how to compute and print the nonlocal reaction field energy of a single Na+ ion (modelled as a sperical vacuum-filled system) in water:

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
@printf "Reaction field energy: %f kJ/mol" val
```
More examples are available in the `doc/examples/` directory.


## Documentation
A detailed documentation of `NESSie.jl` is available [online](https://tkemmer.github.io/NESSie.jl/latest/).

You can also build the same documentation manually using the following command in the `doc/` directory:
```sh
julia make.jl
```
Please note that this additionally requires the `Documenter.jl` package to be installed on your machine. After the build process has finished successfully, the documentation can be found in the `doc/build/` directory.


## Testing
`NESSie.jl` provides tests for most of its functions (requires `FactCheck.jl`). You can run the test suite with the following command in the `test/` directory:
```sh
julia runtests.jl
```
The output can optionally be colorized when using Julia with the `--color=yes` argument. The current version of `FactCheck.jl` unfortunately produces some deprecation warnings with Julia 0.6. These warnings can be suppressed with the `--depwarn=no` argument:

```sh
julia --color=yes --depwarn=no runtests.jl
```

## Tools
`NESSie.jl` ships with the following extra tools:

 * **Converter.jl**, a simple Julia script to convert surface and volume meshes into different formats
 * **Mesher**, a wrapper application for the [GAMer](http://www.fetk.org/codes/gamer/) library to create surface and volume meshes from [PDB](https://www.rcsb.org/pdb) files

Please refer to the respective tool's README file for more details.
