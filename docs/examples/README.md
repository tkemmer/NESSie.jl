# Example code
This directory contains example code for electrostatics computations using the local and
nonlocal BEM solvers. The potentials are usually computed for simple test models and
compared to their analytical solutions (or approximations of it) as provided by the
`NESSie.TestModel` module.


## Prerequisites
In addition to NESSie.jl, all code examples presented here additionally require
[Plots.jl](https://github.com/JuliaPlots/Plots.jl). Make sure to have this package installed:
```sh
pkg> add Plots
```

## `bornbem.jl`
Comparison of local and nonlocal potentials of Born ions.

### Usage
```bash
julia -i bornbem.jl BORN_ION
```
where `ION_NAME` is a valid Born ion name (see
[Documentation](
    https://tkemmer.github.io/NESSie.jl/dev/lib/models/#NESSie.TestModel.bornion
)).


## `xiebem.jl`
Comparison of nonlocal electrostatic and reaction field potentials for a Poisson test model
generated from a given PQR file.

### Usage
```bash
julia -i xiebem.jl [PQR_FILE]
```
