# Example code
This directory contains example code for electrostatics computations using the local and
nonlocal BEM solvers. The potentials are usually computed for simple test models and
compared to their analytical solutions (or approximations of it) as provided by the
`NESSie.TestModel` module. The example scripts can be used with the data provided in
the `data/` directory of this repository.


## Prerequisites
All code examples presented here additionally require
[PyPlot.jl](https://github.com/JuliaPy/PyPlot.jl). Make sure to have this package installed:
```sh
pkg> add PyPlot
```

## `bornbem.jl`
Comparison of local and nonlocal exterior potentials of Born ions.

### Usage
```bash
julia bornbem.jl BORN_ION
```
where `ION_NAME` is a valid Born ion name (see
[Documentation](
    https://tkemmer.github.io/NESSie.jl/latest/lib/util.html#NESSie.Born.bornion
)).


## `xiebem.jl`
Comparison of nonlocal reaction field potentials (or electrostatic potentials) for a
Poisson test model generated from a given PQR file.

### Usage
```bash
julia xiebem.jl [PQR_FILE]
```
