# `NESSie.jl` mesh generator

`NESSie.jl` provides two simple wrapper applications for the
[GAMer](http://fetk.org/codes/gamer/) library to either generate surface and volume meshes
for molecular systems or to generate surface meshes for spheres.

## Requirements
 * [GCC](https://gcc.gnu.org/) 5 or above
 * [GAMer](http://fetk.org/codes/gamer/)

## Installation
In this directory, simply call
```bash
make
```
to compile both applications `mesher` and `sphere`.

## `mesher`: molecule mesher
Generates both surface and volume meshs for the given molecular system.

### Usage
```bash
mesher <file in> <dir out> [<option>...]
```

### Supported input formats
| Name | File extension |
|------|----------------|
| PDB  | `.pdb`         |
| PQR  | `.pqr`         |

### Supported options

| Option                   | Default value | Description                            |
|--------------------------|------:|------------------------------------------------|
| `--coarsen-rate=<float>` |   0.3 | Coarsening rate                                |
| `--iso-value=<float>`    |   1.4 | Iso value for gaussian surface (marching cube) |
| `--max-iterations=<int>` |     6 | Maximum number of surface smoothing steps      |
| `--mesh-nodes-max=<int>` | 80000 | Maximum number of surface mesh vertices        |
| `--mesh-nodes-min=<int>` | 10000 | Minimum number of surface mesh vertices        |
| `--sphere-quality=<int>` |     6 | Number of bounding sphere vertices: ``4^n+2``  |
| `--sphere-ratio=<float>` |   2.0 | Sphere-to-molecule size ratio                  |


## `sphere`: sphere generator
Generates a surface mesh for a sphere with the given radius.

### Usage
```bash
sphere <file out> <sphere radius> <sphere quality>
```
where sphere quality ``n`` determines the number of vertices via ``4^n+2``.
