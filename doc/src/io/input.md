# Input file support

```@meta
    CurrentModule = ProteinES.IO
```

Currently supported input file formats with different models:

| File type          | Surface model | Volume model | Charges included |
|--------------------|:-------------:|:------------:|:----------------:|
| [HMO](@ref)        | ✓             |              | ✓                |
| [Mcsf](@ref)       |               | ✓            |                  |

## HMO
```@docs
    readhmo
```

## Mcsf
```@docs
    readmcsf
```

## Internal
```@meta
    DocTestSetup = quote
        using ProteinES.IO: readhmo_nodes, readhmo_elements, readhmo_charges,
                            readmcsf_nodes, readmcsf_elements
    end
```

```@docs
    readhmo_nodes
    readhmo_elements
    readhmo_charges
    readmcsf_nodes
    readmcsf_elements
```

```@meta
    DocTestSetup = nothing
```
