# Input file support

```@meta
    CurrentModule = ProteinES.IO
```

Currently supported input file formats with different models:

| File type          | Surface model | Volume model | Charges included |
|--------------------|:-------------:|:------------:|:----------------:|
| [HMO](@ref)        | ✓             |              | ✓                |
| [Mcsf](@ref)       |               | ✓            |                  |
| [MSMS](@ref)       | ✓             |              |                  |

## HMO
```@docs
    readhmo
```

## Mcsf
```@docs
    readmcsf
```

## MSMS
```@docs
    readmsms
```

## Internal
```@meta
    DocTestSetup = quote
        using ProteinES.IO: readhmo_nodes, readhmo_elements, readhmo_charges,
                            readmcsf_nodes, readmcsf_elements, readmsms_nodes,
                            readmsms_elements
    end
```

```@docs
    readhmo_nodes
    readhmo_elements
    readhmo_charges
    readmcsf_nodes
    readmcsf_elements
    readmsms_nodes
    readmsms_elements
```

```@meta
    DocTestSetup = nothing
```
