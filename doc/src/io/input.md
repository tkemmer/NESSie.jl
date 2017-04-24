# Input file support

```@meta
    CurrentModule = ProteinES.IO
```

Currently supported input file formats with different models:

| File type          | Surface model | Volume model | Charges included |
|--------------------|:-------------:|:------------:|:----------------:|
| [HMO](@ref)        | ✓             |              | ✓                |

## HMO
```@docs
    readhmo
```

## Internal
```@meta
    DocTestSetup = quote
        using ProteinES.IO: readhmo_nodes, readhmo_elements, readhmo_charges
    end
```

```@docs
    readhmo_nodes
    readhmo_elements
    readhmo_charges
```

```@meta
    DocTestSetup = nothing
```
