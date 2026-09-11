# Output formats

```@meta
    CurrentModule = NESSie.Format
```

Currently supported output file formats with different models:

| File type          | Surface model | Volume model | Charges included |
|--------------------|:-------------:|:------------:|:----------------:|
| [HMO](@ref hmoout) | ✓             |              | ✓                |
| [OBJ](@ref)        | ✓             |              |                  |
| [SKEL](@ref)       | ✓             | ✓            |                  |
| [STL](@ref stlout) | ✓             |              |                  |
| [VTK](@ref)        | ✓             | ✓            |                  |

## [HMO](@id hmoout)
```@docs
    writehmo
```

## OBJ
```@docs
    writeobj
```

## SKEL
```@docs
    writeskel
```

## [STL](@id stlout)
```@docs
    writestl
```

## VTK
```@docs
    writevtk
```
