# Input formats

```@meta
CurrentModule = NESSie.Format
```

Currently supported input file formats with different models:

| File type          | Surface model | Charges included |
|--------------------|:-------------:|:----------------:|
| [HMO](@ref)        | ✓            | ✓                |
| [MSMS](@ref)       | ✓            |                  |
| [OBJ](@ref)        | ✓            |                  |
| [OFF](@ref)        | ✓            |                  |
| [PQR](@ref)        |               | ✓ (charges only) |
| [STL](@ref)        | ✓            |                  |

## HMO
```@docs
readhmo
```

## MSMS
```@docs
readmsms
```

## OBJ
```@docs
readobj
```

## OFF
```@docs
readoff
```

## PQR
```@docs
readpqr
```

## STL
```@docs
readstl
```
