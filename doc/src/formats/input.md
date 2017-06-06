# Input formats

```@meta
    CurrentModule = NESSie.Format
```

Currently supported input file formats with different models:

| File type          | Surface model | Volume model | Charges included |
|--------------------|:-------------:|:------------:|:----------------:|
| [HMO](@ref)        | ✓             |              | ✓                |
| [Mcsf](@ref)       |               | ✓            |                  |
| [MSMS](@ref)       | ✓             |              |                  |
| [OFF](@ref)        | ✓             |              |                  |
| [PQR](@ref)        |               |              | ✓ (charges only) |

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

## OFF
```@docs
    readoff
```

## PQR
```@docs
    readpqr
```
