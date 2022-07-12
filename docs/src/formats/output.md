# Output formats

```@meta
    CurrentModule = NESSie.Format
```

Currently supported output file formats with different models:

| File type          | Point cloud | Surface model | Volume model | Charges included |
|--------------------|:-----------:|:-------------:|:------------:|:----------------:|
| [HMO](@ref hmoout) |             | ✓             |              | ✓                |
| [OBJ](@ref)        |             | ✓             |              |                  |
| [SKEL](@ref)       |             | ✓             | ✓            |                  |
| [STL](@ref stlout) |             | ✓             |              |                  |
| [VTK](@ref)        |             | ✓             | ✓            |                  |
| [XML3D](@ref)/JSON | ✓           | ✓             |              |                  |
| [XML3D](@ref)/XML  | ✓           |               |              |                  |

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

## XML3D
```@docs
    writexml3d_json
    writexml3d_xml
```
