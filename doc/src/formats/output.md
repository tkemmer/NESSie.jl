# Output formats

```@meta
    CurrentModule = NESSie.Format
```

Currently supported output file formats with different models:

| File type          | Point cloud | Surface model | Volume model |
|--------------------|:-----------:|:-------------:|:------------:|
| [SKEL](@ref)       |             | ✓             | ✓            |
| [STL](@ref)        |             | ✓             |              |
| [VTK](@ref)        |             | ✓             | ✓            |
| [XML3D](@ref)/JSON | ✓           | ✓             |              |
| [XML3D](@ref)/XML  | ✓           |               |              |

## SKEL
```@docs
    writeskel
```

## STL
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
