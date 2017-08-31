# `NESSie.jl` file format converter

## Usage:

```sh
julia converter.jl [<format in>] <file in> <format out> <file out>
```

If not specified, the input format will be guessed by file extension.
Multi-file input, as provided by some file formats (e.g., MSMS), has
to be specified without file extension. For instance,

```sh
julia converter.jl msms mymesh surface.json mymesh.json
```

will read the MSMS-generated surface mesh from both files mymesh.face
and mymesh.vert and convert it into a XML3D-compatible JSON format.

### Available input formats:

| Argument   | Mesh type | File type/extension |
|------------|-----------|---------------------|
| **`off`**  | Surface   | GAMer/`.off`        |
| **`mcsf`** | Volume    | GAMer/`.m`          |
| **`msms`** | Surface   | MSMS/`.{face,vert}` |
| **`hmo`**  | Surface   | HyperMesh/`.hmo`    |

### Available output formats:

| Argument           | Mesh type | File type/extension |
|--------------------|-----------|---------------------|
| **`nodes.xml`**    | Nodes     | XML3D/`.xml`        |
| **`surface.skel`** | Surface   | SKEL/`.skel`        |
| **`volume.vtu`**   | Volume    | VTK/`.vtu`          |
| **`nodes.json`**   | Nodes     | XML3D/`.json`       |
| **`volume.skel`**  | Volume    | SKEL/`.skel`        |
| **`surface.vtp`**  | Surface   | VTK/`.vtp`          |
| **`surface.json`** | Surface   | XML3D/`.json`       |
