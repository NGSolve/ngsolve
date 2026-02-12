# Netgen .vol to Gmsh/VTK Converter

## Overview

Convert Netgen `.vol` mesh files to Gmsh `.msh` or VTK `.vtk` format for visualization.

## Requirements

- NGSolve / Netgen (installed)
- Gmsh (optional, for viewing .msh files)
- ParaView (optional, for viewing .vtk files)

## Usage

### Basic Gmsh Conversion

```bash
# Convert .vol to .msh (Gmsh format)
python vol_to_gmsh.py input.vol

# Output: input.msh
```

### VTK Conversion

```bash
# Convert .vol to .vtk (VTK format)
python vol_to_gmsh.py input.vol --vtk

# Output: input.vtk
```

### Custom Output Name

```bash
# Specify output filename
python vol_to_gmsh.py input.vol output.msh
python vol_to_gmsh.py input.vol output.vtk --vtk
```

### View in Gmsh

```bash
# Convert and open in Gmsh
python vol_to_gmsh.py input.vol --view
```

## Examples

### Example 1: Convert cube mesh

```bash
python vol_to_gmsh.py install_ksugahar/share/ngsolve/cube.vol
```

Output:
```
Loading: install_ksugahar/share/ngsolve/cube.vol
  Vertices: 228
  Volume elements: 756
  Surface elements: 338
Exporting to Gmsh format: install_ksugahar\share\ngsolve\cube.msh
[OK] Conversion complete: install_ksugahar\share\ngsolve\cube.msh

To view in Gmsh:
  gmsh install_ksugahar\share\ngsolve\cube.msh
```

### Example 2: Convert coil mesh

```bash
python vol_to_gmsh.py install_ksugahar/share/ngsolve/coil.vol
```

Output:
```
Loading: install_ksugahar/share/ngsolve/coil.vol
  Vertices: 331
  Volume elements: 1709
  Surface elements: 320
Exporting to Gmsh format: install_ksugahar\share\ngsolve\coil.msh
[OK] Conversion complete: install_ksugahar\share\ngsolve\coil.msh
```

### Example 3: VTK format for ParaView

```bash
python vol_to_gmsh.py install_ksugahar/share/ngsolve/square.vol --vtk
```

Output:
```
Loading NGSolve mesh: install_ksugahar/share/ngsolve/square.vol
Exporting to VTK: install_ksugahar\share\ngsolve\square.vtk
[OK] VTK export complete: install_ksugahar\share\ngsolve\square.vtk

To view in ParaView:
  paraview install_ksugahar\share\ngsolve\square.vtk
```

## Gmsh Format Details

The script exports to **Gmsh Format 2.x** (ASCII format):

- Section `$Nodes`: Vertex coordinates (1-indexed)
- Section `$Elements`: Elements with type IDs
- Supports up to 2nd order elements

### Element Type IDs (Gmsh 2.x)

| Type ID | Element | Nodes |
|---------|---------|-------|
| 2 | Triangle | 3 |
| 4 | Tetrahedron | 4 |
| 5 | Hexahedron | 8 |
| 6 | Wedge/Prism | 6 |

For details, see `gmsh.pdf` (Gmsh reference manual).

## VTK Format Details

VTK format is suitable for:
- ParaView visualization
- Python analysis with PyVista
- Legacy VTK readers

## Tested Files

Successfully converted:
- `cube.vol` (228 vertices, 756 elements)
- `coil.vol` (331 vertices, 1709 elements)
- `shaft.vol` (558 vertices, 1622 elements)
- `square.vol` (2D mesh)

## Troubleshooting

### Error: File not found

Check that the input `.vol` file exists:
```bash
ls -l input.vol
```

### Error: Gmsh not found in PATH

Install Gmsh from https://gmsh.info/ or disable `--view` option.

### Unicode encoding errors (Windows)

Script now uses ASCII characters `[OK]` instead of Unicode checkmarks for Windows cp932 compatibility.

## Related Tools

- **Gmsh**: https://gmsh.info/
- **ParaView**: https://www.paraview.org/
- **PyVista**: https://docs.pyvista.org/

## Notes

- Gmsh format version can be selected with `--format gmsh` or `--format gmsh2`
- Default: Gmsh 2.x format (widely supported)
- VTK format includes mesh geometry only (no field data)
