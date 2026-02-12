#!/usr/bin/env python
"""
Convert Netgen .vol mesh to Gmsh .msh format for visualization.

Usage:
    python vol_to_gmsh.py input.vol [output.msh]

Features:
- Converts Netgen volume mesh (.vol) to Gmsh format (.msh)
- Preserves element types (hex, tet, wedge)
- Supports NGSolve mesh with materials
- Can open in Gmsh for visualization

Requirements:
- NGSolve / Netgen
- Gmsh (for viewing, optional)

License:
    MIT License
    Copyright (c) 2026 NGSolve Contributors

Author:
    Created for NGSolve project mesh conversion utilities
"""

import sys
import os
from pathlib import Path


def convert_vol_to_gmsh(vol_file, gmsh_file=None, format='gmsh2'):
    """
    Convert Netgen .vol to Gmsh .msh format.

    Parameters:
    -----------
    vol_file : str
        Input .vol file path
    gmsh_file : str, optional
        Output .msh file path (default: same as input with .msh extension)
    format : str
        'gmsh' or 'gmsh2' (default: gmsh2 for modern Gmsh)

    Returns:
    --------
    str
        Path to the generated Gmsh file

    Raises:
    -------
    RuntimeError
        If mesh loading or export fails
    """
    from netgen.meshing import Mesh as NetgenMesh

    # Load Netgen mesh
    print(f"Loading: {vol_file}")
    ngmesh = NetgenMesh()

    try:
        ngmesh.Load(vol_file)
    except Exception as e:
        raise RuntimeError(f"Failed to load mesh file: {e}")

    # Get mesh info (updated API)
    nv = len(ngmesh.Points())
    ne = len(ngmesh.Elements3D())
    nse = len(ngmesh.Elements2D())

    print(f"  Vertices: {nv}")
    print(f"  Volume elements: {ne}")
    print(f"  Surface elements: {nse}")

    # Check for empty mesh
    if nv == 0:
        raise RuntimeError("Mesh has no vertices")

    # Warn if no volume elements (surface-only mesh)
    if ne == 0 and nse > 0:
        print("  Note: Surface-only mesh (no volume elements)")

    # Determine output filename
    if gmsh_file is None:
        vol_path = Path(vol_file)
        gmsh_file = str(vol_path.with_suffix('.msh'))

    # Export to Gmsh format
    print(f"Exporting to Gmsh format: {gmsh_file}")

    # Note: Netgen Export() uses format string to determine output format
    # For Gmsh: use 'Gmsh Format' or 'Gmsh2 Format'
    if format == 'gmsh2':
        format_str = 'Gmsh2 Format'
    else:
        format_str = 'Gmsh Format'

    try:
        ngmesh.Export(gmsh_file, format_str)
    except Exception as e:
        raise RuntimeError(f"Failed to export to Gmsh format: {e}")

    print(f"[OK] Conversion complete: {gmsh_file}")
    print(f"\nTo view in Gmsh:")
    print(f"  gmsh {gmsh_file}")

    return gmsh_file


def convert_vol_to_vtk(vol_file, vtk_file=None):
    """
    Convert Netgen .vol to VTK format (alternative to Gmsh).

    Uses NGSolve VTKOutput for conversion.

    Parameters:
    -----------
    vol_file : str
        Input .vol file path
    vtk_file : str, optional
        Output .vtk file path (default: same as input with .vtk extension)

    Returns:
    --------
    str
        Path to the generated VTK file

    Raises:
    -------
    RuntimeError
        If mesh loading or export fails
    """
    from ngsolve import Mesh, VTKOutput

    print(f"Loading NGSolve mesh: {vol_file}")

    try:
        mesh = Mesh(vol_file)
    except Exception as e:
        raise RuntimeError(f"Failed to load mesh with NGSolve: {e}")

    # Determine output filename
    if vtk_file is None:
        vol_path = Path(vol_file)
        vtk_file = str(vol_path.with_suffix(''))  # VTKOutput adds .vtk

    print(f"Exporting to VTK: {vtk_file}.vtu")

    # Export using VTKOutput (creates .vtu file - VTK Unstructured Grid)
    try:
        vtk = VTKOutput(mesh, coefs=[], names=[], filename=vtk_file)
        vtk.Do()
    except Exception as e:
        raise RuntimeError(f"Failed to export to VTK format: {e}")

    print(f"[OK] VTK export complete: {vtk_file}.vtu")
    print(f"\nTo view in ParaView:")
    print(f"  paraview {vtk_file}.vtu")

    return vtk_file + '.vtu'


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Convert Netgen .vol to Gmsh .msh or VTK format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert to Gmsh (default)
  python vol_to_gmsh.py mesh.vol

  # Convert to Gmsh with custom output name
  python vol_to_gmsh.py mesh.vol output.msh

  # Convert to VTK instead
  python vol_to_gmsh.py mesh.vol --vtk

  # View in Gmsh after conversion
  python vol_to_gmsh.py mesh.vol --view
        """
    )

    parser.add_argument('vol_file', help='Input .vol file')
    parser.add_argument('output_file', nargs='?', help='Output .msh or .vtk file (optional)')
    parser.add_argument('--vtk', action='store_true', help='Export to VTK instead of Gmsh')
    parser.add_argument('--format', choices=['gmsh', 'gmsh2'], default='gmsh2',
                        help='Gmsh format version (default: gmsh2)')
    parser.add_argument('--view', action='store_true', help='Open in Gmsh after conversion')

    args = parser.parse_args()

    # Check input file exists
    if not os.path.exists(args.vol_file):
        print(f"Error: File not found: {args.vol_file}")
        sys.exit(1)

    # Convert
    try:
        if args.vtk:
            output_file = convert_vol_to_vtk(args.vol_file, args.output_file)
        else:
            output_file = convert_vol_to_gmsh(args.vol_file, args.output_file, args.format)
    except RuntimeError as e:
        print(f"Error during conversion: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)

    # View in Gmsh if requested
    if args.view and not args.vtk:
        print(f"\nLaunching Gmsh...")
        import subprocess
        try:
            subprocess.run(['gmsh', output_file])
        except FileNotFoundError:
            print("Error: Gmsh not found in PATH")
            print("Install Gmsh: https://gmsh.info/")
        except Exception as e:
            print(f"Error launching Gmsh: {e}")


if __name__ == '__main__':
    main()
