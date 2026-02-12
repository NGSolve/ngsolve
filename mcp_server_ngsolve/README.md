# NGSolve MCP Server

Model Context Protocol server for NGSolve FEM with Radia coupling.

## Overview

This server provides NGSolve mesh generation and Radia field integration through MCP. Works in conjunction with `mcp_server_radia` for complete electromagnetic simulation workflow.

## Features

- **Mesh Generation**: Netgen-based mesh creation (box, cylinder)
- **Mesh Import**: GMSH file import
- **Radia Coupling**: Import Radia fields from shared workspace
- **Interpolation**: Create interpolated fields from Radia data
- **Kelvin Transformation**: Unbounded domain simulation using Kelvin transformation

## Installation

```bash
# Install NGSolve
pip install ngsolve==6.2.2405

# Install scipy for interpolation
pip install scipy

# Install MCP SDK
pip install mcp
```

## Tools (15 total)

### Mesh Generation (4)
- `ngsolve_mesh_create_box` - Create box mesh
- `ngsolve_mesh_create_cylinder` - Create cylinder mesh
- `ngsolve_mesh_import_file` - Import GMSH mesh
- `ngsolve_mesh_get_info` - Get mesh statistics

### Radia Coupling (4)
- `ngsolve_radia_import_object` - Import Radia object from workspace
- `ngsolve_radia_get_field_data` - Get pre-computed field data
- `ngsolve_radia_create_interpolated_field` - Create interpolated CF
- `ngsolve_workspace_list_radia_objects` - List available Radia objects

### Kelvin Transformation (7)
- `kelvin_create_mesh_with_transform` - Create mesh with Kelvin transformation (sphere, cylinder, box)
- `kelvin_omega_reduced_omega_solve` - Solve using Ω-Reduced Ω method
- `kelvin_compute_perturbation_energy` - Compute perturbation field energy
- `kelvin_export_vtk` - Export solution to VTK format
- `kelvin_compare_analytical` - Compare with analytical solution (sphere)
- `kelvin_adaptive_mesh_refinement` - Adaptive mesh refinement (planned)
- `kelvin_check_availability` - Check NGSolve availability

## Usage

### Claude Desktop Configuration

```json
{
  "mcpServers": {
    "ngsolve": {
      "command": "python",
      "args": ["-m", "mcp_server_ngsolve.server"],
      "cwd": "S:\\NGSolve\\01_GitHub\\ngsolve_ksugahar",
      "env": {
        "PYTHONPATH": "S:\\NGSolve\\01_GitHub\\ngsolve_ksugahar;S:\\Radia\\01_Github"
      }
    }
  }
}
```

### Mesh Generation Example

```
User: "Create a 20cm cubic iron mesh with 5mm elements"

Claude calls:
  ngsolve_mesh_create_box(
    pmin=[-0.1, -0.1, -0.1],
    pmax=[0.1, 0.1, 0.1],
    maxh=0.005,
    material_name="iron"
  )
```

### Radia Coupling Example

```
User: "Import the magnet from Radia session abc123 and create interpolated field"

Claude calls:
  ngsolve_radia_import_object(
    session_id="abc123...",
    radia_object_name="magnet"
  )
  ngsolve_radia_create_interpolated_field(
    session_id="abc123...",
    radia_object_name="magnet",
    mesh_name="iron_mesh"
  )
```

## Shared Workspace

This server requires access to `mcp_shared` module for workspace communication.

**Setup via Symbolic Link** (Already configured):
```powershell
# Symbolic link created at:
# S:\NGSolve\01_GitHub\ngsolve_ksugahar\mcp_shared
#   -> S:\Radia\01_Github\mcp_shared
```

**Or via PYTHONPATH**:
Set `PYTHONPATH` to include `S:\Radia\01_Github\mcp_shared`

## See Also

- [Dual Server Deployment](../docs/DUAL_SERVER_DEPLOYMENT.md)
- [MCP Server Architecture](../docs/MCP_SERVER_ARCHITECTURE.md)
- [Radia MCP Server](../mcp_server_radia/README.md)
