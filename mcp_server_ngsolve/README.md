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

### Kelvin Transformation Example

```
User: "Solve magnetostatics with Kelvin transformation for unbounded domain"

Claude calls:
  kelvin_create_mesh_with_transform(
    geometry_type="cylinder",
    inner_radius=0.15,
    outer_radius=0.30,
    height=0.40,
    maxh=0.015,
    kelvin_radius=0.30
  )
  kelvin_omega_reduced_omega_solve(
    mesh_name="kelvin_mesh",
    permeability_inner=1000.0,
    permeability_outer=1.0,
    magnetization=[0, 0, 954930],
    solver="direct"
  )
  kelvin_export_vtk(
    solution_name="H_solution",
    output_file="kelvin_result.vtk"
  )

Result: VTK file with H-field solution including transformed outer region
```

### Complete Workflow Example

```
User: "Create 10cm NdFeB magnet in Radia, export to workspace, import to NGSolve, and solve with Kelvin transform"

Step 1 - Radia MCP Server:
  radia_geometry_set_units(units="m")
  radia_geometry_create_recmag(
    center=[0, 0, 0],
    dimensions=[0.1, 0.1, 0.1],
    magnetization=[0, 0, 954930]
  )
  radia_workspace_export_object(
    object_name="magnet",
    export_geometry=True,
    export_fields=True
  )

Step 2 - NGSolve MCP Server:
  ngsolve_radia_import_object(
    session_id="<from_step1>",
    radia_object_name="magnet"
  )
  kelvin_create_mesh_with_transform(
    geometry_type="box",
    inner_radius=0.10,
    outer_radius=0.25,
    kelvin_radius=0.25,
    maxh=0.015
  )
  kelvin_omega_reduced_omega_solve(
    mesh_name="kelvin_mesh",
    permeability_inner=1.0,
    permeability_outer=1.0,
    radia_session_id="<from_step1>",
    solver="iterative"
  )

Result: Complete unbounded domain simulation using Radia+NGSolve
```

## Best Practices

### Units and Coordinate Systems

**Always use meters (SI units):**

```python
# In Radia
import radia as rad
rad.FldUnits('m')  # Required for NGSolve integration

# NGSolve uses meters by default
```

**Coordinate consistency:**
- Radia origin → NGSolve mesh center
- Both use right-handed Cartesian coordinates (x, y, z)
- Ensure geometric alignment when importing Radia objects

### Kelvin Transformation Guidelines

**Choosing Kelvin radius (R):**

| Inner domain radius | Kelvin radius (R) | Outer domain radius | Transformation quality |
|-------------------|------------------|-------------------|----------------------|
| 0.10 m | 0.20 m | 0.30 m | Good (R = 2× inner) |
| 0.10 m | 0.25 m | 0.40 m | Better (R = 2.5× inner) |
| 0.10 m | 0.30 m | 0.50 m | Best (R = 3× inner) |

**Rule:** Set Kelvin radius R = 2-3× inner domain radius for accurate far-field representation.

**Permeability transformation:**
- Inner region (Ω): μ(r) = material permeability
- Outer region (Ω'): μ'(r') = (R/r')² × μ₀
- Automatic transformation applied by solver

### Mesh Resolution

**For accurate field evaluation:**

| Region | Recommended maxh | Purpose |
|--------|-----------------|---------|
| Inner domain (Ω) | h < L/10 | Field accuracy in physical region |
| Kelvin boundary | h < R/20 | Smooth transformation |
| Outer domain (Ω') | h < 2×(R/10) | Far-field representation |

Where L = characteristic length of inner geometry, R = Kelvin radius.

### Solver Selection

**Direct solver:**
- Use for: Small problems (< 50k DOFs)
- Pros: Exact solution, no convergence issues
- Cons: O(N³) complexity, high memory

**Iterative solver (BiCGSTAB):**
- Use for: Large problems (> 50k DOFs)
- Pros: O(N) memory, faster for large systems
- Cons: Requires good preconditioner
- Set tolerance: 1e-8 for high accuracy

### Performance Optimization

**Radia H-matrix acceleration:**
- Enable for batch evaluation (> 100 points)
- Expected speedup: 10-100× for large datasets
- Precision: Use eps=1e-6 for field calculations

**NGSolve mesh refinement:**
- Start with coarse mesh (maxh=0.020)
- Refine near boundaries and material interfaces
- Use adaptive refinement for critical regions

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

## Development Policy

### Branch Management
- **Feature branches**: Create feature branches for Pull Request purposes
  - Feature branches are temporary and should be deleted after PR approval and merge
  - Naming convention: `feature/mcp-tools-enhancement`, `fix/kelvin-mesh-bug`
- **Master branch**: Always kept in sync with the latest version
  - Research lab internal policy: master branch is continuously updated with latest stable code
  - Direct commits to master are allowed for internal development
  - Master branch always reflects the current working state

### Pull Request Workflow (for upstream contributions)
1. Create a feature branch from master
2. Make your changes and commit with descriptive messages
3. Push to your fork and create a Pull Request
4. After PR approval and merge, delete the feature branch
5. Keep master branch clean and stable

### Internal Development (Research Lab)
- Master branch is the primary development branch
- Always keep master in sync with the latest version
- Feature branches are used only for PR submissions to upstream
- Internal changes can be committed directly to master

## Technical Background

### Kelvin Transformation (Ω-Reduced Ω Method)

The Kelvin transformation maps an unbounded exterior domain to a bounded domain, enabling FEM simulation of infinite-extent problems.

**Mathematical formulation:**
- Transformation: r' = R²/r (inversion about sphere of radius R)
- Physical domain (Ω): r < R (bounded)
- Transformed domain (Ω'): R < r' < ∞ → R < R²/r < R (bounded)

**Permeability transformation:**
```
μ'(r') = (R/r')² μ₀
```

**Field relationships:**
```
H(r) in Ω  →  H'(r') in Ω' with transformed permeability
B = μH      →  B' = μ'H' = (R/r')² μ₀ H'
```

**Advantages:**
- Converts unbounded problem to bounded FEM domain
- Exact representation of far-field behavior (r → ∞)
- No artificial boundary conditions needed
- Compatible with standard FEM solvers

**Implementation:**
```python
# Create mesh with Kelvin transformation
kelvin_create_mesh_with_transform(
    inner_radius=0.15,    # Physical domain radius
    outer_radius=0.30,    # Computational domain radius
    kelvin_radius=0.25    # Transformation radius R
)

# Solve with automatic permeability transformation
kelvin_omega_reduced_omega_solve(
    permeability_inner=μᵣ,     # Physical domain
    permeability_outer=1.0      # Air (μ₀ in outer domain)
)
```

**Verification:**
- Analytical solutions available for sphere geometry
- Use `kelvin_compare_analytical()` to verify implementation
- Expected error: < 5% for proper mesh resolution

### Radia-NGSolve Coupling Architecture

**Data flow:**
```
Radia (MMM) → Workspace → NGSolve (FEM)
    ↓                         ↓
  Geometry                 Mesh + BC
  Materials                Interpolation
  Field data               FEM solution
```

**Workspace format:**
- JSON metadata: geometry, materials, units
- NPZ field data: structured grids for scipy.interpolate
- Session-based isolation

**Integration methods:**
1. **Field interpolation:** Import pre-computed Radia field as NGSolve CoefficientFunction
2. **Magnetization import:** Use Radia-computed M(r) as FEM input
3. **Hybrid solver:** Radia for PM regions, NGSolve for μᵣ(H) nonlinear regions

## See Also

- [Radia MCP Server](https://github.com/ksugahar/Radia) - Companion server for Radia
- [NGSolve Official](https://github.com/NGSolve/ngsolve) - Upstream NGSolve repository
- [Kelvin Transformation Examples](../../../NGSolve/2025_12_14_Kelvin変換/)
- [Dual Server Deployment](../docs/DUAL_SERVER_DEPLOYMENT.md)
