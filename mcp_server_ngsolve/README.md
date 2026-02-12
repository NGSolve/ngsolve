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

## Tools (31 total)

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

### Diagnostic & Debugging (4)
- `ngsolve_server_info` - Get server version and status
- `ngsolve_list_objects` - List all objects in server state
- `ngsolve_get_object_info` - Get detailed information about an object
- `ngsolve_clear_state` - Clear all objects from server state (reset)

### Eddy Current Analysis (T-Omega Method) (5)
- `ngsolve_compute_genus` - Compute genus (number of holes) of multiply-connected domain
- `ngsolve_compute_loop_fields` - Compute loop fields for domains with holes
- `ngsolve_t_omega_setup` - Set up T-Omega FE spaces (HCurl + H1)
- `ngsolve_t_omega_solve_coupled` - Solve T-Omega system with loop current coupling
- `ngsolve_loop_current_analysis` - Analyze loop currents (resistance, inductance)

### Magnetostatic Analysis (2-Scalar Method) (7)
- `ngsolve_two_scalar_setup` - Set up Reduced-Total scalar potential FE spaces
- `ngsolve_compute_h0_coil` - Compute H₀ source field from current-carrying coil
- `ngsolve_compute_h0_pm` - Compute H₀ source field from permanent magnet
- `ngsolve_two_scalar_solve` - Solve magnetostatic problem with 2-scalar method
- `ngsolve_h_to_omega` - Convert H field to Omega scalar potential on boundary
- `ngsolve_compute_theta_field` - Compute Theta (Θ) field for coil jump (MMF = N×I)
- `ngsolve_two_scalar_solve_with_jump` - Solve with coil current jump using Theta field

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

## Validation Examples

### Rotating Magnet Eddy Current Analysis (Time Domain)

Complete validation examples for **transient eddy current analysis** using Radia-NGSolve coupling are available in the Radia repository at [`examples/NGSolve_Integration/rotating_magnets/`](https://github.com/ksugahar/Radia/tree/master/examples/NGSolve_Integration/rotating_magnets).

**Physical Model:**
- Rotating 1mm³ permanent magnet (Br = 0.2 T) moving over 0.5mm copper plate (σ = 5.8×10⁷ S/m)
- 180 timesteps with 4°/step rotation (total 720°, 2 full rotations)
- Time-dependent analysis: Backward Euler method for time discretization

**Two Time-Domain Formulation Comparison:**

#### 1. A-Φ Method (Vector-Scalar Potential) - Not yet in MCP

**File:** [`comparison_A_Phi_method.py`](https://github.com/ksugahar/Radia/blob/master/examples/NGSolve_Integration/rotating_magnets/comparison_A_Phi_method.py)

**Formulation:**
```
Vector potential: A_total = A_ext + A_r
Electric potential: Φ
Governing equations:
  (1) ∇×(1/μ ∇×A_r) + σ(∂A_r/∂t + ∇Φ) = -σ∂A_ext/∂t  (Ampère + Faraday)
  (2) ∇·[σ(∂A_r/∂t + ∇Φ)] = -∇·[σ∂A_ext/∂t]          (Current continuity)

Field reconstruction:
  B = curl(A_total) = curl(A_ext + A_r)
  E = -∂A_total/∂t - grad(Φ)
  J = σE (eddy current density)
```

**Implementation features:**
- Radia provides A_ext via `'a'` field type
- HCurl(nograds=True) for A_r + H1 for Φ (mixed formulation)
- Tree-cotree gauge automatically applied via `nograds=True`
- Direct computation of Joule loss: P = ∫ J·J/σ dV

**Current MCP status:** ❌ **Not implemented** - Candidate for future v1.4.0 release

#### 2. T-Ω Method (Current-Magnetic Scalar Potential) - Available in MCP v1.2.0+

**File:** [`comparison_T_Omega_method.py`](https://github.com/ksugahar/Radia/blob/master/examples/NGSolve_Integration/rotating_magnets/comparison_T_Omega_method.py)

**Formulation:**
```
Current potential: T (J = curl(T))
Magnetic scalar potential: Ω (H = H_ext - grad(Ω))
Governing equations:
  (conductor) ∇×(ν∇×T) + σ∂T/∂t = -σ∇(∂Ω/∂t)
  (all domain) ∇·μ(H_ext - ∇Ω) = 0

Field reconstruction:
  J = curl(T) (eddy current density, conductor only)
  H = H_ext - grad(Ω)
  B = μH
```

**Implementation features:**
- Radia provides H_ext via `'h'` field type
- HCurl(nograds=True) for T (conductor only) + H1 for Ω (global)
- T defined only in conductor using `definedon` → DOF reduction
- Loop fields handled for multiply-connected conductors

**Current MCP status:** ✅ **Implemented in v1.2.0** - Use these tools:
```python
# 1. Topology analysis
ngsolve_compute_genus(mesh_name="mesh")

# 2. Loop fields (if genus > 0)
ngsolve_compute_loop_fields(mesh_name="mesh", domain="conductor", order=2)

# 3. T-Omega setup
ngsolve_t_omega_setup(mesh_name="mesh", conductor_domain="conductor", order=2)

# 4. Solve T-Omega system
ngsolve_t_omega_solve_coupled(
    fespace_name="t_omega_space",
    sigma=5.8e7,
    mu=1.257e-6,
    conductor_domain="conductor"
)
```

**Validation Results:**

| Aspect | Result | Notes |
|--------|--------|-------|
| Maxwell relation | curl(A_ext) ≈ B_ext/μ₀ | Relative error < 0.1% |
| Eddy current pattern | Both methods agree qualitatively | Peak under moving magnet |
| Magnetic energy | Consistent between methods | W_mag = (1/2μ) ∫ B·B dV |
| Joule loss | A-Φ: Direct calculation | P = ∫ J·J/σ dV |
| Time evolution | 180 steps successfully | Backward Euler stable |

**Key Insights for MCP Implementation:**

1. **External field handling:**
   - A-Φ: Requires `A_ext` from Radia → Use `radia_ngsolve_create_field(field_type='a')`
   - T-Ω: Requires `H_ext` from Radia → Use `radia_ngsolve_create_field(field_type='h')`

2. **Time discretization:**
   - Both use Backward Euler: (u^(n+1) - u^n)/Δt for stability
   - Requires storing previous timestep solution
   - Implicit scheme ensures unconditional stability

3. **Gauge fixing:**
   - Both methods use `nograds=True` for automatic tree-cotree gauge
   - Essential for uniqueness of curl-based potentials

4. **DOF efficiency:**
   - A-Φ: A_r defined in all domains (higher DOF)
   - T-Ω: T defined only in conductor (lower DOF, more efficient)

**Future MCP Enhancement (A-Φ Method):**

To fully support the A-Φ validation example, the following tools would be needed:

```python
# Proposed for v1.4.0
ngsolve_a_phi_setup(mesh_name, conductor_domain, order=2)
ngsolve_a_phi_solve_transient(
    fespace_name,
    a_ext_field,  # From Radia
    sigma,
    mu,
    dt,
    num_steps
)
ngsolve_compute_eddy_current_a_phi(solution_name)
```

For detailed documentation of the validation setup, formulations, and results, see the [README.md](https://github.com/ksugahar/Radia/blob/master/examples/NGSolve_Integration/rotating_magnets/README.md) in the validation directory.

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

### 2-Scalar Method with Coil Jump

**When to use coil jump handling:**
- Current-carrying coils (electromagnets, motors)
- Situations where scalar potential Ω must be multi-valued
- When Ampère's law ∮H·dl = N×I needs explicit representation

**Theta (Θ) field setup:**

```python
# Step 1: Define cut surface (where Ω jumps)
# Cut surface must intersect the current loop completely
cut_minus = "coil_cut_minus"  # Θ = 0 side
cut_plus = "coil_cut_plus"    # Θ = N×I side

# Step 2: Compute Theta field
NI = 1000  # 100 turns × 10 A = 1000 A·turns (MMF)
ngsolve_compute_theta_field(
    mesh_name="mesh",
    coil_domain="coil",
    cut_surface_minus=cut_minus,
    cut_surface_plus=cut_plus,
    magnetomotive_force=NI
)

# Step 3: Solve with jump
ngsolve_two_scalar_solve_with_jump(
    fespace_name="two_scalar_space",
    theta_field_name="theta_field",
    coil_domain="coil",
    air_domain="air"
)
```

**Cut surface placement:**
- Choose a surface that crosses the coil current path once
- Typically: planar surface perpendicular to coil axis
- Must be consistent with mesh boundaries
- Example for z-axis coil: use x=0 plane with y≥0

**Verification:**
```python
# Check Ampère's law: ∮H·dl = N×I
# Integrate H around a closed path enclosing the coil
# Result should equal magnetomotive force
```

**Common pitfalls:**
- **Double counting**: Don't include coil current in both H₀ and Θ
- **Cut surface**: Must be a complete surface, not just edges
- **Sign convention**: Θ jump direction follows right-hand rule with current

## Troubleshooting

### Kelvin Transformation Issues

**Problem: Large errors at Kelvin boundary (r = R)**

**原因 (Causes):**
- Kelvin radius R too small (R < 2× inner domain radius)
- Insufficient mesh resolution at boundary (maxh > R/20)
- Permeability discontinuity not properly handled

**解決策 (Solutions):**
```python
# 1. Increase Kelvin radius
kelvin_create_mesh_with_transform(
    inner_radius=0.10,
    kelvin_radius=0.25,  # Was 0.15 → increase to 2.5× inner
    outer_radius=0.40
)

# 2. Refine mesh at boundary
kelvin_create_mesh_with_transform(
    inner_radius=0.10,
    kelvin_radius=0.25,
    outer_radius=0.40,
    maxh=0.012  # Was 0.020 → refine to < R/20
)

# 3. Check permeability values
kelvin_omega_reduced_omega_solve(
    permeability_inner=μᵣ,
    permeability_outer=1.0,  # Must be 1.0 for air in outer domain
    ...
)
```

**Problem: Solution diverges or solver fails to converge**

**原因 (Causes):**
- Inner/outer domain overlap (inner_radius ≥ kelvin_radius)
- Improper boundary conditions at r = R
- Iterative solver without preconditioner

**解決策 (Solutions):**
```python
# 1. Check domain sizes
assert inner_radius < kelvin_radius < outer_radius
# Example: 0.10 < 0.25 < 0.40 ✓

# 2. Use direct solver for debugging
kelvin_omega_reduced_omega_solve(
    solver="direct",  # Not "iterative"
    ...
)

# 3. If using iterative, check tolerance
kelvin_omega_reduced_omega_solve(
    solver="iterative",
    tolerance=1e-8,  # Default may be too loose
    ...
)
```

**Problem: Field values incorrect far from magnet (r → R)**

**原因 (Causes):**
- Transformation not properly applied
- Outer domain radius too small (outer_radius < 1.5× kelvin_radius)
- Field type mismatch (evaluating wrong field)

**解決策 (Solutions):**
```python
# 1. Verify transformation parameters
# Rule: outer_radius ≥ 1.5× kelvin_radius for smooth decay
kelvin_create_mesh_with_transform(
    kelvin_radius=0.25,
    outer_radius=0.40,  # = 1.6× kelvin_radius ✓
    ...
)

# 2. Check field type
# Use H-field (not B-field) in outer domain for better accuracy
kelvin_omega_reduced_omega_solve(
    field_type="h",  # Recommended for Kelvin transform
    ...
)
```

### Radia-NGSolve Coupling Issues

**Problem: Field interpolation errors at mesh boundaries**

**原因 (Causes):**
- Coordinate system mismatch between Radia and NGSolve
- Units inconsistency (Radia in mm, NGSolve in m)
- Mesh extends into Radia magnet geometry

**解決策 (Solutions):**
```python
# 1. Always use meters in Radia
import radia as rad
rad.FldUnits('m')  # CRITICAL before creating geometry

# 2. Verify coordinate alignment
# Radia magnet center should match NGSolve mesh reference point

# 3. Keep mesh away from magnet surfaces
# Rule: mesh boundaries > 0.01m from Radia object surfaces
```

**Problem: NGSolve cannot import Radia object**

**原因 (Causes):**
- Workspace session not found
- mcp_shared module not accessible
- Session expired or cleared

**解決策 (Solutions):**
```python
# 1. List available sessions
ngsolve_workspace_list_radia_objects()

# 2. Check symbolic link
# Verify: S:\NGSolve\01_GitHub\ngsolve_ksugahar\mcp_shared
#      → S:\Radia\01_Github\mcp_shared

# 3. Re-export from Radia
# In Radia MCP:
radia_workspace_export_object(
    object_name="magnet",
    export_geometry=True,
    export_fields=True
)
```

### Performance Issues

**Problem: Solver too slow for large problems**

**解決策 (Solutions):**
```python
# 1. Use iterative solver for large problems (> 50k DOFs)
kelvin_omega_reduced_omega_solve(
    solver="iterative",  # Not "direct"
    tolerance=1e-8
)

# 2. Enable Radia H-matrix for field evaluation
# In Radia MCP:
radia_ngsolve_enable_hmatrix(enable=True, precision=1e-6)

# 3. Reduce mesh resolution in outer domain
# Inner: maxh = 0.010
# Outer: maxh = 0.020  # Can be coarser
```

**Problem: Memory exhausted during solve**

**解決策 (Solutions):**
```python
# 1. Reduce mesh resolution
kelvin_create_mesh_with_transform(
    maxh=0.020,  # Was 0.010 → reduce by 2×
    ...
)

# 2. Use iterative solver (lower memory)
kelvin_omega_reduced_omega_solve(
    solver="iterative",
    ...
)

# 3. Reduce outer domain extent if possible
# Smaller outer_radius → fewer elements
```

### Verification and Validation

**ケルビン変換の正しさを確認する手順 (Steps to verify Kelvin transformation):**

```python
# 1. Use analytical comparison for sphere geometry
kelvin_compare_analytical(
    solution_name="H_solution",
    geometry_type="sphere",
    radius=0.10
)
# Expected error: < 5%

# 2. Check field continuity at Kelvin boundary
# Sample H-field at r = R from both sides
# |H(R-ε) - H(R+ε)| should be small

# 3. Verify far-field decay
# H(r) should decay as 1/r³ for dipole field
# Check at multiple radii: 1.1R, 1.3R, 1.5R

# 4. Energy conservation
kelvin_compute_perturbation_energy(
    solution_name="H_solution"
)
# Total energy should match analytical value within 10%
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
