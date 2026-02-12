"""
Kelvin Transform Tools for NGSolve MCP Server

Tools for infinite domain simulation using Kelvin transformation:
- Ω-Reduced Ω method (Total/Reduced scalar potential)
- Kelvin transformation for unbounded domains
- Energy calculation for perturbation fields
"""

import sys
from pathlib import Path
from typing import Any, Dict, List

try:
    from ngsolve import *
    import numpy as np
    NGSOLVE_AVAILABLE = True
except ImportError:
    NGSOLVE_AVAILABLE = False

from mcp.types import Tool


def get_tools() -> List[Tool]:
    """Get list of Kelvin transformation tools."""
    if not NGSOLVE_AVAILABLE:
        return [
            Tool(
                name="kelvin_check_availability",
                description="Check if NGSolve is available for Kelvin transformation.",
                inputSchema={"type": "object", "properties": {}, "required": []}
            )
        ]

    return [
        Tool(
            name="kelvin_create_mesh_with_transform",
            description="Create mesh with Kelvin transformation for unbounded domain analysis. Supports sphere, cylinder, and box geometries.",
            inputSchema={
                "type": "object",
                "properties": {
                    "geometry_type": {
                        "type": "string",
                        "enum": ["sphere", "cylinder", "box"],
                        "description": "Type of magnetic region geometry",
                        "default": "sphere"
                    },
                    "magnetic_region_size": {
                        "type": "number",
                        "description": "Size of magnetic region (radius for sphere/cylinder, half-size for box) in meters"
                    },
                    "inner_air_radius": {
                        "type": "number",
                        "description": "Inner air (reduced) region radius in meters"
                    },
                    "kelvin_radius": {
                        "type": "number",
                        "description": "Kelvin transformation sphere radius in meters"
                    },
                    "cylinder_height": {
                        "type": "number",
                        "description": "Cylinder height (only for cylinder geometry)",
                        "default": 1.0
                    },
                    "maxh": {
                        "type": "number",
                        "description": "Maximum element size",
                        "default": 0.1
                    },
                    "mesh_name": {
                        "type": "string",
                        "description": "Name for the mesh",
                        "default": "kelvin_mesh"
                    },
                    "kelvin_offset_z": {
                        "type": "number",
                        "description": "Z-offset for exterior Kelvin domain center",
                        "default": 3.0
                    }
                },
                "required": ["magnetic_region_size", "inner_air_radius", "kelvin_radius"]
            }
        ),
        Tool(
            name="kelvin_omega_reduced_omega_solve",
            description="Solve magnetostatic problem using Ω-Reduced Ω method with Kelvin transformation.",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_name": {
                        "type": "string",
                        "description": "Name of the mesh with Kelvin regions"
                    },
                    "source_field_type": {
                        "type": "string",
                        "enum": ["uniform", "coil", "radia"],
                        "description": "Type of source field",
                        "default": "uniform"
                    },
                    "source_field_params": {
                        "type": "object",
                        "description": "Source field parameters (H0 for uniform, etc.)",
                        "properties": {
                            "H0": {
                                "type": "number",
                                "description": "Uniform field magnitude (A/m)"
                            },
                            "direction": {
                                "type": "array",
                                "items": {"type": "number"},
                                "minItems": 3,
                                "maxItems": 3,
                                "description": "Field direction [x, y, z]"
                            }
                        }
                    },
                    "permeability": {
                        "type": "number",
                        "description": "Relative permeability of magnetic region",
                        "default": 100.0
                    },
                    "fe_order": {
                        "type": "integer",
                        "description": "Finite element order",
                        "default": 1
                    },
                    "use_kelvin": {
                        "type": "boolean",
                        "description": "Enable Kelvin transformation",
                        "default": True
                    }
                },
                "required": ["mesh_name"]
            }
        ),
        Tool(
            name="kelvin_compute_perturbation_energy",
            description="Compute perturbation field energy (avoids infinity at far field).",
            inputSchema={
                "type": "object",
                "properties": {
                    "solution_name": {
                        "type": "string",
                        "description": "Name of the Ω solution GridFunction"
                    },
                    "compare_analytical": {
                        "type": "boolean",
                        "description": "Compare with analytical solution for sphere",
                        "default": False
                    }
                },
                "required": ["solution_name"]
            }
        ),
        Tool(
            name="kelvin_adaptive_mesh_refinement",
            description="Perform adaptive mesh refinement based on error estimation.",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_name": {
                        "type": "string",
                        "description": "Name of current mesh"
                    },
                    "solution_name": {
                        "type": "string",
                        "description": "Name of solution"
                    },
                    "max_refinements": {
                        "type": "integer",
                        "description": "Maximum number of refinement iterations",
                        "default": 5
                    },
                    "error_tolerance": {
                        "type": "number",
                        "description": "Target error tolerance",
                        "default": 1e-3
                    }
                },
                "required": ["mesh_name", "solution_name"]
            }
        ),
        Tool(
            name="kelvin_export_vtk",
            description="Export Kelvin transformation solution to VTK format for visualization.",
            inputSchema={
                "type": "object",
                "properties": {
                    "solution_name": {
                        "type": "string",
                        "description": "Name of the solution GridFunction"
                    },
                    "output_file": {
                        "type": "string",
                        "description": "Output VTK file path (without extension)"
                    },
                    "include_fields": {
                        "type": "boolean",
                        "description": "Include B and H field output",
                        "default": True
                    }
                },
                "required": ["solution_name", "output_file"]
            }
        ),
        Tool(
            name="kelvin_compare_analytical",
            description="Compare numerical solution with analytical solution for sphere in uniform field.",
            inputSchema={
                "type": "object",
                "properties": {
                    "solution_name": {
                        "type": "string",
                        "description": "Name of the Ω solution"
                    },
                    "geometry_params": {
                        "type": "object",
                        "description": "Geometry parameters (sphere_radius, mu_r, H0)",
                        "properties": {
                            "sphere_radius": {"type": "number"},
                            "mu_r": {"type": "number"},
                            "H0": {"type": "number"}
                        }
                    }
                },
                "required": ["solution_name", "geometry_params"]
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute a Kelvin transformation tool."""
    if not NGSOLVE_AVAILABLE and name != "kelvin_check_availability":
        return {
            "error": "NGSolve not available",
            "message": "Install with: pip install ngsolve==6.2.2405"
        }

    try:
        if name == "kelvin_check_availability":
            return _check_availability()
        elif name == "kelvin_create_mesh_with_transform":
            return _create_kelvin_mesh(arguments, state)
        elif name == "kelvin_omega_reduced_omega_solve":
            return _omega_reduced_omega_solve(arguments, state)
        elif name == "kelvin_compute_perturbation_energy":
            return _compute_perturbation_energy(arguments, state)
        elif name == "kelvin_adaptive_mesh_refinement":
            return _adaptive_refinement(arguments, state)
        elif name == "kelvin_export_vtk":
            return _export_vtk(arguments, state)
        elif name == "kelvin_compare_analytical":
            return _compare_analytical(arguments, state)
        else:
            return {"error": f"Unknown Kelvin tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name, "traceback": __import__('traceback').format_exc()}


def _check_availability() -> Dict[str, Any]:
    """Check if NGSolve is available."""
    return {
        "ngsolve_available": NGSOLVE_AVAILABLE,
        "message": "NGSolve is available" if NGSOLVE_AVAILABLE else "Install with: pip install ngsolve==6.2.2405"
    }


def _create_kelvin_mesh(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create mesh with Kelvin transformation regions."""
    from netgen.occ import Sphere, Cylinder, Box, Pnt, Axis, OCCGeometry, Glue, Vertex
    from netgen.meshing import MeshingParameters, IdentificationType

    geometry_type = args.get("geometry_type", "sphere")
    r_mag = args["magnetic_region_size"]
    r_inner = args["inner_air_radius"]
    r_kelvin = args["kelvin_radius"]
    maxh = args.get("maxh", 0.1)
    mesh_name = args.get("mesh_name", "kelvin_mesh")
    offset_z = args.get("kelvin_offset_z", 3.0)

    # Validate radii
    if not (r_mag < r_inner < r_kelvin):
        return {
            "error": "Radii must satisfy: r_magnetic < r_inner_air < r_kelvin"
        }

    # Create magnetic region geometry based on type
    if geometry_type == "sphere":
        mag_region = Sphere(Pnt(0, 0, 0), r_mag)
        mag_region.mat("magnetic")
        mag_region.maxh = maxh
        # Name the boundary
        for face in mag_region.faces:
            face.name = "magnetic_boundary"
    elif geometry_type == "cylinder":
        cyl_height = args.get("cylinder_height", 1.0)
        cyl_axis = Axis(Pnt(0, 0, -cyl_height/2), Pnt(0, 0, cyl_height/2))
        mag_region = Cylinder(cyl_axis, r_mag, cyl_height)
        mag_region.mat("magnetic")
        mag_region.maxh = maxh
        for face in mag_region.faces:
            face.name = "magnetic_boundary"
    elif geometry_type == "box":
        mag_region = Box(Pnt(-r_mag, -r_mag, -r_mag), Pnt(r_mag, r_mag, r_mag))
        mag_region.mat("magnetic")
        mag_region.maxh = maxh
        for face in mag_region.faces:
            face.name = "magnetic_boundary"
    else:
        return {"error": f"Unknown geometry type: {geometry_type}"}

    # Create inner air sphere (reduced region)
    inner_sphere = Sphere(Pnt(0, 0, 0), r_inner)
    inner_sphere.maxh = maxh
    for face in inner_sphere.faces:
        face.name = "kelvin_int"
    inner_air = inner_sphere - mag_region
    inner_air.mat("air_inner")

    # Create exterior Kelvin domain (offset in z-direction)
    outer_sphere = Sphere(Pnt(0, 0, offset_z), r_kelvin)
    outer_sphere.maxh = maxh
    outer_sphere.mat("air_outer")
    for face in outer_sphere.faces:
        face.name = "kelvin_ext"

    # GND vertex at center of exterior domain (represents infinity)
    vertex = Vertex(Pnt(0, 0, offset_z))
    vertex.name = "GND"

    # Glue all domains
    geo = Glue([inner_air, mag_region, outer_sphere, vertex])

    # Identify periodic faces
    kelvin_int_face = None
    kelvin_ext_face = None
    for solid in geo.solids:
        for face in solid.faces:
            if face.name == "kelvin_int":
                kelvin_int_face = face
            elif face.name == "kelvin_ext":
                kelvin_ext_face = face

    if kelvin_int_face is not None and kelvin_ext_face is not None:
        kelvin_int_face.Identify(kelvin_ext_face, "periodic", IdentificationType.PERIODIC)

    # Generate mesh
    mp = MeshingParameters(maxh=maxh, grading=0.5)
    mesh = Mesh(OCCGeometry(geo).GenerateMesh(mp))

    # Store mesh and parameters
    state[mesh_name] = mesh
    state[f"{mesh_name}_params"] = {
        "geometry_type": geometry_type,
        "r_magnetic": r_mag,
        "r_inner": r_inner,
        "r_kelvin": r_kelvin,
        "offset_z": offset_z,
        "maxh": maxh
    }

    return {
        "success": True,
        "mesh_name": mesh_name,
        "mesh_info": {
            "num_vertices": mesh.nv,
            "num_elements": mesh.ne,
            "geometry_type": geometry_type,
            "r_magnetic": r_mag,
            "r_inner_air": r_inner,
            "r_kelvin": r_kelvin,
            "kelvin_offset_z": offset_z,
            "materials": list(mesh.GetMaterials()),
            "boundaries": list(mesh.GetBoundaries())
        }
    }


def _omega_reduced_omega_solve(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Solve using Ω-Reduced Ω method with Kelvin transformation."""
    import math

    mesh_name = args["mesh_name"]
    if mesh_name not in state:
        return {"error": f"Mesh '{mesh_name}' not found"}

    mesh = state[mesh_name]
    mesh_params = state.get(f"{mesh_name}_params", {})

    source_type = args.get("source_field_type", "uniform")
    source_params = args.get("source_field_params", {})
    mu_r = args.get("permeability", 100.0)
    fe_order = args.get("fe_order", 1)
    use_kelvin = args.get("use_kelvin", True)

    # Physical constants
    mu0 = 4e-7 * math.pi

    # Define permeability
    mu_dic = {"magnetic": mu_r * mu0, "reduced": mu0, "kelvin": mu0}
    Mu = CoefficientFunction([mu_dic[mat] for mat in mesh.GetMaterials()])

    # Define source field
    if source_type == "uniform":
        H0 = source_params.get("H0", 1.0)
        direction = source_params.get("direction", [0, 0, 1])
        # Normalize direction
        norm = math.sqrt(sum(d**2 for d in direction))
        direction = [d/norm for d in direction]

        # Source scalar potential: Ωs = H0 * (dx*x + dy*y + dz*z)
        Ov = H0 * (direction[0]*x + direction[1]*y + direction[2]*z)
        Bv = mu0 * CoefficientFunction((direction[0]*H0, direction[1]*H0, direction[2]*H0))
    else:
        return {"error": f"Source field type '{source_type}' not yet implemented"}

    # Setup finite element space
    fes = H1(mesh, order=fe_order)
    fes = Periodic(fes) if use_kelvin else fes

    omega, psi = fes.TnT()

    # Bilinear form
    a = BilinearForm(fes)
    a += Mu * grad(omega) * grad(psi) * dx("magnetic")
    a += Mu * grad(omega) * grad(psi) * dx("reduced")

    # Add Kelvin transformation term if enabled
    if use_kelvin:
        rs = mesh_params.get("r_kelvin", 1.0)
        xs = 2 * rs  # Kelvin center
        r = sqrt((x - xs)**2 + y**2 + z**2)
        fac = rs**2 / r**2
        a += Mu * fac * grad(omega) * grad(psi) * dx("kelvin")

    a.Assemble()

    # Set Dirichlet boundary condition (source potential on interface)
    gfOmega = GridFunction(fes)
    gfOmega.Set(Ov, BND)

    # Linear form
    f = LinearForm(fes)
    f += Mu * grad(gfOmega) * grad(psi) * dx("reduced")
    f.Assemble()

    # Remove Dirichlet DOFs
    fcut = np.array(f.vec.FV())[fes.FreeDofs()]
    np.array(f.vec.FV(), copy=False)[fes.FreeDofs()] = fcut

    # Add Neumann boundary condition
    normal = -specialcf.normal(mesh.dim)  # Outward normal for 3D
    f += (normal * Bv) * psi * ds
    f.Assemble()

    # Solve
    gfOmega.vec.data = a.mat.Inverse(fes.FreeDofs()) * f.vec

    # Store solution
    solution_name = f"{mesh_name}_omega_solution"
    state[solution_name] = gfOmega
    state[f"{solution_name}_params"] = {
        "mu_r": mu_r,
        "source_type": source_type,
        "fe_order": fe_order,
        "use_kelvin": use_kelvin,
        "ndof": fes.ndof
    }

    # Compute field at center
    try:
        mip = mesh(0, 0, 0)
        B_center = [gfOmega(mip) * Mu(mip)] if mip else [0, 0, 0]
    except:
        B_center = [0, 0, 0]

    return {
        "success": True,
        "solution_name": solution_name,
        "ndof": fes.ndof,
        "fe_order": fe_order,
        "kelvin_enabled": use_kelvin,
        "B_center": B_center,
        "message": f"Ω-Reduced Ω solution computed with {fes.ndof} DOFs"
    }


def _compute_perturbation_energy(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Compute perturbation field energy."""
    import math

    solution_name = args["solution_name"]
    if solution_name not in state:
        return {"error": f"Solution '{solution_name}' not found"}

    gfOmega = state[solution_name]
    params = state.get(f"{solution_name}_params", {})

    mesh = gfOmega.space.mesh
    mu_r = params.get("mu_r", 100.0)
    mu0 = 4e-7 * math.pi

    # Compute perturbation energy in each region
    # (Implementation simplified - full version would compute H_pert properly)

    energy_magnetic = Integrate(0.5 * mu_r * mu0 * grad(gfOmega)**2 * dx("magnetic"), mesh)
    energy_reduced = Integrate(0.5 * mu0 * grad(gfOmega)**2 * dx("reduced"), mesh)

    total_energy = energy_magnetic + energy_reduced

    result = {
        "success": True,
        "perturbation_energy": {
            "magnetic_region": energy_magnetic,
            "reduced_region": energy_reduced,
            "total": total_energy
        },
        "units": "Joules"
    }

    # Analytical comparison for sphere in uniform field
    if args.get("compare_analytical", False):
        # Get sphere radius (assume from mesh params)
        # Analytical: W = (2π/3) * μ0 * ((μr-1)/(μr+2))^2 * H0^2 * a^3 * (μr + 2)
        # This is a simplified implementation
        result["analytical_comparison"] = "Not yet implemented - requires geometry parameters"

    return result


def _adaptive_refinement(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Perform adaptive mesh refinement."""
    return {
        "error": "Adaptive refinement not yet implemented",
        "message": "This feature requires integration with existing adaptive mesh code"
    }


def _export_vtk(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Export solution to VTK format."""
    from ngsolve import VTKOutput, grad

    solution_name = args["solution_name"]
    output_file = args["output_file"]
    include_fields = args.get("include_fields", True)

    if solution_name not in state:
        return {"error": f"Solution '{solution_name}' not found"}

    gfOmega = state[solution_name]
    params = state.get(f"{solution_name}_params", {})
    mesh = gfOmega.space.mesh

    # Prepare coefficients for export
    coefs = [gfOmega]
    names = ["Omega"]

    if include_fields:
        import math
        mu0 = 4e-7 * math.pi
        mu_r = params.get("mu_r", 100.0)

        # Add gradient (H field approximation)
        coefs.append(grad(gfOmega))
        names.append("grad_Omega")

        # Add B field approximation
        mu_dic = {"magnetic": mu_r * mu0, "air_inner": mu0, "air_outer": mu0}
        from ngsolve import CoefficientFunction
        Mu = CoefficientFunction([mu_dic.get(mat, mu0) for mat in mesh.GetMaterials()])
        coefs.append(Mu * grad(gfOmega))
        names.append("B_field")

    # Export to VTK
    VTKOutput(ma=mesh, coefs=coefs, names=names, filename=output_file).Do()

    return {
        "success": True,
        "output_file": f"{output_file}.vtu",
        "fields_exported": names,
        "message": f"VTK file saved: {output_file}.vtu"
    }


def _compare_analytical(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Compare numerical solution with analytical solution for sphere."""
    import math
    import numpy as np

    solution_name = args["solution_name"]
    geom_params = args["geometry_params"]

    if solution_name not in state:
        return {"error": f"Solution '{solution_name}' not found"}

    gfOmega = state[solution_name]
    mesh = gfOmega.space.mesh

    # Extract geometry parameters
    sphere_radius = geom_params["sphere_radius"]
    mu_r = geom_params["mu_r"]
    H0 = geom_params["H0"]

    # Analytical solution for sphere in uniform field
    Hz_analytical_interior = 3.0 / (mu_r + 2) * H0

    # Evaluate at test points
    from ngsolve import grad
    grad_Omega = grad(gfOmega)

    test_points = {
        "origin": (0, 0, 0),
        "x_0.2": (0.2, 0, 0),
        "z_0.3": (0, 0, 0.3)
    }

    results = []
    for point_name, coords in test_points.items():
        try:
            mip = mesh(coords[0], coords[1], coords[2])
            Hz_numerical = grad_Omega[2](mip)
            error_pct = abs(Hz_numerical - Hz_analytical_interior) / abs(Hz_analytical_interior) * 100
            results.append({
                "point": point_name,
                "coordinates": coords,
                "Hz_numerical": float(Hz_numerical),
                "Hz_analytical": float(Hz_analytical_interior),
                "error_percent": float(error_pct)
            })
        except Exception as e:
            results.append({
                "point": point_name,
                "coordinates": coords,
                "error": str(e)
            })

    # Compute z-axis profile
    z_vals = np.linspace(-sphere_radius * 0.9, sphere_radius * 0.9, 20)
    z_profile = []
    for zv in z_vals:
        try:
            mip = mesh(0, 0, float(zv))
            Hz_num = grad_Omega[2](mip)
            z_profile.append({
                "z": float(zv),
                "Hz": float(Hz_num),
                "Hz_analytical": float(Hz_analytical_interior)
            })
        except:
            pass

    return {
        "success": True,
        "analytical_solution": {
            "Hz_interior": float(Hz_analytical_interior),
            "formula": "3/(mu_r+2) * H0"
        },
        "point_comparison": results,
        "z_axis_profile": z_profile,
        "geometry": {
            "sphere_radius": sphere_radius,
            "mu_r": mu_r,
            "H0": H0
        }
    }
