"""
Two-Scalar Potential Method Tools for NGSolve MCP Server

Reduced-Total scalar potential method for magnetostatic problems with
permanent magnets or current sources.
"""

from typing import Any, Dict, List
import numpy as np

try:
    from ngsolve import *
    import ngsolve
    NGSOLVE_AVAILABLE = True
except ImportError:
    NGSOLVE_AVAILABLE = False

from mcp.types import Tool


def get_tools() -> List[Tool]:
    """Get list of two-scalar method tools."""
    if not NGSOLVE_AVAILABLE:
        return []

    return [
        Tool(
            name="ngsolve_two_scalar_setup",
            description="Set up Reduced-Total scalar potential FE spaces for 2-scalar method",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_name": {
                        "type": "string",
                        "description": "Name of the mesh object"
                    },
                    "source_domain": {
                        "type": "string",
                        "description": "Source domain name (permanent magnet or coil region)"
                    },
                    "air_domain": {
                        "type": "string",
                        "description": "Air/non-source domain name"
                    },
                    "outer_boundary": {
                        "type": "string",
                        "description": "Outer boundary name for Dirichlet BC",
                        "default": "outer"
                    },
                    "order": {
                        "type": "integer",
                        "description": "Finite element order (default: 2)",
                        "default": 2
                    }
                },
                "required": ["mesh_name", "source_domain", "air_domain"]
            }
        ),
        Tool(
            name="ngsolve_compute_h0_coil",
            description="Compute H₀ source field from current-carrying coil using simplified model",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_name": {
                        "type": "string",
                        "description": "Name of the mesh object"
                    },
                    "coil_center": {
                        "type": "array",
                        "description": "Coil center position [x, y, z]",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3
                    },
                    "coil_radius": {
                        "type": "number",
                        "description": "Coil radius (m)"
                    },
                    "coil_current": {
                        "type": "number",
                        "description": "Coil current (A)"
                    },
                    "coil_turns": {
                        "type": "integer",
                        "description": "Number of turns (default: 1)",
                        "default": 1
                    },
                    "coil_axis": {
                        "type": "array",
                        "description": "Coil axis direction [x, y, z] (default: [0,0,1])",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "default": [0, 0, 1]
                    }
                },
                "required": ["mesh_name", "coil_center", "coil_radius", "coil_current"]
            }
        ),
        Tool(
            name="ngsolve_compute_h0_pm",
            description="Compute H₀ source field from permanent magnet with uniform magnetization",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_name": {
                        "type": "string",
                        "description": "Name of the mesh object"
                    },
                    "magnetization": {
                        "type": "array",
                        "description": "Magnetization vector [Mx, My, Mz] (A/m)",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3
                    },
                    "pm_domain": {
                        "type": "string",
                        "description": "Permanent magnet domain name"
                    }
                },
                "required": ["mesh_name", "magnetization", "pm_domain"]
            }
        ),
        Tool(
            name="ngsolve_two_scalar_solve",
            description="Solve magnetostatic problem using 2-scalar potential method",
            inputSchema={
                "type": "object",
                "properties": {
                    "fespace_name": {
                        "type": "string",
                        "description": "Name of the 2-scalar FE space"
                    },
                    "h0_field_name": {
                        "type": "string",
                        "description": "Name of the H₀ source field CoefficientFunction"
                    },
                    "mu_source": {
                        "type": "number",
                        "description": "Permeability in source region (H/m)",
                        "default": 1.2566370614e-6
                    },
                    "mu_air": {
                        "type": "number",
                        "description": "Permeability in air region (H/m)",
                        "default": 1.2566370614e-6
                    },
                    "source_domain": {
                        "type": "string",
                        "description": "Source domain name"
                    },
                    "air_domain": {
                        "type": "string",
                        "description": "Air domain name"
                    },
                    "tolerance": {
                        "type": "number",
                        "description": "Solver tolerance",
                        "default": 1e-10
                    },
                    "max_iterations": {
                        "type": "integer",
                        "description": "Maximum iterations",
                        "default": 1000
                    }
                },
                "required": ["fespace_name", "h0_field_name", "source_domain", "air_domain"]
            }
        ),
        Tool(
            name="ngsolve_h_to_omega",
            description="Convert H field to Omega scalar potential on boundary",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_name": {
                        "type": "string",
                        "description": "Name of the mesh object"
                    },
                    "h_field_name": {
                        "type": "string",
                        "description": "Name of the H field CoefficientFunction"
                    },
                    "boundary_name": {
                        "type": "string",
                        "description": "Boundary name for conversion"
                    },
                    "order": {
                        "type": "integer",
                        "description": "FE order for Omega (default: 2)",
                        "default": 2
                    },
                    "tolerance": {
                        "type": "number",
                        "description": "Solver tolerance",
                        "default": 1e-12
                    }
                },
                "required": ["mesh_name", "h_field_name", "boundary_name"]
            }
        ),
        Tool(
            name="ngsolve_compute_theta_field",
            description="Compute Theta (Θ) scalar potential field for coil jump in 2-scalar method. Theta represents magnetomotive force (MMF) = N×I.",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_name": {
                        "type": "string",
                        "description": "Name of the mesh object"
                    },
                    "coil_domain": {
                        "type": "string",
                        "description": "Coil domain name where current flows"
                    },
                    "cut_surface_minus": {
                        "type": "string",
                        "description": "Cut surface boundary name (Θ=0 side)"
                    },
                    "cut_surface_plus": {
                        "type": "string",
                        "description": "Cut surface boundary name (Θ=NI side)"
                    },
                    "magnetomotive_force": {
                        "type": "number",
                        "description": "N×I (turns × current) in Ampere-turns"
                    },
                    "order": {
                        "type": "integer",
                        "description": "FE order (default: 2)",
                        "default": 2
                    },
                    "tolerance": {
                        "type": "number",
                        "description": "Solver tolerance",
                        "default": 1e-10
                    }
                },
                "required": ["mesh_name", "coil_domain", "cut_surface_minus", "cut_surface_plus", "magnetomotive_force"]
            }
        ),
        Tool(
            name="ngsolve_two_scalar_solve_with_jump",
            description="Solve magnetostatic problem with coil current jump using 2-scalar method + Theta field",
            inputSchema={
                "type": "object",
                "properties": {
                    "fespace_name": {
                        "type": "string",
                        "description": "Name of the 2-scalar FE space"
                    },
                    "theta_field_name": {
                        "type": "string",
                        "description": "Name of the Theta field (from ngsolve_compute_theta_field)"
                    },
                    "h0_field_name": {
                        "type": "string",
                        "description": "Name of the H₀ source field (optional)",
                        "default": None
                    },
                    "mu_coil": {
                        "type": "number",
                        "description": "Permeability in coil region (H/m)",
                        "default": 1.2566370614e-6
                    },
                    "mu_air": {
                        "type": "number",
                        "description": "Permeability in air region (H/m)",
                        "default": 1.2566370614e-6
                    },
                    "coil_domain": {
                        "type": "string",
                        "description": "Coil domain name"
                    },
                    "air_domain": {
                        "type": "string",
                        "description": "Air domain name"
                    },
                    "tolerance": {
                        "type": "number",
                        "description": "Solver tolerance",
                        "default": 1e-10
                    },
                    "max_iterations": {
                        "type": "integer",
                        "description": "Maximum iterations",
                        "default": 1000
                    }
                },
                "required": ["fespace_name", "theta_field_name", "coil_domain", "air_domain"]
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute a two-scalar method tool."""
    if not NGSOLVE_AVAILABLE:
        return {
            "error": "NGSolve not available",
            "message": "Install NGSolve: pip install ngsolve"
        }

    try:
        if name == "ngsolve_two_scalar_setup":
            return _two_scalar_setup(arguments, state)
        elif name == "ngsolve_compute_h0_coil":
            return _compute_h0_coil(arguments, state)
        elif name == "ngsolve_compute_h0_pm":
            return _compute_h0_pm(arguments, state)
        elif name == "ngsolve_two_scalar_solve":
            return await _two_scalar_solve(arguments, state)
        elif name == "ngsolve_h_to_omega":
            return await _h_to_omega(arguments, state)
        elif name == "ngsolve_compute_theta_field":
            return await _compute_theta_field(arguments, state)
        elif name == "ngsolve_two_scalar_solve_with_jump":
            return await _two_scalar_solve_with_jump(arguments, state)
        else:
            return {"error": f"Unknown two-scalar tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name, "traceback": __import__('traceback').format_exc()}


def _two_scalar_setup(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Set up Reduced-Total scalar potential FE spaces."""
    mesh_name = args["mesh_name"]
    source_domain = args["source_domain"]
    air_domain = args["air_domain"]
    outer_boundary = args.get("outer_boundary", "outer")
    order = args.get("order", 2)

    if mesh_name not in state:
        return {"error": f"Mesh '{mesh_name}' not found"}

    mesh = state[mesh_name]

    try:
        # Ωᵣ (reduced scalar potential) in source region
        fesR = H1(mesh, order=order, definedon=source_domain)

        # Ωₜ (total scalar potential) in air region
        fesT = H1(mesh, order=order, definedon=air_domain,
                 dirichlet=outer_boundary)

        # Mixed space
        fespace = fesR * fesT

        fespace_name = f"{mesh_name}_two_scalar_space"
        state[fespace_name] = fespace

        return {
            "success": True,
            "fespace_name": fespace_name,
            "ndof_reduced": fesR.ndof,
            "ndof_total": fesT.ndof,
            "ndof_combined": fespace.ndof,
            "order": order,
            "source_domain": source_domain,
            "air_domain": air_domain,
            "message": f"2-scalar FE space created with {fespace.ndof} DOFs"
        }

    except Exception as e:
        return {
            "error": str(e),
            "traceback": __import__('traceback').format_exc()
        }


def _compute_h0_coil(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Compute H₀ from current-carrying coil (circular loop approximation)."""
    mesh_name = args["mesh_name"]
    center = np.array(args["coil_center"])
    radius = args["coil_radius"]
    current = args["coil_current"]
    turns = args.get("coil_turns", 1)
    axis = np.array(args.get("coil_axis", [0, 0, 1]))

    if mesh_name not in state:
        return {"error": f"Mesh '{mesh_name}' not found"}

    mesh = state[mesh_name]

    try:
        # Normalize axis
        axis = axis / np.linalg.norm(axis)

        # Simplified H₀ for circular coil (on-axis approximation)
        # H(z) = (N*I*R²) / (2*(R² + z²)^(3/2)) along axis
        # This is a simplified model

        total_current = current * turns
        cx, cy, cz = center
        ax, ay, az = axis

        # Create CoefficientFunction for H₀
        # On-axis field (simplified)
        x, y, z = ngsolve.x, ngsolve.y, ngsolve.z

        # Distance from coil center
        dx = x - cx
        dy = y - cy
        dz = z - cz

        # Axial distance
        z_ax = dx*ax + dy*ay + dz*az

        # Radial distance from axis
        r_perp_sq = dx**2 + dy**2 + dz**2 - z_ax**2

        # On-axis approximation
        denom = (radius**2 + z_ax**2)**(3/2)
        H_mag = (total_current * radius**2) / (2 * denom)

        # Direction along axis
        H0_cf = CoefficientFunction((H_mag * ax, H_mag * ay, H_mag * az))

        h0_name = f"{mesh_name}_H0_coil"
        state[h0_name] = H0_cf

        return {
            "success": True,
            "h0_field_name": h0_name,
            "coil_center": center.tolist(),
            "coil_radius": radius,
            "total_current": total_current,
            "coil_axis": axis.tolist(),
            "message": f"H₀ field from coil computed (simplified on-axis model)"
        }

    except Exception as e:
        return {
            "error": str(e),
            "traceback": __import__('traceback').format_exc()
        }


def _compute_h0_pm(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Compute H₀ from permanent magnet."""
    mesh_name = args["mesh_name"]
    magnetization = np.array(args["magnetization"])
    pm_domain = args["pm_domain"]

    if mesh_name not in state:
        return {"error": f"Mesh '{mesh_name}' not found"}

    mesh = state[mesh_name]

    try:
        # For permanent magnet with uniform magnetization:
        # H₀ = M / μ₀ in the PM region
        # (simplified model - actual field is more complex)

        mu0 = 4e-7 * np.pi
        Mx, My, Mz = magnetization

        # H₀ = M / μ₀
        Hx = Mx / mu0
        Hy = My / mu0
        Hz = Mz / mu0

        H0_cf = CoefficientFunction((Hx, Hy, Hz))

        h0_name = f"{mesh_name}_H0_pm"
        state[h0_name] = H0_cf

        return {
            "success": True,
            "h0_field_name": h0_name,
            "magnetization": magnetization.tolist(),
            "h0_magnitude": float(np.linalg.norm([Hx, Hy, Hz])),
            "pm_domain": pm_domain,
            "message": f"H₀ field from PM computed (uniform magnetization model)"
        }

    except Exception as e:
        return {
            "error": str(e),
            "traceback": __import__('traceback').format_exc()
        }


async def _two_scalar_solve(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Solve 2-scalar potential problem."""
    fespace_name = args["fespace_name"]
    h0_field_name = args["h0_field_name"]
    mu_source = args.get("mu_source", 4e-7 * np.pi)
    mu_air = args.get("mu_air", 4e-7 * np.pi)
    source_domain = args["source_domain"]
    air_domain = args["air_domain"]
    tol = args.get("tolerance", 1e-10)
    max_iter = args.get("max_iterations", 1000)

    if fespace_name not in state:
        return {"error": f"FE space '{fespace_name}' not found"}
    if h0_field_name not in state:
        return {"error": f"H₀ field '{h0_field_name}' not found"}

    fespace = state[fespace_name]
    H0 = state[h0_field_name]

    try:
        mesh = fespace.mesh
        (Omega_r, Omega_t), (psi_r, psi_t) = fespace.TnT()

        # System matrix
        a = BilinearForm(fespace)
        # Source region: ∇·μ(H₀ - ∇Ωᵣ) = 0
        a += mu_source * grad(Omega_r) * grad(psi_r) * dx(source_domain)
        # Air region: ∇·μ(-∇Ωₜ) = 0
        a += mu_air * grad(Omega_t) * grad(psi_t) * dx(air_domain)

        # Source term from H₀
        f = LinearForm(fespace)
        f += mu_source * H0 * grad(psi_r) * dx(source_domain)

        # Assemble
        with TaskManager():
            a.Assemble()
            f.Assemble()

        # Solve
        gf = GridFunction(fespace)
        inv = CGSolver(a.mat, fespace.FreeDofs(), maxsteps=max_iter, tol=tol)
        gf.vec.data = inv * f.vec

        # Extract components
        gf_Omega_r, gf_Omega_t = gf.components

        # Reconstruct H and B fields
        H_source = H0 - grad(gf_Omega_r)
        H_air = -grad(gf_Omega_t)

        # Combined H field (using material indicator)
        # Note: This is simplified - proper implementation uses IfPos
        H_total_cf = H0 - grad(gf_Omega_r)  # Simplified

        B_source = mu_source * H_source
        B_air = mu_air * H_air

        # Store results
        solution_name = f"{fespace_name}_solution"
        state[solution_name] = gf
        state[f"{solution_name}_Omega_r"] = gf_Omega_r
        state[f"{solution_name}_Omega_t"] = gf_Omega_t
        state[f"{solution_name}_H_source"] = H_source
        state[f"{solution_name}_H_air"] = H_air
        state[f"{solution_name}_B_source"] = B_source
        state[f"{solution_name}_B_air"] = B_air

        # Compute some statistics
        H_norm_source = sqrt(Integrate(H_source**2 * dx(source_domain), mesh))
        H_norm_air = sqrt(Integrate(H_air**2 * dx(air_domain), mesh))

        return {
            "success": True,
            "solution_name": solution_name,
            "omega_r_name": f"{solution_name}_Omega_r",
            "omega_t_name": f"{solution_name}_Omega_t",
            "h_source_name": f"{solution_name}_H_source",
            "h_air_name": f"{solution_name}_H_air",
            "b_source_name": f"{solution_name}_B_source",
            "b_air_name": f"{solution_name}_B_air",
            "h_norm_source": float(H_norm_source),
            "h_norm_air": float(H_norm_air),
            "message": "2-scalar problem solved successfully"
        }

    except Exception as e:
        return {
            "error": str(e),
            "traceback": __import__('traceback').format_exc()
        }


async def _h_to_omega(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Convert H field to Omega on boundary."""
    mesh_name = args["mesh_name"]
    h_field_name = args["h_field_name"]
    boundary_name = args["boundary_name"]
    order = args.get("order", 2)
    tol = args.get("tolerance", 1e-12)

    if mesh_name not in state:
        return {"error": f"Mesh '{mesh_name}' not found"}
    if h_field_name not in state:
        return {"error": f"H field '{h_field_name}' not found"}

    mesh = state[mesh_name]
    H = state[h_field_name]

    try:
        # H1 space on boundary
        fesOmega = H1(mesh, order=order,
                     definedon=mesh.Boundaries(boundary_name),
                     complex=False)
        omega, psi = fesOmega.TnT()

        # Minimize ||∇Ω - H||² on boundary
        a = BilinearForm(fesOmega)
        a += grad(omega).Trace() * grad(psi).Trace() * ds

        f = LinearForm(fesOmega)
        f += (grad(psi).Trace() * H) * ds

        with TaskManager():
            a.Assemble()
            f.Assemble()

        gfOmega = GridFunction(fesOmega)
        inv = CGSolver(a.mat, fesOmega.FreeDofs(), maxsteps=200, tol=tol)
        gfOmega.vec.data = inv * f.vec

        omega_name = f"{mesh_name}_{boundary_name}_omega"
        state[omega_name] = gfOmega

        # Compute norm for verification
        norm = sqrt(Integrate(gfOmega**2 * ds(boundary_name), mesh))

        return {
            "success": True,
            "omega_name": omega_name,
            "boundary_name": boundary_name,
            "ndof": fesOmega.ndof,
            "omega_norm": float(norm),
            "message": f"H field converted to Omega on boundary '{boundary_name}'"
        }

    except Exception as e:
        return {
            "error": str(e),
            "traceback": __import__('traceback').format_exc()
        }


async def _compute_theta_field(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compute Theta scalar potential field for coil jump.

    Theta represents the magnetomotive force (MMF) jump across the coil:
    [Ω]⁺₋ = Θ = N×I (turns × current)

    Solves Laplace equation ∇²Θ = 0 in coil domain with boundary conditions:
    - Θ = 0 on cut_surface_minus
    - Θ = NI on cut_surface_plus
    """
    mesh_name = args["mesh_name"]
    coil_domain = args["coil_domain"]
    cut_minus = args["cut_surface_minus"]
    cut_plus = args["cut_surface_plus"]
    NI = args["magnetomotive_force"]
    order = args.get("order", 2)
    tol = args.get("tolerance", 1e-10)

    if mesh_name not in state:
        return {"error": f"Mesh '{mesh_name}' not found"}

    mesh = state[mesh_name]

    try:
        # H1 space on coil domain
        fesTheta = H1(mesh, order=order, definedon=coil_domain)
        theta, psi = fesTheta.TnT()

        # Laplace equation for harmonic Theta field
        a = BilinearForm(fesTheta)
        a += grad(theta) * grad(psi) * dx(coil_domain)

        # Grid function for Theta
        gfTheta = GridFunction(fesTheta)

        # Set boundary conditions on cut surfaces
        gfTheta.Set(0, BND, mesh.Boundaries(cut_minus))
        gfTheta.Set(NI, BND, mesh.Boundaries(cut_plus))

        # Right-hand side with Dirichlet BC incorporated
        f = LinearForm(fesTheta)
        f += -grad(gfTheta) * grad(psi) * dx(coil_domain)

        with TaskManager():
            a.Assemble()
            f.Assemble()

        # Remove Dirichlet boundary components
        fcut = np.array(f.vec.FV())[fesTheta.FreeDofs()]
        np.array(f.vec.FV(), copy=False)[fesTheta.FreeDofs()] = fcut

        # Solve for interior DOFs
        inv = CGSolver(a.mat, fesTheta.FreeDofs(), maxsteps=1000, tol=tol)
        gfTheta.vec.data += inv * f.vec

        # Store in state
        theta_name = f"{mesh_name}_theta_field"
        state[theta_name] = gfTheta

        # Compute gradient magnitude for verification
        grad_theta = grad(gfTheta)
        grad_norm = sqrt(Integrate(grad_theta * grad_theta * dx(coil_domain), mesh))

        return {
            "success": True,
            "theta_field_name": theta_name,
            "coil_domain": coil_domain,
            "magnetomotive_force_NI": NI,
            "ndof": fesTheta.ndof,
            "grad_theta_norm": float(grad_norm),
            "message": f"Theta field computed with MMF={NI} A·turns, ||∇Θ||={grad_norm:.6e}"
        }

    except Exception as e:
        return {
            "error": str(e),
            "traceback": __import__('traceback').format_exc()
        }


async def _two_scalar_solve_with_jump(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Solve magnetostatic problem with coil current jump using 2-scalar method.

    Incorporates Theta field to handle the Ω jump across coil:
    - True Ωᵣ = Ω_continuous + Θ
    - Solves for Ω_continuous and Ωₜ
    - Reconstructs H = H₀ - ∇(Ω_continuous + Θ) in coil
    - H = -∇Ωₜ in air
    """
    fespace_name = args["fespace_name"]
    theta_field_name = args["theta_field_name"]
    h0_field_name = args.get("h0_field_name", None)
    mu_coil = args.get("mu_coil", 1.2566370614e-6)
    mu_air = args.get("mu_air", 1.2566370614e-6)
    coil_domain = args["coil_domain"]
    air_domain = args["air_domain"]
    tol = args.get("tolerance", 1e-10)
    max_iter = args.get("max_iterations", 1000)

    if fespace_name not in state:
        return {"error": f"FE space '{fespace_name}' not found"}
    if theta_field_name not in state:
        return {"error": f"Theta field '{theta_field_name}' not found"}

    fespace = state[fespace_name]
    gfTheta = state[theta_field_name]
    mesh = gfTheta.space.mesh

    # Get H₀ if provided
    H0 = None
    if h0_field_name and h0_field_name in state:
        H0 = state[h0_field_name]

    try:
        (Omega_r, Omega_t), (psi_r, psi_t) = fespace.TnT()

        # Bilinear form
        a = BilinearForm(fespace)
        # Coil region: ∇·μ(H₀ - ∇Ωᵣ - ∇Θ) = 0
        # Since Θ is known, this becomes: ∇·μ(H₀ - ∇Ωᵣ) = ∇·μ(∇Θ)
        a += mu_coil * grad(Omega_r) * grad(psi_r) * dx(coil_domain)

        # Air region: ∇·μ(-∇Ωₜ) = 0
        a += mu_air * grad(Omega_t) * grad(psi_t) * dx(air_domain)

        # Linear form
        f = LinearForm(fespace)

        # Source term from H₀ (if exists)
        if H0 is not None:
            f += mu_coil * H0 * grad(psi_r) * dx(coil_domain)

        # Source term from Theta gradient
        # ∇·μ(∇Θ) in weak form: -∫ μ∇Θ·∇ψ dx
        f += -mu_coil * grad(gfTheta) * grad(psi_r) * dx(coil_domain)

        with TaskManager():
            a.Assemble()
            f.Assemble()

        # Solve
        gf = GridFunction(fespace)
        inv = CGSolver(a.mat, fespace.FreeDofs(), maxsteps=max_iter, tol=tol)
        gf.vec.data = inv * f.vec

        # Extract components
        gf_Omega_r_cont, gf_Omega_t = gf.components

        # Reconstruct full Ωᵣ = Ω_continuous + Θ
        # Create GridFunction on coil domain
        fesR_full = H1(mesh, order=gfTheta.space.globalorder, definedon=coil_domain)
        gf_Omega_r_full = GridFunction(fesR_full)
        gf_Omega_r_full.Set(gf_Omega_r_cont + gfTheta, VOL, definedon=coil_domain)

        # Reconstruct H and B fields
        if H0 is not None:
            H_coil = H0 - grad(gf_Omega_r_full)
        else:
            H_coil = -grad(gf_Omega_r_full)

        H_air = -grad(gf_Omega_t)

        B_coil = mu_coil * H_coil
        B_air = mu_air * H_air

        # Store results
        solution_name = f"{fespace_name}_solution"
        state[solution_name] = gf
        state[f"{solution_name}_Omega_r_full"] = gf_Omega_r_full
        state[f"{solution_name}_Omega_t"] = gf_Omega_t
        state[f"{solution_name}_H_coil"] = H_coil
        state[f"{solution_name}_H_air"] = H_air
        state[f"{solution_name}_B_coil"] = B_coil
        state[f"{solution_name}_B_air"] = B_air

        # Compute norms for verification
        B_coil_norm = sqrt(Integrate(B_coil * B_coil * dx(coil_domain), mesh))
        B_air_norm = sqrt(Integrate(B_air * B_air * dx(air_domain), mesh))

        # Compute energy
        energy_coil = 0.5 / mu_coil * Integrate(B_coil * B_coil * dx(coil_domain), mesh)
        energy_air = 0.5 / mu_air * Integrate(B_air * B_air * dx(air_domain), mesh)
        total_energy = energy_coil + energy_air

        return {
            "success": True,
            "solution_name": solution_name,
            "ndof": fespace.ndof,
            "B_coil_norm": float(B_coil_norm),
            "B_air_norm": float(B_air_norm),
            "magnetic_energy_coil": float(energy_coil),
            "magnetic_energy_air": float(energy_air),
            "total_magnetic_energy": float(total_energy),
            "stored_fields": {
                "Omega_r_full": f"{solution_name}_Omega_r_full",
                "Omega_t": f"{solution_name}_Omega_t",
                "H_coil": f"{solution_name}_H_coil",
                "H_air": f"{solution_name}_H_air",
                "B_coil": f"{solution_name}_B_coil",
                "B_air": f"{solution_name}_B_air"
            },
            "message": f"2-scalar with coil jump solved: ||B_coil||={B_coil_norm:.6e}, ||B_air||={B_air_norm:.6e}, Energy={total_energy:.6e} J"
        }

    except Exception as e:
        return {
            "error": str(e),
            "traceback": __import__('traceback').format_exc()
        }
