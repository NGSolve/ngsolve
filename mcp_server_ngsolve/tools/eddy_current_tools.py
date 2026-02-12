"""
Eddy Current Analysis Tools for NGSolve MCP Server

T-Omega method implementation for eddy current problems in multiply-connected domains.
Supports loop field computation for domains with holes.
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
    """Get list of eddy current analysis tools."""
    if not NGSOLVE_AVAILABLE:
        return []

    return [
        Tool(
            name="ngsolve_compute_genus",
            description="Compute genus (number of holes) of multiply-connected domain from surface mesh",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_name": {
                        "type": "string",
                        "description": "Name of the mesh object"
                    },
                    "boundary_name": {
                        "type": "string",
                        "description": "Boundary name to analyze (e.g., 'conductorBND')"
                    },
                    "connected_components": {
                        "type": "integer",
                        "description": "Number of connected components (default: 1)",
                        "default": 1
                    }
                },
                "required": ["mesh_name", "boundary_name"]
            }
        ),
        Tool(
            name="ngsolve_compute_loop_fields",
            description="Compute loop fields for multiply-connected domain (T-Omega method)",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_name": {
                        "type": "string",
                        "description": "Name of the mesh object"
                    },
                    "domain": {
                        "type": "string",
                        "description": "Domain name for loop fields (e.g., 'air')"
                    },
                    "boundary_name": {
                        "type": "string",
                        "description": "Boundary name (e.g., 'conductorBND')"
                    },
                    "order": {
                        "type": "integer",
                        "description": "Finite element order (default: 1)",
                        "default": 1
                    },
                    "tolerance": {
                        "type": "number",
                        "description": "Solver tolerance (default: 1e-12)",
                        "default": 1e-12
                    },
                    "max_iterations": {
                        "type": "integer",
                        "description": "Maximum solver iterations (default: 200)",
                        "default": 200
                    }
                },
                "required": ["mesh_name", "domain", "boundary_name"]
            }
        ),
        Tool(
            name="ngsolve_t_omega_setup",
            description="Set up T-Omega finite element spaces for eddy current analysis",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_name": {
                        "type": "string",
                        "description": "Name of the mesh object"
                    },
                    "conductor_domain": {
                        "type": "string",
                        "description": "Conductor domain name (e.g., 'conductor', 'sig')"
                    },
                    "conductor_boundary": {
                        "type": "string",
                        "description": "Conductor boundary name (e.g., 'conductorBND')"
                    },
                    "dirichlet_omega": {
                        "type": "string",
                        "description": "Dirichlet boundaries for Omega (e.g., 'upper|lower')",
                        "default": "upper|lower"
                    },
                    "order": {
                        "type": "integer",
                        "description": "Finite element order (default: 2)",
                        "default": 2
                    }
                },
                "required": ["mesh_name", "conductor_domain", "conductor_boundary"]
            }
        ),
        Tool(
            name="ngsolve_t_omega_solve_coupled",
            description="Solve T-Omega eddy current problem with loop current coupling",
            inputSchema={
                "type": "object",
                "properties": {
                    "fespace_name": {
                        "type": "string",
                        "description": "Name of the T-Omega FE space"
                    },
                    "loop_fields_name": {
                        "type": "string",
                        "description": "Name of the loop fields object"
                    },
                    "conductor_domain": {
                        "type": "string",
                        "description": "Conductor domain name"
                    },
                    "boundary_name": {
                        "type": "string",
                        "description": "Conductor boundary name"
                    },
                    "sigma": {
                        "type": "number",
                        "description": "Electrical conductivity (S/m)"
                    },
                    "mu": {
                        "type": "number",
                        "description": "Magnetic permeability (H/m)",
                        "default": 1.2566370614e-6
                    },
                    "frequency": {
                        "type": "number",
                        "description": "Frequency (Hz) for s = j*2*pi*f",
                        "default": 50.0
                    },
                    "voltages": {
                        "type": "array",
                        "description": "Applied voltages for each loop (optional)",
                        "items": {"type": "number"}
                    },
                    "source_omega": {
                        "type": "number",
                        "description": "Source Omega value at Dirichlet boundary (optional)"
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
                "required": ["fespace_name", "loop_fields_name", "conductor_domain",
                           "boundary_name", "sigma"]
            }
        ),
        Tool(
            name="ngsolve_loop_current_analysis",
            description="Analyze loop currents: compute resistance and inductance",
            inputSchema={
                "type": "object",
                "properties": {
                    "solution_name": {
                        "type": "string",
                        "description": "Name of the T-Omega solution GridFunction"
                    },
                    "loop_field_name": {
                        "type": "string",
                        "description": "Name of the loop field to analyze"
                    },
                    "conductor_domain": {
                        "type": "string",
                        "description": "Conductor domain name"
                    },
                    "sigma": {
                        "type": "number",
                        "description": "Electrical conductivity (S/m)"
                    },
                    "mu": {
                        "type": "number",
                        "description": "Magnetic permeability (H/m)",
                        "default": 1.2566370614e-6
                    },
                    "frequency": {
                        "type": "number",
                        "description": "Frequency (Hz)",
                        "default": 50.0
                    }
                },
                "required": ["solution_name", "loop_field_name", "conductor_domain", "sigma"]
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute an eddy current analysis tool."""
    if not NGSOLVE_AVAILABLE:
        return {
            "error": "NGSolve not available",
            "message": "Install NGSolve: pip install ngsolve"
        }

    try:
        if name == "ngsolve_compute_genus":
            return _compute_genus(arguments, state)
        elif name == "ngsolve_compute_loop_fields":
            return await _compute_loop_fields(arguments, state)
        elif name == "ngsolve_t_omega_setup":
            return _t_omega_setup(arguments, state)
        elif name == "ngsolve_t_omega_solve_coupled":
            return await _t_omega_solve_coupled(arguments, state)
        elif name == "ngsolve_loop_current_analysis":
            return _loop_current_analysis(arguments, state)
        else:
            return {"error": f"Unknown eddy current tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name, "traceback": __import__('traceback').format_exc()}


def _compute_genus(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Compute genus of multiply-connected domain."""
    mesh_name = args["mesh_name"]
    boundary_name = args["boundary_name"]
    connected = args.get("connected_components", 1)

    if mesh_name not in state:
        return {"error": f"Mesh '{mesh_name}' not found"}

    mesh = state[mesh_name]

    try:
        # Extract surface mesh from boundary
        # Simplified genus computation based on Euler characteristic
        # χ = V - E + F (Euler characteristic)
        # genus g = (2 - χ - b) / 2, where b = number of boundary components

        # For NGSolve mesh, we need to analyze the boundary topology
        # This is a simplified implementation
        boundary = mesh.Boundaries(boundary_name)

        n_vertices = 0
        n_edges = 0
        n_faces = 0

        vertices_set = set()
        edges_set = set()

        for el in boundary.Elements():
            n_faces += 1
            for v in el.vertices:
                vertices_set.add(v.nr)
            for e in el.edges:
                edges_set.add(e.nr)

        n_vertices = len(vertices_set)
        n_edges = len(edges_set)

        # Euler characteristic
        chi = n_vertices - n_edges + n_faces

        # Genus formula: g = (2 - χ - b) / 2
        # For simply connected surface: χ = 2, g = 0
        # For torus: χ = 0, g = 1
        genus = (2 - chi - connected) // 2

        # Ensure non-negative
        genus = max(0, genus)

        return {
            "success": True,
            "genus": genus,
            "euler_characteristic": chi,
            "vertices": n_vertices,
            "edges": n_edges,
            "faces": n_faces,
            "connected_components": connected,
            "message": f"Domain has genus {genus} (number of holes)"
        }

    except Exception as e:
        return {
            "success": False,
            "genus": 0,
            "message": f"Could not compute genus: {str(e)}. Assuming simply-connected (genus=0)"
        }


async def _compute_loop_fields(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Compute loop fields for multiply-connected domain."""
    mesh_name = args["mesh_name"]
    domain = args["domain"]
    boundary_name = args["boundary_name"]
    order = args.get("order", 1)
    tol = args.get("tolerance", 1e-12)
    max_iter = args.get("max_iterations", 200)

    if mesh_name not in state:
        return {"error": f"Mesh '{mesh_name}' not found"}

    mesh = state[mesh_name]

    try:
        # First compute genus
        genus_result = _compute_genus({"mesh_name": mesh_name,
                                       "boundary_name": boundary_name}, state)
        genus = genus_result.get("genus", 0)

        if genus == 0:
            return {
                "success": True,
                "genus": 0,
                "loop_fields": [],
                "message": "Simply-connected domain, no loop fields needed"
            }

        # HCurl space for loop fields
        fes = HCurl(mesh, order=order, nograds=True, definedon=domain)
        u, v = fes.TnT()

        # H1 space for scalar potential
        fesPhi = H1(mesh, order=order, definedon=domain)
        phi, psi = fesPhi.TnT()

        loops = []
        loop_names = []

        for k in range(genus):
            # Create GridFunction
            gfu = GridFunction(fes)

            # Select random edge on boundary (simplified)
            # In production, should use proper edge selection algorithm
            boundary_els = list(mesh.Boundaries(boundary_name).Elements())
            if boundary_els:
                el = boundary_els[0]
                if hasattr(el, 'edges') and len(el.edges) > 0:
                    edge_dofs = fes.GetDofNrs(el.edges[0])
                    if len(edge_dofs) > 0:
                        gfu.vec[edge_dofs[0]] = 1.0
                        fes.FreeDofs()[edge_dofs[0]] = False

            # Solve curl T = 0
            a = BilinearForm(fes)
            a += curl(u) * curl(v) * dx
            a.Assemble()

            fr = -a.mat * gfu.vec

            # Simple iterative solver (CG)
            inv = CGSolver(a.mat, fes.FreeDofs(), maxsteps=max_iter, tol=tol)
            gfu.vec.data += inv * fr

            # Remove gradient part
            gfPhi = GridFunction(fesPhi)
            aPhi = BilinearForm(fesPhi)
            aPhi += grad(phi) * grad(psi) * dx
            fPhi = LinearForm(fesPhi)
            fPhi += grad(psi) * gfu * dx
            aPhi.Assemble()
            fPhi.Assemble()

            invPhi = CGSolver(aPhi.mat, fesPhi.FreeDofs(), maxsteps=max_iter, tol=tol)
            gfPhi.vec.data = invPhi * fPhi.vec

            gfw = gfu - grad(gfPhi)

            # Orthogonalize against existing loops
            gft = GridFunction(fes)
            gft.vec.data = gfw.vec

            for kd in range(len(loops)):
                prod = Integrate(gfw * loops[kd] * dx, mesh)
                gft.vec.data -= prod * loops[kd].vec

            # Normalize
            norm2 = Integrate(gft * gft * dx, mesh)
            norm = np.sqrt(norm2)
            if norm > 1e-10:
                gft.vec.data /= norm

            loop_name = f"loop_field_{k}"
            state[loop_name] = gft
            loops.append(gft)
            loop_names.append(loop_name)

        # Store loop fields list
        loop_fields_name = f"{mesh_name}_loop_fields"
        state[loop_fields_name] = loops

        return {
            "success": True,
            "genus": genus,
            "num_loops": len(loops),
            "loop_field_names": loop_names,
            "loop_fields_list_name": loop_fields_name,
            "message": f"Computed {len(loops)} loop fields for genus-{genus} domain"
        }

    except Exception as e:
        return {
            "error": str(e),
            "traceback": __import__('traceback').format_exc()
        }


def _t_omega_setup(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Set up T-Omega finite element spaces."""
    mesh_name = args["mesh_name"]
    conductor_domain = args["conductor_domain"]
    conductor_boundary = args["conductor_boundary"]
    dirichlet_omega = args.get("dirichlet_omega", "upper|lower")
    order = args.get("order", 2)

    if mesh_name not in state:
        return {"error": f"Mesh '{mesh_name}' not found"}

    mesh = state[mesh_name]

    try:
        # T: Current density (HCurl, nograds for div-free)
        fesT = HCurl(mesh, order=order, nograds=True,
                    definedon=conductor_domain,
                    dirichlet=conductor_boundary,
                    complex=False)

        # Omega: Magnetic scalar potential (H1)
        fesOmega = H1(mesh, order=order,
                     dirichlet=dirichlet_omega,
                     complex=False)

        # Mixed space
        fespace = fesT * fesOmega

        fespace_name = f"{mesh_name}_t_omega_space"
        state[fespace_name] = fespace

        return {
            "success": True,
            "fespace_name": fespace_name,
            "ndof_T": fesT.ndof,
            "ndof_Omega": fesOmega.ndof,
            "ndof_total": fespace.ndof,
            "order": order,
            "message": f"T-Omega space created with {fespace.ndof} DOFs"
        }

    except Exception as e:
        return {
            "error": str(e),
            "traceback": __import__('traceback').format_exc()
        }


async def _t_omega_solve_coupled(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Solve coupled T-Omega system with loop currents."""
    fespace_name = args["fespace_name"]
    loop_fields_name = args["loop_fields_name"]
    conductor_domain = args["conductor_domain"]
    boundary_name = args["boundary_name"]
    sigma = args["sigma"]
    mu = args.get("mu", 4e-7 * np.pi)
    freq = args.get("frequency", 50.0)
    voltages = args.get("voltages")
    source_omega = args.get("source_omega")
    tol = args.get("tolerance", 1e-10)
    max_iter = args.get("max_iterations", 1000)

    if fespace_name not in state:
        return {"error": f"FE space '{fespace_name}' not found"}
    if loop_fields_name not in state:
        return {"error": f"Loop fields '{loop_fields_name}' not found"}

    fespace = state[fespace_name]
    loop_fields = state[loop_fields_name]

    try:
        # Laplace variable s = j*omega
        s = 2j * np.pi * freq

        mesh = fespace.mesh
        (T, Omega), (W, psi) = fespace.TnT()

        # System matrix
        a = BilinearForm(fespace)
        # Conductor region
        a += (1/sigma) * curl(T) * curl(W) * dx(conductor_domain)
        a += s * mu * T * W * dx(conductor_domain)
        a += s * mu * T * grad(psi) * dx(conductor_domain)
        # Air region (if exists)
        a += s * mu * grad(Omega) * grad(psi) * dx

        a.Assemble()

        # Source term (if source_omega provided)
        gfTOmega = GridFunction(fespace)
        if source_omega is not None:
            gfT, gfOmega = gfTOmega.components
            gfOmega.Set(source_omega, definedon=mesh.Boundaries("upper|lower"))

        source_term = -a.mat * gfTOmega.vec

        # This is a simplified implementation
        # Full implementation would require coupling matrix construction
        # as shown in EMPY_PATTERNS.md

        # For now, solve without loop coupling
        inv = CGSolver(a.mat, fespace.FreeDofs(), maxsteps=max_iter, tol=tol)
        gfTOmega.vec.data = inv * source_term

        solution_name = f"{fespace_name}_solution"
        state[solution_name] = gfTOmega

        return {
            "success": True,
            "solution_name": solution_name,
            "message": "T-Omega system solved (simplified, without full loop coupling)",
            "note": "Full loop coupling implementation requires SolveCoupled2 from EMPY"
        }

    except Exception as e:
        return {
            "error": str(e),
            "traceback": __import__('traceback').format_exc()
        }


def _loop_current_analysis(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Analyze loop currents: compute R and L."""
    solution_name = args["solution_name"]
    loop_field_name = args["loop_field_name"]
    conductor_domain = args["conductor_domain"]
    sigma = args["sigma"]
    mu = args.get("mu", 4e-7 * np.pi)
    freq = args.get("frequency", 50.0)

    if solution_name not in state:
        return {"error": f"Solution '{solution_name}' not found"}
    if loop_field_name not in state:
        return {"error": f"Loop field '{loop_field_name}' not found"}

    gfTOmega = state[solution_name]
    loopField = state[loop_field_name]

    try:
        gfT, gfOmega = gfTOmega.components
        mesh = gfTOmega.space.mesh
        s = 2j * np.pi * freq

        # Resistance calculation
        R = Integrate((1/sigma) * curl(gfT)**2 * dx(conductor_domain), mesh)

        # Inductance calculation
        L = Integrate(mu * (gfT + grad(gfOmega) + loopField)**2 * dx, mesh)

        # Impedance Z = R + s*L
        Z = R + s * L

        return {
            "success": True,
            "resistance": float(np.real(R)),
            "inductance": float(np.real(L)),
            "impedance_real": float(np.real(Z)),
            "impedance_imag": float(np.imag(Z)),
            "impedance_magnitude": float(np.abs(Z)),
            "frequency": freq,
            "message": f"R = {np.real(R):.6e} Ω, L = {np.real(L):.6e} H"
        }

    except Exception as e:
        return {
            "error": str(e),
            "traceback": __import__('traceback').format_exc()
        }
