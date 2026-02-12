"""
Mesh Generation Tools for NGSolve MCP Server

Tools for mesh generation and import:
- Netgen mesh generation (box, cylinder, etc.)
- GMSH mesh import
- Mesh statistics and info
"""

import sys
from pathlib import Path
from typing import Any, Dict, List

try:
    from ngsolve import Mesh
    from netgen.occ import Box, Cylinder, Pnt, OCCGeometry, Axis
    from netgen.meshing import MeshingParameters
    NGSOLVE_AVAILABLE = True
except ImportError:
    NGSOLVE_AVAILABLE = False

from mcp.types import Tool


def get_tools() -> List[Tool]:
    """Get list of mesh generation tools."""
    if not NGSOLVE_AVAILABLE:
        return [
            Tool(
                name="ngsolve_check_availability",
                description="Check if NGSolve is available. Install with: pip install ngsolve==6.2.2405",
                inputSchema={"type": "object", "properties": {}, "required": []}
            )
        ]

    return [
        Tool(
            name="ngsolve_mesh_create_box",
            description="Create a box mesh using Netgen.",
            inputSchema={
                "type": "object",
                "properties": {
                    "pmin": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Minimum point [x, y, z] in meters"
                    },
                    "pmax": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Maximum point [x, y, z] in meters"
                    },
                    "maxh": {
                        "type": "number",
                        "description": "Maximum element size",
                        "default": 0.05
                    },
                    "material_name": {
                        "type": "string",
                        "description": "Material name for the volume",
                        "default": "iron"
                    },
                    "mesh_name": {
                        "type": "string",
                        "description": "Name for the mesh object",
                        "default": "box_mesh"
                    }
                },
                "required": ["pmin", "pmax"]
            }
        ),
        Tool(
            name="ngsolve_mesh_create_cylinder",
            description="Create a cylindrical mesh using Netgen.",
            inputSchema={
                "type": "object",
                "properties": {
                    "center": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Center point [x, y, z]"
                    },
                    "axis": {
                        "type": "array",
                        "items": {"type": "number"},
                        "minItems": 3,
                        "maxItems": 3,
                        "description": "Axis direction [dx, dy, dz]",
                        "default": [0, 0, 1]
                    },
                    "radius": {
                        "type": "number",
                        "description": "Cylinder radius"
                    },
                    "height": {
                        "type": "number",
                        "description": "Cylinder height"
                    },
                    "maxh": {
                        "type": "number",
                        "description": "Maximum element size",
                        "default": 0.05
                    },
                    "material_name": {
                        "type": "string",
                        "description": "Material name",
                        "default": "iron"
                    },
                    "mesh_name": {
                        "type": "string",
                        "description": "Name for the mesh",
                        "default": "cylinder_mesh"
                    }
                },
                "required": ["center", "radius", "height"]
            }
        ),
        Tool(
            name="ngsolve_mesh_import_file",
            description="Import a mesh from GMSH file (.msh format).",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_file": {
                        "type": "string",
                        "description": "Path to mesh file (.msh, .vol, .vol.gz)"
                    },
                    "mesh_name": {
                        "type": "string",
                        "description": "Name for the imported mesh",
                        "default": "imported_mesh"
                    }
                },
                "required": ["mesh_file"]
            }
        ),
        Tool(
            name="ngsolve_mesh_get_info",
            description="Get mesh statistics and information.",
            inputSchema={
                "type": "object",
                "properties": {
                    "mesh_name": {
                        "type": "string",
                        "description": "Name of the mesh object"
                    }
                },
                "required": ["mesh_name"]
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute a mesh tool."""
    if not NGSOLVE_AVAILABLE and name != "ngsolve_check_availability":
        return {
            "error": "NGSolve not available",
            "message": "Install with: pip install ngsolve==6.2.2405"
        }

    try:
        if name == "ngsolve_check_availability":
            return _check_availability()
        elif name == "ngsolve_mesh_create_box":
            return _create_box_mesh(arguments, state)
        elif name == "ngsolve_mesh_create_cylinder":
            return _create_cylinder_mesh(arguments, state)
        elif name == "ngsolve_mesh_import_file":
            return _import_mesh(arguments, state)
        elif name == "ngsolve_mesh_get_info":
            return _get_mesh_info(arguments, state)
        else:
            return {"error": f"Unknown mesh tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name}


def _check_availability() -> Dict[str, Any]:
    """Check if NGSolve is available."""
    return {
        "ngsolve_available": NGSOLVE_AVAILABLE,
        "message": "NGSolve is available" if NGSOLVE_AVAILABLE else "Install with: pip install ngsolve==6.2.2405"
    }


def _create_box_mesh(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a box mesh using Netgen."""
    pmin = args["pmin"]
    pmax = args["pmax"]
    maxh = args.get("maxh", 0.05)
    material_name = args.get("material_name", "iron")
    mesh_name = args.get("mesh_name", "box_mesh")

    # Create box geometry
    p1 = Pnt(pmin[0], pmin[1], pmin[2])
    p2 = Pnt(pmax[0], pmax[1], pmax[2])
    box = Box(p1, p2)
    box.mat(material_name)

    # Generate mesh
    geo = OCCGeometry(box)
    mp = MeshingParameters(maxh=maxh)
    mesh = Mesh(geo.GenerateMesh(mp))

    # Store in state
    state[mesh_name] = mesh

    return {
        "success": True,
        "mesh_name": mesh_name,
        "mesh_info": {
            "num_vertices": mesh.nv,
            "num_elements": mesh.ne,
            "dimension": mesh.dim,
            "maxh": maxh,
            "material": material_name,
            "bounds": {"pmin": pmin, "pmax": pmax}
        }
    }


def _create_cylinder_mesh(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create a cylinder mesh using Netgen."""
    center = args["center"]
    axis_dir = args.get("axis", [0, 0, 1])
    radius = args["radius"]
    height = args["height"]
    maxh = args.get("maxh", 0.05)
    material_name = args.get("material_name", "iron")
    mesh_name = args.get("mesh_name", "cylinder_mesh")

    # Create cylinder geometry
    cyl_axis = Axis(Pnt(*center), Pnt(center[0] + axis_dir[0],
                                       center[1] + axis_dir[1],
                                       center[2] + axis_dir[2]))
    cylinder = Cylinder(cyl_axis, radius, height)
    cylinder.mat(material_name)

    # Generate mesh
    geo = OCCGeometry(cylinder)
    mp = MeshingParameters(maxh=maxh)
    mesh = Mesh(geo.GenerateMesh(mp))

    # Store in state
    state[mesh_name] = mesh

    return {
        "success": True,
        "mesh_name": mesh_name,
        "mesh_info": {
            "num_vertices": mesh.nv,
            "num_elements": mesh.ne,
            "dimension": mesh.dim,
            "maxh": maxh,
            "material": material_name,
            "geometry": {
                "center": center,
                "radius": radius,
                "height": height,
                "axis": axis_dir
            }
        }
    }


def _import_mesh(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Import a mesh from file."""
    mesh_file = args["mesh_file"]
    mesh_name = args.get("mesh_name", "imported_mesh")

    # Load mesh
    mesh = Mesh(mesh_file)

    # Store in state
    state[mesh_name] = mesh

    return {
        "success": True,
        "mesh_name": mesh_name,
        "mesh_file": mesh_file,
        "mesh_info": {
            "num_vertices": mesh.nv,
            "num_elements": mesh.ne,
            "dimension": mesh.dim
        }
    }


def _get_mesh_info(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Get mesh information."""
    mesh_name = args["mesh_name"]

    if mesh_name not in state:
        return {"error": f"Mesh '{mesh_name}' not found in state"}

    mesh = state[mesh_name]

    return {
        "success": True,
        "mesh_name": mesh_name,
        "mesh_info": {
            "num_vertices": mesh.nv,
            "num_elements": mesh.ne,
            "num_faces": mesh.nface if hasattr(mesh, 'nface') else None,
            "num_edges": mesh.nedge if hasattr(mesh, 'nedge') else None,
            "dimension": mesh.dim,
            "materials": list(mesh.GetMaterials()) if hasattr(mesh, 'GetMaterials') else []
        }
    }
