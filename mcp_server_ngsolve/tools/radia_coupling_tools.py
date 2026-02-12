"""
Radia Coupling Tools for NGSolve MCP Server

Tools for importing and coupling Radia fields with NGSolve FEM:
- Import Radia objects from shared workspace
- Create RadiaField CoefficientFunction
- Couple FEM solutions with Radia source fields
"""

import sys
import json
from pathlib import Path
from typing import Any, Dict, List

# Add paths
shared_path = Path(__file__).parent.parent.parent / "mcp_shared"
if shared_path.exists() and str(shared_path) not in sys.path:
    sys.path.insert(0, str(shared_path))

radia_src_path = Path(__file__).parent.parent.parent / "src" / "radia"
if radia_src_path.exists() and str(radia_src_path) not in sys.path:
    sys.path.insert(0, str(radia_src_path))

try:
    from workspace import SharedWorkspace
    WORKSPACE_AVAILABLE = True
except ImportError:
    WORKSPACE_AVAILABLE = False

try:
    from ngsolve import Mesh, HDiv, GridFunction, CoefficientFunction
    NGSOLVE_AVAILABLE = True
except ImportError:
    NGSOLVE_AVAILABLE = False

from mcp.types import Tool


def get_tools() -> List[Tool]:
    """Get list of Radia coupling tools."""
    if not WORKSPACE_AVAILABLE or not NGSOLVE_AVAILABLE:
        return [
            Tool(
                name="ngsolve_radia_check_availability",
                description="Check if Radia coupling is available.",
                inputSchema={"type": "object", "properties": {}, "required": []}
            )
        ]

    return [
        Tool(
            name="ngsolve_radia_import_object",
            description="Import Radia object from shared workspace.",
            inputSchema={
                "type": "object",
                "properties": {
                    "session_id": {
                        "type": "string",
                        "description": "Shared workspace session ID"
                    },
                    "radia_object_name": {
                        "type": "string",
                        "description": "Name of Radia object to import"
                    }
                },
                "required": ["session_id", "radia_object_name"]
            }
        ),
        Tool(
            name="ngsolve_radia_get_field_data",
            description="Get pre-computed Radia field data from workspace.",
            inputSchema={
                "type": "object",
                "properties": {
                    "session_id": {
                        "type": "string",
                        "description": "Session ID"
                    },
                    "radia_object_name": {
                        "type": "string",
                        "description": "Radia object name"
                    }
                },
                "required": ["session_id", "radia_object_name"]
            }
        ),
        Tool(
            name="ngsolve_radia_create_interpolated_field",
            description="Create interpolated CoefficientFunction from Radia field data.",
            inputSchema={
                "type": "object",
                "properties": {
                    "session_id": {
                        "type": "string",
                        "description": "Session ID"
                    },
                    "radia_object_name": {
                        "type": "string",
                        "description": "Radia object name"
                    },
                    "mesh_name": {
                        "type": "string",
                        "description": "NGSolve mesh name for interpolation"
                    }
                },
                "required": ["session_id", "radia_object_name", "mesh_name"]
            }
        ),
        Tool(
            name="ngsolve_workspace_list_radia_objects",
            description="List all available Radia objects in a workspace session.",
            inputSchema={
                "type": "object",
                "properties": {
                    "session_id": {
                        "type": "string",
                        "description": "Session ID"
                    }
                },
                "required": ["session_id"]
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute a Radia coupling tool."""
    if not WORKSPACE_AVAILABLE:
        return {"error": "Shared workspace not available"}

    if not NGSOLVE_AVAILABLE:
        return {"error": "NGSolve not available. Install with: pip install ngsolve==6.2.2405"}

    try:
        if name == "ngsolve_radia_check_availability":
            return _check_availability()
        elif name == "ngsolve_radia_import_object":
            return _import_radia_object(arguments, state)
        elif name == "ngsolve_radia_get_field_data":
            return _get_field_data(arguments, state)
        elif name == "ngsolve_radia_create_interpolated_field":
            return _create_interpolated_field(arguments, state)
        elif name == "ngsolve_workspace_list_radia_objects":
            return _list_radia_objects(arguments, state)
        else:
            return {"error": f"Unknown Radia coupling tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name}


def _check_availability() -> Dict[str, Any]:
    """Check availability of Radia coupling."""
    return {
        "workspace_available": WORKSPACE_AVAILABLE,
        "ngsolve_available": NGSOLVE_AVAILABLE,
        "radia_coupling_ready": WORKSPACE_AVAILABLE and NGSOLVE_AVAILABLE,
        "message": "Radia coupling is ready" if (WORKSPACE_AVAILABLE and NGSOLVE_AVAILABLE) else "Missing dependencies"
    }


def _import_radia_object(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Import Radia object from shared workspace."""
    workspace = SharedWorkspace()
    session_id = args["session_id"]
    radia_obj_name = args["radia_object_name"]

    # Import object metadata
    obj_data = workspace.import_radia_object(session_id, radia_obj_name)

    # Store reference in NGSolve server state
    state_key = f"radia_{radia_obj_name}"
    state[state_key] = obj_data

    return {
        "success": True,
        "session_id": session_id,
        "radia_object_name": radia_obj_name,
        "available_data": {
            "metadata": obj_data["metadata"] is not None,
            "geometry": obj_data["geometry_file"] is not None,
            "fields": obj_data["field_file"] is not None
        },
        "geometry_file": obj_data["geometry_file"],
        "field_file": obj_data["field_file"],
        "message": f"Radia object '{radia_obj_name}' imported from session {session_id}"
    }


def _get_field_data(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Get pre-computed Radia field data."""
    session_id = args["session_id"]
    radia_obj_name = args["radia_object_name"]

    state_key = f"radia_{radia_obj_name}"
    if state_key not in state:
        return {
            "error": f"Radia object '{radia_obj_name}' not imported. Use ngsolve_radia_import_object first."
        }

    obj_data = state[state_key]

    if not obj_data["field_file"]:
        return {
            "error": "No pre-computed field data available. Export fields from Radia server first.",
            "suggestion": "Use radia_workspace_export_object with export_fields=True"
        }

    # Load field data
    import numpy as np
    field_data = np.load(obj_data["field_file"])

    points = field_data["points"].tolist()
    field_values = field_data["field_values"].tolist()
    field_type = str(field_data["field_type"])

    return {
        "success": True,
        "radia_object_name": radia_obj_name,
        "field_type": field_type,
        "num_points": len(points),
        "field_statistics": {
            "min": [float(f) for f in field_data["field_values"].min(axis=0)],
            "max": [float(f) for f in field_data["field_values"].max(axis=0)],
            "mean": [float(f) for f in field_data["field_values"].mean(axis=0)]
        },
        "message": f"Field data loaded: {len(points)} points"
    }


def _create_interpolated_field(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Create interpolated CoefficientFunction from Radia field data."""
    import numpy as np
    from scipy.interpolate import LinearNDInterpolator

    session_id = args["session_id"]
    radia_obj_name = args["radia_object_name"]
    mesh_name = args["mesh_name"]

    # Check if Radia object is imported
    state_key = f"radia_{radia_obj_name}"
    if state_key not in state:
        return {
            "error": f"Radia object '{radia_obj_name}' not imported. Use ngsolve_radia_import_object first."
        }

    # Check if mesh exists
    if mesh_name not in state:
        return {
            "error": f"Mesh '{mesh_name}' not found in state"
        }

    obj_data = state[state_key]

    if not obj_data["field_file"]:
        return {
            "error": "No pre-computed field data available."
        }

    # Load field data
    field_data = np.load(obj_data["field_file"])
    points = field_data["points"]
    field_values = field_data["field_values"]

    # Create interpolator
    interp_x = LinearNDInterpolator(points, field_values[:, 0])
    interp_y = LinearNDInterpolator(points, field_values[:, 1])
    interp_z = LinearNDInterpolator(points, field_values[:, 2])

    # Create CoefficientFunction (simplified - stores interpolator info)
    cf_name = f"{radia_obj_name}_field_cf"
    state[cf_name] = {
        "type": "interpolated_radia_field",
        "radia_object": radia_obj_name,
        "interpolators": {
            "x": interp_x,
            "y": interp_y,
            "z": interp_z
        },
        "num_points": len(points)
    }

    return {
        "success": True,
        "coefficient_function_name": cf_name,
        "radia_object_name": radia_obj_name,
        "mesh_name": mesh_name,
        "num_interpolation_points": len(points),
        "message": f"Interpolated CoefficientFunction created: {cf_name}"
    }


def _list_radia_objects(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """List all Radia objects in a session."""
    workspace = SharedWorkspace()
    session_id = args["session_id"]

    try:
        session_info = workspace.get_session_info(session_id)
        radia_objects = session_info.get("radia_objects", [])

        return {
            "success": True,
            "session_id": session_id,
            "num_objects": len(radia_objects),
            "radia_objects": radia_objects
        }
    except Exception as e:
        return {
            "error": f"Session '{session_id}' not found: {str(e)}"
        }
