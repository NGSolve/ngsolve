"""
Diagnostic and Debugging Tools for NGSolve MCP Server

Tools for server health check, state inspection, and debugging.
"""

from typing import Any, Dict, List
from mcp.types import Tool

__version__ = "1.0.0"


def get_tools() -> List[Tool]:
    """Get list of diagnostic tools."""
    return [
        Tool(
            name="ngsolve_server_info",
            description="Get server version and status information",
            inputSchema={
                "type": "object",
                "properties": {},
                "required": []
            }
        ),
        Tool(
            name="ngsolve_list_objects",
            description="List all NGSolve objects currently in server state (meshes, GridFunctions, etc.)",
            inputSchema={
                "type": "object",
                "properties": {},
                "required": []
            }
        ),
        Tool(
            name="ngsolve_get_object_info",
            description="Get detailed information about a specific NGSolve object",
            inputSchema={
                "type": "object",
                "properties": {
                    "object_name": {
                        "type": "string",
                        "description": "Name of the NGSolve object"
                    }
                },
                "required": ["object_name"]
            }
        ),
        Tool(
            name="ngsolve_clear_state",
            description="Clear all objects from server state (reset)",
            inputSchema={
                "type": "object",
                "properties": {
                    "confirm": {
                        "type": "boolean",
                        "description": "Must be true to confirm clearing all state",
                        "default": False
                    }
                },
                "required": ["confirm"]
            }
        ),
    ]


async def execute(name: str, arguments: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Execute a diagnostic tool."""
    try:
        if name == "ngsolve_server_info":
            return _server_info(state)
        elif name == "ngsolve_list_objects":
            return _list_objects(state)
        elif name == "ngsolve_get_object_info":
            return _get_object_info(arguments, state)
        elif name == "ngsolve_clear_state":
            return _clear_state(arguments, state)
        else:
            return {"error": f"Unknown diagnostic tool: {name}"}
    except Exception as e:
        return {"error": str(e), "tool": name, "traceback": __import__('traceback').format_exc()}


def _server_info(state: Dict[str, Any]) -> Dict[str, Any]:
    """Get server information."""
    try:
        import ngsolve
        ngsolve_available = True
        ngsolve_version = ngsolve.__version__
    except ImportError:
        ngsolve_available = False
        ngsolve_version = None

    try:
        import netgen
        netgen_available = True
        netgen_version = netgen.__version__
    except ImportError:
        netgen_available = False
        netgen_version = None

    return {
        "success": True,
        "server": "NGSolve MCP Server",
        "version": __version__,
        "ngsolve_available": ngsolve_available,
        "ngsolve_version": ngsolve_version,
        "netgen_available": netgen_available,
        "netgen_version": netgen_version,
        "state_objects": len(state),
        "object_names": list(state.keys())
    }


def _list_objects(state: Dict[str, Any]) -> Dict[str, Any]:
    """List all objects in state."""
    objects = []

    for name, obj in state.items():
        obj_info = {
            "name": name,
            "type": type(obj).__name__
        }

        # Add NGSolve-specific info
        try:
            from ngsolve import Mesh, GridFunction, FESpace

            if isinstance(obj, Mesh):
                obj_info["mesh_info"] = {
                    "nv": obj.nv,
                    "ne": obj.ne,
                    "dim": obj.dim
                }
            elif isinstance(obj, GridFunction):
                obj_info["gf_info"] = {
                    "space": type(obj.space).__name__,
                    "ndof": obj.space.ndof
                }
            elif isinstance(obj, FESpace):
                obj_info["space_info"] = {
                    "ndof": obj.ndof,
                    "order": getattr(obj, 'order', 'N/A')
                }
        except ImportError:
            pass

        objects.append(obj_info)

    return {
        "success": True,
        "total_objects": len(objects),
        "objects": objects
    }


def _get_object_info(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Get detailed information about a specific object."""
    obj_name = args["object_name"]

    if obj_name not in state:
        return {"error": f"Object '{obj_name}' not found in state"}

    obj = state[obj_name]

    info = {
        "success": True,
        "name": obj_name,
        "type": type(obj).__name__,
        "value": str(obj)[:200]  # Truncate long strings
    }

    # Add detailed NGSolve-specific information
    try:
        from ngsolve import Mesh, GridFunction, FESpace, CoefficientFunction

        if isinstance(obj, Mesh):
            info["details"] = {
                "vertices": obj.nv,
                "elements": obj.ne,
                "dimension": obj.dim,
                "boundaries": len(obj.GetBoundaries()) if hasattr(obj, 'GetBoundaries') else 'N/A'
            }
        elif isinstance(obj, GridFunction):
            info["details"] = {
                "space_type": type(obj.space).__name__,
                "ndof": obj.space.ndof,
                "order": getattr(obj.space, 'order', 'N/A'),
                "components": obj.components
            }
        elif isinstance(obj, FESpace):
            info["details"] = {
                "ndof": obj.ndof,
                "order": getattr(obj, 'order', 'N/A'),
                "type": type(obj).__name__
            }
        elif isinstance(obj, CoefficientFunction):
            info["details"] = {
                "dim": obj.dim
            }
    except ImportError:
        pass

    return info


def _clear_state(args: Dict[str, Any], state: Dict[str, Any]) -> Dict[str, Any]:
    """Clear all objects from state."""
    if not args.get("confirm", False):
        return {
            "error": "Must set confirm=true to clear state",
            "warning": "This will delete all objects from server memory"
        }

    # Store count before clearing
    object_count = len(state)
    object_names = list(state.keys())

    # Clear the state dictionary
    state.clear()

    return {
        "success": True,
        "message": "State cleared successfully",
        "objects_removed": object_count,
        "removed_names": object_names
    }
