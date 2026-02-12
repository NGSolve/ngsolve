#!/usr/bin/env python3
"""
NGSolve MCP Server

Model Context Protocol server for NGSolve FEM with Radia coupling.
Provides mesh generation, FEM analysis, and Radia field integration.

Features:
- Mesh generation (Netgen)
- Mesh import (GMSH)
- Radia field coupling via shared workspace
- FEM analysis with Radia source fields
"""

import sys
import json
import logging
from typing import Any, Dict, List, Optional, Sequence
from pathlib import Path

try:
    from mcp.server import Server
    from mcp.server.stdio import stdio_server
    from mcp.types import Tool, TextContent
except ImportError:
    print("Error: MCP SDK not installed. Install with: pip install mcp", file=sys.stderr)
    sys.exit(1)

# Import tool modules
from .tools import (
    mesh_tools,
    radia_coupling_tools,
    kelvin_transform_tools,
    diagnostic_tools,
    eddy_current_tools,
    two_scalar_tools,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("ngsolve-mcp-server")


class NGSolveMCPServer:
    """NGSolve MCP Server implementation with Radia coupling."""

    def __init__(self):
        """Initialize the NGSolve MCP server."""
        self.server = Server("ngsolve-mcp-server")
        self.ngsolve_state: Dict[str, Any] = {}

        # Register tool handlers
        self._register_handlers()

        logger.info("NGSolve MCP Server initialized")

    def _register_handlers(self):
        """Register all tool handlers."""

        @self.server.list_tools()
        async def list_tools() -> List[Tool]:
            """List all available NGSolve tools."""
            tools = []

            # Mesh generation tools
            tools.extend(mesh_tools.get_tools())

            # Radia coupling tools
            tools.extend(radia_coupling_tools.get_tools())

            # Kelvin transformation tools
            tools.extend(kelvin_transform_tools.get_tools())

            # Diagnostic tools
            tools.extend(diagnostic_tools.get_tools())

            # Eddy current analysis tools
            tools.extend(eddy_current_tools.get_tools())

            # Two-scalar potential method tools
            tools.extend(two_scalar_tools.get_tools())

            logger.info(f"Listing {len(tools)} available tools")
            return tools

        @self.server.call_tool()
        async def call_tool(name: str, arguments: Dict[str, Any]) -> Sequence[TextContent]:
            """Execute an NGSolve tool."""
            logger.info(f"Calling tool: {name} with arguments: {arguments}")

            try:
                # Dispatch to appropriate tool handler
                if name.startswith("ngsolve_mesh_") or name == "ngsolve_check_availability":
                    result = await mesh_tools.execute(name, arguments, self.ngsolve_state)
                elif name.startswith("ngsolve_radia_") or name.startswith("ngsolve_workspace_"):
                    result = await radia_coupling_tools.execute(name, arguments, self.ngsolve_state)
                elif name.startswith("kelvin_"):
                    result = await kelvin_transform_tools.execute(name, arguments, self.ngsolve_state)
                elif name.startswith("ngsolve_server_") or name.startswith("ngsolve_list_") or name.startswith("ngsolve_get_") or name.startswith("ngsolve_clear_"):
                    result = await diagnostic_tools.execute(name, arguments, self.ngsolve_state)
                # Two-scalar tools (check before generic ngsolve_compute_)
                elif name.startswith("ngsolve_two_scalar_") or name.startswith("ngsolve_h_to_") or name.startswith("ngsolve_compute_h0_") or name == "ngsolve_compute_theta_field":
                    result = await two_scalar_tools.execute(name, arguments, self.ngsolve_state)
                # Eddy current tools (generic ngsolve_compute_ after two-scalar specific ones)
                elif name.startswith("ngsolve_compute_") or name.startswith("ngsolve_t_omega_") or name.startswith("ngsolve_loop_"):
                    result = await eddy_current_tools.execute(name, arguments, self.ngsolve_state)
                else:
                    raise ValueError(f"Unknown tool: {name}")

                # Format result
                return [TextContent(
                    type="text",
                    text=json.dumps(result, indent=2)
                )]

            except Exception as e:
                logger.error(f"Error executing tool {name}: {str(e)}", exc_info=True)
                return [TextContent(
                    type="text",
                    text=json.dumps({
                        "error": str(e),
                        "tool": name,
                        "arguments": arguments
                    }, indent=2)
                )]

    async def run(self):
        """Run the MCP server."""
        logger.info("Starting NGSolve MCP Server")
        async with stdio_server() as (read_stream, write_stream):
            await self.server.run(
                read_stream,
                write_stream,
                self.server.create_initialization_options()
            )


async def main():
    """Main entry point for the NGSolve MCP server."""
    server = NGSolveMCPServer()
    await server.run()


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())
