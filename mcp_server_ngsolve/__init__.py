"""
NGSolve MCP Server

Model Context Protocol server for NGSolve FEM with Radia coupling.
Provides mesh generation, FEM analysis, and Radia field integration.
"""

__version__ = "1.0.0"

from .server import NGSolveMCPServer, main

__all__ = ["NGSolveMCPServer", "main"]
