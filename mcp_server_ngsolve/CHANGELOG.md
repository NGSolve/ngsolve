# Changelog

All notable changes to the NGSolve MCP Server will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-02-12

### Added
- **Diagnostic and Debugging Tools** (4 new tools)
  - `ngsolve_server_info`: Get server version, NGSolve/Netgen availability, and state summary
  - `ngsolve_list_objects`: List all objects (meshes, GridFunctions, FE spaces) in server state
  - `ngsolve_get_object_info`: Get detailed information about specific objects (mesh statistics, DOF counts)
  - `ngsolve_clear_state`: Reset server state (clear all objects)

- **Documentation Enhancements**
  - Best Practices section:
    * Units and coordinate systems guidance
    * Kelvin transformation guidelines with radius selection table
    * Mesh resolution recommendations by region
    * Solver selection guidance (direct vs iterative)
    * Performance optimization strategies
  - Technical Background section:
    * Kelvin transformation mathematics (Ω-Reduced Ω method)
    * Permeability transformation formula: μ'(r') = (R/r')² μ₀
    * Radia-NGSolve coupling architecture
    * Integration methods (interpolation, magnetization import, hybrid solver)
  - Complete workflow examples:
    * Kelvin transformation usage with complete setup
    * End-to-end Radia→NGSolve integration pipeline
  - Development policy and branch management guidelines

### Improved
- README structure with comprehensive examples
- Server architecture documentation
- Error handling with detailed object information
- State inspection capabilities with FEM-specific details

### Technical Details
- Total tools: 19 (up from 15)
- Enhanced diagnostic tools with NGSolve-specific information (mesh stats, FE space details)
- Comprehensive documentation for Kelvin transformation parameters
- Performance guidance for large-scale FEM problems

## [1.0.0] - 2026-02-11

### Added
- Initial release of NGSolve MCP Server
- **Mesh Generation Tools** (4 tools)
  - Box and cylinder mesh creation via Netgen
  - GMSH mesh file import
  - Mesh statistics and information

- **Radia Coupling Tools** (4 tools)
  - Import Radia objects from shared workspace
  - Access pre-computed field data
  - Create interpolated CoefficientFunction from Radia fields
  - List available Radia objects in workspace

- **Kelvin Transformation Tools** (7 tools)
  - Mesh creation with Kelvin transformation (sphere, cylinder, box)
  - Ω-Reduced Ω method solver for unbounded domains
  - Perturbation field energy computation
  - VTK export for visualization
  - Analytical solution comparison (sphere geometry)
  - Adaptive mesh refinement (planned)
  - NGSolve availability check

- **Infrastructure**
  - MCP protocol implementation with stdio transport
  - Asynchronous tool execution
  - State management for NGSolve objects (meshes, GridFunctions, FE spaces)
  - Shared workspace integration via symbolic link
  - Comprehensive error handling and logging

### Documentation
- Installation instructions
- Tool reference documentation
- Usage examples for Claude Desktop integration
- Mesh generation and Radia coupling workflows
- Shared workspace configuration

## Release Notes

### Version 1.1.0 Highlights

This release significantly enhances the NGSolve MCP server with:

1. **Comprehensive Documentation**: Added detailed best practices, technical background, and complete workflow examples. Documentation now covers Kelvin transformation theory, mesh resolution guidelines, and performance optimization strategies.

2. **Diagnostic Tools**: New debugging and monitoring capabilities with NGSolve-specific information (mesh statistics, FE space DOF counts, object type inspection).

3. **Best Practices**: Guidelines for:
   - Choosing optimal Kelvin radius (R = 2-3× inner domain radius)
   - Mesh resolution by region (inner domain, Kelvin boundary, outer domain)
   - Solver selection based on problem size (direct for <50k DOFs, iterative for larger)
   - Radia H-matrix acceleration for batch field evaluation

4. **Technical Documentation**: Mathematical formulation of Kelvin transformation, permeability transformation, and Radia-NGSolve coupling architecture.

### Migration Notes

No breaking changes. All existing tools maintain backward compatibility.

### Known Issues

None reported.

### Future Enhancements

Planned for version 1.2.0:
- Adaptive mesh refinement implementation (currently planned)
- H-formulation comparison workflows
- Magnetization import from NGSolve solutions

### Contributors

- Research Lab Team
- Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
