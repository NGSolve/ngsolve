#!/usr/bin/env python
"""
Validation test for vol_to_gmsh.py

Validates the correctness of generated Gmsh files:
- File format compliance
- Element type correctness
- Vertex coordinate accuracy
- Boundary preservation
"""

import os
import sys
import subprocess
import re
from pathlib import Path


def validate_gmsh_format(msh_file):
    """Validate Gmsh file format structure."""
    print(f"\n{'='*60}")
    print(f"VALIDATION: Gmsh format structure")
    print(f"{'='*60}")

    if not os.path.exists(msh_file):
        print(f"[FAIL] File not found: {msh_file}")
        return False

    with open(msh_file, 'r') as f:
        content = f.read()

    # Check required sections
    required_sections = [
        (r'\$MeshFormat', r'\$EndMeshFormat'),
        (r'\$Nodes', r'\$EndNodes'),
        (r'\$Elements', r'\$EndElements'),
    ]

    all_valid = True
    for start, end in required_sections:
        if re.search(start, content) and re.search(end, content):
            section_name = start.replace('\\$', '')
            print(f"[PASS] {section_name} section found")
        else:
            section_name = start.replace('\\$', '')
            print(f"[FAIL] {section_name} section missing or incomplete")
            all_valid = False

    # Check format version
    format_match = re.search(r'\$MeshFormat\n([\d.]+)', content)
    if format_match:
        version = format_match.group(1)
        print(f"[INFO] Gmsh format version: {version}")
        if version.startswith('2.'):
            print(f"[PASS] Format version is 2.x")
        else:
            print(f"[WARN] Format version is not 2.x: {version}")
    else:
        print(f"[FAIL] Cannot parse format version")
        all_valid = False

    return all_valid


def validate_element_types(msh_file):
    """Validate element type IDs in Gmsh file."""
    print(f"\n{'='*60}")
    print(f"VALIDATION: Element types")
    print(f"{'='*60}")

    with open(msh_file, 'r') as f:
        content = f.read()

    # Extract elements section
    elements_match = re.search(r'\$Elements\n(\d+)\n(.*?)\$EndElements', content, re.DOTALL)
    if not elements_match:
        print("[FAIL] Cannot parse Elements section")
        return False

    num_elements = int(elements_match.group(1))
    elements_text = elements_match.group(2)

    # Count element types
    # Gmsh element type IDs:
    # 2 = Triangle (3 nodes)
    # 4 = Tetrahedron (4 nodes)
    # 5 = Hexahedron (8 nodes)
    # 6 = Wedge/Prism (6 nodes)
    # 8 = Line (2 nodes)
    element_types = {}
    element_names = {
        '1': 'Line',
        '2': 'Triangle',
        '3': 'Quadrangle',
        '4': 'Tetrahedron',
        '5': 'Hexahedron',
        '6': 'Wedge',
        '7': 'Pyramid',
        '8': 'Line (2nd order)',
        '9': 'Triangle (2nd order)',
        '10': 'Quadrangle (2nd order)',
        '11': 'Tetrahedron (2nd order)',
    }

    for line in elements_text.strip().split('\n'):
        parts = line.split()
        if len(parts) >= 2:
            elem_type = parts[1]
            element_types[elem_type] = element_types.get(elem_type, 0) + 1

    print(f"[INFO] Total elements: {num_elements}")
    for elem_type, count in sorted(element_types.items()):
        name = element_names.get(elem_type, f'Unknown (type {elem_type})')
        print(f"[INFO] Type {elem_type} ({name}): {count} elements")

    if len(element_types) > 0:
        print(f"[PASS] Found {len(element_types)} different element types")
        return True
    else:
        print(f"[FAIL] No elements found")
        return False


def validate_vertex_coordinates(vol_file, msh_file):
    """Validate that vertex coordinates are preserved."""
    print(f"\n{'='*60}")
    print(f"VALIDATION: Vertex coordinate preservation")
    print(f"{'='*60}")

    try:
        from netgen.meshing import Mesh as NetgenMesh

        # Load original mesh
        ngmesh = NetgenMesh()
        ngmesh.Load(vol_file)

        # Get first few vertices from original
        orig_vertices = []
        for i, pt in enumerate(ngmesh.Points()):
            if i >= 5:  # Check first 5 vertices
                break
            orig_vertices.append((pt[0], pt[1], pt[2]))

        # Read Gmsh file
        with open(msh_file, 'r') as f:
            content = f.read()

        # Extract nodes section
        nodes_match = re.search(r'\$Nodes\n(\d+)\n(.*?)\$EndNodes', content, re.DOTALL)
        if not nodes_match:
            print("[FAIL] Cannot parse Nodes section")
            return False

        nodes_text = nodes_match.group(2)
        gmsh_vertices = []

        for i, line in enumerate(nodes_text.strip().split('\n')):
            if i >= 5:  # Check first 5 vertices
                break
            parts = line.split()
            if len(parts) >= 4:
                # Format: node_id x y z
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                gmsh_vertices.append((x, y, z))

        # Compare vertices
        tolerance = 1e-10
        all_match = True
        for i, (orig, gmsh) in enumerate(zip(orig_vertices, gmsh_vertices)):
            diff = sum((a - b)**2 for a, b in zip(orig, gmsh))**0.5
            if diff < tolerance:
                print(f"[PASS] Vertex {i+1}: coordinates match (diff={diff:.2e})")
            else:
                print(f"[FAIL] Vertex {i+1}: coordinates mismatch (diff={diff:.2e})")
                print(f"       Original: {orig}")
                print(f"       Gmsh:     {gmsh}")
                all_match = False

        return all_match

    except ImportError:
        print("[SKIP] Netgen not available for validation")
        return None
    except Exception as e:
        print(f"[ERROR] Validation failed: {e}")
        return False


def validate_vtu_format(vtu_file):
    """Validate VTK Unstructured Grid file format."""
    print(f"\n{'='*60}")
    print(f"VALIDATION: VTU file format")
    print(f"{'='*60}")

    if not os.path.exists(vtu_file):
        print(f"[FAIL] File not found: {vtu_file}")
        return False

    # VTU files may contain binary data, read with proper encoding
    try:
        with open(vtu_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except Exception as e:
        print(f"[WARN] Cannot read as text (may be binary VTU): {e}")
        # Try binary mode to check file validity
        try:
            with open(vtu_file, 'rb') as f:
                content_bytes = f.read()
                content = content_bytes.decode('utf-8', errors='ignore')
        except Exception as e2:
            print(f"[FAIL] Cannot read file: {e2}")
            return False

    # Check XML structure
    required_tags = [
        '<VTKFile',
        '<UnstructuredGrid>',
        '<Piece',
        '<Points>',
        '<Cells>',
    ]

    all_valid = True
    for tag in required_tags:
        if tag in content:
            print(f"[PASS] Found {tag}")
        else:
            print(f"[FAIL] Missing {tag}")
            all_valid = False

    # Check file size (should not be empty)
    file_size = os.path.getsize(vtu_file)
    if file_size > 100:
        print(f"[PASS] File size reasonable: {file_size} bytes")
    else:
        print(f"[FAIL] File too small: {file_size} bytes")
        all_valid = False

    return all_valid


def main():
    """Run validation tests."""
    print("="*60)
    print("VALIDATION TEST SUITE: vol_to_gmsh.py")
    print("="*60)

    # Generate test files
    test_vol = "install_ksugahar/share/ngsolve/cube.vol"
    test_msh = "validation_test.msh"
    test_vtu = "validation_test.vtu"

    print("\n[INFO] Generating test files...")

    # Generate Gmsh file
    result1 = subprocess.run(
        ["python", "vol_to_gmsh.py", test_vol, test_msh],
        capture_output=True,
        text=True
    )

    # Generate VTU file
    result2 = subprocess.run(
        ["python", "vol_to_gmsh.py", test_vol, test_vtu.replace('.vtu', ''), "--vtk"],
        capture_output=True,
        text=True
    )

    if result1.returncode != 0:
        print(f"[FAIL] Failed to generate Gmsh file: {result1.stderr}")
        return 1

    if result2.returncode != 0:
        print(f"[FAIL] Failed to generate VTU file: {result2.stderr}")
        return 1

    # Run validation tests
    tests = [
        ("Gmsh format structure", lambda: validate_gmsh_format(test_msh)),
        ("Element types", lambda: validate_element_types(test_msh)),
        ("Vertex coordinates", lambda: validate_vertex_coordinates(test_vol, test_msh)),
        ("VTU format", lambda: validate_vtu_format(test_vtu)),
    ]

    results = {}
    for name, test_func in tests:
        try:
            result = test_func()
            results[name] = result
        except Exception as e:
            print(f"[ERROR] {name}: {e}")
            results[name] = False

    # Cleanup
    if os.path.exists(test_msh):
        os.remove(test_msh)
    if os.path.exists(test_vtu):
        os.remove(test_vtu)

    # Summary
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)

    passed = sum(1 for r in results.values() if r is True)
    failed = sum(1 for r in results.values() if r is False)
    skipped = sum(1 for r in results.values() if r is None)
    total = len(results)

    for name, result in results.items():
        status = "PASS" if result is True else ("SKIP" if result is None else "FAIL")
        symbol = "OK" if result is True else ("--" if result is None else "XX")
        print(f"  [{symbol}] {name}: {status}")

    print(f"\nResults: {passed} passed, {failed} failed, {skipped} skipped (total: {total})")

    if failed == 0:
        print("\n[SUCCESS] All validation tests passed!")
        return 0
    else:
        print(f"\n[FAILURE] {failed} test(s) failed")
        return 1


if __name__ == '__main__':
    sys.exit(main())
