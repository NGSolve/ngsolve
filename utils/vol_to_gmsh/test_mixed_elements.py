#!/usr/bin/env python
"""
Mixed element mesh validation test

Specifically tests coilshield.vol which contains:
- Wedge elements (6 nodes)
- Tetrahedron elements (4 nodes)
- Pyramid elements (5 nodes)

Validates that all element types are correctly converted.
"""

import os
import sys
import subprocess
import re


def validate_mixed_element_conversion():
    """Validate coilshield.vol mixed element conversion."""
    print("="*60)
    print("MIXED ELEMENT VALIDATION TEST")
    print("="*60)

    vol_file = "install_ksugahar/share/ngsolve/coilshield.vol"
    msh_file = "test_mixed_coilshield.msh"

    if not os.path.exists(vol_file):
        print(f"[SKIP] Test file not found: {vol_file}")
        return None

    # Analyze original mesh
    print(f"\n[1] Analyzing original mesh...")
    try:
        from netgen.meshing import Mesh as NetgenMesh

        ngmesh = NetgenMesh()
        ngmesh.Load(vol_file)

        # Count element types in original
        elem_counts = {}
        elem_type_names = {4: 'Tet', 5: 'Pyramid', 6: 'Wedge', 8: 'Hex'}

        for el in ngmesh.Elements3D():
            el_type = len(el.vertices)
            elem_counts[el_type] = elem_counts.get(el_type, 0) + 1

        print("    Original mesh element types:")
        for el_type, count in sorted(elem_counts.items()):
            name = elem_type_names.get(el_type, f'Unknown({el_type})')
            print(f"      {name} ({el_type} nodes): {count}")

        expected_elem_types = elem_counts.copy()

    except ImportError:
        print("[SKIP] Netgen not available")
        return None
    except Exception as e:
        print(f"[FAIL] Error analyzing mesh: {e}")
        return False

    # Convert to Gmsh
    print(f"\n[2] Converting to Gmsh format...")
    result = subprocess.run(
        ["python", "vol_to_gmsh.py", vol_file, msh_file],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"[FAIL] Conversion failed: {result.stderr}")
        return False

    # Analyze Gmsh file
    print(f"\n[3] Analyzing Gmsh file element types...")

    with open(msh_file, 'r') as f:
        content = f.read()

    # Extract elements section
    elements_match = re.search(r'\$Elements\n(\d+)\n(.*?)\$EndElements', content, re.DOTALL)
    if not elements_match:
        print("[FAIL] Cannot parse Elements section")
        os.remove(msh_file)
        return False

    elements_text = elements_match.group(2)

    # Gmsh element type IDs:
    # 1 = Line (2 nodes)
    # 2 = Triangle (3 nodes)
    # 3 = Quadrangle (4 nodes)
    # 4 = Tetrahedron (4 nodes)
    # 5 = Hexahedron (8 nodes)
    # 6 = Wedge/Prism (6 nodes)
    # 7 = Pyramid (5 nodes)
    # 8 = Line (2 nodes, 2nd order)

    gmsh_elem_types = {}
    gmsh_type_names = {
        1: 'Line', 2: 'Triangle', 3: 'Quad',
        4: 'Tet', 5: 'Hex', 6: 'Wedge', 7: 'Pyramid'
    }

    for line in elements_text.strip().split('\n'):
        parts = line.split()
        if len(parts) >= 2:
            elem_type = int(parts[1])
            gmsh_elem_types[elem_type] = gmsh_elem_types.get(elem_type, 0) + 1

    print("    Gmsh file element types:")
    for elem_type, count in sorted(gmsh_elem_types.items()):
        name = gmsh_type_names.get(elem_type, f'Unknown({elem_type})')
        print(f"      Type {elem_type} ({name}): {count}")

    # Validate conversion
    print(f"\n[4] Validating element type conversion...")

    # Map Netgen element types to Gmsh types
    # Netgen: 4 nodes → Gmsh: type 4 (Tet)
    # Netgen: 5 nodes → Gmsh: type 7 (Pyramid)
    # Netgen: 6 nodes → Gmsh: type 6 (Wedge)
    # Netgen: 8 nodes → Gmsh: type 5 (Hex)

    netgen_to_gmsh = {4: 4, 5: 7, 6: 6, 8: 5}

    validation_passed = True

    for netgen_type, expected_count in expected_elem_types.items():
        gmsh_type = netgen_to_gmsh[netgen_type]
        actual_count = gmsh_elem_types.get(gmsh_type, 0)

        netgen_name = elem_type_names[netgen_type]
        gmsh_name = gmsh_type_names[gmsh_type]

        if actual_count == expected_count:
            print(f"[PASS] {netgen_name}: {expected_count} elements correctly converted to Gmsh type {gmsh_type} ({gmsh_name})")
        else:
            print(f"[FAIL] {netgen_name}: Expected {expected_count}, got {actual_count} in Gmsh")
            validation_passed = False

    # Check for surface elements (should also be present)
    surface_types = {2, 3}  # Triangle, Quad
    surface_count = sum(gmsh_elem_types.get(t, 0) for t in surface_types)

    if surface_count > 0:
        print(f"[INFO] Surface elements: {surface_count} (Triangle + Quad)")
    else:
        print(f"[WARN] No surface elements found")

    # Cleanup
    os.remove(msh_file)

    if validation_passed:
        print(f"\n[PASS] All element types correctly converted!")
        return True
    else:
        print(f"\n[FAIL] Some element types not converted correctly")
        return False


def test_element_type_preservation():
    """Test that all element types are preserved in order."""
    print("\n" + "="*60)
    print("ELEMENT ORDER PRESERVATION TEST")
    print("="*60)

    vol_file = "install_ksugahar/share/ngsolve/coilshield.vol"
    msh_file = "test_order_coilshield.msh"

    if not os.path.exists(vol_file):
        print(f"[SKIP] Test file not found")
        return None

    # Convert
    result = subprocess.run(
        ["python", "vol_to_gmsh.py", vol_file, msh_file],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"[FAIL] Conversion failed")
        return False

    # Read Gmsh file and check element ordering
    with open(msh_file, 'r') as f:
        content = f.read()

    # Extract elements
    elements_match = re.search(r'\$Elements\n(\d+)\n(.*?)\$EndElements', content, re.DOTALL)
    elements_text = elements_match.group(2)

    # Parse elements
    elements = []
    for line in elements_text.strip().split('\n'):
        parts = line.split()
        if len(parts) >= 2:
            elem_id = int(parts[0])
            elem_type = int(parts[1])
            elements.append((elem_id, elem_type))

    # Check sequential numbering
    expected_id = 1
    all_sequential = True
    for elem_id, elem_type in elements:
        if elem_id != expected_id:
            print(f"[WARN] Non-sequential element ID: expected {expected_id}, got {elem_id}")
            all_sequential = False
            break
        expected_id += 1

    if all_sequential:
        print(f"[PASS] Element IDs are sequential (1 to {len(elements)})")
    else:
        print(f"[FAIL] Element IDs not sequential")

    # Check element type distribution
    type_sequence = [t for _, t in elements[:50]]  # First 50 elements
    print(f"[INFO] First 50 element types: {type_sequence[:20]}...")

    os.remove(msh_file)
    return all_sequential


def main():
    """Run mixed element validation tests."""
    results = {}

    results['Mixed element conversion'] = validate_mixed_element_conversion()
    results['Element order preservation'] = test_element_type_preservation()

    # Summary
    print("\n" + "="*60)
    print("MIXED ELEMENT TEST SUMMARY")
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
        print("\n[SUCCESS] All mixed element tests passed!")
        return 0
    else:
        print(f"\n[FAILURE] {failed} test(s) failed")
        return 1


if __name__ == '__main__':
    sys.exit(main())
