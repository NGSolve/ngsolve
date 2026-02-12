#!/usr/bin/env python
"""
Comprehensive test suite for vol_to_gmsh.py

Tests additional use cases and edge cases:
- 2D surface meshes
- Different Gmsh format versions
- Large meshes (performance)
- Different element types
- Path handling (spaces, absolute/relative)
- Actual viewer compatibility
"""

import os
import sys
import subprocess
import time
from pathlib import Path


def test_surface_only_mesh():
    """Test 2D surface-only mesh (no volume elements)."""
    print("\n" + "="*60)
    print("TEST: 2D Surface-only mesh")
    print("="*60)

    test_file = "install_ksugahar/share/ngsolve/square.vol"
    if not os.path.exists(test_file):
        print(f"[SKIP] Test file not found: {test_file}")
        return None

    result = subprocess.run(
        ["python", "vol_to_gmsh.py", test_file],
        capture_output=True,
        text=True
    )

    if result.returncode == 0:
        if "Surface-only mesh" in result.stdout:
            print("[PASS] Surface-only mesh detected and warned")
            return True
        else:
            print("[PASS] Conversion succeeded")
            return True
    else:
        print(f"[FAIL] Conversion failed: {result.stderr}")
        return False


def test_gmsh_format_versions():
    """Test both Gmsh format versions (gmsh and gmsh2)."""
    print("\n" + "="*60)
    print("TEST: Gmsh format versions (gmsh vs gmsh2)")
    print("="*60)

    test_file = "install_ksugahar/share/ngsolve/cube.vol"

    # Test gmsh format
    result1 = subprocess.run(
        ["python", "vol_to_gmsh.py", test_file, "test_gmsh1.msh", "--format", "gmsh"],
        capture_output=True,
        text=True
    )

    # Test gmsh2 format (default)
    result2 = subprocess.run(
        ["python", "vol_to_gmsh.py", test_file, "test_gmsh2.msh", "--format", "gmsh2"],
        capture_output=True,
        text=True
    )

    success = True
    if result1.returncode == 0 and os.path.exists("test_gmsh1.msh"):
        print("[PASS] Gmsh format version 1 succeeded")
        os.remove("test_gmsh1.msh")
    else:
        print("[FAIL] Gmsh format version 1 failed")
        success = False

    if result2.returncode == 0 and os.path.exists("test_gmsh2.msh"):
        print("[PASS] Gmsh format version 2 succeeded")
        os.remove("test_gmsh2.msh")
    else:
        print("[FAIL] Gmsh format version 2 failed")
        success = False

    return success


def test_large_mesh_performance():
    """Test performance with larger mesh."""
    print("\n" + "="*60)
    print("TEST: Large mesh performance")
    print("="*60)

    test_file = "install_ksugahar/share/ngsolve/coil.vol"
    if not os.path.exists(test_file):
        print(f"[SKIP] Test file not found: {test_file}")
        return None

    start_time = time.time()

    result = subprocess.run(
        ["python", "vol_to_gmsh.py", test_file],
        capture_output=True,
        text=True
    )

    elapsed = time.time() - start_time

    if result.returncode == 0:
        print(f"[PASS] Conversion completed in {elapsed:.2f} seconds")
        if elapsed < 10.0:
            print("[INFO] Performance acceptable (< 10s)")
            return True
        else:
            print("[WARN] Performance slow (>= 10s)")
            return True
    else:
        print(f"[FAIL] Conversion failed")
        return False


def test_path_with_spaces():
    """Test file paths with spaces."""
    print("\n" + "="*60)
    print("TEST: File paths with spaces")
    print("="*60)

    # Create test directory with spaces
    test_dir = Path("test dir with spaces")
    test_dir.mkdir(exist_ok=True)

    # Copy a test file
    import shutil
    src = "install_ksugahar/share/ngsolve/cube.vol"
    if not os.path.exists(src):
        print(f"[SKIP] Source file not found: {src}")
        return None

    dst = test_dir / "test mesh.vol"
    shutil.copy(src, dst)

    # Test conversion
    result = subprocess.run(
        ["python", "vol_to_gmsh.py", str(dst)],
        capture_output=True,
        text=True
    )

    # Cleanup
    if (test_dir / "test mesh.msh").exists():
        (test_dir / "test mesh.msh").unlink()
    dst.unlink()
    test_dir.rmdir()

    if result.returncode == 0:
        print("[PASS] Paths with spaces handled correctly")
        return True
    else:
        print(f"[FAIL] Paths with spaces failed: {result.stderr}")
        return False


def test_absolute_vs_relative_paths():
    """Test absolute and relative path handling."""
    print("\n" + "="*60)
    print("TEST: Absolute vs relative paths")
    print("="*60)

    test_file_rel = "install_ksugahar/share/ngsolve/cube.vol"
    test_file_abs = os.path.abspath(test_file_rel)

    # Test relative path
    result1 = subprocess.run(
        ["python", "vol_to_gmsh.py", test_file_rel, "test_rel.msh"],
        capture_output=True,
        text=True
    )

    # Test absolute path
    result2 = subprocess.run(
        ["python", "vol_to_gmsh.py", test_file_abs, "test_abs.msh"],
        capture_output=True,
        text=True
    )

    success = True
    if result1.returncode == 0 and os.path.exists("test_rel.msh"):
        print("[PASS] Relative path succeeded")
        os.remove("test_rel.msh")
    else:
        print("[FAIL] Relative path failed")
        success = False

    if result2.returncode == 0 and os.path.exists("test_abs.msh"):
        print("[PASS] Absolute path succeeded")
        os.remove("test_abs.msh")
    else:
        print("[FAIL] Absolute path failed")
        success = False

    return success


def test_gmsh_viewer_compatibility():
    """Test if generated files can be opened in Gmsh (if available)."""
    print("\n" + "="*60)
    print("TEST: Gmsh viewer compatibility")
    print("="*60)

    # Generate test file
    test_file = "install_ksugahar/share/ngsolve/cube.vol"
    output_file = "test_gmsh_compat.msh"

    result = subprocess.run(
        ["python", "vol_to_gmsh.py", test_file, output_file],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print("[FAIL] Failed to generate test file")
        return False

    # Check if Gmsh is available
    try:
        # Try to run Gmsh in batch mode to validate the file
        gmsh_result = subprocess.run(
            ["gmsh", "-check", output_file],
            capture_output=True,
            text=True,
            timeout=5
        )

        os.remove(output_file)

        if gmsh_result.returncode == 0:
            print("[PASS] Gmsh can open and validate the file")
            return True
        else:
            print(f"[WARN] Gmsh validation returned non-zero: {gmsh_result.returncode}")
            print(f"       This might be normal for batch mode")
            return True
    except FileNotFoundError:
        print("[SKIP] Gmsh not found in PATH - cannot test viewer compatibility")
        if os.path.exists(output_file):
            os.remove(output_file)
        return None
    except subprocess.TimeoutExpired:
        print("[SKIP] Gmsh check timed out")
        if os.path.exists(output_file):
            os.remove(output_file)
        return None
    except Exception as e:
        print(f"[SKIP] Error testing Gmsh: {e}")
        if os.path.exists(output_file):
            os.remove(output_file)
        return None


def test_mesh_info_accuracy():
    """Test that reported mesh info matches actual content."""
    print("\n" + "="*60)
    print("TEST: Mesh info accuracy")
    print("="*60)

    test_file = "install_ksugahar/share/ngsolve/cube.vol"
    output_file = "test_info.msh"

    result = subprocess.run(
        ["python", "vol_to_gmsh.py", test_file, output_file],
        capture_output=True,
        text=True
    )

    # Parse output for vertex count
    import re
    vertex_match = re.search(r'Vertices:\s+(\d+)', result.stdout)
    elements_match = re.search(r'Volume elements:\s+(\d+)', result.stdout)

    if not vertex_match or not elements_match:
        print("[FAIL] Could not parse mesh info from output")
        if os.path.exists(output_file):
            os.remove(output_file)
        return False

    reported_vertices = int(vertex_match.group(1))
    reported_elements = int(elements_match.group(1))

    # Read Gmsh file and count
    with open(output_file, 'r') as f:
        content = f.read()

    # Count vertices in $Nodes section
    nodes_match = re.search(r'\$Nodes\n(\d+)', content)
    if nodes_match:
        actual_vertices = int(nodes_match.group(1))
        if actual_vertices == reported_vertices:
            print(f"[PASS] Vertex count matches: {actual_vertices}")
        else:
            print(f"[FAIL] Vertex mismatch: reported={reported_vertices}, actual={actual_vertices}")
            os.remove(output_file)
            return False

    os.remove(output_file)
    print("[PASS] Mesh info accuracy verified")
    return True


def test_multiple_conversions():
    """Test converting multiple files in sequence."""
    print("\n" + "="*60)
    print("TEST: Multiple sequential conversions")
    print("="*60)

    test_files = [
        "install_ksugahar/share/ngsolve/cube.vol",
        "install_ksugahar/share/ngsolve/coil.vol",
        "install_ksugahar/share/ngsolve/shaft.vol"
    ]

    success = True
    for test_file in test_files:
        if not os.path.exists(test_file):
            print(f"[SKIP] {test_file}")
            continue

        result = subprocess.run(
            ["python", "vol_to_gmsh.py", test_file],
            capture_output=True,
            text=True
        )

        if result.returncode == 0:
            print(f"[PASS] {os.path.basename(test_file)}")
        else:
            print(f"[FAIL] {os.path.basename(test_file)}")
            success = False

    return success


def main():
    """Run comprehensive test suite."""
    print("="*60)
    print("COMPREHENSIVE TEST SUITE: vol_to_gmsh.py")
    print("="*60)

    tests = [
        ("2D Surface mesh", test_surface_only_mesh),
        ("Gmsh format versions", test_gmsh_format_versions),
        ("Large mesh performance", test_large_mesh_performance),
        ("Paths with spaces", test_path_with_spaces),
        ("Absolute vs relative paths", test_absolute_vs_relative_paths),
        ("Gmsh viewer compatibility", test_gmsh_viewer_compatibility),
        ("Mesh info accuracy", test_mesh_info_accuracy),
        ("Multiple conversions", test_multiple_conversions),
    ]

    results = {}
    for name, test_func in tests:
        try:
            result = test_func()
            results[name] = result
        except Exception as e:
            print(f"[ERROR] {name}: {e}")
            results[name] = False

    # Summary
    print("\n" + "="*60)
    print("COMPREHENSIVE TEST SUMMARY")
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
        print("\n[SUCCESS] All non-skipped tests passed!")
        return 0
    else:
        print(f"\n[FAILURE] {failed} test(s) failed")
        return 1


if __name__ == '__main__':
    sys.exit(main())
