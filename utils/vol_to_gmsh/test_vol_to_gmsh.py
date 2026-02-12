#!/usr/bin/env python
"""
Test script for vol_to_gmsh.py

Run automated tests for Netgen .vol to Gmsh/VTK conversion.
"""

import os
import sys
import subprocess
from pathlib import Path


def run_test(description, command, expect_success=True):
    """Run a single test command."""
    print(f"\n{'='*60}")
    print(f"TEST: {description}")
    print(f"{'='*60}")
    print(f"Command: {' '.join(command)}")

    result = subprocess.run(command, capture_output=True, text=True)

    if expect_success:
        if result.returncode == 0:
            print("[PASS] Command succeeded as expected")
            return True
        else:
            print("[FAIL] Command failed unexpectedly")
            print(f"STDERR: {result.stderr}")
            return False
    else:
        if result.returncode != 0:
            print("[PASS] Command failed as expected")
            return True
        else:
            print("[FAIL] Command succeeded unexpectedly")
            return False


def main():
    """Run all tests."""
    print("vol_to_gmsh.py Test Suite")
    print("=" * 60)

    test_dir = Path("install_ksugahar/share/ngsolve")
    tests_passed = 0
    tests_total = 0

    # Test 1: Convert cube.vol to Gmsh
    tests_total += 1
    if run_test(
        "Convert cube.vol to Gmsh format",
        ["python", "vol_to_gmsh.py", str(test_dir / "cube.vol")],
        expect_success=True
    ):
        tests_passed += 1
        # Check output file exists
        if (test_dir / "cube.msh").exists():
            print("[PASS] Output file cube.msh exists")
            tests_passed += 1
        else:
            print("[FAIL] Output file cube.msh not found")
        tests_total += 1

    # Test 2: Convert coil.vol to VTK
    tests_total += 1
    if run_test(
        "Convert coil.vol to VTK format",
        ["python", "vol_to_gmsh.py", str(test_dir / "coil.vol"), "--vtk"],
        expect_success=True
    ):
        tests_passed += 1
        # Check output file exists (VTKOutput creates .vtu files)
        if (test_dir / "coil.vtu").exists():
            print("[PASS] Output file coil.vtu exists")
            tests_passed += 1
        else:
            print("[FAIL] Output file coil.vtu not found")
        tests_total += 1

    # Test 3: Error handling - missing file
    tests_total += 1
    if run_test(
        "Error handling: missing input file",
        ["python", "vol_to_gmsh.py", "nonexistent.vol"],
        expect_success=False
    ):
        tests_passed += 1

    # Test 4: Custom output filename
    tests_total += 1
    output_file = test_dir / "test_output.msh"
    if output_file.exists():
        output_file.unlink()

    if run_test(
        "Custom output filename",
        ["python", "vol_to_gmsh.py", str(test_dir / "shaft.vol"), str(output_file)],
        expect_success=True
    ):
        tests_passed += 1
        if output_file.exists():
            print("[PASS] Custom output file exists")
            tests_passed += 1
            # Cleanup
            output_file.unlink()
        else:
            print("[FAIL] Custom output file not found")
        tests_total += 1

    # Summary
    print(f"\n{'='*60}")
    print(f"TEST SUMMARY")
    print(f"{'='*60}")
    print(f"Tests passed: {tests_passed}/{tests_total}")

    if tests_passed == tests_total:
        print("\n[SUCCESS] All tests passed!")
        return 0
    else:
        print(f"\n[FAILURE] {tests_total - tests_passed} test(s) failed")
        return 1


if __name__ == '__main__':
    sys.exit(main())
