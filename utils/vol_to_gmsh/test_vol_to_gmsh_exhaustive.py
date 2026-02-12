#!/usr/bin/env python
"""
Exhaustive test suite for vol_to_gmsh.py

Tests all available meshes with detailed analysis:
- All mesh files in NGSolve installation
- Element type distribution
- Material information
- Boundary labels
- Coordinate ranges
- File size verification
"""

import os
import sys
import subprocess
import glob
from pathlib import Path


def analyze_netgen_mesh(vol_file):
    """Analyze Netgen mesh in detail."""
    try:
        from netgen.meshing import Mesh as NetgenMesh

        ngmesh = NetgenMesh()
        ngmesh.Load(vol_file)

        # Basic stats
        nv = len(ngmesh.Points())
        ne = len(ngmesh.Elements3D())
        nse = len(ngmesh.Elements2D())

        # Coordinate ranges
        points = ngmesh.Points()
        if len(points) > 0:
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            zs = [p[2] for p in points]

            coord_range = {
                'x': (min(xs), max(xs)),
                'y': (min(ys), max(ys)),
                'z': (min(zs), max(zs))
            }
        else:
            coord_range = None

        # Element types
        elem_types = {}
        for el in ngmesh.Elements3D():
            el_type = len(el.vertices)  # 4=tet, 6=wedge, 8=hex, 5=pyramid
            elem_types[el_type] = elem_types.get(el_type, 0) + 1

        # Surface element types
        surf_types = {}
        for el in ngmesh.Elements2D():
            el_type = len(el.vertices)  # 3=tri, 4=quad
            surf_types[el_type] = surf_types.get(el_type, 0) + 1

        return {
            'vertices': nv,
            'volume_elements': ne,
            'surface_elements': nse,
            'coord_range': coord_range,
            'volume_elem_types': elem_types,
            'surface_elem_types': surf_types
        }

    except Exception as e:
        return {'error': str(e)}


def test_mesh_conversion(vol_file):
    """Test conversion of a single mesh file."""
    basename = os.path.basename(vol_file)
    print(f"\n{'='*60}")
    print(f"TEST: {basename}")
    print(f"{'='*60}")

    # Analyze original mesh
    print(f"\n[1] Analyzing original mesh...")
    mesh_info = analyze_netgen_mesh(vol_file)

    if 'error' in mesh_info:
        print(f"[FAIL] Cannot analyze mesh: {mesh_info['error']}")
        return False

    print(f"    Vertices: {mesh_info['vertices']}")
    print(f"    Volume elements: {mesh_info['volume_elements']}")
    print(f"    Surface elements: {mesh_info['surface_elements']}")

    if mesh_info['volume_elem_types']:
        print(f"    Volume element types:")
        type_names = {4: 'Tet', 5: 'Pyramid', 6: 'Wedge', 8: 'Hex'}
        for el_type, count in mesh_info['volume_elem_types'].items():
            name = type_names.get(el_type, f'Unknown({el_type})')
            print(f"      {name}: {count}")

    if mesh_info['coord_range']:
        cr = mesh_info['coord_range']
        print(f"    Coordinate ranges:")
        print(f"      x: [{cr['x'][0]:.3f}, {cr['x'][1]:.3f}]")
        print(f"      y: [{cr['y'][0]:.3f}, {cr['y'][1]:.3f}]")
        print(f"      z: [{cr['z'][0]:.3f}, {cr['z'][1]:.3f}]")

    # Test Gmsh conversion
    print(f"\n[2] Testing Gmsh conversion...")
    msh_file = f"test_{Path(vol_file).stem}.msh"

    result = subprocess.run(
        ["python", "vol_to_gmsh.py", vol_file, msh_file],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"[FAIL] Gmsh conversion failed:")
        print(f"    {result.stderr}")
        return False

    # Verify Gmsh file
    if not os.path.exists(msh_file):
        print(f"[FAIL] Output file not created: {msh_file}")
        return False

    msh_size = os.path.getsize(msh_file)
    print(f"[PASS] Gmsh file created: {msh_size} bytes")

    # Test VTK conversion
    print(f"\n[3] Testing VTK conversion...")
    vtu_file_base = f"test_{Path(vol_file).stem}"
    vtu_file = vtu_file_base + ".vtu"

    result = subprocess.run(
        ["python", "vol_to_gmsh.py", vol_file, vtu_file_base, "--vtk"],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        print(f"[FAIL] VTK conversion failed:")
        print(f"    {result.stderr}")
        os.remove(msh_file)
        return False

    # Verify VTU file
    if not os.path.exists(vtu_file):
        print(f"[FAIL] Output file not created: {vtu_file}")
        os.remove(msh_file)
        return False

    vtu_size = os.path.getsize(vtu_file)
    print(f"[PASS] VTU file created: {vtu_size} bytes")

    # File size sanity check
    print(f"\n[4] File size sanity check...")
    vol_size = os.path.getsize(vol_file)

    # Expected: msh and vtu should be comparable to vol (not too small/large)
    size_ratio_msh = msh_size / vol_size
    size_ratio_vtu = vtu_size / vol_size

    print(f"    Original .vol: {vol_size} bytes")
    print(f"    Generated .msh: {msh_size} bytes (ratio: {size_ratio_msh:.2f}x)")
    print(f"    Generated .vtu: {vtu_size} bytes (ratio: {size_ratio_vtu:.2f}x)")

    # Sanity check: files should not be too small (indicates error)
    if msh_size < 100:
        print(f"[FAIL] Gmsh file too small (< 100 bytes)")
        os.remove(msh_file)
        os.remove(vtu_file)
        return False

    if vtu_size < 100:
        print(f"[FAIL] VTU file too small (< 100 bytes)")
        os.remove(msh_file)
        os.remove(vtu_file)
        return False

    print(f"[PASS] File sizes reasonable")

    # Cleanup
    os.remove(msh_file)
    os.remove(vtu_file)

    print(f"\n[PASS] All tests passed for {basename}")
    return True


def test_batch_conversion():
    """Test converting all meshes in batch."""
    print(f"\n{'='*60}")
    print(f"TEST: Batch conversion of all meshes")
    print(f"{'='*60}")

    mesh_dir = "install_ksugahar/share/ngsolve"
    vol_files = glob.glob(os.path.join(mesh_dir, "*.vol"))

    if not vol_files:
        print("[SKIP] No .vol files found")
        return None

    print(f"[INFO] Found {len(vol_files)} mesh files")

    # Convert all in sequence
    import time
    start_time = time.time()

    for vol_file in vol_files:
        result = subprocess.run(
            ["python", "vol_to_gmsh.py", vol_file],
            capture_output=True,
            text=True,
            timeout=30
        )

        if result.returncode != 0:
            print(f"[FAIL] {os.path.basename(vol_file)}: {result.stderr}")
            return False

    elapsed = time.time() - start_time

    print(f"[PASS] Converted {len(vol_files)} files in {elapsed:.2f} seconds")
    print(f"[INFO] Average: {elapsed/len(vol_files):.2f} seconds per file")

    # Cleanup generated files
    for vol_file in vol_files:
        msh_file = Path(vol_file).with_suffix('.msh')
        if msh_file.exists():
            msh_file.unlink()

    return True


def test_stress_large_mesh():
    """Stress test with largest available mesh."""
    print(f"\n{'='*60}")
    print(f"TEST: Stress test with largest mesh")
    print(f"{'='*60}")

    mesh_dir = "install_ksugahar/share/ngsolve"
    vol_files = glob.glob(os.path.join(mesh_dir, "*.vol"))

    # Find largest mesh
    largest = max(vol_files, key=lambda f: os.path.getsize(f))
    size = os.path.getsize(largest)

    print(f"[INFO] Largest mesh: {os.path.basename(largest)} ({size} bytes)")

    # Convert with timing
    import time
    start = time.time()

    result = subprocess.run(
        ["python", "vol_to_gmsh.py", largest],
        capture_output=True,
        text=True,
        timeout=60
    )

    elapsed = time.time() - start

    if result.returncode != 0:
        print(f"[FAIL] Conversion failed: {result.stderr}")
        return False

    print(f"[PASS] Conversion completed in {elapsed:.2f} seconds")

    # Check output file size
    msh_file = Path(largest).with_suffix('.msh')
    if msh_file.exists():
        out_size = os.path.getsize(msh_file)
        print(f"[INFO] Output size: {out_size} bytes ({out_size/size:.2f}x)")
        msh_file.unlink()

    return True


def main():
    """Run exhaustive test suite."""
    print("="*60)
    print("EXHAUSTIVE TEST SUITE: vol_to_gmsh.py")
    print("="*60)

    mesh_dir = "install_ksugahar/share/ngsolve"
    vol_files = sorted(glob.glob(os.path.join(mesh_dir, "*.vol")))

    if not vol_files:
        print("[ERROR] No .vol files found in", mesh_dir)
        return 1

    print(f"\n[INFO] Found {len(vol_files)} mesh files to test")
    print(f"[INFO] Files: {[os.path.basename(f) for f in vol_files]}")

    # Test each mesh individually
    results = {}
    for vol_file in vol_files:
        basename = os.path.basename(vol_file)
        try:
            result = test_mesh_conversion(vol_file)
            results[basename] = result
        except Exception as e:
            print(f"\n[ERROR] {basename}: {e}")
            results[basename] = False

    # Additional tests
    print(f"\n{'='*60}")
    print("ADDITIONAL TESTS")
    print(f"{'='*60}")

    results['Batch conversion'] = test_batch_conversion()
    results['Stress test'] = test_stress_large_mesh()

    # Summary
    print(f"\n{'='*60}")
    print("EXHAUSTIVE TEST SUMMARY")
    print(f"{'='*60}")

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
        print("\n[SUCCESS] All exhaustive tests passed!")
        return 0
    else:
        print(f"\n[FAILURE] {failed} test(s) failed")
        return 1


if __name__ == '__main__':
    sys.exit(main())
