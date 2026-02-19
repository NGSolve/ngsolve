"""Tests for vol_to_gmsh converter."""

import os
import pytest
from pathlib import Path
from vol_to_gmsh import convert_vol_to_gmsh, convert_vol_to_vtk


@pytest.fixture
def sample_mesh(tmp_path):
    """Create a simple test mesh."""
    from netgen.csg import unit_cube
    mesh = unit_cube.GenerateMesh(maxh=0.5)
    vol_file = str(tmp_path / "test.vol")
    mesh.Save(vol_file)
    return vol_file


def test_convert_to_gmsh(sample_mesh, tmp_path):
    """Test basic Gmsh conversion."""
    output = str(tmp_path / "test.msh")
    result = convert_vol_to_gmsh(sample_mesh, output)
    assert os.path.exists(result)
    assert result == output

    with open(result) as f:
        content = f.read()
    assert "$MeshFormat" in content
    assert "$Nodes" in content
    assert "$Elements" in content


def test_convert_to_gmsh_default_name(sample_mesh):
    """Test Gmsh conversion with auto-generated filename."""
    result = convert_vol_to_gmsh(sample_mesh)
    expected = str(Path(sample_mesh).with_suffix('.msh'))
    assert result == expected
    assert os.path.exists(result)


def test_convert_to_vtk(sample_mesh, tmp_path):
    """Test VTK conversion."""
    output = str(tmp_path / "test")
    result = convert_vol_to_vtk(sample_mesh, output)
    assert os.path.exists(result)
    assert result.endswith('.vtu')


def test_missing_file():
    """Test error handling for missing file."""
    with pytest.raises(RuntimeError):
        convert_vol_to_gmsh("nonexistent.vol")


def test_mixed_elements(tmp_path):
    """Test conversion of mesh with multiple element types."""
    from netgen.occ import Box, Cylinder, Glue, OCCGeometry
    from netgen.occ import Pnt, Ax

    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    geo = OCCGeometry(box)
    mesh = geo.GenerateMesh(maxh=0.5)

    vol_file = str(tmp_path / "mixed.vol")
    mesh.Save(vol_file)

    output = str(tmp_path / "mixed.msh")
    result = convert_vol_to_gmsh(vol_file, output)
    assert os.path.exists(result)
