import pytest
from ngsolve import VectorialVoxelCoefficient, Mesh
from netgen.csg import unit_cube
from netgen.geom2d import unit_square
import numpy as np

unit_cube_mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))
unit_square_mesh = Mesh(unit_square.GenerateMesh(maxh=0.2))


def test_vector_voxelcf_2d():
    vec_data = np.array([1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0]).reshape(2, 2, 2)
    vec_vcf = VectorialVoxelCoefficient(start=(0, 0), end=(1, 1), values=vec_data, linear=False, shape=(2,))
    mip = unit_square_mesh(0.1, 0.1)
    assert np.all(np.isclose(np.array([1.0, 0.0]), vec_vcf(mip)))
    return


def test_vector_voxelcf_2d_complex():
    vec_data = np.array([1j, 0.0, 0.0, 1j, -1j, 0.0, 0.0, -1j]).reshape(2, 2, 2)
    vec_vcf = VectorialVoxelCoefficient(start=(0, 0), end=(1, 1), values=vec_data, linear=False, shape=(2,))
    mip = unit_square_mesh(0.1, 0.1)
    assert np.all(np.isclose(np.array([1j, 0.0]), vec_vcf(mip)))
    return


def test_vector_voxelcf_3d():
    # 8 quadrants, 1 voxel per quadrant
    vec_data = np.array([1.0, 0.0, 0.0,  # x-direction
                         0.0, 1.0, 0.0,  # y-direction
                         -1.0, 0.0, 0.0,  # -x-direction
                         0.0, -1.0, 0.0,  # -y-direction
                         1.0, 0.0, 1.0,  # xz-direction
                         0.0, 1.0, 1.0,  # yz-direction
                         -1.0, 0.0, -1.0,  # -xz-direction
                         0.0, -1.0, -1.0,  # -yz-direction
                         ]).reshape(2, 2, 2, 3)
    vec_vcf = VectorialVoxelCoefficient(start=(0, 0, 0), end=(1, 1, 1), values=vec_data, linear=False, shape=(3,))
    mip = unit_cube_mesh(0.1, 0.1, 0.1)
    assert np.all(np.isclose(np.array([1.0, 0.0, 0.0]), vec_vcf(mip)))
    return


def test_vector_voxelcf_3d_complex():
    # 8 quadrants, 1 voxel per quadrant
    vec_data = np.array([1j, 0.0, 0.0,  # x-direction
                         0.0, 1j, 0.0,  # y-direction
                         -1j, 0.0, 0.0,  # -x-direction
                         0.0, -1j, 0.0,  # -y-direction
                         1j, 0.0, 1j,  # xz-direction
                         0.0, 1j, 1j,  # yz-direction
                         -1j, 0.0, -1j,  # -xz-direction
                         0.0, -1j, -1j,  # -yz-direction
                         ]).reshape(2, 2, 2, 3)
    vec_vcf = VectorialVoxelCoefficient(start=(0, 0, 0), end=(1, 1, 1), values=vec_data, linear=False, shape=(3,))
    mip = unit_cube_mesh(0.1, 0.1, 0.1)
    assert np.all(np.isclose(np.array([1j, 0.0, 0.0]), vec_vcf(mip)))
    return


def test_matrix_voxelcf_2d():
    mat_data = np.array([1.0, 2.0, 3.0, 4.0,
                         1.0, 2.0, 3.0, 4.0,
                         1.0, 2.0, 3.0, 4.0,
                         1.0, 2.0, 3.0, 4.0]).reshape(2, 2, 2, 2)

    mat_vcf = VectorialVoxelCoefficient(start=(0, 0), end=(1, 1), values=mat_data, linear=False, shape=(2, 2))
    mip = unit_square_mesh(0.1, 0.1)
    assert np.all(np.isclose(np.array([1.0, 2.0, 3.0, 4.0]), mat_vcf(mip)))
    return


def test_matrix_voxelcf_2d_complex():
    mat_data = np.array([1.0j, 2.0j, 3.0j, 4.0j,
                         1.0j, 2.0j, 3.0j, 4.0j,
                         1.0j, 2.0j, 3.0j, 4.0j,
                         1.0j, 2.0j, 3.0j, 4.0j]).reshape(2, 2, 2, 2)

    mat_vcf = VectorialVoxelCoefficient(start=(0, 0), end=(1, 1), values=mat_data, linear=False, shape=(2, 2))
    mip = unit_square_mesh(0.1, 0.1)
    assert np.all(np.isclose(np.array([1.0j, 2.0j, 3.0j, 4.0j]), mat_vcf(mip)))
    return


def test_matrix_voxelcf_3d():
    mat_data = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                         1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                         1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                         1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                         1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                         1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                         1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                         1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]).reshape(2, 2, 2, 3, 3)
    mat_vcf = VectorialVoxelCoefficient(start=(0, 0, 0), end=(1, 1, 1), values=mat_data, linear=False, shape=(3, 3))
    mip = unit_cube_mesh(0.1, 0.1, 0.1)
    assert np.all(np.isclose(np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]), mat_vcf(mip)))
    return


def test_matrix_voxelcf_3d_complex():
    mat_data = np.array([1.0j, 2.0j, 3.0j, 4.0j, 5.0j, 6.0j, 7.0j, 8.0j, 9.0j,
                         1.0j, 2.0j, 3.0j, 4.0j, 5.0j, 6.0j, 7.0j, 8.0j, 9.0j,
                         1.0j, 2.0j, 3.0j, 4.0j, 5.0j, 6.0j, 7.0j, 8.0j, 9.0j,
                         1.0j, 2.0j, 3.0j, 4.0j, 5.0j, 6.0j, 7.0j, 8.0j, 9.0j,
                         1.0j, 2.0j, 3.0j, 4.0j, 5.0j, 6.0j, 7.0j, 8.0j, 9.0j,
                         1.0j, 2.0j, 3.0j, 4.0j, 5.0j, 6.0j, 7.0j, 8.0j, 9.0j,
                         1.0j, 2.0j, 3.0j, 4.0j, 5.0j, 6.0j, 7.0j, 8.0j, 9.0j,
                         1.0j, 2.0j, 3.0j, 4.0j, 5.0j, 6.0j, 7.0j, 8.0j, 9.0j]).reshape(2, 2, 2, 3, 3)
    mat_vcf = VectorialVoxelCoefficient(start=(0, 0, 0), end=(1, 1, 1), values=mat_data, linear=False, shape=(3, 3))
    mip = unit_cube_mesh(0.1, 0.1, 0.1)
    assert np.all(np.isclose(np.array([1.0j, 2.0j, 3.0j, 4.0j, 5.0j, 6.0j, 7.0j, 8.0j, 9.0j]), mat_vcf(mip)))
    return


def test_matrix_vec_product_voxelcf_2d():
    vec_data = np.array([1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0]).reshape(2, 2, 2)
    vec_vcf = VectorialVoxelCoefficient(start=(0, 0), end=(1, 1), values=vec_data, linear=False, shape=(2,))
    voxel_identity_mat = np.array([1.0, 0.0,
                                   0.0, 1.0,
                                   1.0, 0.0,
                                   0.0, 1.0,
                                   1.0, 0.0,
                                   0.0, 1.0,
                                   1.0, 0.0,
                                   0.0, 1.0]).reshape(2, 2, 2, 2)
    voxel_identity_mat_cf = VectorialVoxelCoefficient(start=(0, 0), end=(1, 1), values=voxel_identity_mat, linear=False, shape=(2, 2))
    prod_cf = voxel_identity_mat_cf * vec_vcf
    mip = unit_square_mesh(0.1, 0.1)
    assert np.all(np.isclose(prod_cf(mip), vec_vcf(mip)))
    return


def test_matrix_vec_product_voxelcf_2d_complex():
    vec_data = np.array([1j, 0.0, 0.0, 1j, -1j, 0.0, 0.0, -1j]).reshape(2, 2, 2)
    vec_vcf = VectorialVoxelCoefficient(start=(0, 0), end=(1, 1), values=vec_data, linear=False, shape=(2,))
    voxel_identity_mat = np.array([1.0, 0.0,
                                   0.0, 1.0,
                                   1.0, 0.0,
                                   0.0, 1.0,
                                   1.0, 0.0,
                                   0.0, 1.0,
                                   1.0, 0.0,
                                   0.0, 1.0], dtype=np.complex128).reshape(2, 2, 2, 2)
    voxel_identity_mat_cf = VectorialVoxelCoefficient(start=(0, 0), end=(1, 1), values=voxel_identity_mat, linear=False, shape=(2, 2))
    prod_cf = voxel_identity_mat_cf * vec_vcf
    mip = unit_square_mesh(0.1, 0.1)
    assert np.all(np.isclose(prod_cf(mip), vec_vcf(mip)))
    return


def test_matrix_vec_product_voxelcf_3d():
    # 8 quadrants, 1 voxel per quadrant
    vec_data = np.array([1.0, 0.0, 0.0,  # x-direction
                         0.0, 1.0, 0.0,  # y-direction
                         -1.0, 0.0, 0.0,  # -x-direction
                         0.0, -1.0, 0.0,  # -y-direction
                         1.0, 0.0, 1.0,  # xz-direction
                         0.0, 1.0, 1.0,  # yz-direction
                         -1.0, 0.0, -1.0,  # -xz-direction
                         0.0, -1.0, -1.0,  # -yz-direction
                         ]).reshape(2, 2, 2, 3)

    vec_vcf = VectorialVoxelCoefficient(start=(0, 0, 0), end=(1, 1, 1), values=vec_data, linear=False, shape=(3,))

    voxel_identity_mat = np.array([1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0]).reshape(2, 2, 2, 3, 3)
    voxel_identity_mat_cf = VectorialVoxelCoefficient(start=(0, 0, 0), end=(1, 1, 1), values=voxel_identity_mat, linear=False, shape=(3, 3))
    prod_cf = voxel_identity_mat_cf * vec_vcf
    mip = unit_cube_mesh(0.1, 0.1, 0.1)
    assert np.all(np.isclose(prod_cf(mip), vec_vcf(mip)))


def test_matrix_vec_product_voxelcf_3d_complex():
    # 8 quadrants, 1 voxel per quadrant
    vec_data = np.array([1j, 0.0, 0.0,  # x-direction
                         0.0, 1j, 0.0,  # y-direction
                         -1j, 0.0, 0.0,  # -x-direction
                         0.0, -1j, 0.0,  # -y-direction
                         1j, 0.0, 1j,  # xz-direction
                         0.0, 1j, 1j,  # yz-direction
                         -1j, 0.0, -1j,  # -xz-direction
                         0.0, -1j, -1j,  # -yz-direction
                         ]).reshape(2, 2, 2, 3)

    vec_vcf = VectorialVoxelCoefficient(start=(0, 0, 0), end=(1, 1, 1), values=vec_data, linear=False, shape=(3,))

    voxel_identity_mat = np.array([1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0,
                                   1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0], dtype=np.complex128).reshape(2, 2, 2, 3, 3)
    voxel_identity_mat_cf = VectorialVoxelCoefficient(start=(0, 0, 0), end=(1, 1, 1), values=voxel_identity_mat, linear=False, shape=(3, 3))
    prod_cf = voxel_identity_mat_cf * vec_vcf
    mip = unit_cube_mesh(0.1, 0.1, 0.1)
    assert np.all(np.isclose(prod_cf(mip), vec_vcf(mip)))
