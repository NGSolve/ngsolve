
from netgen.csg import *
from ngsolve import *


def test_node_access():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))
    cell_done = [False] * mesh.ne
    face_done = [False] * mesh.nface
    edge_done = [False] * mesh.nedge
    vert_done = [False] * mesh.nv

    for cell in mesh.nodes(CELL):
        cell_done[cell.nr] = True
        assert len(cell.vertices) == 4
        assert len(cell.edges) == 6
        assert len(cell.faces) == 4
        for vert in cell.vertices:
            vert_done[vert.nr] = True
        for edge in cell.edges:
            edge_done[edge.nr] = True
            assert len(mesh[edge].vertices) == 2
            for v in mesh[edge].vertices:
                assert v in cell.vertices
        for face in cell.faces:
            face_done[face.nr] = True
            assert len(mesh[face].vertices) == 3
            for v in mesh[face].vertices:
                assert v in cell.vertices
                
            print (mesh[face].vertices)
            print (mesh[face].edges)
            
            assert len(mesh[face].edges) == 3
            for e in mesh[face].edges:
                assert e in cell.edges

    assert all(cell_done)
    assert all(face_done)
    assert all(edge_done)
    assert all(vert_done)

def test_dofs():
    mesh = Mesh(unit_cube.GenerateMesh(maxh=0.2))
    for order in range(1,6):
        fes = H1(mesh,order=order)

        dof_done = [False] * fes.ndof

        for vert in mesh.vertices:
            dofs = fes.GetDofNrs(vert)
            assert len(dofs) == 1
            dof_done[dofs[0]] = True

        for edge in mesh.edges:
            dofs = fes.GetDofNrs(edge)
            assert len(dofs) == (order-1 if order > 1 else 0)
            for dof in dofs:
                dof_done[dof] = True

        for face in mesh.faces:
            dofs = fes.GetDofNrs(face)
            assert len(dofs) == (((order-1) * (order-2))/2 if order > 2 else 0)
            for dof in dofs:
                dof_done[dof] = True

        for cell in mesh.nodes(CELL):
            dofs = fes.GetDofNrs(cell)
            assert len(dofs) == (((order-3)*(order-2)*(order-1))/6 if order > 3 else 0)
            for dof in dofs:
                dof_done[dof] = True

        assert all(dof_done)

if __name__ == "__main__":
    test_node_access()
    test_dofs()
