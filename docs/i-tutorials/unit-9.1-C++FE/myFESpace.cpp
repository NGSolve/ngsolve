/*********************************************************************/
/* File:   myFESpace.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/


/*

My own FESpace for linear and quadratic triangular elements.

A fe-space provides the connection between the local reference
element, and the global mesh.

*/


#include <comp.hpp>    // provides FESpace, ...
// #include <h1lofe.hpp>
// #include <regex>
#include <python_comp.hpp>
#include "myElement.hpp"
#include "myFESpace.hpp"
#include "myDiffOp.hpp"


namespace ngcomp
{

  MyFESpace :: MyFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    cout << "Constructor of MyFESpace" << endl;
    cout << "Flags = " << flags << endl;

    // this is needed for pickling and needs to be the same as used in
    // RegisterFESpace later
    type = "myfespace";

    secondorder = flags.GetDefineFlag ("secondorder");

    if (!secondorder)
      cout << "You have chosen first order elements" << endl;
    else
      cout << "You have chosen second order elements" << endl;

    // needed for symbolic integrators and to draw solution
    evaluator[VOL] = make_shared<T_DifferentialOperator<MyDiffOpId>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<MyDiffOpGradient>>();
  }

  DocInfo MyFESpace :: GetDocu()
  {
    auto docu = FESpace::GetDocu();
    docu.Arg("secondorder") = "bool = False\n"
      "  Use second order basis functions";
    return docu;
  }

  void MyFESpace :: Update()
  {
    // some global update:
    cout << "Update MyFESpace, #vert = " << ma->GetNV()
         << ", #edge = " << ma->GetNEdges() << endl;

    // number of vertices
    nvert = ma->GetNV();

    // number of dofs:
    size_t ndof = nvert;
    if (secondorder)
      ndof += ma->GetNEdges();  // num vertics + num edges

    SetNDof (ndof);
  }

  void MyFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    // returns dofs of element ei
    // may be a volume triangle or boundary segment

    dnums.SetSize(0);

    // first dofs are vertex numbers:
    for (auto v : ma->GetElVertices(ei))
      dnums.Append (v);

    if (secondorder)
      {
        // more dofs on edges:
        for (auto e : ma->GetElEdges(ei))
          dnums.Append (nvert+e);
      }
  }

  FiniteElement & MyFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    if (ei.IsVolume())
      {
        if (!secondorder)
          return * new (alloc) MyLinearTrig;
        else
          return * new (alloc) MyQuadraticTrig;
      }
    else
      throw Exception("Boundary elements not implemented yet!");
    // else
    //   {
    //     if (!secondorder)
    //       return * new (alloc) MyLinearSegm;
    //     else
    //       return * new (alloc) MyQuadraticSegm;
    //   }
  }

  /*
    register fe-spaces
    Object of type MyFESpace can be defined in the pde-file via
    "define fespace v -type=myfespace"
  */

  static RegisterFESpace<MyFESpace> initifes ("myfespace");
}

void ExportMyFESpace(py::module m)
{
  cout << "called ExportMyFESpace" << endl;
  using namespace ngcomp;

  ExportFESpace<MyFESpace>(m, "MyFESpace")
    .def("GetNVert", &MyFESpace::GetNVert)
    ;
}
