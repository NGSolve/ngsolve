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
#include <diffop_impl.hpp>

#include "myElement.hpp"
#include "myFESpace.hpp"


namespace ngcomp
{

  MyFESpace :: MyFESpace (const MeshAccess & ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    cout << "Constructor of MyFESpace" << endl;
    cout << "Flags = " << flags << endl;

    secondorder = flags.GetDefineFlag ("secondorder");

    if (!secondorder)
      cout << "You have chosen first order elements" << endl;
    else
      cout << "You have chosen second order elements" << endl;


    if (!secondorder)
      reference_element = new MyLinearTrig;
    else
      reference_element = new MyQuadraticTrig;


    if (!secondorder)
      reference_surface_element = new FE_Segm1;
    else
      reference_surface_element = new FE_Segm2;

    
    // needed to draw solution function
    evaluator = new T_DifferentialOperator<DiffOpId<2> >;

    static ConstantCoefficientFunction one(1);
    integrator = GetIntegrators().CreateBFI("mass", ma.GetDimension(), &one);
  }
    
  
  MyFESpace :: ~MyFESpace ()
  {
    // nothing to do
  }

  
  void MyFESpace :: Update(LocalHeap & lh)
  {
    // some global update:

    cout << "Update MyFESpace, #vert = " << ma.GetNV() 
         << ", #edge = " << ma.GetNEdges() << endl;


    nvert = ma.GetNV();
    // number of dofs:
    if (!secondorder)
      ndof = nvert;  // number of vertices
    else
      ndof = nvert + ma.GetNEdges();  // num vertics + num edges
  }

  
  void MyFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    // returns dofs of element number elnr

    dnums.SetSize(0);

    Array<int> vert_nums;
    ma.GetElVertices (elnr, vert_nums);  // global vertex numbers
    
    // first 3 dofs are vertex numbers:
    for (int i = 0; i < 3; i++)
      dnums.Append (vert_nums[i]);

    if (secondorder)
      {
        // 3 more dofs on edges:

        Array<int> edge_nums;
        ma.GetElEdges (elnr, edge_nums);    // global edge numbers

        for (int i = 0; i < 3; i++)
          dnums.Append (nvert+edge_nums[i]);
      }
  }


  void MyFESpace :: GetSDofNrs (int elnr, Array<int> & dnums) const
  {
    // the same for the surface elements

    dnums.SetSize(0);

    Array<int> vert_nums;
    ma.GetSElVertices (elnr, vert_nums);  
    
    // first 2 dofs are vertex numbers:
    for (int i = 0; i < 2; i++)
      dnums.Append (vert_nums[i]);

    if (secondorder)
      {
        // 1 more dof on the edge:

        Array<int> edge_nums;
        ma.GetSElEdges (elnr, edge_nums);    // global edge number

        dnums.Append (nvert+edge_nums[0]);
      }
  }

  /// returns the reference-element 
  const FiniteElement & MyFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    if (ma.GetElType(elnr) == ET_TRIG)
      return *reference_element;

    // should return different elements for mixed-element type meshes
    throw Exception ("Sorry, only triangular elements are supported");
  }

  /// the same for the surface elements
  const FiniteElement & MyFESpace :: GetSFE (int selnr, LocalHeap & lh) const
  {
    return *reference_surface_element;
  }




  /*
    register fe-spaces
    Object of type MyFESpace can be defined in the pde-file via
    "define fespace v -type=myfespace"
  */

  static RegisterFESpace<MyFESpace> initifes ("myfespace");
}
