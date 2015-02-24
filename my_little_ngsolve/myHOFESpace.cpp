/*********************************************************************/
/* File:   myFESpace.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*

My own FESpace for high order finite elements

*/


#include <comp.hpp>    // provides FESpace, ...

#include "myHOElement.hpp"
#include "myHOFESpace.hpp"


namespace ngcomp
{

  MyHighOrderFESpace :: MyHighOrderFESpace (const MeshAccess & ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    cout << "Constructor of MyHighOrderFESpace" << endl;

    order = int(flags.GetNumFlag ("order", 2));

    // needed to draw solution function
    evaluator = new T_DifferentialOperator<DiffOpId<2> >;

    static ConstantCoefficientFunction one(1);
    integrator = GetIntegrators().CreateBFI("mass", ma.GetDimension(), &one);
  }
    
  
  MyHighOrderFESpace :: ~MyHighOrderFESpace ()
  {
    // nothing to do
  }

  
  void MyHighOrderFESpace :: Update(LocalHeap & lh)
  {
    // some global update:

    int n_vert = ma.GetNV();  
    int n_edge = ma.GetNEdges(); 
    int n_cell = ma.GetNE();  

    first_edge_dof.SetSize (n_edge+1);
    int ii = n_vert;
    for (int i = 0; i < n_edge; i++, ii+=order-1)
      first_edge_dof[i] = ii;
    first_edge_dof[n_edge] = ii;
      
    first_cell_dof.SetSize (n_cell+1);
    for (int i = 0; i < n_cell; i++, ii+=(order-1)*(order-2)/2)
      first_cell_dof[i] = ii;
    first_cell_dof[n_cell] = ii;

    // cout << "first_edge_dof = " << endl << first_edge_dof << endl;
    // cout << "first_cell_dof = " << endl << first_cell_dof << endl;

    ndof = ii;
  }

  
  void MyHighOrderFESpace :: GetDofNrs (int elnr, Array<int> & dnums) const
  {
    // returns dofs of element number elnr

    dnums.SetSize(0);

    Ngs_Element ngel = ma.GetElement (elnr);

    // vertex dofs
    for (int i = 0; i < 3; i++)
      dnums.Append (ngel.vertices[i]);

    // edge dofs
    for (int i = 0; i < 3; i++)
      {
        int first = first_edge_dof[ngel.edges[i]];
        int next  = first_edge_dof[ngel.edges[i]+1];
        for (int j = first; j < next; j++)
          dnums.Append (j);
      }

    int first = first_cell_dof[elnr];
    int next  = first_cell_dof[elnr+1];
    for (int j = first; j < next; j++)
      dnums.Append (j);

    // cout << "dnums = " << dnums << endl;
  }


  void MyHighOrderFESpace :: GetSDofNrs (int elnr, Array<int> & dnums) const
  {
    // the same for the surface elements

    dnums.SetSize(0);

    Ngs_Element ngel = ma.GetSElement (elnr);

    // vertex dofs
    for (int i = 0; i < 2; i++)
      dnums.Append (ngel.vertices[i]);

    // edge dofs
    int first = first_edge_dof[ngel.edges[0]];
    int next  = first_edge_dof[ngel.edges[0]+1];
    for (int j = first; j < next; j++)
      dnums.Append (j);
  }


  /// returns the reference-element 
  const FiniteElement & MyHighOrderFESpace :: GetFE (int elnr, LocalHeap & lh) const
  {
    MyHighOrderTrig * trig = new (lh) MyHighOrderTrig(order);

    Ngs_Element ngel = ma.GetElement (elnr);

    for (int i = 0; i < 3; i++)
      trig->SetVertexNumber (i, ngel.vertices[i]);

    return *trig;
  }


  /// the same for the surface elements
  const FiniteElement & MyHighOrderFESpace :: GetSFE (int elnr, LocalHeap & lh) const
  {
    MyHighOrderSegm * segm = new (lh) MyHighOrderSegm(order);

    Ngs_Element ngel = ma.GetSElement (elnr);

    for (int i = 0; i < 2; i++)
      segm->SetVertexNumber (i, ngel.vertices[i]);

    return *segm;
  }



  static RegisterFESpace<MyHighOrderFESpace> initifes ("myhofespace");
}
