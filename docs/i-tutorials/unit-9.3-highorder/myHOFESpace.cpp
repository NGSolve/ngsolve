/*

My own FESpace for high order finite elements

*/


#include <comp.hpp>  
#include "myHOElement.hpp"
#include "myHOFESpace.hpp"


namespace ngcomp
{

  MyHighOrderFESpace :: MyHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    type = "myhofespace";

    // Get the order from the flags, default is 2
    order = int(flags.GetNumFlag ("order", 2));

    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
    evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
  }
    
  void MyHighOrderFESpace :: Update()
  {
    // some global update:

    int n_vert = ma->GetNV();  
    int n_edge = ma->GetNEdges(); 
    int n_cell = ma->GetNE();  

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

    SetNDof (ii);
  }

  void MyHighOrderFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    // returns dofs of element number elnr
    dnums.SetSize(0);
    auto ngel = ma->GetElement (ei);

    // vertex dofs
    for (auto v : ngel.Vertices())
      dnums.Append(v);

    // edge dofs
    for (auto e : ngel.Edges())
      for(auto j : Range(first_edge_dof[e], first_edge_dof[e+1]))
        dnums.Append (j);

    // inner dofs
    if (ei.IsVolume())
      for(auto j : Range(first_cell_dof[ei.Nr()], first_cell_dof[ei.Nr()+1]))
        dnums.Append (j);
  }

  
  FiniteElement & MyHighOrderFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    auto ngel = ma->GetElement (ei);
    switch (ngel.GetType())
      {
      case ET_TRIG:
        {
          auto trig = new (alloc) MyHighOrderTrig(order);
          trig->SetVertexNumbers (ngel.vertices);
          return *trig;
        }
      case ET_SEGM:
        {
          auto segm = new (alloc) MyHighOrderSegm(order);
          segm->SetVertexNumbers (ngel.vertices);
          return *segm;
        }
      default:
        throw Exception (string("Element type ")+ToString(ngel.GetType())+" not supported");
      }
  }
  
}
