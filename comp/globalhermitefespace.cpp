#include "globalhermitefespace.hpp"
//using namespace ngmg;

namespace ngcomp
{
  GlobalHermiteFESpace::GlobalHermiteFESpace (int aspacial_dim, Array<int> & aorders,const Flags & flags, bool parseflags)
    : FESpace (NULL, flags), spacial_dim(aspacial_dim)
  {
    orders.SetSize(aorders.Size());
    orders = aorders;
    T_Ansatz.SetSize(orders.Size());
    V_Ansatz.SetSize(orders.Size());
    bool usehmfuncs = flags.GetDefineFlag("usefunctions");
    T_Ansatz = usehmfuncs? 0.5 : 1.0;
    T_Ansatz[orders.Size()/2]*=1.0;
    T_Ansatz[orders.Size()/3]*=1.0;
    for(int i=0;i<V_Ansatz.Size();i++)
    {
        V_Ansatz[i] = new Vector<>(spacial_dim);
        *V_Ansatz[i] = 0.0;
    }
    
    //V_Ansatz[orders.Size()/5]=Vec<1>(0.3);
    if(spacial_dim == 1)
    {
      hermitefel = new Distribution<1>(orders[0]);
      ndof = dynamic_cast<Distribution<1> *>(hermitefel)->GetNDof<NODAL>();
      evaluator = make_shared<T_DifferentialOperator<DiffOpId<1>>>();
      flux_evaluator = make_shared<T_DifferentialOperator<DiffOpGradient<1>>>();
      boundary_evaluator = make_shared<T_DifferentialOperator<DiffOpIdBoundary<1>>>();
      dynamic_cast<Distribution<1> *>(hermitefel)->SetUseHermiteFunctions(usehmfuncs);
    }
    if(spacial_dim == 2)
    {
      hermitefel = new Distribution<2>(orders[0]);
      ndof = dynamic_cast<Distribution<2> *>(hermitefel)->GetNDof<NODAL>();
      evaluator = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
      flux_evaluator = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
      boundary_evaluator = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
      dynamic_cast<Distribution<2> *>(hermitefel)->SetUseHermiteFunctions(usehmfuncs);
    }
    if(spacial_dim == 3)
    {
      hermitefel = new Distribution<3>(orders[0]);
      ndof = dynamic_cast<Distribution<3> *>(hermitefel)->GetNDof<NODAL>();
      evaluator = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
      flux_evaluator = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
      boundary_evaluator = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
      dynamic_cast<Distribution<3> *>(hermitefel)->SetUseHermiteFunctions(usehmfuncs);
    }
  }
  void GlobalHermiteFESpace::Update(LocalHeap & lh)
  {
  }
}