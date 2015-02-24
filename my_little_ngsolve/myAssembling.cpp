/*********************************************************************/
/* File:   myAssembling.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   3. May. 2010                                              */
/*********************************************************************/


/*

Assembling the matrix

*/

#include <solve.hpp>

using namespace ngsolve;

namespace myAssembling
{
  class NumProcMyAssembling : public NumProc
  {
  protected:
    shared_ptr<GridFunction> gfu;

  public:
    
    NumProcMyAssembling (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    { 
      cout << "We assemble matrix and rhs vector" << endl;

      gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));
    }
  
    virtual string GetClassName () const
    {
      return "MyAssembling";
    }


    virtual void Do (LocalHeap & lh)
    {
      shared_ptr<FESpace> fes = gfu -> GetFESpace();

      int ndof = fes->GetNDof();
      int ne = GetMeshAccess()->GetNE();
    

      // setup element->dof table:


      // first we get the number of dofs per element ...
      Array<int> dnums;
      Array<int> cnt(ne);
      for (int i = 0; i < ne; i++)
	{
	  fes->GetDofNrs (i, dnums);
	  cnt[i] = dnums.Size();
	}	  
      
      // allocate the table in compressed form ...
      Table<int> el2dof(cnt);

      // and fill it
      for (int i = 0; i < ne; i++)
	{
	  fes->GetDofNrs (i, dnums);
          el2dof[i] = dnums;
	}
      cout << "el2dof - table: " << el2dof << endl;

      // generate sparse matrix from element-to-dof table
      SparseMatrixSymmetric<double> & mat = *new SparseMatrixSymmetric<double> (ndof, el2dof);

      VVector<double> vecf (fes->GetNDof());

      LaplaceIntegrator<2> laplace (make_shared<ConstantCoefficientFunction> (1));
      SourceIntegrator<2> source (make_shared<ConstantCoefficientFunction> (1));

      mat = 0.0;
      vecf = 0.0;

      for (int i = 0; i < ne; i++)   // loop over elements
	{  
	  HeapReset hr(lh); 

	  const ElementTransformation & eltrans = ma->GetTrafo (i, 0, lh);
	  
	  fes->GetDofNrs (i, dnums);
	  const FiniteElement & fel =  fes->GetFE (i, lh);
	  
	  FlatMatrix<> elmat (dnums.Size(), lh);
	  laplace.CalcElementMatrix (fel, eltrans, elmat, lh);
	  mat.AddElementMatrix (dnums, elmat);

	  FlatVector<> elvec (dnums.Size(), lh);
	  source.CalcElementVector (fel, eltrans, elvec, lh);
	  vecf.AddIndirect (dnums, elvec);
	} 

      *testout << "mat = " << mat << endl;
      *testout << "vecf = " << vecf << endl;

      shared_ptr<BaseMatrix> inv = mat.InverseMatrix (fes->GetFreeDofs());

      gfu -> GetVector() = (*inv) * vecf;
    }
  };

  static RegisterNumProc<NumProcMyAssembling> npinit1("myassembling");
}
