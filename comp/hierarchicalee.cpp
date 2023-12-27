
remove file


/*********************************************************************/
/* File:   hierarchicalee.cpp                                        */
/* Author: Joachim Schoeberl                                         */
/* Date:   12. May. 2003                                             */
/*********************************************************************/

/* 
   Hierarchical error estimator
*/

#include <comp.hpp>

namespace ngcomp
{
  using namespace ngcomp;

  template <class SCAL>
  void CalcErrorHierarchical (const S_BilinearForm<SCAL> & bfa,
			      const S_BilinearForm<SCAL> & bfa2,
			      const S_LinearForm<SCAL> & lff,
			      S_GridFunction<SCAL> & gfu,
			      const FESpace & festest,
			      FlatVector<double> & err,
			      LocalHeap & lh)
  {
#ifdef ISTSCHONALT
    int i, j, k;

    int ne = ma.GetNE();
    int nse = ma.GetNSE();

    const FESpace & fes = gfu.GetFESpace();
    const S_BaseVector<SCAL> & vecu = 
      dynamic_cast<const S_BaseVector<SCAL>&> (gfu.GetVector());
    ElementTransformation eltrans;

    int dim     = fes.GetDimension();
    Array<int> dnums, dnumstest;

    VVector<SCAL> res(festest.GetNDof());
    VVector<SCAL> diag(festest.GetNDof());
    res = 0;
    diag = 0;
    for (i = 0; i < ne; i++)
      {
	lh.CleanUp();

	const FiniteElement & fel = fes.GetFE(i, lh);
	const FiniteElement & feltest = festest.GetFE(i, lh);

	ma.GetElementTransformation (i, eltrans, lh);
	fes.GetDofNrs (i, dnums);
	festest.GetDofNrs (i, dnumstest);

	FlatVector<SCAL> elu(dnums.Size() * dim, lh);
	FlatVector<SCAL> elres(dnumstest.Size() * dim, lh);
	FlatVector<SCAL> eldiag(dnumstest.Size() * dim, lh);
	vecu.GetIndirect (dnums, elu);
	fes.TransformVec (i, 0, elu, TRANSFORM_SOL);

	elres = 0;
	eldiag = 0;

	for (j = 0; j < lff.NumIntegrators(); j++)
	  {
	    if (lff.GetIntegrator(j) -> BoundaryForm()) continue;
	    FlatVector<SCAL> elvec;
	    lff.GetIntegrator(j) -> CalcElementVector (feltest, eltrans, elvec, lh);
	    elres += elvec;
	  }
	
	for (j = 0; j < bfa.NumIntegrators(); j++)
	  {
	    if (bfa.GetIntegrator(j) -> BoundaryForm()) continue;
	    FlatVector<SCAL> elvec (dnumstest.Size()*dim, lh);
	    bfa.GetIntegrator(j) -> ApplyMixedElementMatrix (fel, feltest, eltrans, elu, elvec, lh);
	    elres -= elvec;
	  }

	for (j = 0; j < bfa2.NumIntegrators(); j++)
	  {
	    if (bfa2.GetIntegrator(j) -> BoundaryForm()) continue;
	    FlatVector<SCAL> elvec (dnumstest.Size()*dim, lh);

	    bfa2.GetIntegrator(j) -> CalcElementMatrixDiag (feltest, eltrans, elvec, lh);
	    eldiag += elvec;
	  }

	festest.TransformVec (i, 0, elres, TRANSFORM_RHS);
	res.AddIndirect (dnumstest, elres);
	diag.AddIndirect (dnumstest, eldiag);
	cout << "\rcompute element " << i << "/" << ne << flush;
      }

    (*testout) << "res = " << endl << res << endl;
    (*testout) << "diag = " << endl << diag << endl;


    for (i = 0; i < ne; i++)
      {
	lh.CleanUp();
	
	festest.GetDofNrs (i, dnumstest);

	double sum = 0;
	for (j = 0; j < dnumstest.Size(); j++)
	  {
	    int di = dnumstest[j];
	    sum += abs (res(di)*res(di)/diag(di));
	  }
	
	err(i) = sum;
      }
#endif
  }
 

  template
  void CalcErrorHierarchical<double> (const S_BilinearForm<double> & bfa,
				      const S_BilinearForm<double> & bfa2,
				      const S_LinearForm<double> & lff,
				      S_GridFunction<double> & bu,
				      const FESpace & festest,
				      FlatVector<double> & err,
				      LocalHeap & lh);

}
