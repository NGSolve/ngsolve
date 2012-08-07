#include <comp.hpp>
#include <solve.hpp>


namespace ngcomp
{

 
  template <class SCAL, class TV>
  class BDDCMatrix : public BaseMatrix
  {
    const BilinearForm & bfa;

    BaseMatrix *harmonicext, *harmonicexttrans, 
      *innersolve, *pwbmat;    

    SparseMatrix<SCAL,TV,TV> *sparse_innersolve, 
      *sparse_harmonicext, *sparse_harmonicexttrans;


    Array<double> weight;
    
    bool block;
    bool hypre;

    BaseMatrix * inv;
    BaseMatrix * inv_coarse;
    string inversetype;
    BitArray * free_dofs;

    BaseVector * tmp;
    BaseVector * tmp2;

  public:

    void SetHypre (bool ah = true) { hypre = ah; }
    
    BDDCMatrix (const BilinearForm & abfa, 
		const string & ainversetype, bool ablock, bool ahypre)
      : bfa(abfa), block(ablock), inversetype(ainversetype)
    {
      static Timer timer ("BDDC Constructor");

      hypre = ahypre;

      pwbmat = NULL;
      inv = NULL;
      inv_coarse = NULL;
      tmp = NULL;
      tmp2 = NULL;

      RegionTimer reg(timer);

      const FESpace & fes = bfa.GetFESpace();
      const MeshAccess & ma = fes.GetMeshAccess();

      Array<int> wbdcnt(ma.GetNE()+ma.GetNSE());
      Array<int> ifcnt(ma.GetNE()+ma.GetNSE());
      wbdcnt = 0;
      ifcnt = 0;
      const BitArray & freedofs = *fes.GetFreeDofs();
      
      for (int bound = 0, ii = 0; bound <= 1; bound++)
	for (int i = 0; i < (bound ? ma.GetNSE() : ma.GetNE()); i++, ii++)
	  {
	    Array<int> dnums;
	    fes.GetDofNrs (i, bound, dnums);
	    for (int j = 0; j < dnums.Size(); j++)
	      {
		if (dnums[j] == -1) continue;
		if (!freedofs.Test(dnums[j])) continue;
		COUPLING_TYPE ct = fes.GetDofCouplingType(dnums[j]);
		if (ct == LOCAL_DOF && bfa.UsesEliminateInternal()) continue;
		
		if (ct == WIREBASKET_DOF)
		  wbdcnt[ii]++;
		else
		  ifcnt[ii]++;
	      }
	  }
      
      Table<int> el2wbdofs(wbdcnt);   // wirebasket dofs on each element
      Table<int> el2ifdofs(ifcnt);    // interface dofs on each element
      
      for (int bound = 0, ii = 0; bound <= 1; bound++)
	for (int i = 0; i < (bound ? ma.GetNSE() : ma.GetNE()); i++, ii++)
	  {
	    Array<int> dnums;
	    fes.GetDofNrs (i, bound, dnums);
	    
	    int lifcnt = 0;
	    int lwbcnt = 0;
	    
	    for (int j = 0; j < dnums.Size(); j++)
	      {
		if (dnums[j] == -1) continue;
		if (!freedofs.Test(dnums[j])) continue;
		COUPLING_TYPE ct = fes.GetDofCouplingType(dnums[j]);
		if (ct == LOCAL_DOF && bfa.UsesEliminateInternal()) continue;
		
		if (ct == WIREBASKET_DOF)
		  el2wbdofs[ii][lwbcnt++] = dnums[j];
		else
		  el2ifdofs[ii][lifcnt++] = dnums[j];
	      } 
	  }
      

      int ndof = fes.GetNDof();      
      
      if (!bfa.IsSymmetric())
	{
	  harmonicexttrans = sparse_harmonicexttrans =
	    new SparseMatrix<SCAL,TV,TV>(ndof, el2wbdofs, el2ifdofs, false);
	  harmonicexttrans -> AsVector() = 0.0;
	}
      else
	harmonicexttrans = sparse_harmonicexttrans = NULL;


      innersolve = sparse_innersolve = bfa.IsSymmetric() 
	? new SparseMatrixSymmetric<SCAL,TV>(ndof, el2ifdofs)
	: new SparseMatrix<SCAL,TV,TV>(ndof, el2ifdofs, el2ifdofs, bfa.IsSymmetric());
      innersolve->AsVector() = 0.0;

      harmonicext = sparse_harmonicext =
	new SparseMatrix<SCAL,TV,TV>(ndof, el2ifdofs, el2wbdofs, false);
      harmonicext->AsVector() = 0.0;

      pwbmat = bfa.IsSymmetric() && !hypre
	? new SparseMatrixSymmetric<SCAL,TV>(ndof, el2wbdofs)
	: new SparseMatrix<SCAL,TV,TV>(ndof, el2wbdofs, el2wbdofs, bfa.IsSymmetric() && !hypre);
      pwbmat -> AsVector() = 0.0;
      pwbmat -> SetInverseType (inversetype);
      
      weight.SetSize (fes.GetNDof());
      weight = 0;
    }

    
    void AddMatrix (FlatMatrix<SCAL> elmat, Array<int> & dnums)

    {
      static Timer timer ("BDDC - Addmatrix");
      RegionTimer reg(timer);
      static Timer timer2("BDDC - Add to sparse");
      static Timer timer3("BDDC - compute");

      const FESpace & fes = bfa.GetFESpace();
      
      ArrayMem<int, 20> localwbdofs, localintdofs;
      
      for (int k = 0; k < dnums.Size(); k++)
	{
	  COUPLING_TYPE ct = fes.GetDofCouplingType(dnums[k]);	      
	  if (ct == WIREBASKET_DOF)
	    localwbdofs.Append (k);
	  else
	    localintdofs.Append (k);
	}
      
      int sizew = localwbdofs.Size();
      int sizei = localintdofs.Size();
      
      ArrayMem<double, 20> el2ifweight(sizei);
      for (int k = 0; k < sizei; k++)
	if (typeid(SCAL) == typeid(double))
	  el2ifweight[k] = fabs (elmat(localintdofs[k],
				       localintdofs[k]));
	else
	  el2ifweight[k] = 1;


      
      Matrix<SCAL> a = elmat.Rows(localwbdofs).Cols(localwbdofs);
      Matrix<SCAL> b = elmat.Rows(localwbdofs).Cols(localintdofs);
      Matrix<SCAL> c = elmat.Rows(localintdofs).Cols(localwbdofs);
      Matrix<SCAL> d = elmat.Rows(localintdofs).Cols(localintdofs);
      Matrix<SCAL> het (sizew, sizei);
      Matrix<SCAL> he (sizei, sizew);
	  
      if (sizei)
	{      
	  RegionTimer reg(timer3);
	  timer3.AddFlops (sizei*sizei*sizei + 2*sizei*sizei*sizew);

	  LapackInverse (d);
	  
	  if (sizew)
	    {
	      he = SCAL(0.0);

	      he -= d*c   | Lapack;
	      a += b*he   | Lapack;
	  
	      //R * E
	      for (int k = 0; k < sizei; k++)
		he.Row(k) *= el2ifweight[k]; 
	  
	      if (!bfa.IsSymmetric())
		{	      
		  het = SCAL(0.0);
		  het -= b*d | Lapack;
		  
		  //E * R^T
		  for (int l = 0; l < sizei; l++)
		    het.Col(l) *= el2ifweight[l];
		}
	    }
	  
	  //R * A_ii^(-1) * R^T
	  for (int k = 0; k < sizei; k++) d.Row(k) *= el2ifweight[k]; 
	  for (int l = 0; l < sizei; l++) d.Col(l) *= el2ifweight[l]; 
	}


      RegionTimer regadd(timer2);

      ArrayMem<int, 20> wbdofs, intdofs;   
      wbdofs = dnums[localwbdofs];
      intdofs = dnums[localintdofs];

#pragma omp critical(bddcaddelmat)
      {
	for (int j = 0; j < intdofs.Size(); j++)
	  weight[intdofs[j]] += el2ifweight[j];

	sparse_harmonicext->AddElementMatrix(intdofs,wbdofs,he);
	
	if (!bfa.IsSymmetric())
	  sparse_harmonicexttrans->AddElementMatrix(wbdofs,intdofs,het);
	
	sparse_innersolve -> AddElementMatrix(intdofs,intdofs,d);
	
	dynamic_cast<SparseMatrix<SCAL,TV,TV>*>(pwbmat)
	  ->AddElementMatrix(wbdofs,wbdofs,a);
      }
    }


    
    void Finalize()
    {
      cout << IM(1) << "Setup BDDC preconditioner" << endl;

      static Timer timer ("BDDC Constructor 3");
      RegionTimer reg(timer);

      const FESpace & fes = bfa.GetFESpace();
      int ndof = fes.GetNDof();      


#ifdef PARALLEL
      AllReduceDofData (weight, MPI_SUM, fes.GetParallelDofs());
#endif

      for (int i = 0; i < sparse_innersolve->Height(); i++)
	{
	  FlatArray<int> cols = sparse_innersolve -> GetRowIndices(i);
	  for (int j = 0; j < cols.Size(); j++)
	    if (weight[i] && weight[cols[j]])
	      sparse_innersolve->GetRowValues(i)(j) /= (weight[i] * weight[cols[j]]);
	}
      
      
      for (int i = 0; i < sparse_harmonicext->Height(); i++)
	if (weight[i])
	  sparse_harmonicext->GetRowValues(i) /= weight[i];
      
      if (!bfa.IsSymmetric())
	for (int i = 0; i < sparse_harmonicexttrans->Height(); i++)
	  if (weight[i])
	    sparse_harmonicexttrans->GetRowValues(i) /= weight[i];


      // now generate wire-basked solver

      free_dofs = new BitArray (ndof);
      free_dofs->Clear();

      // *free_dofs = wbdof;
      for (int i = 0; i < ndof; i++)
	if (fes.GetDofCouplingType(i) == WIREBASKET_DOF)
	  free_dofs -> Set(i);


      if (fes.GetFreeDofs())
	free_dofs -> And (*fes.GetFreeDofs());
      int cntfreedofs = 0;
      for (int i = 0; i < free_dofs->Size(); i++)
	if (free_dofs->Test(i)) cntfreedofs++;


      if (block)
	{
	  //Smoothing Blocks
	  Flags flags;
	  flags.SetFlag("eliminate_internal");
	  flags.SetFlag("subassembled");
	  cout << "call Create Smoothing Blocks of " << bfa.GetFESpace().GetName() << endl;
	  Table<int> & blocks = *(bfa.GetFESpace().CreateSmoothingBlocks(flags));
	  cout << "has blocks" << endl << endl;
	  // *testout << "blocks = " << endl << blocks << endl;
	  // *testout << "pwbmat = " << endl << *pwbmat << endl;
	  cout << "call block-jacobi inverse" << endl;
	  inv = dynamic_cast<BaseSparseMatrix*> (pwbmat)->CreateBlockJacobiPrecond(blocks, 0, 0, 0);      
	  cout << "has inverse" << endl << endl;
	  // *testout << "blockjacobi = " << endl << *inv << endl;
	  
	  //Coarse Grid of Wirebasket
	  cout << "call directsolverclusters inverse" << endl;
	  Array<int> & clusters = *(bfa.GetFESpace().CreateDirectSolverClusters(flags));
	  cout << "has clusters" << endl << endl;
	  
	  cout << "call coarse wirebasket grid inverse" << endl;
	  inv_coarse = pwbmat->InverseMatrix(&clusters);
	  cout << "has inverse" << endl << endl;
	  
	  tmp = new VVector<>(ndof);
	  tmp2 = new VVector<>(ndof);
	}
      else
	{
#ifdef PARALLEL
	  if (bfa.GetFESpace().IsParallel())
	    {
	      ParallelDofs * pardofs = &bfa.GetFESpace().GetParallelDofs();

	      pwbmat = new ParallelMatrix (pwbmat, pardofs);
	      pwbmat -> SetInverseType (inversetype);

#ifdef HYPRE
	      if (hypre)
		inv = new HyprePreconditioner (*pwbmat, free_dofs);
	      else
#endif
		inv = pwbmat -> InverseMatrix (free_dofs);

	      tmp = new ParallelVVector<TV>(ndof, pardofs);
	      innersolve = new ParallelMatrix (innersolve, pardofs);
	      harmonicext = new ParallelMatrix (harmonicext, pardofs);
	      harmonicexttrans = new ParallelMatrix (harmonicexttrans, pardofs);
	    }
	  else
#endif
	    {
	      cout << "call wirebasket inverse ( with " << cntfreedofs 
		   << " free dofs out of " << pwbmat->Height() << " )" << endl;

	      inv = pwbmat->InverseMatrix(free_dofs);
	      cout << "has inverse" << endl;
	      tmp = new VVector<TV>(ndof);
	    }
	}
    }

    ~BDDCMatrix()
    {
      delete inv;
      delete pwbmat;
      delete inv_coarse;
      delete harmonicext;
      delete harmonicexttrans;
      delete innersolve;

      delete tmp;
      delete tmp2;
    }
    
    virtual BaseVector * CreateVector () const
    {
      return bfa.GetMatrix().CreateVector();
    }

    
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      static Timer timer ("Apply BDDC preconditioner");
      static Timer timerifs ("Apply BDDC preconditioner - apply ifs");
      static Timer timerwb ("Apply BDDC preconditioner - wb solve");
      static Timer timerharmonicext ("Apply BDDC preconditioner - harmonic extension");
      static Timer timerharmonicexttrans ("Apply BDDC preconditioner - harmonic extension trans");
      

      RegionTimer reg (timer);

      x.Cumulate();
      y = x;

      timerharmonicexttrans.Start();

      if (bfa.IsSymmetric())
	y += Transpose(*harmonicext) * x; 
      else
	y += *harmonicexttrans * x;

      timerharmonicexttrans.Stop();

      timerwb.Start();
      *tmp = 0;
      if (block)
	{
	  if (true) //GS
	    {
	      dynamic_cast<BaseBlockJacobiPrecond*>(inv)->GSSmoothResiduum (*tmp, y, *tmp2 ,1);
	      
	      if (inv_coarse)
		*tmp += (*inv_coarse) * *tmp2; 
	      dynamic_cast<BaseBlockJacobiPrecond*>(inv)->GSSmoothBack (*tmp, y);
	    }
	  else
	    { //jacobi only (old)
	      *tmp = (*inv) * y;
	      *tmp += (*inv_coarse) * y; 
	    }
	}
      else
	{
	  *tmp = (*inv) * y;
	}
      timerwb.Stop();

      timerifs.Start();
      *tmp += *innersolve * x;
      timerifs.Stop();

      timerharmonicext.Start();
      
      y = *tmp;
      y += *harmonicext * *tmp;

      timerharmonicext.Stop();

      y.Cumulate();
    }
  };









  template <class SCAL, class TV = SCAL>
  class NGS_DLL_HEADER BDDCPreconditioner : public Preconditioner
  {
    const S_BilinearForm<SCAL> * bfa;
    BDDCMatrix<SCAL,TV> * pre;
    string inversetype;
    bool block, hypre;
  public:
    BDDCPreconditioner (const PDE & pde, const Flags & aflags, const string & aname)
      : Preconditioner (&pde, aflags, aname)
    {
      bfa = dynamic_cast<const S_BilinearForm<SCAL>*>(pde.GetBilinearForm (aflags.GetStringFlag ("bilinearform", NULL)));
      const_cast<S_BilinearForm<SCAL>*> (bfa) -> SetPreconditioner (this);
      inversetype = flags.GetStringFlag("inverse", "sparsecholesky");
      if (flags.GetDefineFlag("refelement")) Exception ("refelement - BDDC not supported");
      block = flags.GetDefineFlag("block");
      hypre = flags.GetDefineFlag("usehypre");
      pre = NULL;
    }
    
    virtual ~BDDCPreconditioner()
    {
      delete pre;
    }
    
    virtual void InitLevel () 
    {
      delete pre;
      pre = new BDDCMatrix<SCAL,TV>(*bfa, inversetype, block, hypre);
      pre -> SetHypre (hypre);
    }

    virtual void FinalizeLevel () 
    {  
      pre -> Finalize();
    }


    virtual void AddElementMatrix (const Array<int> & dnums,
				   const FlatMatrix<SCAL> & elmat,
				   bool inner_element, int elnr,
				   LocalHeap & lh);

    virtual void Update ()
    {
      if (test) Test();
    }  

    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa->GetMatrix();
    }

    virtual const BaseMatrix & GetMatrix() const
    {
      return *pre;
    }

    virtual void CleanUpLevel ()
    {
      delete pre;
      pre = NULL;
    }


    virtual void Mult (const BaseVector & x, BaseVector & y) const
    {
      y = 0.0;
      pre -> MultAdd (1, x, y);
    }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      pre -> MultAdd (s, x, y);
    }


    virtual const char * ClassName() const
    { return "BDDC Preconditioner"; }
  };




  template <class SCAL, class TV>
  void BDDCPreconditioner<SCAL, TV> ::
  AddElementMatrix (const Array<int> & dnums,
		    const FlatMatrix<SCAL> & elmat,
		    bool inner_element, int elnr,
		    LocalHeap & lh)
  {
    const FESpace & fes = bfa->GetFESpace();

    int used = 0;
    for (int i = 0; i < dnums.Size(); i++)
      if (dnums[i] != -1 && fes.GetFreeDofs()->Test(dnums[i])) used++;
    
    ArrayMem<int,20> hdnums(used);
    Matrix<SCAL> helmat(used);
    
    for (int i = 0, ii = 0; i < dnums.Size(); i++)
      if (dnums[i] != -1 && fes.GetFreeDofs()->Test(dnums[i]))
	{
	  hdnums[ii] = dnums[i];
	  
	  for (int j = 0, jj = 0; j < dnums.Size(); j++)
	    if (dnums[j] != -1 && fes.GetFreeDofs()->Test(dnums[j]))
	      {
		helmat(ii,jj) = elmat(i,j);
		jj++;
	      }
	  ii++;
	}
    
    if (L2Norm (helmat) != 0)
      pre -> AddMatrix(helmat, hdnums);
  }
  






  template class BDDCPreconditioner<double>;
  template class BDDCPreconditioner<double, Complex>;

  static RegisterPreconditioner<BDDCPreconditioner<double> > initpre ("bddc");
  static RegisterPreconditioner<BDDCPreconditioner<Complex> > initpre2 ("bddcc");
  static RegisterPreconditioner<BDDCPreconditioner<double,Complex> > initpre3 ("bddcrc");
}
