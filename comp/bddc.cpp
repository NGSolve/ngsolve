#include <comp.hpp>
#include <solve.hpp>


namespace ngcomp
{

 
  template <class SCAL, class TV>
  class BDDCMatrix : public BaseMatrix
  {
    const BilinearForm & bfa;

    Array<Matrix<SCAL>*> & elmats;
    Array<Array<int>*> & eldnums;

    BaseMatrix * harmonicext;
    BaseMatrix * harmonicexttrans;
    BaseMatrix * innersolve;    
    BaseMatrix * pwbmat;    

    bool block;

    BaseMatrix * inv;
    BaseMatrix * inv_coarse;
    string inversetype;
    BitArray * free_dofs;

    BaseVector * tmp;
    BaseVector * tmp2;

  public:
    
    BDDCMatrix (const BilinearForm & abfa, 
		Array<Matrix<SCAL>*> & aelmats, Array<Array<int>*> & aeldnums,
		const string & ainversetype, bool ablock)
      : bfa(abfa), elmats(aelmats), eldnums(aeldnums), 
	block(ablock), inversetype(ainversetype)
    {
      static Timer timer ("BDDC Constructor");
      static Timer timer1 ("BDDC Constructor 1");
      static Timer timer2 ("BDDC Constructor 2");
      static Timer timer3 ("BDDC Constructor 3");

      pwbmat = NULL;
      inv = NULL;
      inv_coarse = NULL;
      tmp = NULL;
      tmp2 = NULL;


      RegionTimer reg(timer);
      timer1.Start();

      MyMPI_Barrier ();
      cout << IM(1) << "BDDC-marix constructor" << endl;
      MyMPI_Barrier ();

      const FESpace & fes = bfa.GetFESpace();
      int ne = elmats.Size();
      int ndof = fes.GetNDof();      

      
      Array<int> wbdcnt(ne);  //count number of wirebasket dofs on each element
      Array<int> ifcnt(ne);   //count number of interface dofs on each element
      wbdcnt = 0;
      ifcnt = 0;

      for (int i = 0; i < ne; i++)
        {
	  if (!eldnums[i]) continue;
	  Array<int> & dnums = *eldnums[i];
	  for (int j = 0; j < dnums.Size(); j++)
	    {
	      if (dnums[j] == -1) continue;
	      COUPLING_TYPE ct = fes.GetDofCouplingType(dnums[j]);
	      if (ct == WIREBASKET_DOF)
		wbdcnt[i]++;
	      else
		ifcnt[i]++;
	    }
        }

      Table<int> el2wbdofs(wbdcnt);   // wirebasket dofs on each element
      Table<int> el2ifdofs(ifcnt);    // interface dofs on each element

      for (int i = 0; i < ne; i++)
        {
	  if (!eldnums[i]) continue;
	  Array<int> & dnums = *eldnums[i];

	  int lifcnt = 0;
	  int lwbcnt = 0;
	  for (int j = 0; j < dnums.Size(); j++)
	    {
	      if (dnums[j] == -1) continue;
	      COUPLING_TYPE ct = fes.GetDofCouplingType(dnums[j]);
	      
	      if (ct == WIREBASKET_DOF)
		el2wbdofs[i][lwbcnt++] = dnums[j];
	      else
		el2ifdofs[i][lifcnt++] = dnums[j];
	    }
        }


      cout << IM(3) << ( bfa.IsSymmetric() ? "symmetric" : "non-symmetric") << endl;

      SparseMatrix<SCAL,TV,TV> * sparse_innersolve;
      SparseMatrix<SCAL,TV,TV> * sparse_harmonicext;
      SparseMatrix<SCAL,TV,TV> * sparse_harmonicexttrans;

      if (!bfa.IsSymmetric())
	{
	  harmonicexttrans = sparse_harmonicexttrans =
	    new SparseMatrix<SCAL,TV,TV>(ndof, el2wbdofs, el2ifdofs, false);
	  harmonicexttrans -> AsVector() = 0.0;
	}

      innersolve = sparse_innersolve = bfa.IsSymmetric() 
	? new SparseMatrixSymmetric<SCAL,TV>(ndof, el2ifdofs)
	: new SparseMatrix<SCAL,TV,TV>(ndof, el2ifdofs, el2ifdofs, bfa.IsSymmetric());
      innersolve->AsVector() = 0.0;

      harmonicext = sparse_harmonicext =
	new SparseMatrix<SCAL,TV,TV>(ndof, el2ifdofs, el2wbdofs, false);
      harmonicext->AsVector() = 0.0;

      pwbmat = bfa.IsSymmetric() 
	? new SparseMatrixSymmetric<SCAL,TV>(ndof, el2wbdofs)
	: new SparseMatrix<SCAL,TV,TV>(ndof, el2wbdofs, el2wbdofs, bfa.IsSymmetric());
      pwbmat -> AsVector() = 0.0;
      pwbmat -> SetInverseType (inversetype);
      
      timer1.Stop();






      Array<double> weight(fes.GetNDof());      
      weight = 0;

      timer2.Start();

      for (int i = 0; i < ne; i++)
	{
	    if (!eldnums[i]) continue;

	    FlatMatrix<SCAL> elmat = *elmats[i]; 
	    Array<int> & dnums = *eldnums[i];

	    FlatArray<int> interfacedofs = el2ifdofs[i];
	    FlatArray<int> wirebasketdofs = el2wbdofs[i];
	    Array<int> localwbdofs;    
	    Array<int> localintdofs;   
	    Array<int> ldnums;
	    Array<double> el2ifweight;

	    for (int k = 0; k < dnums.Size(); k++)
	      {
		COUPLING_TYPE ct = fes.GetDofCouplingType(dnums[k]);	      
		if (ct == WIREBASKET_DOF)
		  localwbdofs.Append (k);
		else
		  {
		    localintdofs.Append (k);

		    if (typeid(SCAL) == typeid(double))
		      el2ifweight.Append (fabs (elmat(k,k)));
		    else
		      el2ifweight.Append (1);
		  }
	      }

	    for (int j = 0; j < interfacedofs.Size(); j++)
	      weight[interfacedofs[j]] += el2ifweight[j];
	    

	    for (int k = 0; k < localwbdofs.Size(); k++)
	      if (dnums[localwbdofs[k]] != wirebasketdofs[k])
		cout << "wbERROR" << endl;
	    for (int k = 0; k < localintdofs.Size(); k++)
	      if (dnums[localintdofs[k]] != interfacedofs[k])
		cout << "intERROR" << endl;

	    int sizew = localwbdofs.Size();
	    int sizei = localintdofs.Size();

	    Matrix<SCAL> a = elmat.Rows(localwbdofs).Cols(localwbdofs);

	    if (sizei)
	      {      
		Matrix<SCAL> b = elmat.Rows(localwbdofs).Cols(localintdofs);
		Matrix<SCAL> c = elmat.Rows(localintdofs).Cols(localwbdofs);
		Matrix<SCAL> d = elmat.Rows(localintdofs).Cols(localintdofs);
	      
		LapackInverse (d);
	      
		Matrix<SCAL> he (sizei, sizew);
		he = SCAL(0.0);
		
		if (sizew)
		  {
		    he -= d*c   | Lapack;
		    a += b*he   | Lapack;
		  }

		//R * E
		for (int k = 0; k < sizei; k++)
		  he.Row(k) *= el2ifweight /* [i] */[k]; 

		sparse_harmonicext->AddElementMatrix(interfacedofs,wirebasketdofs,he);
	    
		if (!bfa.IsSymmetric())
		  {
		    Matrix<SCAL> het (sizew, sizei);
		
		    het = SCAL(0.0);
		    LapackMultAddAB(b,d,-1.0,het);
		
		    //E * R^T
		    for (int l = 0; l < sizei; l++)
		      het.Col(l) *= el2ifweight /* [i] */[l];

		    sparse_harmonicexttrans->AddElementMatrix(wirebasketdofs,interfacedofs,het);
		  }
	    
		//R * A_ii^(-1) * R^T
		for (int k = 0; k < sizei; k++) d.Row(k) *= el2ifweight[k]; 
		for (int l = 0; l < sizei; l++) d.Col(l) *= el2ifweight[l]; 
		
		sparse_innersolve -> AddElementMatrix(interfacedofs,interfacedofs,d);
	      }
	  
	    dynamic_cast<SparseMatrix<SCAL,TV,TV>*>(pwbmat)
	      ->AddElementMatrix(wirebasketdofs,wirebasketdofs,a);
	  }





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

      timer2.Stop();
      timer3.Start();

      for (int i = 0; i < elmats.Size(); i++)
	{
	  delete elmats[i];
	  delete eldnums[i];
	}

      elmats.SetSize(0);
      eldnums.SetSize(0);
      
      


      
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
	  if (MyMPI_GetNTasks() > 1)
	    {
	      MyMPI_Barrier();

	      ParallelMatrix parwb(*pwbmat, bfa.GetFESpace().GetParallelDofs());
	      parwb.SetInverseType (inversetype);
	      inv = parwb.InverseMatrix (free_dofs);

	      tmp = new ParallelVVector<TV>(ndof, &bfa.GetFESpace().GetParallelDofs());
	      innersolve = new ParallelMatrix (*innersolve, bfa.GetFESpace().GetParallelDofs());
	      harmonicext = new ParallelMatrix (*harmonicext, bfa.GetFESpace().GetParallelDofs());
	      harmonicexttrans = new ParallelMatrix (*harmonicexttrans, bfa.GetFESpace().GetParallelDofs());
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
      timer3.Stop();
    }

    ~BDDCMatrix()
    {
      delete inv;
      delete pwbmat;
      delete inv_coarse;
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
    BaseMatrix * pre;
    string inversetype;
    bool block;

    Array<Matrix<SCAL>*> elmats;
    Array<Array<int>*> eldnums;
  public:
    BDDCPreconditioner (const PDE & pde, const Flags & aflags, const string & aname)
      : Preconditioner (&pde, aflags, aname)
    {
      bfa = dynamic_cast<const S_BilinearForm<SCAL>*>(pde.GetBilinearForm (aflags.GetStringFlag ("bilinearform", NULL)));
      const_cast<S_BilinearForm<SCAL>*> (bfa) -> SetPreconditioner (this);
      inversetype = flags.GetStringFlag("inverse", "sparsecholesky");
      if (flags.GetDefineFlag("refelement")) Exception ("refelement - BDDC not supported");
      block = flags.GetDefineFlag("block");
      pre = NULL;
    }
    
    virtual ~BDDCPreconditioner()
    {
      delete pre;
    }
    


    virtual void InitLevel () { ; }

    virtual void FinishLevel () { ; }


    virtual void AddElementMatrix (const Array<int> & dnums,
				   const FlatMatrix<SCAL> & elmat,
				   bool inner_element, int elnr,
				   LocalHeap & lh);

    virtual void Update ();

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
    int hnr = elnr;
    if (!inner_element)   // is boundary element
      hnr += ma.GetNE();   

#pragma omp critical(bddcaddelmat)
    {
      if (hnr >= elmats.Size())
	{
	  int oldsize = elmats.Size();
	  elmats.SetSize (hnr+1);
	  eldnums.SetSize (hnr+1);
	  for (int j = oldsize; j < hnr+1; j++)
	    elmats[j] = NULL;
	  for (int j = oldsize; j < hnr+1; j++)
	    eldnums[j] = NULL;
	}

      int used = 0;
      for (int i = 0; i < dnums.Size(); i++)
	if (dnums[i] != -1 && fes.GetFreeDofs()->Test(dnums[i])) used++;

      delete eldnums[hnr];
      eldnums[hnr] = new Array<int> (used);

      for (int i = 0, ii = 0; i < dnums.Size(); i++)
	if (dnums[i] != -1 && fes.GetFreeDofs()->Test(dnums[i])) 
	  {
	    (*eldnums[hnr])[ii] = dnums[i];
	    ii++;
	  }

      delete elmats[hnr];
      elmats[hnr] = new Matrix<SCAL> (used);

      for (int i = 0, ii = 0; i < dnums.Size(); i++)
	if (dnums[i] != -1 && fes.GetFreeDofs()->Test(dnums[i]))
	  {
	    for (int j = 0, jj = 0; j < dnums.Size(); j++)
	      if (dnums[j] != -1 && fes.GetFreeDofs()->Test(dnums[j]))
		{
		  (*elmats[hnr])(ii,jj) = elmat(i,j);
		  jj++;
		}
	    ii++;
	  }

      if (L2Norm (*elmats[hnr]) == 0)
	{
	  delete elmats[hnr];
	  delete eldnums[hnr];
	  elmats[hnr] = NULL;
	  eldnums[hnr] = NULL;
	}
    }
  }
  


  template <class SCAL, class TV>
  void BDDCPreconditioner<SCAL, TV> ::
  Update ()
  {
    cout << IM(3) << "update bddc, inversetype = " << inversetype << endl;

    if (elmats.Size() || (MyMPI_GetId() == 0 && MyMPI_GetNTasks() > 1))
      {
	delete pre;
	pre = new BDDCMatrix<SCAL,TV>(*bfa, elmats, eldnums, inversetype, block);
      }
    else
      {
	cerr << "don't update precond, since I don't have matrices" << endl;
      }
    
    for (int i = 0; i < elmats.Size(); i++)
      {
	delete elmats[i];
	delete eldnums[i];
      }
    elmats.SetSize(0);
    eldnums.SetSize(0);

    if (test) Test();
  }  

  template class BDDCPreconditioner<double>;
  template class BDDCPreconditioner<double, Complex>;
    


  static RegisterPreconditioner<BDDCPreconditioner<double> > initpre ("bddc");
  static RegisterPreconditioner<BDDCPreconditioner<Complex> > initpre2 ("bddcc");
  static RegisterPreconditioner<BDDCPreconditioner<double,Complex> > initpre3 ("bddcrc");
}
