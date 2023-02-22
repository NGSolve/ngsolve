#include <comp.hpp>
// #include <solve.hpp>
#include "hypre_precond.hpp"


namespace ngcomp
{

 
  template <class SCAL, class TV>
  class BDDCMatrix : public BaseMatrix
  {
    shared_ptr<BilinearForm> bfa;
    shared_ptr<FESpace> fes;

    shared_ptr<BaseMatrix> harmonicext, harmonicexttrans, 
      innersolve;

    shared_ptr<BaseMatrix> pwbmat;
    shared_ptr<BaseSparseMatrix> sparse_pwbmat;

    shared_ptr<SparseMatrix<SCAL,TV,TV>> sparse_innersolve, 
      sparse_harmonicext, sparse_harmonicexttrans;


    Array<double> weight;
    
    bool block;
    bool hypre;
    bool coarse;
    bool local; // act as bddc for the local matrix
    
    shared_ptr<BaseMatrix> inv;
    shared_ptr<BaseMatrix> inv_coarse;
    string inversetype;   //sparsecholesky or pardiso or ....
    string coarsetype;    //general precond.. (e.g. AMG)

    shared_ptr<BaseVector> tmp;
    shared_ptr<BaseVector> tmp2;

    shared_ptr<BitArray> wb_free_dofs;

  public:

    void SetHypre (bool ah = true) { hypre = ah; }
    void SetCoarseGridPreconditioner (bool ah = true) { coarse = ah; }
    
    BDDCMatrix (shared_ptr<BilinearForm> abfa, 
                Flags flags,
		const string & ainversetype, 
		const string & acoarsetype, 
                bool ablock, 
                bool ahypre
      )
      : bfa(abfa), block(ablock), inversetype(ainversetype), coarsetype(acoarsetype)
    {
      static Timer timer ("BDDC Constructor");

      fes = bfa->GetFESpace();
      
      coarse = (coarsetype != "none");

      hypre = ahypre;

      local = flags.GetDefineFlag("local");
      
      // pwbmat = NULL;
      inv = NULL;

      inv_coarse = NULL;
      // tmp = NULL;
      // tmp2 = NULL;
      RegionTimer reg(timer);

      // auto fes = bfa -> GetFESpace();
      shared_ptr<MeshAccess> ma = fes->GetMeshAccess();

      Array<int> wbdcnt(ma->GetNE()+ma->GetNSE()+ma->GetNCD2E()+bfa->GetSpecialElements().Size());
      Array<int> ifcnt(ma->GetNE()+ma->GetNSE()+ma->GetNCD2E()+bfa->GetSpecialElements().Size());
      wbdcnt = 0;
      ifcnt = 0;
      const BitArray & freedofs = *fes->GetFreeDofs();
      

      LocalHeap lh(10000, "BDDC-constr, dummy heap");
      
      for (auto vb : { VOL, BND, BBND })
        IterateElements 
          (*fes, vb, lh, 
           [&] (FESpace::Element el, LocalHeap & lh)
           {
             int base = (vb == VOL) ? 0 : ((vb == BND) ? ma->GetNE() : ma->GetNE() + ma->GetNSE());
             for (auto d : el.GetDofs())
               {
                 if (!IsRegularDof(d)) continue;
                 if (!freedofs.Test(d)) continue;
                 COUPLING_TYPE ct = fes->GetDofCouplingType(d);
                 if ( ((ct & CONDENSABLE_DOF) != 0) && bfa->UsesEliminateInternal()) continue;
		
                 int ii = base + el.Nr();
                 if (ct == WIREBASKET_DOF)
                   wbdcnt[ii]++;
                 else
                   ifcnt[ii]++;
               }
           });
      Array<int> dnums;
      int base = ma->GetNE() + ma->GetNSE() + ma->GetNCD2E();
      for(auto i : Range(bfa->GetSpecialElements()))
        {
          auto& el = bfa->GetSpecialElements()[i];
          dnums.SetSize0();
          el->GetDofNrs(dnums);
          for (auto d : dnums) {
            if (!IsRegularDof(d)) continue;
            if (!freedofs.Test(d)) continue;
            COUPLING_TYPE ct = fes->GetDofCouplingType(d);
            if (((ct & CONDENSABLE_DOF) != 0) && bfa->UsesEliminateInternal())
              continue;

            int ii = base + i;
            if (ct == WIREBASKET_DOF)
              wbdcnt[ii]++;
            else
              ifcnt[ii]++;
          }
        }

      Table<int> el2wbdofs(wbdcnt);   // wirebasket dofs on each element
      Table<int> el2ifdofs(ifcnt);    // interface dofs on each element
      
      for (auto vb : { VOL, BND, BBND })
        IterateElements 
          (*fes, vb, lh, 
           [&] (FESpace::Element el, LocalHeap & lh)
           {
             int base = (vb == VOL) ? 0 : ((vb == BND) ? ma->GetNE() : ma->GetNE() + ma->GetNSE());
             int lifcnt = 0;
             int lwbcnt = 0;

             for (auto d : el.GetDofs())
               {
                 if (!IsRegularDof(d)) continue;
                 if (!freedofs.Test(d)) continue;
                 COUPLING_TYPE ct = fes->GetDofCouplingType(d);
                 if ( ((ct & CONDENSABLE_DOF) != 0) && bfa->UsesEliminateInternal()) continue;
		
                 int ii = base + el.Nr();
                 if (ct == WIREBASKET_DOF)
                   el2wbdofs[ii][lwbcnt++] = d;
                 else
                   el2ifdofs[ii][lifcnt++] = d;
               }
           });

      base = ma->GetNE() + ma->GetNSE() + ma->GetNCD2E();
      for(auto i : Range(bfa->GetSpecialElements()))
        {
          auto& el = bfa->GetSpecialElements()[i];
          int lifcnt = 0;
          int lwbcnt = 0;
          dnums.SetSize0();
          el->GetDofNrs(dnums);
          for (auto d : dnums) {
            if (!IsRegularDof(d))
              continue;
            if (!freedofs.Test(d))
              continue;
            COUPLING_TYPE ct = fes->GetDofCouplingType(d);
            if (((ct & CONDENSABLE_DOF) != 0) && bfa->UsesEliminateInternal())
              continue;

            int ii = base + i;
            if (ct == WIREBASKET_DOF)
              el2wbdofs[ii][lwbcnt++] = d;
            else
              el2ifdofs[ii][lifcnt++] = d;
          }
        }

      auto ndof = fes->GetNDof();      

      wb_free_dofs = make_shared<BitArray> (ndof);
      wb_free_dofs->Clear();

      // *wb_free_dofs = wbdof;
      for (auto i : Range(ndof))
	if (fes->GetDofCouplingType(i) == WIREBASKET_DOF)
	  wb_free_dofs -> SetBit(i);


      if (fes->GetFreeDofs())
	wb_free_dofs -> And (*fes->GetFreeDofs());
      
      if (!bfa->SymmetricStorage()) 
	{
	  harmonicexttrans = sparse_harmonicexttrans =
	    make_shared<SparseMatrix<SCAL,TV,TV>>(ndof, ndof, el2wbdofs, el2ifdofs, false);
	  harmonicexttrans -> AsVector() = 0.0;
	}
      else
	harmonicexttrans = sparse_harmonicexttrans = nullptr;


      innersolve = sparse_innersolve = bfa->SymmetricStorage() 
	? make_shared<SparseMatrixSymmetric<SCAL,TV>>(ndof, el2ifdofs)
	: make_shared<SparseMatrix<SCAL,TV,TV>>(ndof, ndof, el2ifdofs, el2ifdofs, false); // bfa.IsSymmetric());
      innersolve->AsVector() = 0.0;

      harmonicext = sparse_harmonicext =
	make_shared<SparseMatrix<SCAL,TV,TV>>(ndof, ndof, el2ifdofs, el2wbdofs, false);
      harmonicext->AsVector() = 0.0;
      if (bfa->SymmetricStorage() && !hypre)
        pwbmat = make_shared<SparseMatrixSymmetric<SCAL,TV>>(ndof, el2wbdofs);
      else
        pwbmat = make_shared<SparseMatrix<SCAL,TV,TV>>(ndof, ndof, el2wbdofs, el2wbdofs, false); // bfa.IsSymmetric() && !hypre);
      pwbmat -> AsVector() = 0.0;
      pwbmat -> SetInverseType (inversetype);
      sparse_pwbmat = dynamic_pointer_cast<BaseSparseMatrix>(pwbmat);
      sparse_pwbmat-> SetSPD ( bfa->IsSPD() );
      
      weight.SetSize (fes->GetNDof());
      weight = 0;

      GetMemoryTracer().SetName("BDDCMatrix");
      GetMemoryTracer().Track(
          weight, "weight",
          *sparse_innersolve, "mat_innersolve",
          *sparse_harmonicext, "mat_harmonicext",
          *sparse_harmonicexttrans, "mat_harmonicexttrans"
          );

      if (coarse)
      {
        flags.SetFlag ("not_register_for_auto_update");
	auto creator = GetPreconditionerClasses().GetPreconditioner(coarsetype);
	if(creator == nullptr)
	  throw Exception("Nothing known about preconditioner " + coarsetype);
        inv = creator->creatorbf (bfa, flags, "wirebasket"+coarsetype);
        dynamic_pointer_cast<Preconditioner>(inv) -> InitLevel(wb_free_dofs);
      }
    }

    bool IsComplex() const override { return pwbmat -> IsComplex(); }
    
    void AddMatrix (FlatMatrix<SCAL> elmat, FlatArray<int> dnums, 
		    ElementId id, LocalHeap & lh)

    {
      static Timer timer ("BDDC - Addmatrix", NoTracing);
      RegionTimer reg (timer);
      static Timer timer2("BDDC - Add to sparse", NoTracing, NoTiming);
      static Timer timer3("BDDC - compute", NoTracing, NoTiming);

      HeapReset hr(lh);

      // auto fes = bfa->GetFESpace();
      
      ArrayMem<int, 100> localwbdofs, localintdofs;
      
      for (int k : Range(dnums))
	{
	  COUPLING_TYPE ct = fes->GetDofCouplingType(dnums[k]);	      
	  if (ct == WIREBASKET_DOF)
	    localwbdofs.Append (k);
	  else
	    localintdofs.Append (k);
	}

      int sizew = localwbdofs.Size();
      int sizei = localintdofs.Size();
      
      FlatArray<double> el2ifweight(sizei, lh);
      for (int k = 0; k < sizei; k++)
	el2ifweight[k] = fabs (elmat(localintdofs[k],
				     localintdofs[k]));

	/*
	if (typeid(SCAL) == typeid(double))
	  el2ifweight[k] = fabs (elmat(localintdofs[k],
				       localintdofs[k]));
	else
	  el2ifweight[k] = 1;
	*/

      
      FlatMatrix<SCAL> a = elmat.Rows(localwbdofs).Cols(localwbdofs) | lh;
      FlatMatrix<SCAL> b = elmat.Rows(localwbdofs).Cols(localintdofs) | lh;
      FlatMatrix<SCAL> c = elmat.Rows(localintdofs).Cols(localwbdofs) | lh;
      FlatMatrix<SCAL> d = elmat.Rows(localintdofs).Cols(localintdofs) | lh;
      FlatMatrix<SCAL> het (sizew, sizei, lh);
      FlatMatrix<SCAL> he (sizei, sizew, lh);
	  

      if (sizei)
	{      
          RegionTimer regcompute (timer3);
          NgProfiler::AddThreadFlops (timer3, TaskManager::GetThreadId(),
                                      sizei*sizei*sizei + 2*sizei*sizei*sizew);

          CalcInverse (d);  // , INVERSE_LIB::INV_NGBLA);
          
	  if (sizew)
	    {
              /*
	      he = SCAL(0.0);

	      he -= d*c   | Lapack;
	      a += b*he   | Lapack;
              */
              he = -d*c;
              a += b*he;
	      //R * E
	      for (size_t k = 0; k < sizei; k++)
		he.Row(k) *= el2ifweight[k]; 

	      if (!bfa->SymmetricStorage())
		{
                  /*
		  het = SCAL(0.0);
		  het -= b*d  | Lapack;
		  */
                  het = -b*d;
                  
		  //E * R^T
		  for (size_t l = 0; l < sizei; l++)
		    het.Col(l) *= el2ifweight[l];
		}
	    }
	  //R * A_ii^(-1) * R^T
	  for (size_t k = 0; k < sizei; k++) d.Row(k) *= el2ifweight[k]; 
	  for (size_t l = 0; l < sizei; l++) d.Col(l) *= el2ifweight[l]; 
	}

      RegionTimer regadd (timer2);
      
      FlatArray<int> wbdofs(localwbdofs.Size(), lh);
      FlatArray<int> intdofs(localintdofs.Size(), lh);   
      wbdofs = dnums[localwbdofs];
      intdofs = dnums[localintdofs];


      for (int j = 0; j < intdofs.Size(); j++)
        weight[intdofs[j]] += el2ifweight[j];
      
      sparse_harmonicext->AddElementMatrix(intdofs,wbdofs,he);
      
      if (!bfa->SymmetricStorage())
        sparse_harmonicexttrans->AddElementMatrix(wbdofs,intdofs,het);
      
      sparse_innersolve -> AddElementMatrix(intdofs,intdofs,d);

      dynamic_pointer_cast<SparseMatrix<SCAL,TV,TV>>(pwbmat)
        ->AddElementMatrix(wbdofs,wbdofs,a);
      if (coarse)
        dynamic_pointer_cast<Preconditioner>(inv)->AddElementMatrix(wbdofs,a,id,lh);
    }


    
    void Finalize()
    {
      static Timer timer ("BDDC Finalize");
      RegionTimer reg(timer);

      // auto fes = bfa->GetFESpace();
      int ndof = fes->GetNDof();      

      if (!local)
	AllReduceDofData (weight, MPI_SUM, fes->GetParallelDofs());
      
      ParallelFor (weight.Size(),
                   [&] (size_t i)
                   {
                     if (weight[i]) weight[i] = 1.0/weight[i];
                   });

      ParallelFor (sparse_innersolve->Height(),
                   [&] (size_t i)
                   {
                     FlatArray<int> cols = sparse_innersolve -> GetRowIndices(i);
                     FlatVector<SCAL> values = sparse_innersolve->GetRowValues(i);
                     double wi = weight[i];
                     for (int j = 0; j < cols.Size(); j++)
                       values(j) *= wi * weight[cols[j]];
                   }, TasksPerThread(5));

      ParallelFor (sparse_harmonicext->Height(),
                   [&] (size_t i)
                   {
                     sparse_harmonicext->GetRowValues(i) *= weight[i];                     
                   }, TasksPerThread(5));
      
      if (!bfa->SymmetricStorage())
        {
          ParallelFor (// sparse_harmonicexttrans->Height(),
                       sparse_harmonicexttrans->GetBalancing(),
                       [&] (size_t i)
                       {
                         FlatArray<int> rowind = sparse_harmonicexttrans->GetRowIndices(i);
                         FlatVector<SCAL> values = sparse_harmonicexttrans->GetRowValues(i);
                         for (int j = 0; j < rowind.Size(); j++)
                           values[j] *= weight[rowind[j]];
                       }, TasksPerThread(5));
        }
      
      // now generate wire-basked solver

      if (block)
	{
          if (coarse)
            throw Exception("combination of coarse and block not implemented! ");

	  //Smoothing Blocks
	  Flags flags;
	  flags.SetFlag("eliminate_internal");
	  flags.SetFlag("subassembled");
	  cout << IM(3) << "call Create Smoothing Blocks of " << bfa->GetFESpace()->GetName() << endl;
	  shared_ptr<Table<int>> blocks = (bfa->GetFESpace()->CreateSmoothingBlocks(flags));
	  cout << IM(3) << "has blocks" << endl << endl;
	  // *testout << "blocks = " << endl << blocks << endl;
	  // *testout << "pwbmat = " << endl << *pwbmat << endl;
	  cout << IM(3) << "call block-jacobi inverse" << endl;

	  inv = sparse_pwbmat->CreateBlockJacobiPrecond(blocks, 0, 0, 0);
	  // inv = dynamic_cast<BaseSparseMatrix*> (pwbmat)->CreateJacobiPrecond(wb_free_dofs);

	  cout << IM(3) << "has inverse" << endl << endl;
	  // *testout << "blockjacobi = " << endl << *inv << endl;
	  
	  //Coarse Grid of Wirebasket
	  cout << IM(3) << "call directsolverclusters inverse" << endl;
	  auto clusters = bfa->GetFESpace()->CreateDirectSolverClusters(flags);
	  cout << IM(3) << "has clusters" << endl << endl;
	  
	  cout << IM(3) << "call coarse wirebasket grid inverse" << endl;
	  inv_coarse = pwbmat->InverseMatrix(clusters);
	  cout << IM(3) << "has inverse" << endl << endl;
	  
	  tmp = make_shared<VVector<>>(ndof);
	  tmp2 = make_shared<VVector<>>(ndof);
	}
      else
	{
	  if (bfa->GetFESpace()->IsParallel() && !local)
	    {
	      shared_ptr<ParallelDofs> pardofs = bfa->GetFESpace()->GetParallelDofs();

	      pwbmat = make_shared<ParallelMatrix> (pwbmat, pardofs);
	      pwbmat -> SetInverseType (inversetype);

#ifdef HYPRE
	      if (hypre)
		inv = make_shared<HyprePreconditioner> (*pwbmat, wb_free_dofs);
	      else
#endif
                if (coarse)
                {
                  dynamic_pointer_cast<Preconditioner>(inv) -> FinalizeLevel(pwbmat.get());
                }
                else
                  inv = pwbmat -> InverseMatrix (wb_free_dofs);

	      tmp = make_unique<ParallelVVector<TV>>(pardofs);
	      innersolve =
                ComposeOperators(make_shared<ParallelMatrix> (innersolve, pardofs, pardofs, C2D),
                                 make_shared<CumulationOperator> (pardofs));
	      innersolve =
                ComposeOperators(make_shared<CumulationOperator> (pardofs), innersolve);
              
	      harmonicext =
                ComposeOperators(make_shared<CumulationOperator> (pardofs),
                                 make_shared<ParallelMatrix> (harmonicext, pardofs, pardofs, C2D));
              
	      if (harmonicexttrans)
		harmonicexttrans = ComposeOperators(make_shared<ParallelMatrix> (harmonicexttrans, pardofs, pardofs, C2D),
                                                    make_shared<CumulationOperator> (pardofs));
	    }
	  else
	    {

              int cntfreedofs = 0;
              for (int i = 0; i < wb_free_dofs->Size(); i++)
                if (wb_free_dofs->Test(i)) cntfreedofs++;

              if (coarse)
              {
                cout << IM(3) << "call wirebasket preconditioner finalize ( with " << cntfreedofs
                     << " free dofs out of " << pwbmat->Height() << " )" << endl;
                // throw Exception("combination of coarse and block not implemented! ");
                dynamic_pointer_cast<Preconditioner>(inv) -> FinalizeLevel(pwbmat.get());
              }
              else
              {
                cout << IM(3) << "call wirebasket inverse ( with " << cntfreedofs
                     << " free dofs out of " << pwbmat->Height() << " )" << endl;
                inv = pwbmat->InverseMatrix(wb_free_dofs);
              }
	      cout << IM(3) << "has inverse" << endl;
	      tmp = make_unique<VVector<TV>>(ndof);
	    }
	}
      if(inv_coarse)
        GetMemoryTracer().Track(*inv_coarse, "CoarseInverse");
      if(inv)
        GetMemoryTracer().Track(*inv, "Inverse");
    }

    ~BDDCMatrix()  { } 


    AutoVector CreateRowVector() const override { return bfa->GetMatrix().CreateColVector(); }
    AutoVector CreateColVector() const override { return bfa->GetMatrix().CreateRowVector(); }

    int VHeight() const override { return bfa->GetMatrix().VHeight(); }
    int VWidth() const override { return bfa->GetMatrix().VHeight(); }

    
    void Mult (const BaseVector & x, BaseVector & y) const override
    {
      static Timer timer ("Apply BDDC preconditioner");
      static Timer timerifs ("Apply BDDC preconditioner - apply ifs");
      static Timer timerwb ("Apply BDDC preconditioner - wb solve");
      static Timer timerharmonicext ("Apply BDDC preconditioner - harmonic extension");
      static Timer timerharmonicexttrans ("Apply BDDC preconditioner - harmonic extension trans");
      

      RegionTimer reg (timer);

      x.Distribute();
      y = x;
      
      timerharmonicexttrans.Start();

      if (bfa->SymmetricStorage())
	y += Transpose(*harmonicext) * x; 
      else
	y += *harmonicexttrans * x;

      timerharmonicexttrans.Stop();

      timerwb.Start();
      *tmp = 0;
      if (block)
	{
          if (coarse)
            throw Exception("combination of coarse and block not implemented! ");
	  if (true) //GS
	    {
	      dynamic_cast<BaseBlockJacobiPrecond*>(inv.get())->GSSmoothResiduum (*tmp, y, *tmp2 ,1);
	      
	      if (inv_coarse)
		*tmp += (*inv_coarse) * *tmp2; 
	      dynamic_cast<BaseBlockJacobiPrecond*>(inv.get())->GSSmoothBack (*tmp, y);
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
      // tmp->Distribute(); 
      *tmp += *innersolve * x;
      timerifs.Stop();

      timerharmonicext.Start();

      // tmp->Cumulate();
      y = *tmp;
      y += *harmonicext * *tmp;

      timerharmonicext.Stop();

      y.Cumulate();
    }

    shared_ptr<BaseMatrix> CreateDeviceMatrix() const override
    {
      shared_ptr<BaseMatrix> op;
      if (bfa->SymmetricStorage())
        op = harmonicext * inv * TransposeOperator (harmonicext) + innersolve;
      else
        op = harmonicext * inv * harmonicexttrans + innersolve;
      
      return op->CreateDeviceMatrix();
    }
    
  };









  template <class SCAL, class TV = SCAL>
  class NGS_DLL_HEADER BDDCPreconditioner : public Preconditioner
  {
    shared_ptr<S_BilinearForm<SCAL>> bfa;
    shared_ptr<FESpace> fes;
    shared_ptr<BDDCMatrix<SCAL,TV>> pre;
    shared_ptr<BitArray> freedofs;
    string inversetype;
    string coarsetype;
    bool block, hypre;
  public:
    BDDCPreconditioner (shared_ptr<BilinearForm> abfa, const Flags & aflags,
                        const string aname = "bddcprecond")
      : Preconditioner (abfa, aflags, aname)
    {
      bfa = dynamic_pointer_cast<S_BilinearForm<SCAL>> (abfa);
      // bfa -> SetPreconditioner (this);
      inversetype = flags.GetStringFlag("inverse", "sparsecholesky");
      coarsetype = flags.GetStringFlag("coarsetype", "none");
      if(coarsetype=="myamg_hcurl")
	(dynamic_pointer_cast<HCurlHighOrderFESpace>(bfa->GetFESpace()))->DoCouplingDofUpgrade(false);
      if (flags.GetDefineFlag("refelement")) Exception ("refelement - BDDC not supported");
      block = flags.GetDefineFlag("block");
      hypre = flags.GetDefineFlag("usehypre");
      // pre = NULL;
      fes = bfa->GetFESpace();
    }


    BDDCPreconditioner (const PDE & pde, const Flags & aflags, const string & aname)
      : BDDCPreconditioner (pde.GetBilinearForm (aflags.GetStringFlag ("bilinearform")), 
                            aflags, aname)
    { ; }
    /*
      : Preconditioner (&pde, aflags, aname)
    {
    bfa = dynamic_pointer_cast<S_BilinearForm<SCAL>>  (pde.GetBilinearForm (aflags.GetStringFlag ("bilinearform")));
      const_cast<S_BilinearForm<SCAL>*> (bfa) -> SetPreconditioner (this);
      inversetype = flags.GetStringFlag("inverse", "sparsecholesky");
      coarsetype = flags.GetStringFlag("coarsetype", "none");
      if (flags.GetDefineFlag("refelement")) Exception ("refelement - BDDC not supported");
      block = flags.GetDefineFlag("block");
      hypre = flags.GetDefineFlag("usehypre");
      pre = NULL;
    }
    */


    virtual ~BDDCPreconditioner()
    {
      ; // delete pre;
    }
    
    virtual void InitLevel (shared_ptr<BitArray> _freedofs) 
    {
      freedofs = _freedofs;
      pre = make_shared<BDDCMatrix<SCAL,TV>>(bfa, flags, inversetype, coarsetype, block, hypre);
      pre -> SetHypre (hypre);
      GetMemoryTracer().Track(*pre, "pre");
    }

    virtual void FinalizeLevel (const BaseMatrix *)
    {
      pre -> Finalize();
      if (test) Test();
      timestamp = bfa->GetTimeStamp();      
    }


    virtual void AddElementMatrix (FlatArray<int> dnums,
				   const FlatMatrix<SCAL> & elmat,
				   ElementId id, 
				   LocalHeap & lh);

    virtual void Update ()
    {
      if (timestamp < bfa->GetTimeStamp())
        throw Exception("A BDDC preconditioner must be defined before assembling");
    }  

    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa->GetMatrix();
    }

    virtual const BaseMatrix & GetMatrix() const
    {
      if (!pre)
        ThrowPreconditionerNotReady();        
      return *pre;
    }

    virtual void CleanUpLevel ()
    {
      /*
      delete pre;
      pre = NULL;
      */
      pre.reset();
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
  AddElementMatrix (FlatArray<int> dnums,
		    const FlatMatrix<SCAL> & elmat,
		    ElementId id, 
		    LocalHeap & lh)
  {
    // auto fes = bfa->GetFESpace();
    int used = 0;
    for (int i : Range(dnums))
      if (IsRegularDof(dnums[i]) && freedofs->Test(dnums[i])) used++;
    
    FlatArray<int> compress(used, lh);
    int cnt = 0;
    for (int i : Range(dnums))
      if (IsRegularDof(dnums[i]) && freedofs->Test(dnums[i]))
        compress[cnt++] = i;

    FlatArray<int> hdnums(used, lh);
    FlatMatrix<SCAL> helmat(used,used, lh);
    hdnums = dnums[compress];
    helmat = elmat.Rows(compress).Cols(compress);
    
    if (L2Norm (helmat) != 0)
      pre -> AddMatrix(helmat, hdnums, id, lh);
  }
  






  template class BDDCPreconditioner<double>;
  template class BDDCPreconditioner<double, Complex>;


  template <>
  shared_ptr<Preconditioner> RegisterPreconditioner<BDDCPreconditioner<double>>::
  CreateBF(shared_ptr<BilinearForm> bfa, const Flags & flags, const string & name)
  {
    // cout << "complex bddc ? " << bfa->GetFESpace()->IsComplex() << endl;
    if (bfa->GetFESpace()->IsComplex())
      return make_shared<BDDCPreconditioner<Complex>> (bfa, flags, name);
    else
      return make_shared<BDDCPreconditioner<double>> (bfa, flags, name);
  }

  static RegisterPreconditioner<BDDCPreconditioner<double> > initpre ("bddc");
  static RegisterPreconditioner<BDDCPreconditioner<Complex> > initpre2 ("bddcc");
  static RegisterPreconditioner<BDDCPreconditioner<double,Complex> > initpre3 ("bddcrc");
}
