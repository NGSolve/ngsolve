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
		const string & ainversetype, bool ablock, bool aebe)
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
      Table<double> el2ifweight(ifcnt);

      BitArray ifdof(ndof);
      BitArray wbdof(ndof);
      ifdof.Clear();
      wbdof.Clear();

      Vector<> weight(fes.GetNDof());
      weight = 0;

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
		{
		  wbdof.Set(dnums[j]);
		  el2wbdofs[i][lwbcnt++] = dnums[j];
		}
	      else
		{
		  if (typeid(SCAL) == typeid(double))
		    el2ifweight[i][lifcnt] = fabs ((*elmats[i])(j,j));
		  else
		    el2ifweight[i][lifcnt] = 0.5;
		  weight[dnums[j]] += el2ifweight[i][lifcnt];

		  ifdof.Set(dnums[j]);
		  el2ifdofs[i][lifcnt++] = dnums[j];
		}
	    }
        }


#ifdef PARALLEL
      // accumulate weight

      if (typeid(SCAL) == typeid(double))
	{
	  ParallelVVector<double> pv_weight (fes.GetNDof(), &fes.GetParallelDofs(), DISTRIBUTED);
	  
	  if (ntasks > 1 && id == 0)
	    {
	      pv_weight.Cumulate();  
	    }
	  else
	    {
	      for (int i = 0; i < weight.Size(); i++)
		pv_weight(i) = weight[i];
	      pv_weight.Cumulate();  
	      for (int i = 0; i < weight.Size(); i++)
		weight[i] = pv_weight(i);
	    }
	}
      else
	{
	  ParallelVVector<Complex> pv_weight (fes.GetNDof(), &fes.GetParallelDofs(), DISTRIBUTED);
	  
	  if (ntasks > 1 && id == 0)
	    {
	      pv_weight.Cumulate();  
	    }
	  else
	    {
	      for (int i = 0; i < weight.Size(); i++)
		pv_weight(i) = weight[i];
	      pv_weight.Cumulate();  
	      for (int i = 0; i < weight.Size(); i++)
		weight[i] = pv_weight(i).real();
	    }
	}
#endif

      if (id == 0)
	cout << ( bfa.IsSymmetric() ? "symmetric" : "non-symmetric") << endl;

      if (!bfa.IsSymmetric())
	{
	  harmonicexttrans = 
	    new SparseMatrix<SCAL,TV,TV>(ndof, el2wbdofs, el2ifdofs, false);
	  harmonicexttrans -> AsVector() = 0.0;
	}

      innersolve = bfa.IsSymmetric() 
	? new SparseMatrixSymmetric<SCAL,TV>(ndof, el2ifdofs)
	: new SparseMatrix<SCAL,TV,TV>(ndof, el2ifdofs, el2ifdofs, bfa.IsSymmetric());
      innersolve->AsVector() = 0.0;

      harmonicext =
	new SparseMatrix<SCAL,TV,TV>(ndof, el2ifdofs, el2wbdofs, false);
      harmonicext->AsVector() = 0.0;

      pwbmat = bfa.IsSymmetric() 
	? new SparseMatrixSymmetric<SCAL,TV>(ndof, el2wbdofs)
	: new SparseMatrix<SCAL,TV,TV>(ndof, el2wbdofs, el2wbdofs, bfa.IsSymmetric());
      pwbmat -> AsVector() = 0.0;
      pwbmat -> SetInverseType (inversetype);
      
      timer1.Stop();
      timer2.Start();

      if (ntasks == 1 || id > 0)
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

	    for (int k = 0; k < dnums.Size(); k++)
	      {
		COUPLING_TYPE ct = fes.GetDofCouplingType(dnums[k]);	      
		if (ct == WIREBASKET_DOF)
		  localwbdofs.Append (k);
		else
		  localintdofs.Append (k);
	      }

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
		
		he -= d*c   | Lapack;
		a += b*he   | Lapack;
		
		//R * E
		for (int k = 0; k < sizei; k++)
		  he.Row(k) *= el2ifweight[i][k] / weight[el2ifdofs[i][k]];

		dynamic_cast<SparseMatrix<SCAL,TV,TV>*> (harmonicext)->AddElementMatrix(interfacedofs,wirebasketdofs,he);
	    
		if (!bfa.IsSymmetric())
		  {
		    Matrix<SCAL> het (sizew, sizei);
		
		    het = SCAL(0.0);
		    LapackMultAddAB(b,d,-1.0,het);
		
		    //E * R^T
		    for (int l = 0; l < sizei; l++)
		      het.Col(l) *= el2ifweight[i][l] / weight[ el2ifdofs[i][l] ];

		    dynamic_cast<SparseMatrix<SCAL,TV,TV> *> (harmonicexttrans)
		      ->AddElementMatrix(wirebasketdofs,interfacedofs,het);
		  }
	    
		//R * A_ii^(-1) * R^T
		for (int k = 0; k < sizei; k++)
		  d.Row(k) *= el2ifweight[i][k] / weight[el2ifdofs[i][k]];
		for (int l = 0; l < sizei; l++)
		  d.Col(l) *= el2ifweight[i][l] / weight[el2ifdofs[i][l]];	    
		
		dynamic_cast<SparseMatrix<SCAL,TV,TV> *> (innersolve)
		  -> AddElementMatrix(interfacedofs,interfacedofs,d);
	      }
	  
	    dynamic_cast<SparseMatrix<SCAL,TV,TV>*>(pwbmat)
	      ->AddElementMatrix(wirebasketdofs,wirebasketdofs,a);
	  }

      timer2.Stop();
      timer3.Start();

      for (int i = 0; i < elmats.Size(); i++)
	{
	  delete elmats[i];
	  delete eldnums[i];
	}

      elmats.SetSize(0);
      eldnums.SetSize(0);
      
      
      // cout << "matrix filed" << endl;


      free_dofs = new BitArray (ndof);
      free_dofs->Clear();

      *free_dofs = wbdof;
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
	  if (ntasks > 1)
	    {
	      MyMPI_Barrier();
	      /*
		inv = new MasterInverse<SCAL> (dynamic_cast<SparseMatrixTM<SCAL>&> (*pwbmat), 
		free_dofs, &bfa.GetFESpace().GetParallelDofs());
	      */

	      /*
		inv = new ParallelMumpsInverse<SCAL> (dynamic_cast<SparseMatrixTM<SCAL>&> (*pwbmat), 
		free_dofs, NULL, 
		&bfa.GetFESpace().GetParallelDofs());
	      */

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
      static Timer timeretc ("Apply BDDC preconditioner - etc");
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












  class BDDCMatrixRefElement : public BaseMatrix
  {
    const S_BilinearForm<double> & bfa;
    Array<int> restrict;
    Array<int> multiple;
    BaseMatrix * inv, *inv_coarse;
    
    Table<int> *wbdofs, *dcdofs, *internaldofs, *externaldofs;
    Table<int> *globwbdofs, *globdcdofs, *globinternaldofs, *globexternaldofs;

    Table<int> *interface_dofs;
    Table<double> *weighting;

    Array<int> elclassnr;
    Table<int> * els_of_class;
    Array<Matrix<>*> invdc_ref, extwb_ref, schurwb_ref;
    Array<Matrix<>*> inv_int_ref, extension_int_ref, extension_int_ref_trans;
    int ne, ndof;
    
    string inversetype;
    bool print;

    bool block;
    BitArray * free_dofs_wb;

  public:
    BDDCMatrixRefElement (const S_BilinearForm<double> & abfa, const string & ainversetype)
      : bfa(abfa), inversetype (ainversetype)
    {
      cout << "BDDC MatrixRefElement" << endl;
      print = true;
      block = true;
      LocalHeap lh(500000000, "BDDC - MatrixRefEleemnt");

      if (print)
	*testout << "BDDC - Constructor" << endl;

      const FESpace & fes = bfa.GetFESpace();
      const MeshAccess & ma = fes.GetMeshAccess();

      ne = ma.GetNE();
      ndof = fes.GetNDof();
      Array<int> cnt(ne), cntwb(ne), cntdc(ne), cntinternal(ne), cntexternal(ne);

      Array<int> gwbdofs(ndof);
      Array<int> wbdnums, dnums, extdnums, dcdnums, vnums;
      gwbdofs = 0;

      for (int i = 0; i < ne; i++)
        {
          fes.GetDofNrs (i, dnums);
	  fes.GetDofNrs (i, extdnums, EXTERNAL_DOF);	  
	  fes.GetDofNrs (i, wbdnums, WIREBASKET_DOF);	  

          cnt[i] = dnums.Size();
          cntwb[i] = wbdnums.Size();
          cntdc[i] = extdnums.Size()-wbdnums.Size();
          cntinternal[i] = dnums.Size()-extdnums.Size();
          cntexternal[i] = extdnums.Size();

          for (int j = 0; j < wbdnums.Size(); j++)
            gwbdofs[wbdnums[j]] = 1;
        }

      if (print)
	{
	  *testout << "Glocal wirebasket - dofs:" << endl;
	  *testout << "ndof = " << ndof << endl;
	  *testout << gwbdofs << endl;
	}

      // Table<int> dcdofnrs(cnt);   // discontinuous dofs
      Table<int> dcdofnrs(cntexternal);   // discontinuous dofs

      wbdofs = new Table<int>(cntwb);
      dcdofs = new Table<int>(cntdc);
      internaldofs = new Table<int>(cntinternal);
      externaldofs = new Table<int>(cntexternal);

      globwbdofs = new Table<int>(cntwb);
      globdcdofs = new Table<int>(cntdc);
      globinternaldofs = new Table<int>(cntinternal);
      globexternaldofs = new Table<int>(cntexternal);      

      int totdofs  = fes.GetNDof();
      for (int i = 0; i < ne; i++)
        {
	  //*testout << "element " << i << endl;

          fes.GetDofNrs (i, dnums);
          fes.GetDofNrs (i, extdnums, EXTERNAL_DOF);


	  *testout << "dnums = " << endl << dnums << endl;
	  *testout << "ext dnums = " << endl << extdnums << endl;

	  int nint = 0, next = 0;

	  for (int j = 0; j < dnums.Size(); j++)
	    {
	      if (extdnums.Contains (dnums[j]))
		{
		  (*externaldofs)[i][next] = j;
		  (*globexternaldofs)[i][next] = dnums[j];
		  next++;
		}
	      else
		{
		  (*internaldofs)[i][nint] = j;
		  (*globinternaldofs)[i][nint] = dnums[j];
		  nint++;
		}
	    }
	  
          dcdnums.SetSize (extdnums.Size());

          for (int j = 0; j < extdnums.Size(); j++)
            {
	      if (extdnums[j] == -1)
		dcdnums[j] = -1;
	      else
		{
		  if (gwbdofs[extdnums[j]])
		    dcdnums[j] = extdnums[j];
		  else
		    {
		      dcdnums[j] = totdofs;
		      totdofs++;
		    }
		}
	      /*
		dcdnums[j] = extdnums[j];
		if (dcdnums[j] == -1) continue;
		if (gwbdofs[dcdnums[j]]) continue;
		dcdnums[j] = totdofs;
		totdofs++;
	      */
            }

	  *testout << "dcdnums = " << endl << dcdnums << endl;

          int wbi = 0, dci = 0;
          for (int j = 0; j < extdnums.Size(); j++)
            {
              if (extdnums[j] == -1) continue;
              if (gwbdofs[extdnums[j]])
                {
                  (*wbdofs)[i][wbi] = j;
                  (*globwbdofs)[i][wbi] = extdnums[j];
                  wbi++;
                }
              else
                {
                  (*dcdofs)[i][dci] = j;
                  (*globdcdofs)[i][dci] = dcdnums[j];
                  dci++;
                }
            }
          
          for (int j = 0; j < extdnums.Size(); j++)
            dcdofnrs[i][j] = dcdnums[j];
        }


      if (print)
	{
	  *testout << "dcdofnrs = " << endl;
	  *testout << dcdofnrs << endl;
	}

      int firstdcdof = fes.GetNDof();      
      
      MatrixGraph graph(firstdcdof, *globwbdofs,  *globwbdofs, 1);
      SparseMatrixSymmetric<double> & dcmat 
	= *new SparseMatrixSymmetric<double> (graph, 1);

      dcmat = 0.0;
      // dcmat.SetInverseType ("sparsecholesky");
      dcmat.SetInverseType (inversetype);
      *testout << "dcmat, init = " << dcmat << endl;

      elclassnr.SetSize(ne);
      int maxclass = ne; // 32;

      invdc_ref.SetSize(maxclass);
      extwb_ref.SetSize(maxclass);
      schurwb_ref.SetSize(maxclass);
      
      inv_int_ref.SetSize(maxclass);
      extension_int_ref.SetSize(maxclass);
      extension_int_ref_trans.SetSize(maxclass);

      invdc_ref = NULL;
      extwb_ref = NULL;
      schurwb_ref = NULL;
      inv_int_ref = NULL;
      extension_int_ref = NULL;

      for (int i = 0; i < ne; i++)
        {
          HeapReset hr(lh);
          ma.GetElVertices(i, vnums);

          int classnr = 0;

          int sort[4] = { 0, 1, 2, 3 };
          if (vnums.Size() == 4)
            {
              if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 1; }
              if (vnums[sort[2]] > vnums[sort[3]]) { Swap (sort[2], sort[3]); classnr += 2; }
              if (vnums[sort[0]] > vnums[sort[2]]) { Swap (sort[0], sort[2]); classnr += 4; }
              if (vnums[sort[1]] > vnums[sort[3]]) { Swap (sort[1], sort[3]); classnr += 8; }
              if (vnums[sort[1]] > vnums[sort[2]]) { Swap (sort[1], sort[2]); classnr += 16; }
            }
          else
            {
              if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 1; }
              if (vnums[sort[1]] > vnums[sort[2]]) { Swap (sort[1], sort[2]); classnr += 2; }
              if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 4; }
            }
	  // classnr = i;

          elclassnr[i] = classnr;
          if (invdc_ref[classnr]) continue;

          cout << "assemble class " << classnr << " " << flush;

	  if (print)
	    *testout << "assemble class " << classnr << endl;

          fes.GetDofNrs (i, dnums);
          fes.GetDofNrs (i, extdnums, EXTERNAL_DOF);
	  
          Matrix<> elmat(dnums.Size());
          Matrix<> partelmat(dnums.Size());
          elmat = 0.0;

          // ElementTransformation eltrans;
          // ma.GetElementTransformation (i, eltrans, lh);
	  ElementTransformation & eltrans = ma.GetTrafo (i, false, lh);

          const FiniteElement & fel = fes.GetFE(i, lh);

          for (int j = 0; j < bfa.NumIntegrators(); j++)
            {
	      if (bfa.GetIntegrator(j) -> BoundaryForm()) continue;
              bfa.GetIntegrator(j)->CalcElementMatrix (fel, eltrans, partelmat, lh);
              elmat += partelmat;
            }

	  

	  if (print)
	    *testout << "elmat = " << endl << elmat << endl;
	  
          int nwb = (*wbdofs)[i].Size();
          int ndc = (*dcdofs)[i].Size();
          int nint = (*internaldofs)[i].Size();
          int next = (*externaldofs)[i].Size();

	  cout << " wb = " << nwb << " nifs = " << ndc << " loc = " << nint << endl;
	  // external/internal extension


	  *testout << "nint = " << nint << ", next = " << next << ", nwb = " << nwb << ", ndc = " << ndc << endl;
	  *testout << "wbdofs = " << (*wbdofs)[i] << endl;

          Matrix<> matee(next);
          Matrix<> matie(nint, next);
          Matrix<> extie(nint, next);
          Matrix<> matii(nint);

          for (int j = 0; j < next; j++)
            for (int k = 0; k < next; k++)
              matee(j,k) = elmat((*externaldofs)[i][j], (*externaldofs)[i][k]);
          for (int j = 0; j < nint; j++)
            for (int k = 0; k < next; k++)
              matie(j,k) = elmat((*internaldofs)[i][j], (*externaldofs)[i][k]);
          for (int j = 0; j < nint; j++)
            for (int k = 0; k < nint; k++)
              matii(j,k) = elmat((*internaldofs)[i][j], (*internaldofs)[i][k]);


	  *testout << "matee = " << endl << matee << endl;
	  *testout << "matie = " << endl << matie << endl;
	  *testout << "matii = " << endl << matii << endl;

          LapackInverse (matii);
          extie = matii * matie;
          matee -= Trans (matie) * extie;
	  

	  *testout << "schur-ext = " << endl << matee << endl;

          Matrix<> matwbwb(nwb);
          Matrix<> matdcwb(ndc, nwb);
          Matrix<> extdcwb(ndc, nwb);
          Matrix<> matdcdc(ndc);

          for (int j = 0; j < nwb; j++)
            for (int k = 0; k < nwb; k++)
              matwbwb(j,k) = matee((*wbdofs)[i][j], (*wbdofs)[i][k]);
          for (int j = 0; j < ndc; j++)
            for (int k = 0; k < nwb; k++)
              matdcwb(j,k) = matee((*dcdofs)[i][j], (*wbdofs)[i][k]);
          for (int j = 0; j < ndc; j++)
            for (int k = 0; k < ndc; k++)
              matdcdc(j,k) = matee((*dcdofs)[i][j], (*dcdofs)[i][k]);

          LapackInverse (matdcdc);
          extdcwb = matdcdc * matdcwb;
          matwbwb -= Trans (matdcwb) * extdcwb;

	  *testout << "schur-wb = " << endl << matwbwb << endl;
          
	  invdc_ref[classnr] = new Matrix<> (ndc, ndc);
	  extwb_ref[classnr] = new Matrix<> (ndc, nwb);
	  schurwb_ref[classnr] = new Matrix<> (nwb, nwb);
	  
	  *(invdc_ref[classnr]) = matdcdc;
	  *(extwb_ref[classnr]) = extdcwb;
	  *(schurwb_ref[classnr]) = matwbwb;

	  inv_int_ref[classnr] = new Matrix<> (nint, nint);
	  *(inv_int_ref[classnr]) = matii;
	  extension_int_ref[classnr] = new Matrix<> (nint, next);
	  *(extension_int_ref[classnr]) = extie;

	  extension_int_ref_trans[classnr] = new Matrix<> (next, nint);
	  *(extension_int_ref_trans[classnr]) = Trans(extie);
        }

      TableCreator<int> creator;
      for ( ; !creator.Done(); creator++)
	{
	  for (int i = 0; i < ne; i++)
	    creator.Add (elclassnr[i], i);
	}
      els_of_class = creator.GetTable();

      for (int i = 0; i < ne; i++)
        {
          FlatArray<int> dofs = (*globwbdofs)[i]; 
          FlatMatrix<> matwbwb = *schurwb_ref[elclassnr[i]];
          for (int k = 0; k < dofs.Size(); k++)
            for (int l = 0; l < dofs.Size(); l++)
              if (dofs[k] != -1 && dofs[l] != -1 && dofs[k] >= dofs[l])
                dcmat(dofs[k], dofs[l]) += matwbwb(k,l);
        }

      if (print)
	*testout << "disconnected matrix: " << endl 
		 << dcmat << endl;


      restrict.SetSize(totdofs);
      restrict = -1;
      for (int i = 0; i < ne; i++)
        {
          // fes.GetDofNrs (i, dnums);
          fes.GetDofNrs (i, extdnums, EXTERNAL_DOF);

          for (int j = 0; j < extdnums.Size(); j++)
            {
              if (extdnums[j] != -1)
                restrict[dcdofnrs[i][j]] = extdnums[j];
            }
        }

      multiple.SetSize (fes.GetNDof());
      multiple = 0;
      for (int i = 0; i < restrict.Size(); i++)
        if (restrict[i] != -1)
          multiple[restrict[i]]++;

      for (int i = 0; i < firstdcdof; i++)
        if (gwbdofs[i] == 0)
          dcmat(i,i) = 1;


      const BitArray & free_dofs = *fes.GetFreeDofs();

      
      // Table<int> *interface_dofs;
      // Table<double> *weighting;
      Array<int> cnt_interface(ne);
      cnt_interface = 0;

      for (int i = 0; i < ne; i++)
	cnt_interface[i] = (*globdcdofs)[i].Size();
      
      interface_dofs = new Table<int> (cnt_interface);
      weighting = new Table<double> (cnt_interface);

      for (int i = 0; i < ne; i++)
	{
	  FlatArray<int> dc = (*globdcdofs)[i];
	  for (int j = 0; j < dc.Size(); j++)
	    {
	      int rest = restrict[dc[j]];

	      if (rest == -1) cout << "rest == -1" << endl;
	      if (multiple[rest] <= 0) cout << "multiplie = " << multiple[rest] << endl;

	      (*interface_dofs)[i][j] = rest;

	      if (free_dofs.Test(rest))
		(*weighting)[i][j] = 1.0 / multiple[rest];
	      else
		(*weighting)[i][j] = 0.0;
	    }
	}

      cout << "tables filled" << endl;




      for (int i = 0; i < restrict.Size(); i++)
	if (restrict[i] != -1 && !free_dofs.Test(restrict[i]))
	  restrict[i] = -1;


      free_dofs_wb = new BitArray (free_dofs.Size());

      free_dofs_wb -> Clear();
      for (int i = 0; i < firstdcdof; i++)
        if (gwbdofs[i] != 0)
	  free_dofs_wb->Set(i);
      free_dofs_wb -> And (free_dofs);

      




      cout << "call inverse, type = " << inversetype << endl;
      cout << "dim = " << dcmat.Height() << endl;



      if (block)
	{
	  Flags flags;
	  flags.SetFlag("eliminate_internal");
	  flags.SetFlag("subassembled");
	  Table<int> & blocks = *(bfa.GetFESpace().CreateSmoothingBlocks(flags));
	  /*
	    for (int i = 0; i < blocks.Size(); i++){
	    for (int j = 0; j < blocks[i].Size(); j++)
	    blocks[i][j] = wbdofs[blocks[i][j]];
	    }
	  */
	  cout << "call block-jacobi inverse" << endl;
	  inv = dcmat.CreateBlockJacobiPrecond(blocks, 0, 0, 0);      

	  cout << "dc.height = " << dcmat.Height() << endl;
	  cout << "inv.height = " << inv->Height() << endl;

	  cout << "has inverse" << endl;
	  cout << "call directsolverclusters inverse" << endl;
	  Array<int> & clusters = *(bfa.GetFESpace().CreateDirectSolverClusters(flags));
	  // 	*testout << " clusters = \n " << clusters << endl;
	  /*
	    Array<int> & condclusters = *new Array<int>(nglobalwbdof);
	    for (int i=0; i< clusters.Size(); i++)
	    condclusters[wbdofs[i]] = clusters[i];
	    inv_coarse = dcmat.InverseMatrix(&condclusters);
	  */
	  // 	*testout << " condclusters = \n " << condclusters << endl;
	  inv_coarse = dcmat.InverseMatrix(&clusters);
	  // tmp = new VVector<>(nglobalwbdof);

	  for (int i = 0; i < clusters.Size(); i++)
	    if (clusters[i] && !free_dofs_wb->Test(i))
	      cout << "wrong dirichet bc" << endl;
	}
      else
	{
	  cout << "call inverse" << endl;
	  inv = dcmat.InverseMatrix(free_dofs_wb);
	}

      // inv = dcmat.InverseMatrix(&free_dofs);
      return;




#define sometestsxx
#ifdef sometests

      Array<int> & directblocks = *new Array<int> (free_dofs.Size());
      directblocks = 0;
      // int nv = ma.GetNV();
      int ned = ma.GetNEdges();
      int nfa = ma.GetNFaces();
      for (int i = 0; i < directblocks.Size(); i++)
	if (free_dofs.Test(i))
	  directblocks[i] = 1;

      
      directblocks = 0;

      /*
	for (int i = 0; i < ne; i++)
	{
	fes.GetDofNrs (i, dnums);
	if (!dnums.Size()) continue;

	if (free_dofs.Test(dnums[0]))
	directblocks[dnums[0]] = 1;

	for (int j = 0; j < dnums.Size(); j++)
	if (free_dofs.Test(dnums[j]))
	directblocks[dnums[j]] = 1;
	}
      */

      for (int i = 0; i < -nfa; i++)
	{
	  fes.GetFaceDofNrs (i, dnums);
	  if (!dnums.Size()) continue;

	  if (free_dofs.Test(dnums[0]))
	    directblocks[dnums[0]] = 1;

	  /*
	    for (int j = 1; j < dnums.Size(); j++)
	    if (free_dofs.Test(dnums[j]))
	    directblocks[dnums[j]] = 1;
	  */
	}



      for (int i = 0; i < ned; i++)
	{
	  fes.GetEdgeDofNrs (i, dnums);
	  if (!dnums.Size()) continue;

	  if (free_dofs.Test(dnums[0]))
	    directblocks[dnums[0]] = 1;
	  for (int j = 1; j < dnums.Size(); j++)
	    if (free_dofs.Test(dnums[j]))
	      directblocks[dnums[j]] = i+2;
	}
      
      cout << "directblocks.Size() = " << directblocks.Size() << endl;
      cout << "freedofs.Size() = " << free_dofs.Size() << endl;
      cout << "mat.size = " << dcmat.Height() << endl;

      /*
       *testout << "directblocks,1 = " << endl << directblocks << endl;
      
       directblocks = 0;
       for (int i = 0; i < directblocks.Size(); i++)
       if (free_dofs.Test(i)) 
       directblocks[i] = 1;

       *testout << "dcmat = " << endl << dcmat << endl;
       *testout << "directblocks = " << endl << directblocks << endl;
       *testout << "freedofs = " << endl << free_dofs << endl;
       */

      inv = dcmat.InverseMatrix(&directblocks);
      // inv = dcmat.InverseMatrix(&free_dofs);
      cout << "complete" << endl;

      EigenSystem eigen(dcmat, *inv);
      eigen.Calc();
      eigen.PrintEigenValues (cout);
      if (print)
	{
	  *testout << "restrict = " << endl  << restrict << endl;
	  *testout << "dcmat = " << endl << dcmat << endl;
	  *testout << "freedofs = " << endl << *fes.GetFreeDofs() << endl;
	}
#endif
    }




    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      // cout << "Multadd, |x| = " << L2Norm(x) << ", |y| = " << L2Norm(y) << endl;
      static Timer timer ("Apply BDDC preconditioner");
      static Timer timertransx ("Apply BDDC preconditioner, transx");
      static Timer timertransx1a ("Apply BDDC preconditioner, transx - lapack");
      static Timer timertransx1b ("Apply BDDC preconditioner, transx - read");
      static Timer timertransx1c ("Apply BDDC preconditioner, transx - write");
      // static Timer timertransx1d = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transx - indices");

      static Timer timertransy ("Apply BDDC preconditioner, transy");
      static Timer timertransyr ("Apply BDDC preconditioner, transy - read");
      static Timer timertransyw ("Apply BDDC preconditioner, transy - write");
      // static Timer timertransyi = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transy - ind");


      static Timer timerrest ("Apply BDDC preconditioner, restrict");
      static Timer timerdc ("Apply BDDC preconditioner, decoupled");
      static Timer timerext ("Apply BDDC preconditioner, extend");
      static Timer timerinner ("Apply BDDC preconditioner, inner");
      static Timer timersolve ("Apply BDDC preconditioner, solve");
      // static Timer timeretc  ("Apply BDDC preconditioner, etc");
      RegionTimer reg (timer);


      FlatVector<> fx = x.FVDouble();
      // FlatVector<> fy = y.FVDouble();

      VVector<> transx(ndof);
      VVector<> transy(ndof);
      
      int maxmem = 0;
      for (int cl = 0; cl < els_of_class->Size(); cl++)
	{
	  FlatArray<int> els = (*els_of_class)[cl];
	  if (!els.Size()) continue;
	  
	  const Matrix<> & ext = (*extension_int_ref_trans[cl]);
	  int mem = els.Size() * max (ext.Height(), ext.Width());

	  maxmem = max(maxmem,mem);
	}
      
      Vector<> hmem1(maxmem);
      Vector<> hmem2(maxmem);
      

      {
	RegionTimer reg(timertransx);
	transx.FV() = fx;  



	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & ext = (*extension_int_ref_trans[cl]);

	    FlatMatrix<> vim(els.Size(), ext.Width(), &hmem1(0));
	    FlatMatrix<> vem(els.Size(), ext.Height(), &hmem2(0));

	    timertransx.AddFlops (ext.Width()*ext.Height()*els.Size());
	    timertransx1a.AddFlops (ext.Width()*ext.Height()*els.Size());

	    timertransx1b.Start();

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> dint = (*globinternaldofs)[els[ii]];
		FlatVector<> vim_row = vim.Row(ii);

		for (int j = 0; j < dint.Size(); j++)
		  vim_row(j) = transx(dint[j]);
	      }

	    timertransx1b.Stop();


	    // NgProfiler::AddLoads (timertransx1b, els.Size()*vim.Width()*12);
	    // NgProfiler::AddStores (timertransx1b, els.Size()*vim.Width()*8);


	    
	    timertransx1a.Start();

	    // LapackMultABt (vim, ext, vem);
	    LapackMult (vim, Trans(ext), vem);

	    timertransx1a.Stop();

	    timertransx1c.Start();
	    for (int i = 0, ii = 0; i < els.Size(); i++)
	      {
		FlatArray<int> dext = (*globexternaldofs)[els[i]];

		for (int j = 0; j < dext.Size(); j++, ii++)
		  transx(dext[j]) -= vem(ii);
	      }

	    timertransx1c.Stop();


	    
	    // NgProfiler::AddLoads (timertransx1c, els.Size()*vem.Width()*20);
	    // NgProfiler::AddStores (timertransx1c, els.Size()*vem.Width()*8);
	  }
      }


      // restrict
      { 
	RegionTimer reg(timerrest);

	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & ext = (*extwb_ref[cl]);
	    Matrix<> v1(els.Size(), ext.Height());
	    Matrix<> v2(els.Size(), ext.Width());
	    
	    timerrest.AddFlops (v1.Width()*v2.Width()*els.Size());
		
	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> ifs = (*interface_dofs)[els[ii]];
		FlatArray<double> w = (*weighting)[els[ii]];

		for (int j = 0; j < ifs.Size(); j++)
		  v1(ii,j) = w[j] * transx(ifs[j]);
	      }

	    LapackMult (v1, ext, v2);

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> wb = (*globwbdofs)[els[ii]];
		
		for (int j = 0; j < wb.Size(); j++)
		  transx(wb[j]) -= v2(ii, j);
	      }
	  }
      }


      transy = 0.0;

  

      {
	// cout << "|lx| = " << L2Norm(lx) << flush;
	RegionTimer reg (timersolve);

	BaseVector & subx = *(transx.Range(0, inv->Height()));
	BaseVector & suby = *(transy.Range(0, inv->Height()));
	VVector<> res(inv->Height());

	if (block)
	  {
	    if (false) //GS
	      {
		// cout << "start GS" << endl;
		dynamic_cast<BaseBlockJacobiPrecond*>(inv)->GSSmoothResiduum (suby, subx, res,1);
		// cout << "GS done" << endl;

		if (inv_coarse)
		  suby += (*inv_coarse) * res; 

		// cout << "inv_coarse done" << endl;
		dynamic_cast<BaseBlockJacobiPrecond*>(inv)->GSSmoothBack (suby, subx);
	      }
	    else{ //jacobi only (old)
	      suby = (*inv) * subx;
	      suby += (*inv_coarse) * subx; 
	    }
	  }
	else
	  {
	    suby = (*inv) * subx;
	  }


	// cout << ", |ly| = " << L2Norm(ly) << endl;
      }
      

      // extend
      {
	RegionTimer reg(timerext);

	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & ext = (*extwb_ref[cl]);
	    Matrix<> v1(els.Size(), ext.Width());
	    Matrix<> v2(els.Size(), ext.Height());
	    
	    timerext.AddFlops (v1.Width()*v2.Width()*els.Size());
		
	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> wb = (*globwbdofs)[els[ii]];
		for (int j = 0; j < wb.Size(); j++)
		  v1(ii, j) = transy(wb[j]);
	      }

	    // LapackMultABt (v1, ext, v2);
	    LapackMult (v1, Trans(ext), v2);

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> ifs = (*interface_dofs)[els[ii]];
		FlatArray<double> w = (*weighting)[els[ii]];

		for (int j = 0; j < ifs.Size(); j++)
		  transy(ifs[j]) -= w[j] * v2(ii,j);
	      }
	  }
      }

 


      // solve decoupled
      {
	RegionTimer reg(timerdc);

	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & invdc = (*invdc_ref[cl]);
	    Matrix<> dc1(els.Size(), invdc.Height());
	    Matrix<> dc2(els.Size(), invdc.Height());

	    timerdc.AddFlops (invdc.Width()*invdc.Height()*els.Size());
		
	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> ifs = (*interface_dofs)[els[ii]];
		FlatArray<double> w = (*weighting)[els[ii]];

		for (int j = 0; j < ifs.Size(); j++)
		  dc1(ii,j) = w[j] * transx(ifs[j]);
	      }

	    LapackMult (dc1, invdc, dc2);

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> ifs = (*interface_dofs)[els[ii]];
		FlatArray<double> w = (*weighting)[els[ii]];

		for (int j = 0; j < ifs.Size(); j++)
		  transy(ifs[j]) += w[j] * dc2(ii,j);
	      }
	  }
      }





      {
	RegionTimer reg(timertransy);

	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & ext = (*extension_int_ref_trans[cl]);
	    Matrix<> vim(els.Size(), ext.Width());
	    Matrix<> vem(els.Size(), ext.Height());

	    timertransy.AddFlops (ext.Width()*ext.Height()*els.Size());
		


	    timertransyr.Start();

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> dext = (*globexternaldofs)[els[ii]];
		for (int j = 0; j < dext.Size(); j++)
		  vem(ii, j) = transy(dext[j]);
	      }

	    timertransyr.Stop();


	    LapackMult (vem, ext, vim);


	    timertransyw.Start();

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> dint = (*globinternaldofs)[els[ii]];
		
		for (int j = 0; j < dint.Size(); j++)
		  transy(dint[j]) -= vim(ii, j);
	      }

	    timertransyw.Stop();
	  }
      }
      

      {
	RegionTimer reg(timerinner);

	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & inv_int = (*inv_int_ref[cl]);
	    Matrix<> v1(els.Size(), inv_int.Width());
	    Matrix<> v2(els.Size(), inv_int.Width());

	    timerinner.AddFlops (v1.Width()*v2.Width()*els.Size());
		
	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> dint = (*globinternaldofs)[els[ii]];
		for (int j = 0; j < dint.Size(); j++)
		  v1(ii, j) = transx(dint[j]);
	      }

	    LapackMult (v1, inv_int, v2);

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> dint = (*globinternaldofs)[els[ii]];
		for (int j = 0; j < dint.Size(); j++)
		  transy(dint[j]) += v2(ii, j);
	      }

	  }
      }

      y += s * transy;
    }

  };













  template <class SCAL, class TV>
  BDDCPreconditioner<SCAL, TV> ::
  BDDCPreconditioner (const PDE & pde, const Flags & aflags)
    : Preconditioner (&pde, aflags)
  {
    bfa = dynamic_cast<const S_BilinearForm<SCAL>*>(pde.GetBilinearForm (aflags.GetStringFlag ("bilinearform", NULL)));
    const_cast<S_BilinearForm<SCAL>*> (bfa) -> SetPreconditioner (this);
    inversetype = flags.GetStringFlag("inverse", "sparsecholesky");
    refelement = flags.GetDefineFlag("refelement");
    block = flags.GetDefineFlag("block");
    ebe = flags.GetDefineFlag("ebe");
    pre = NULL;
  }

  template <class SCAL, class TV>
  BDDCPreconditioner<SCAL, TV> ::
  ~BDDCPreconditioner()
  {
    delete pre;
  }


  template <class SCAL, class TV>
  void BDDCPreconditioner<SCAL, TV> :: CleanUpLevel ()
  {
    delete pre;
    pre = NULL;
  }


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

      /*
       *testout << "entry " << hnr << endl;
       *testout << "dnums = " << *eldnums[hnr] << endl;
       *testout << "elmat = " << *elmats[hnr] << endl;
       */

      if (L2Norm (*elmats[hnr]) == 0)
	{
	  /*
	    static int cnt_warn = 0;
	    cnt_warn++;
	    if (cnt_warn < 5)
	    cout << "msg: BDDC::got 0 matrix, delete it" << endl;
	    if (cnt_warn == 5) cout << "(suppress more msgs)" << endl;
	  */
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

    if (refelement)
      {
	cerr << "ref element bddc currently not supported" << endl;
	// pre = new BDDCMatrixRefElement(*bfa, inversetype);
      }
    else
      {
	if (elmats.Size() || (id == 0 && ntasks > 1))
	  {
	    delete pre;
	    pre = new BDDCMatrix<SCAL,TV>(*bfa, elmats, eldnums, inversetype, block,ebe);
	  }
	else
	  {
	    cerr << "don't update precond, since I don't have matrices" << endl;
	  }
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
