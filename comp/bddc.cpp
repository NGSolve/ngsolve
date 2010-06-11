#include <comp.hpp>
#include <solve.hpp>

#include <bddc.hpp>
namespace ngcomp
{

  class BDDCMatrix : public BaseMatrix
  {
    const ElementByElement_BilinearForm<double> & bfa;
    Array<int> restrict;
    Array<int> multiple;
    BaseMatrix * inv;
    string inversetype;
    BitArray * free_dofs;
    const MeshAccess & ma;
  public:
    BDDCMatrix (const ElementByElement_BilinearForm<double> & abfa, const string & inversetype)
      : ma(abfa.GetFESpace().GetMeshAccess()),bfa(abfa) 
    {
      // LocalHeap lh(10000000);
      const FESpace & fes = bfa.GetFESpace();
      const MeshAccess & ma = fes.GetMeshAccess();

      int ne = ma.GetNE();
      Array<int> cnt(ne);

      Array<int> wbdofs(fes.GetNDof()), lwbdofs, dnums;
      wbdofs = 0;

      for (int i = 0; i < ne; i++)
        {
          fes.GetDofNrs (i, dnums, EXTERNAL_DOF);	  
          cnt[i] = dnums.Size();

	  fes.GetDofNrs(i,lwbdofs,WIREBASKET_DOF);
          for (int j = 0; j < lwbdofs.Size(); j++)
            wbdofs[lwbdofs[j]] = 1;
	  *testout << "dnums = " << endl << dnums << endl;
	  *testout << "lwbdofs = " << endl << lwbdofs << endl;
	  
        }

      Table<int> dcdofs(cnt);   // discontinuous dofs

      

      int firstdcdof = fes.GetNDof();
      // for (int i = 0; i < wbdofs.Size(); i++)
      // if (wbdofs[i]) firstdcdof = i+1;
      int ndcdof = fes.GetNDof();

      for (int i = 0; i < ne; i++)
        {
          fes.GetDofNrs (i, dnums, EXTERNAL_DOF);	  

          for (int j = 0; j < dnums.Size(); j++)
            {
              if (dnums[j] == -1) continue;
              if (wbdofs[dnums[j]]) continue;
              dnums[j] = ndcdof;
              firstdcdof++;
	      ndcdof++;
            }

          for (int j = 0; j < dnums.Size(); j++)
            dcdofs[i][j] = dnums[j];
        }

      // *testout << "dcdofs = " << endl << dcdofs << endl;

      restrict.SetSize(ndcdof);
      restrict = -1;
      for (int i = 0; i < ne; i++)
        {
          fes.GetDofNrs (i, dnums, EXTERNAL_DOF);	  

          for (int j = 0; j < dnums.Size(); j++)
            {
              if (dnums[j] != -1)
                restrict[dcdofs[i][j]] = dnums[j];
            }
        }

      // *testout << "restrict = " << endl << restrict << endl;

      multiple.SetSize (fes.GetNDof());
      multiple = 0;
      for (int i = 0; i < restrict.Size(); i++)
        if (restrict[i] != -1)
          multiple[restrict[i]]++;


      cout << "now build graph" << endl;

      MatrixGraph graph(firstdcdof, dcdofs, bfa.IsSymmetric());

      cout << "now allocate matrix" << endl;

      SparseMatrix<double> * pdcmat;
      if (bfa.IsSymmetric()){
	pdcmat = new SparseMatrixSymmetric<double>(graph,1);
	cout << "symmetric" << endl;
      }
      else{
	pdcmat = new SparseMatrix<double>(graph,1);
	cout << "nonsymmetric" << endl;
      }
      SparseMatrix<double>& dcmat=*pdcmat;

      dcmat.SetInverseType (inversetype);
      
      cout << "have matrix" << endl;

      for (int i = 0; i < ne; i++)
        {
          FlatMatrix<> elmat = 
            dynamic_cast<const ElementByElementMatrix<double>&> (bfa.GetMatrix()) . GetElementMatrix (i);

          FlatArray<int> dofs = dcdofs[i];
          for (int k = 0; k < dofs.Size(); k++)
            for (int l = 0; l < dofs.Size(); l++)
              if (dofs[k] != -1 && dofs[l] != -1 && (!bfa.IsSymmetric() || dofs[k] >= dofs[l]))
                dcmat(dofs[k], dofs[l]) += elmat(k,l);
        }
      
      cout << "matrix filled" << endl;

      for (int i = 0; i < restrict.Size(); i++)
	if (restrict[i] == -1)
	  dcmat(i,i) = 1;

      // *testout << "dcmat = " << endl << dcmat << endl;
      *testout << "freedofs = " << *fes.GetFreeDofs() << endl;


      free_dofs = new BitArray (ndcdof);
      free_dofs->Clear();
      
      for (int i = 0; i < ndcdof; i++)
	{
	  if (restrict[i] != -1)
	    if (fes.GetFreeDofs()->Test(restrict[i]))
	      free_dofs->Set(i);
	}
      

      // *testout << "dcmat = " << endl << dcmat << endl;
      // *testout << "restrict = " << endl << restrict << endl;
      cout << "call inverse" << endl;
      inv = dcmat.InverseMatrix(free_dofs);
      cout << "has inverse" << endl;
      *testout << "inverse = " << (*inv) << endl;
      delete pdcmat;
    }

    ~BDDCMatrix(){delete inv;}
    
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      FlatVector<> fx = x.FVDouble();
      FlatVector<> fy = y.FVDouble();

      VVector<> lx(restrict.Size());
      VVector<> ly(restrict.Size());

      for (int i = 0; i < restrict.Size(); i++)
        if (restrict[i] != -1)
          lx(i) = fx(restrict[i]) / multiple[restrict[i]];

      ly = (*inv) * lx;
            
      fy = 0.0;
      for (int i = 0; i < restrict.Size(); i++)
        if (restrict[i] != -1)
          fy(restrict[i]) += ly(i) / multiple[restrict[i]];
    }

  };












  class BDDCMatrixRefElement : public BaseMatrix
  {
    const S_BilinearForm<double> & bfa;
    Array<int> restrict;
    Array<int> multiple;
    BaseMatrix * inv;
    
    Table<int> *wbdofs, *dcdofs, *internaldofs, *externaldofs;
    Table<int> *globwbdofs, *globdcdofs, *globinternaldofs, *globexternaldofs;

    Array<int> elclassnr;
    Array<Matrix<>*> invdc_ref, extwb_ref, schurwb_ref;
    Array<Matrix<>*> inv_int_ref, extension_int_ref, extension_int_ref_trans;
    int ne, ndof;
    
    string inversetype;
    bool print;

  public:
    BDDCMatrixRefElement (const S_BilinearForm<double> & abfa, const string & ainversetype)
      : bfa(abfa), inversetype (ainversetype)
    {
      cout << "BDDC MatrixRefElement" << endl;
      print = true;
      LocalHeap lh(500000000);

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
	  *testout << "element " << i << endl;

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
      
      MatrixGraph graph(firstdcdof, *globwbdofs, 1);
      SparseMatrixSymmetric<double> dcmat(graph, 1);
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

          cout << "assemble class " << classnr << endl;

	  if (print)
	    *testout << "assemble class " << classnr << endl;

          fes.GetDofNrs (i, dnums);
          fes.GetDofNrs (i, extdnums, EXTERNAL_DOF);
	  
          Matrix<> elmat(dnums.Size());
          Matrix<> partelmat(dnums.Size());
          elmat = 0.0;

          ElementTransformation eltrans;
          ma.GetElementTransformation (i, eltrans, lh);
          const FiniteElement & fel = fes.GetFE(i, lh);

          for (int j = 0; j < bfa.NumIntegrators(); j++)
            {
	      if (bfa.GetIntegrator(j) -> BoundaryForm()) continue;
              bfa.GetIntegrator(j)->AssembleElementMatrix (fel, eltrans, partelmat, lh);
              elmat += partelmat;
            }

	  

	  if (print)
	    *testout << "elmat = " << endl << elmat << endl;
	  
          int nwb = (*wbdofs)[i].Size();
          int ndc = (*dcdofs)[i].Size();
          int nint = (*internaldofs)[i].Size();
          int next = (*externaldofs)[i].Size();

	  
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
      for (int i = 0; i < restrict.Size(); i++)
	if (restrict[i] != -1 && !free_dofs.Test(restrict[i]))
	  restrict[i] = -1;


      cout << "call inverse, type = " << inversetype << endl;
      cout << "dim = " << dcmat.Height() << endl;
      /*
      inv = dcmat.InverseMatrix(&free_dofs);
      return;
      */
#define sometests
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

      for (int i = 0; i < nfa; i++)
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
	      directblocks[dnums[j]] = i;
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



    /*
      virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
      {
      static int timer = NgProfiler::CreateTimer ("Apply BDDC preconditioner");
      static int timerrest = NgProfiler::CreateTimer ("Apply BDDC preconditioner, restrict");
      static int timerdc = NgProfiler::CreateTimer ("Apply BDDC preconditioner, decoupled");
      static int timerext = NgProfiler::CreateTimer ("Apply BDDC preconditioner, extend");
      NgProfiler::RegionTimer reg (timer);

      LocalHeap lh(100000);
      FlatVector<> fx = x.FVDouble();
      FlatVector<> fy = y.FVDouble();

      VVector<> lx(restrict.Size());
      VVector<> ly(restrict.Size());

      lx = 0.0;
      for (int i = 0; i < restrict.Size(); i++)
      if (restrict[i] != -1)
      lx(i) = fx(restrict[i]) / multiple[restrict[i]];


      // restrict
      NgProfiler::StartTimer (timerrest);
      #pragma omp parallel
      {
      LocalHeap lh(100000);
      #pragma omp for
      for (int i = 0; i < ne; i++)
      {
      HeapReset hr(lh);
      FlatArray<int> wbd = (*globwbdofs)[i];
      FlatArray<int> dcd = (*globdcdofs)[i];

      FlatVector<> wb(wbd.Size(), lh);
      FlatVector<> dc(dcd.Size(), lh);

      for (int j = 0; j < dcd.Size(); j++)
      dc(j) = lx(dcd[j]);
          
      // wb = Trans(*extwb[i]) * dc;
      wb = Trans(*extwb_ref[elclassnr[i]]) * dc;
          
      for (int j = 0; j < wbd.Size(); j++)
      lx(wbd[j]) -= wb(j);

      NgProfiler::AddFlops (timer, 2*wbd.Size()*dcd.Size());
      NgProfiler::AddFlops (timerrest, wbd.Size()*dcd.Size());
      }
      }
      NgProfiler::StopTimer (timerrest);

      ly = 0.0;
      ly = (*inv) * lx;


      // solve decoupled
      NgProfiler::StartTimer (timerdc);
      #pragma omp parallel
      {
      LocalHeap lh(100000);
      #pragma omp for
      for (int i = 0; i < ne; i++)
      {
      HeapReset hr(lh);
      FlatArray<int> dcd = (*globdcdofs)[i];
            
      FlatVector<> dc1(dcd.Size(), lh);
      FlatVector<> dc2(dcd.Size(), lh);
            
      for (int j = 0; j < dcd.Size(); j++)
      dc1(j) = lx(dcd[j]);
            
      // dc2 = Trans(*invdc[i]) * dc1;
      dc2 = (*invdc_ref[elclassnr[i]]) * dc1;
            
      for (int j = 0; j < dcd.Size(); j++)
      ly(dcd[j]) = dc2(j);
            
      NgProfiler::AddFlops (timer, sqr(dcd.Size()));
      NgProfiler::AddFlops (timerdc, sqr(dcd.Size()));
      }
      }
      NgProfiler::StopTimer (timerdc);

      // extend
      NgProfiler::StartTimer (timerext);
      #pragma omp parallel
      {
      LocalHeap lh(100000);
      #pragma omp for
      for (int i = 0; i < ne; i++)
      {
      HeapReset hr(lh);
      FlatArray<int> wbd = (*globwbdofs)[i];
      FlatArray<int> dcd = (*globdcdofs)[i];

      FlatVector<> wb(wbd.Size(), lh);
      FlatVector<> dc(dcd.Size(), lh);
      for (int j = 0; j < wbd.Size(); j++)
      wb(j) = ly(wbd[j]);

      // dc = (*extwb[i]) * wb;
      dc = (*extwb_ref[elclassnr[i]]) * wb;

      for (int j = 0; j < dcd.Size(); j++)
      ly(dcd[j]) -= dc(j);
      NgProfiler::AddFlops (timerext, wbd.Size()*dcd.Size());
      }
      }
      NgProfiler::StopTimer (timerext);

      fy = 0.0;
      for (int i = 0; i < restrict.Size(); i++)
      if (restrict[i] != -1)
      fy(restrict[i]) += ly(i) / multiple[restrict[i]];
      }
    */



    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      static int timer = NgProfiler::CreateTimer ("Apply BDDC preconditioner");
      static int timertransx = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transx");
      static int timertransy = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transy");
      static int timerrest = NgProfiler::CreateTimer ("Apply BDDC preconditioner, restrict");
      static int timerdc = NgProfiler::CreateTimer ("Apply BDDC preconditioner, decoupled");
      static int timerext = NgProfiler::CreateTimer ("Apply BDDC preconditioner, extend");
      static int timerinner = NgProfiler::CreateTimer ("Apply BDDC preconditioner, inner");
      static int timersolve = NgProfiler::CreateTimer ("Apply BDDC preconditioner, solve");
      NgProfiler::RegionTimer reg (timer);



      LocalHeap lh(1000014);
      FlatVector<> fx = x.FVDouble();
      FlatVector<> fy = y.FVDouble();

      VVector<> transx(ndof);
      VVector<> transy(ndof);

      VVector<> lx(restrict.Size());
      VVector<> ly(restrict.Size());
      
      {
	NgProfiler::RegionTimer reg(timertransx);
	transx = 1.0 * x;

	/*

	#pragma omp parallel
	{
	LocalHeap slh = lh.Split();
	#pragma omp for
	for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	for (int i = 0; i < ne; i++)
	{
	if (elclassnr[i] != cl) continue;

	HeapReset hr(slh);
	FlatArray<int> dint = (*globinternaldofs)[i];
	FlatArray<int> dext = (*globexternaldofs)[i];
	      
	FlatVector<> vi(dint.Size(), slh);
	FlatVector<> ve(dext.Size(), slh);
	      
	for (int j = 0; j < dint.Size(); j++)
	vi(j) = transx(dint[j]);
	      
	ve = (*extension_int_ref_trans[elclassnr[i]]) * vi;

	#pragma omp critical (bddc_transx)
	{
	for (int j = 0; j < ve.Size(); j++)
	transx(dext[j]) -= ve(j);
	NgProfiler::AddFlops (timertransx, dint.Size()*dext.Size());
	}
	}
	}
	*/

	for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	  {
	    int cnt = 0;
	    for (int i = 0; i < ne; i++)
	      if (elclassnr[i] == cl) cnt++;
	    if (cnt == 0) continue;

	    const Matrix<> & ext = (*extension_int_ref_trans[cl]);
	    Matrix<> vim(cnt, ext.Width());
	    Matrix<> vem(cnt, ext.Height());

	    NgProfiler::AddFlops (timertransx, ext.Width()*ext.Height()*cnt);
		
	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dint = (*globinternaldofs)[i];
		
		for (int j = 0; j < dint.Size(); j++)
		  vim(cnt, j) = transx(dint[j]);
		cnt++;
	      }

	    // vem = vim * Trans(ext);
	    LapackMultABt (vim, ext, vem);

	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dext = (*globexternaldofs)[i];

		for (int j = 0; j < dext.Size(); j++)
		  transx(dext[j]) -= vem(cnt,j);
		cnt++;
	      }
	  }
      }

      lx = 0.0;
      for (int i = 0; i < restrict.Size(); i++)
        if (restrict[i] != -1)
          lx(i) = transx(restrict[i]) / multiple[restrict[i]];


      // restrict
      {
	NgProfiler::RegionTimer reg(timerrest);

	// #pragma omp parallel
	{
	  LocalHeap slh = lh.Split();
	  // #pragma omp for
	  for (int i = 0; i < ne; i++)
	    {
	      HeapReset hr(slh);
	      FlatArray<int> wbd = (*globwbdofs)[i];
	      FlatArray<int> dcd = (*globdcdofs)[i];
	      
	      FlatVector<> wb(wbd.Size(), slh);
	      FlatVector<> dc(dcd.Size(), slh);
	    
	      for (int j = 0; j < dcd.Size(); j++)
		dc(j) = lx(dcd[j]);
	      
	      wb = Trans(*extwb_ref[elclassnr[i]]) * dc;
	      
	      // #pragma omp critical (bddc_restrict)
	      {
		for (int j = 0; j < wbd.Size(); j++)
		  lx(wbd[j]) -= wb(j);
	      }
	    }
	}
      }
      
      {
	NgProfiler::RegionTimer reg (timersolve);
	ly = 0.0;
	ly = (*inv) * lx;
      }
	

      // solve decoupled
      {
	NgProfiler::RegionTimer reg(timerdc);

	/*


	#pragma omp parallel
	{
	LocalHeap slh = lh.Split();
	#pragma omp for
	for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	for (int i = 0; i < ne; i++)
	{
	if (elclassnr[i] != cl) continue;
	HeapReset hr(slh);
	FlatArray<int> dcd = (*globdcdofs)[i];
	      
	FlatVector<> dc1(dcd.Size(), slh);
	FlatVector<> dc2(dcd.Size(), slh);
	      
	for (int j = 0; j < dcd.Size(); j++)
	dc1(j) = lx(dcd[j]);
	      
	// dc2 = Trans(*invdc[i]) * dc1;
	dc2 = (*invdc_ref[elclassnr[i]]) * dc1;
	      
	for (int j = 0; j < dcd.Size(); j++)
	ly(dcd[j]) += s * dc2(j);
	      
	NgProfiler::AddFlops (timerdc, sqr (dcd.Size()));
	}
	}
	*/

	for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	  {
	    int cnt = 0;
	    for (int i = 0; i < ne; i++)
	      if (elclassnr[i] == cl) cnt++;
	    if (cnt == 0) continue;

	    const Matrix<> & invdc = (*invdc_ref[cl]);
	    Matrix<> dc1(cnt, invdc.Height());
	    Matrix<> dc2(cnt, invdc.Height());

	    NgProfiler::AddFlops (timerdc, invdc.Width()*invdc.Height()*cnt);
		
	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dcd = (*globdcdofs)[i];

		for (int j = 0; j < dcd.Size(); j++)
		  dc1(cnt, j) = lx(dcd[j]);
		cnt++;
	      }

	    // vem = vim * Trans(ext);
	    LapackMultAB (dc1, invdc, dc2);
	    // dc2 = dc1 * invdc;

	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dcd = (*globdcdofs)[i];
		
		for (int j = 0; j < dcd.Size(); j++)
		  ly(dcd[j]) += s*dc2(cnt, j);
		cnt++;
	      }

	  }




      }

      // extend
      {
	NgProfiler::RegionTimer reg(timerext);

	// #pragma omp parallel
	{
	  LocalHeap slh = lh.Split();
	  // #pragma omp for
	  for (int i = 0; i < ne; i++)
	    {
	      HeapReset hr(slh);
	      FlatArray<int> wbd = (*globwbdofs)[i];
	      FlatArray<int> dcd = (*globdcdofs)[i];
	      
	      FlatVector<> wb(wbd.Size(), slh);
	      FlatVector<> dc(dcd.Size(), slh);
	      for (int j = 0; j < wbd.Size(); j++)
		wb(j) = ly(wbd[j]);
	      
	      // dc = (*extwb[i]) * wb;
	      dc = (*extwb_ref[elclassnr[i]]) * wb;
	      
	      for (int j = 0; j < dcd.Size(); j++)
		ly(dcd[j]) -= dc(j);
	    }
	}
      }

      transy = 0.0;
      for (int i = 0; i < restrict.Size(); i++)
	if (restrict[i] != -1)
	  transy(restrict[i]) += s * ly(i) / multiple[restrict[i]];

      
      {
	NgProfiler::RegionTimer reg(timertransy);
	/*
	  #pragma omp parallel
	  {
	  LocalHeap slh = lh.Split();
	  #pragma omp for
	  for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	  for (int i = 0; i < ne; i++)
	  {
	  if (elclassnr[i] != cl) continue;
	  HeapReset hr(slh);
	  FlatArray<int> dint = (*globinternaldofs)[i];
	  FlatArray<int> dext = (*globexternaldofs)[i];
	      
	  FlatVector<> vi(dint.Size(), slh);
	  FlatVector<> ve(dext.Size(), slh);
	      
	  for (int j = 0; j < dext.Size(); j++)
	  ve(j) = transy(dext[j]);
	      
	  vi = (*extension_int_ref[elclassnr[i]]) * ve;
	      
	  for (int j = 0; j < dint.Size(); j++)
	  transy(dint[j]) -= vi(j);
	  }
	  }
	*/



	for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	  {
	    int cnt = 0;
	    for (int i = 0; i < ne; i++)
	      if (elclassnr[i] == cl) cnt++;
	    if (cnt == 0) continue;

	    const Matrix<> & ext = (*extension_int_ref_trans[cl]);
	    Matrix<> vim(cnt, ext.Width());
	    Matrix<> vem(cnt, ext.Height());

	    NgProfiler::AddFlops (timertransy, ext.Width()*ext.Height()*cnt);
		
	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dext = (*globexternaldofs)[i];

		for (int j = 0; j < dext.Size(); j++)
		  vem(cnt, j) = transy(dext[j]);
		cnt++;
	      }

	    // vem = vim * Trans(ext);
	    LapackMultAB (vem, ext, vim);

	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dint = (*globinternaldofs)[i];
		
		for (int j = 0; j < dint.Size(); j++)
		  transy(dint[j]) -= vim(cnt, j);
		cnt++;
	      }

	  }




      }
      
      y += s * transy;

      {
	NgProfiler::RegionTimer reg(timerinner);
	// #pragma omp parallel
	{
	  LocalHeap slh = lh.Split();
	  // #pragma omp for
	  for (int i = 0; i < ne; i++)
	    {
	      HeapReset hr(slh);
	      FlatArray<int> dint = (*globinternaldofs)[i];
	      
	      FlatVector<> v1(dint.Size(), slh);
	      FlatVector<> v2(dint.Size(), slh);
            
	      for (int j = 0; j < dint.Size(); j++)
		v1(j) = fx(dint[j]);
            
	      v2 = (*inv_int_ref[elclassnr[i]]) * v1;
            
	      for (int j = 0; j < dint.Size(); j++)
		fy(dint[j]) += s * v2(j);
	    }
	}
      }
    }

  };













  template <class SCAL>
  BDDCPreconditioner<SCAL> ::
  BDDCPreconditioner (const PDE * pde, const Flags & aflags, const string aname)
    : Preconditioner (pde, aflags, aname)
  {
    bfa = dynamic_cast<const S_BilinearForm<SCAL>*>(pde->GetBilinearForm (aflags.GetStringFlag ("bilinearform", NULL)));
    inversetype = flags.GetStringFlag("inverse", "sparsecholesky");
    refelement = flags.GetDefineFlag("refelement");
  }


  template <class SCAL>
  void BDDCPreconditioner<SCAL> ::
  Update ()
  {
    cout << "update bddc" << endl;
    if (refelement)
      pre = new BDDCMatrixRefElement(*bfa, inversetype);
    else
      pre = new BDDCMatrix(dynamic_cast<const ElementByElement_BilinearForm<double>&> (*bfa), inversetype);

    if (test) Test();
    
  }  

  template class BDDCPreconditioner<double>;
    
  namespace bddc_cpp
  {
   
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetPreconditionerClasses().AddPreconditioner ("bddc", BDDCPreconditioner<double>::Create);
    }
    
    Init init;
  }


}
