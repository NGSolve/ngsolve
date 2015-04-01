/* *************************************************************************/
/* File:   sparseldl.cc                                                    */
/* Author: Joachim Schoeberl                                               */
/* Date:   18. Nov. 2001                                                   */
/* *************************************************************************/

#include <la.hpp>

#include "concurrentqueue.h" 


typedef moodycamel::ConcurrentQueue<int> TQueue; 
typedef moodycamel::ProducerToken TPToken; 
typedef moodycamel::ConsumerToken TCToken; 
 




namespace ngla
{



  template <class TM>
  void SetIdentity( TM &identity )
  {
    for (int i=0; i<identity.Width(); i++ ) identity(i,i) = 1;
  }
  void SetIdentity( double &identity ) { identity = 1; }
  void SetIdentity( Complex &identity ) { identity = 1; }




  template <class TM, class TV_ROW, class TV_COL>
  SparseCholesky<TM, TV_ROW, TV_COL> :: 
  SparseCholesky (const SparseMatrix<TM, TV_ROW, TV_COL> & a, 
		  const BitArray * ainner,
		  const Array<int> * acluster,
		  bool allow_refactor)
    : SparseFactorization (a, ainner, acluster), mat(a)
  { 
    static Timer t("SparseCholesky - total");
    static Timer ta("SparseCholesky - allocate");
    static Timer tf("SparseCholesky - fill factor");
    RegionTimer reg(t);
    // (*testout) << "matrix = " << a << endl;
    // (*testout) << "diag a = ";
    // for ( int i=0; i<a.Height(); i++ ) (*testout) << i << ", " << a(i,i) << endl;

    int n = a.Height();
    height = n;

    int printstat = 0;
    
    if (printstat)
      cout << IM(4) << "Minimal degree ordering: N = " << n << endl;
    
    clock_t starttime, endtime;
    starttime = clock();
    
    mdo = new MinimumDegreeOrdering (n);
    
    if (!inner && !cluster)
      for (int i = 0; i < n; i++)
	for (int j = 0; j < a.GetRowIndices(i).Size(); j++)
	  {
	    int col = a.GetRowIndices(i)[j];
	    if (col <= i)
	      mdo->AddEdge (i, col);
	  }

    else if (inner)
      for (int i = 0; i < n; i++)
	for (int j = 0; j < a.GetRowIndices(i).Size(); j++)
	  {
	    int col = a.GetRowIndices(i)[j];
	    if (col <= i)
	      if ( (inner->Test(i) && inner->Test(col)) || i==col)
		mdo->AddEdge (i, col);
	  }

    else 
      for (int i = 0; i < n; i++)
	{
	  FlatArray<int> row = a.GetRowIndices(i);
	  for (int j = 0; j < row.Size(); j++)
	    {
	      int col = row[j];
	      if (col <= i)
		if ( ( ((*cluster)[i] == (*cluster)[col]) && (*cluster)[i]) ||
		     i == col )
		  mdo->AddEdge (i, col);
	    }
	}
    
    for (int i = 0; i < n; i++)
      if (a.GetPositionTest (i,i) == numeric_limits<size_t>::max())
	{
	  mdo->AddEdge (i, i);
	  *testout << "add unsused position " << i << endl;
	}

    if (printstat)
      cout << IM(4) << "start ordering" << endl;
    
    // mdo -> PrintCliques ();
    mdo->Order();
    
    endtime = clock();
    if (printstat)
      cout << IM(4) << "ordering time = "
	   << double (endtime - starttime) / CLOCKS_PER_SEC 
	   << " secs" << endl;
    
    starttime = endtime;
    
    if (printstat)
      cout << IM(4) << "," << flush;
    ta.Start();
    Allocate (mdo->order,  mdo->vertices, &mdo->blocknr[0]);
    ta.Stop();

    tf.Start();
    delete mdo;
    mdo = 0;

    diag.SetSize(n);
    lfact.SetSize (nze);
    lfact = TM(0.0);     // first touch
    
    endtime = clock();
    if (printstat)
      (cout) << "allocation time = "
	     << double (endtime - starttime) / CLOCKS_PER_SEC << " secs" << endl;
    
    starttime = endtime;

    TM id;
    id = 0.0;
    SetIdentity(id);
    
    for (int i = 0; i < n; i++)
      if (a.GetPositionTest (i,i) == numeric_limits<size_t>::max())
	SetOrig (i, i, id);

    
    if (!inner && !cluster)
      for (int i = 0; i < n; i++)
	for (int j = 0; j < a.GetRowIndices(i).Size(); j++)
	  {
	    int col = a.GetRowIndices(i)[j];
	    if (col <= i)
	      SetOrig (i, col, a.GetRowValues(i)[j]);
	  }
    
    else if (inner)
      // for (int i = 0; i < n; i++)
      ParallelFor 
        (Range(n), [&](int i)
         {
           for (int j = 0; j < a.GetRowIndices(i).Size(); j++)
             {
               int col = a.GetRowIndices(i)[j];
               if (col <= i)
                 {
                   if ( (inner->Test(i) && inner->Test(col)) )
                     SetOrig (i, col, a.GetRowValues(i)[j]);
                   else
                     if (i==col)
                       SetOrig (i, col, id);
                 }
             }
         });
    else
      for (int i = 0; i < n; i++)
	{
	  FlatArray<int> row = a.GetRowIndices(i);
	  for (int j = 0; j < row.Size(); j++)
	    {
	      int col = row[j];
	      if (col <= i)
		if ( ( (*cluster)[i] == (*cluster)[col] && (*cluster)[i]) ||
		     i == col )
                  SetOrig (i, col, a.GetRowValues(i)[j]);
	    }
	}

    tf.Stop();
    
    if (printstat)
      cout << IM(4) << "do factor " << flush;
    
    
#ifdef LAPACK
    /*
    auto aspd = dynamic_cast<const SparseMatrixSymmetricTM<double>*> (&a);
    if (aspd && aspd -> IsSPD())
    */
    if (a.IsSPD())
      FactorSPD();
    else
#endif
      Factor(); 



    for (int i = 0; i < n; i++)
      if (a.GetPositionTest (i,i) == numeric_limits<size_t>::max())
	diag[order[i]] = TM(0.0);

    if (inner)
      {
	for (int i = 0; i < n; i++)
	  if (!inner->Test(i))
	    diag[order[i]] = TM(0.0);
      }

    if (cluster)
      {
	for (int i = 0; i < n; i++)
	  if (!(*cluster)[i])
	    diag[order[i]] = TM(0.0);
      }


    if (printstat)
      cout << IM(4) << "done" << endl;
    
    endtime = clock();

    if (printstat)
      (cout) << " factoring time = " << double(endtime - starttime) / CLOCKS_PER_SEC << " sec" << endl;
  }
  

  
  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  Allocate (const Array<int> & aorder, 
	    // const Array<CliqueEl*> & cliques,
	    const Array<MDOVertex> & vertices,
	    const int * in_blocknr)
  {
    int n = aorder.Size();

    order.SetSize (n);
    blocknrs.SetSize (n);
    
    // order: now inverse map 
    for (int i = 0; i < n; i++)
      order[aorder[i]] = i;

    for (int i = 0; i < n; i++)
      blocknrs[i] = in_blocknr[i];

    *testout << "blocknrs = " << endl << blocknrs << endl;

    int cnt = 0;
    int cnt_master = 0;

    for (int i = 0; i < n; i++)
      {
	cnt += vertices[aorder[blocknrs[i]]].nconnected - (i-blocknrs[i]);
	if (blocknrs[i] == i)
	  cnt_master += vertices[aorder[i]].nconnected;
      }

    nze = cnt;

    if (n > 2000)
      cout << IM(4) << " " << cnt*sizeof(TM)+cnt_master*sizeof(int) << " Bytes " << flush;
    //cout << IM(4) <<"(cnt="<<cnt<<", sizeof(TM)="<<sizeof(TM)<< ", cnt_master=" << cnt_master << ", sizeof(int)=" << sizeof(int) <<") " << flush;


    /* 
     *testout << " Sparse Cholesky mem needed " << double(cnt*sizeof(TM)+cnt_master*sizeof(int))*1e-6 << " MBytes " << endl; 
     */  

    firstinrow.SetSize(n+1);
    firstinrow_ri.SetSize(n+1);
    rowindex2.SetSize (cnt_master);


    cnt = 0;
    cnt_master = 0;
    maxrow = 0;
    for (int i = 0; i < n; i++)
      {
	firstinrow[i] = cnt;
	int ii = aorder[i];
	int ncon = vertices[ii].nconnected;

	if (blocknrs[i] == i)
	  {
	    firstinrow_ri[i] = cnt_master;

	    for (int j = 0; j < ncon; j++)
	      rowindex2[firstinrow_ri[i]+j] = order[vertices[ii].connected[j]];

	    QuickSort (FlatArray<int> (ncon, &rowindex2[firstinrow_ri[i]]));

	    cnt_master += ncon;
	    cnt += ncon;
	    maxrow = max2 (maxrow, ncon+1);
	  }
	else
	  {
	    firstinrow_ri[i] = firstinrow_ri[i-1]+1;
	    cnt += firstinrow[i]-firstinrow[i-1]-1;
	  }
      }
    firstinrow[n] = cnt;
    firstinrow_ri[n] = cnt_master;


    for (int i = 1; i < blocknrs.Size(); i++)
      {
        if (blocknrs[i] < blocknrs[i-1]) blocknrs[i] = blocknrs[i-1];
        if (blocknrs[i] <= i-256) blocknrs[i] = i;
      }


    // cout << "finding block-dependeny ... " << endl;

    for (int i = 0; i < n; i++)
      if(blocknrs[i] == i) blocks.Append(i);
    blocks.Append(n);

    // find block dependency
    Array<int> block_of_dof(n);
    for (int i = 0; i < blocks.Size()-1; i++)
      block_of_dof[Range(blocks[i], blocks[i+1])] = i;

    DynamicTable<int> dep(blocks.Size()-1);
    for (int i = 0; i < n; i++)
      {
        auto cols = rowindex2.Range(firstinrow_ri[i], firstinrow_ri[i+1]);
        for (int j : cols)
          if (block_of_dof[i] != block_of_dof[j])
            dep.AddUnique (block_of_dof[i], block_of_dof[j]);
      }

    // generate compressed table
    TableCreator<int> creator(dep.Size());
    for ( ; !creator.Done(); creator++)
      for (int i = 0; i < dep.Size(); i++)
        for (int j : dep[i])
          creator.Add(i, j);
    block_dependency = creator.MoveTable();



    cout << "generate micro-tasks" << endl;
    // genare micro-tasks:
    Array<int> first_microtask;
    for (int i = 0; i < blocks.Size()-1; i++)
      {
        auto extdofs = BlockExtDofs (i);
        first_microtask.Append (microtasks.Size());

        if ( (extdofs.Size() == 0) && (BlockDofs(i).Size() == 1) ) continue;

        int nb = extdofs.Size() / 256 + 1;

        MicroTask mt;
        mt.blocknr = i;
        mt.solveL = true;
        mt.bblock = 0;
        mt.nbblocks = 0;
        microtasks.Append (mt);
        
        for (int j = 0; j < nb; j++)
          {
            MicroTask mt;
            mt.blocknr = i;
            mt.solveL = false;
            mt.bblock = j;
            mt.nbblocks = nb;
            microtasks.Append (mt);
          }
      }
    first_microtask.Append (microtasks.Size());
    cout << "have micro-tasks, now build tables" << endl;
    {

      TableCreator<int> creator(microtasks.Size());
      TableCreator<int> creator_trans(microtasks.Size());

      for ( ; !creator.Done(); creator++, creator_trans++)
        {
          for (int i = 0; i < first_microtask.Size()-1; i++)
            {
              for (int b = first_microtask[i]+1; b < first_microtask[i+1]; b++)
                {
                  // L to B dependency
                  creator.Add (first_microtask[i], b);
                  creator_trans.Add (b, first_microtask[i]);

                  // B to L dependency
                  for (int o : block_dependency[i])
                    {
                      creator.Add (b, first_microtask[o]);
                      creator_trans.Add (first_microtask[o], b);
                    }
                }
            }
        }

      cout << "have creator" << endl;

      micro_dependency = creator.MoveTable();
      micro_dependency_trans = creator_trans.MoveTable();

      cout << "have tables" << endl;
      cout << "dag.size = " << micro_dependency.Size() << endl;
    }
  }
  




  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  FactorNew (const SparseMatrix<TM, TV_ROW, TV_COL> & a)
  {
    if ( height != a.Height() )
      {
	cout << IM(4) << "SparseCholesky::FactorNew called with matrix of different size." << endl;
	return;
      }

    TM id;
    id = 0.0;
    SetIdentity(id);

    for (int i = 0; i < nze; i++) lfact[i] = 0.0;

    for (int i = 0; i < height; i++)
      for (int j = 0; j < a.GetRowIndices(i).Size(); j++)
	{
	  int col = a.GetRowIndices(i)[j];
	
	  if ((!inner && !cluster) || 
	      (inner && inner->Test(i) && inner->Test(col)) ||
	      (!inner && cluster && (*cluster)[i] == (*cluster)[col] && (*cluster)[i]) )
	    {
	      if ( col <= i ) SetOrig (i, col, a.GetRowValues(i)[j]);
	    }
	  else if (i == col)
	    SetOrig (i, i, id);
	}
    
    Factor(); 
  }
 



  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: Factor () 
  {
    static Timer factor_timer("SparseCholesky::Factor");

    static Timer timerb("SparseCholesky::Factor - B");
    static Timer timerc("SparseCholesky::Factor - C");

    RegionTimer reg (factor_timer);

    
    int n = Height();
    if (n > 2000)
      cout << IM(4) << " factor " << flush;

    // to avoid aliasing:
    int * hfirstinrow = firstinrow.Addr(0);
    int * hfirstinrow_ri = firstinrow_ri.Addr(0);
    int * hrowindex2 = rowindex2.Addr(0);
    TM * hlfact = lfact.Addr(0);
    
    enum { BS = 4 };

    Array<TM> sum(BS*maxrow);

    double flops1 = 0;
    double flops2 = 0;
    // starttime1 = clock();

    for (int i1 = 0; i1 < n;  )
      {
	int last_same = i1;
	while (last_same < n && blocknrs[last_same] == blocknrs[i1])
	  last_same++;

        timerb.Start();

	// same rows
	int mi = last_same - i1;
        int miBS = (mi / BS) * BS;

        for (int jj = 0; jj < miBS; jj+=4)
          {
            for (int j2 = jj; j2 < jj+4; j2++)
              if (n > 2000 && (i1+j2) % 1000 == 999)
                {
                  if ((i1+j2) % 10000 == 9999)
                    cout << IM(4) << "+" << flush;
                  else
                    cout << IM(4) << "." << flush;
                }
            
            int nk = hfirstinrow[i1+jj+1]-hfirstinrow[i1+jj];
            
            for (int k = 0; k < nk*BS; k++) 
              sum[k] = 0;

            for (int i2 = 0; i2 < jj; i2++)
              {
                int firsti = hfirstinrow[i1+i2] + jj-i2;
                TM * hli = hlfact+firsti;
                
                TM q1 = - diag[i1+i2] * hli[-1];
                TM qtrans1 = Trans (q1);  
                TM q2 = - diag[i1+i2] * hli[0];
                TM qtrans2 = Trans (q2);  
                TM q3 = - diag[i1+i2] * hli[1];
                TM qtrans3 = Trans (q3);  
                TM q4 = - diag[i1+i2] * hli[2];
                TM qtrans4 = Trans (q4);  
                
                sum[   0  ] += qtrans1 * hli[0];
                sum[  BS  ] += qtrans1 * hli[1];
                sum[  BS+1] += qtrans2 * hli[1];
                sum[2*BS  ] += qtrans1 * hli[2];
                sum[2*BS+1] += qtrans2 * hli[2];
                sum[2*BS+2] += qtrans3 * hli[2];
                
                for (int k = 3; k < nk; k++)
                  {
                    TM hv = hli[k];
                    sum[k*BS]    += qtrans1 * hv;
                    sum[k*BS+1]  += qtrans2 * hv;
                    sum[k*BS+2]  += qtrans3 * hv;
                    sum[k*BS+3]  += qtrans4 * hv;
                  }
                
                diag[i1+jj]   += Trans (hli[-1]) * q1;
                diag[i1+jj+1] += Trans (hli[0]) * q2;
                diag[i1+jj+2] += Trans (hli[1]) * q3;
                diag[i1+jj+3] += Trans (hli[2]) * q4;
              }
            
            for (int l = 0; l < BS; l++)
              {
                int firstj = hfirstinrow[i1+jj+l]-l;
                for (int k = l; k < nk; k++)
                  hlfact[firstj+k] += sum[k*BS+l];
              }
            
            
            for (int j2 = jj; j2 < jj+4; j2++)
              {
                int firstj = hfirstinrow[i1+j2];
                for (int i2 = jj; i2 < j2; i2++)
                  {
                    int firsti = hfirstinrow[i1+i2] + j2-i2;
		    
                    TM q = - diag[i1+i2] * hlfact[firsti-1];
                    TM qtrans = Trans (q);  
                    
                    for (int k = 0; k < nk-(j2-jj); k++)
                      hlfact[firstj+k] += qtrans * hlfact[firsti+k];
		    
                    diag[i1+j2] += Trans (hlfact[firsti-1]) * q;
                  }
                
                flops1 += (nk+1)*j2;
                TM aiinv;
                CalcInverse (diag[i1+j2], aiinv);
                diag[i1+j2] = aiinv;
              }
          }

	for (int jj = miBS ; jj < mi; jj++)
	  {
	    if (n > 2000 && (i1+jj) % 1000 == 999)
	      {
		if ((i1+jj) % 10000 == 9999)
		  cout << IM(4) << "+" << flush;
		else
		  cout << IM(4) << "." << flush;
	      }

	    int firstj = hfirstinrow[i1+jj];
	    int nk = hfirstinrow[i1+jj+1]-firstj;

	    for (int i2 = 0; i2 < jj; i2++)
	      {
		int firsti = hfirstinrow[i1+i2] + jj-i2;
		  
		TM q = - diag[i1+i2] * hlfact[firsti-1];
		TM qtrans = Trans (q);  
		  
	        for (int k = 0; k < nk; k++)
		  hlfact[firstj+k] += qtrans * hlfact[firsti+k];

		diag[i1+jj] += Trans (hlfact[firsti-1]) * q;
	      }

            TM aiinv;
	    CalcInverse (diag[i1+jj], aiinv);
	    diag[i1+jj] = aiinv;

	    flops1 += (nk+1)*jj;
	  }


        timerb.Stop();
        timerc.Start();


	// merge rows

	int firsti_ri = hfirstinrow_ri[i1] + last_same-i1-1;
	int firsti = hfirstinrow[i1] + last_same-i1-1;
	int lasti = hfirstinrow[i1+1]-1;
	mi = lasti-firsti+1;

	// loop unrolling for cache
	// for (j = 0; j < mi-BS+1; j+=BS)
        miBS = (mi / BS) * BS;


	// #pragma omp parallel
        {
          
	  // #pragma omp for
          for (int j = 0; j < miBS; j+=BS)
            {
              for (int k = BS*(j+1); k < BS*mi; k++)
                sum[k] = TSCAL_MAT(0.0);
              
              for (int i2 = i1; i2 < last_same; i2++)
                {
                  int first = hfirstinrow[i2] + last_same-i2-1;
                  TM qtrans  = Trans (diag[i2] * hlfact[first+j]);
                  TM qtrans2 = Trans (diag[i2] * hlfact[first+j+1]);
                  TM qtrans3 = Trans (diag[i2] * hlfact[first+j+2]);
                  TM qtrans4 = Trans (diag[i2] * hlfact[first+j+3]);


                  sum[BS*(j+1)] += qtrans * hlfact[first+j+1];
                  sum[BS*(j+2)] += qtrans * hlfact[first+j+2];
                  sum[BS*(j+3)] += qtrans * hlfact[first+j+3];

                  sum[BS*(j+2)+1] += qtrans2 * hlfact[first+j+2];
                  sum[BS*(j+3)+1] += qtrans2 * hlfact[first+j+3];

                  sum[BS*(j+3)+2] += qtrans3 * hlfact[first+j+3];

                  for (int k = j+BS; k < mi; k++)
                    {
                      TM hv = hlfact[first+k];
                      sum[BS*k]   += qtrans  * hv;
                      sum[BS*k+1] += qtrans2 * hv;
                      sum[BS*k+2] += qtrans3 * hv;
                      sum[BS*k+3] += qtrans4 * hv;
                    }
                }

              flops2 += (BS*(mi-j)-6)*(last_same-i1);

              // merge together
              for (int l = 0; l < BS; l++)
                {
                  int firstj = hfirstinrow[hrowindex2[firsti_ri+j+l]];
                  int firstj_ri = hfirstinrow_ri[hrowindex2[firsti_ri+j+l]];

                  for (int k = j+1+l; k < mi; k++)
                    {
                      int kk = hrowindex2[firsti_ri+k];

                      while (hrowindex2[firstj_ri] != kk)
                        {
                          firstj++;
                          firstj_ri++;
                        }
		    
                      lfact[firstj] -= sum[BS*k+l];
                      firstj++;
                      firstj_ri++;
                    }
                }
            }
        }
        
	// for (  ; j < mi; j++)
	for (int j = miBS; j < mi; j++)
	  {
	    for (int k = j+1; k < mi; k++)
	      sum[k] = TSCAL_MAT(0.0);

	    for (int i2 = i1; i2 < last_same; i2++)
	      {
		int first = hfirstinrow[i2] + last_same-i2-1;

		TM qtrans = Trans (diag[i2] * hlfact[first+j]);
		for (int k = j+1; k < mi; k++)
		  sum[k] += qtrans * hlfact[first+k];
	      }

	    flops2 += (mi - j)*(last_same-i1);

	    // merge together
	    int firstj = hfirstinrow[hrowindex2[firsti_ri+j]];
	    int firstj_ri = hfirstinrow_ri[hrowindex2[firsti_ri+j]];

	    for (int k = j+1; k < mi; k++)
	      {
		int kk = hrowindex2[firsti_ri+k];
		while (hrowindex2[firstj_ri] != kk)
		  {
		    firstj++;
		    firstj_ri++;
		  }
		    
		lfact[firstj] -= sum[k];
		firstj++;
		firstj_ri++;
	      }
	  }

	for (int i2 = i1; i2 < last_same; i2++)
	  {
	    int first = hfirstinrow[i2] + last_same-i2-1;
	    int last = hfirstinrow[i2+1];
	    int j_ri = hfirstinrow_ri[i2] + last_same-i2-1;

	    for (int j = first; j < last; j++, j_ri++)
	      {
		TM q = diag[i2] * lfact[j];
		diag[rowindex2[j_ri]] -= Trans (lfact[j]) * q;
	      }
	  }

	i1 = last_same;

        timerc.Stop();
      }

    for (int i = 0, j = 0; i < n; i++)
      {
	TM ai = diag[i];

	int last = hfirstinrow[i+1];
	for ( ; j < last; j++)
	  {
	    TM hm = ai * lfact[j];
	    lfact[j] = hm;
	  }	
      }

    if (n > 2000)
      cout << IM(4) << endl;
  }











#ifdef LAPACK


  void CalcAtA (FlatMatrix<> b1, FlatMatrix<> btb)
  {
    
    // btb = Trans(b1) * b1 | Lapack;

    /*
    if (b1.Width() < 200)
      {
        if (b1.Width() < 10 || b1.Height() < 10)
          btb = Trans(b1) * b1;
        else
          btb = Trans(b1) * b1 | Lapack;
      }

    else
      {
        auto r = Range(b1.Width());
        int parts = r.Size() / 256;
        if (parts == 0) parts = 1;

        // lower half
        for (int i = 0; i < parts; i++)
          for (int j = 0; j <= i; j++)
            {
              auto rowr = r.Split (i,parts);
              auto colr = r.Split (j,parts);
#pragma omp task
              {
                btb.Rows(rowr).Cols(colr) = Trans(b1.Cols(rowr)) * b1.Cols(colr) | Lapack;
                if (i != j)
                  btb.Rows(colr).Cols(rowr) = Trans(btb.Rows(rowr).Cols(colr));                  
              }
            }
#pragma omp taskwait

      }
    */

    int w = b1.Width();
    int nbl = w / 64 + 1;

    if (nbl <= 1)
      {
        btb = Trans(b1) * b1;
      }
    else
    task_manager -> CreateJob ( [&] (const TaskInfo & ti)
                                {
                                  int br = ti.task_nr / nbl;
                                  int bc = ti.task_nr % nbl;
                                  
                                  if (br > bc) return;

                                  auto rowr = Range(w).Split (br, nbl);
                                  auto colr = Range(w).Split (bc, nbl);
                                  
                                  btb.Rows(rowr).Cols(colr) = Trans(b1.Cols(rowr)) * b1.Cols(colr) | Lapack;
                                  
                                  if (br != bc)
                                    btb.Rows(colr).Cols(rowr) = Trans (btb.Rows(rowr).Cols(colr));

                                }, nbl*nbl);

  }




  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: FactorSPD () 
  {
    throw Exception ("FactorSPD called for non-<double,double> matrix");
  }


  template <>
  void SparseCholesky<double,double,double> :: FactorSPD () 
  {
    // task_manager -> StopWorkers();

    static Timer factor_timer("SparseCholesky::Factor SPD");

    static Timer timerb("SparseCholesky::Factor - B");
    static Timer timerc("SparseCholesky::Factor - C");
    static Timer timercla("SparseCholesky::Factor - C(lapack)");
    static Timer timerc1("SparseCholesky::Factor - C1");
    static Timer timerc2("SparseCholesky::Factor - C2");

    RegionTimer reg (factor_timer);
    
    int n = Height();
    if (n > 2000)
      cout << IM(4) << " factor SPD " << flush;

    // to avoid aliasing:
    int * hfirstinrow = firstinrow.Addr(0);
    int * hfirstinrow_ri = firstinrow_ri.Addr(0);
    int * hrowindex2 = rowindex2.Addr(0);
    // double * hlfact = lfact.Addr(0);
    
    // enum { BS = 4 };
    // Array<double> sum(BS*maxrow);


    for (int i1 = 0; i1 < n;  )
      {
	int last_same = i1;
	while (last_same < n && blocknrs[last_same] == blocknrs[i1])
	  last_same++;

        
        timerb.Start();

	// same rows
	int mi = last_same - i1;


        int nk = hfirstinrow[i1+1] - hfirstinrow[i1] + 1;
        Matrix<> a(mi, nk);
        a = 0.0;
	for (int j = 0; j < mi; j++)
	  {
            a(j,j) = diag[i1+j];
            a.Row(j).Range(j+1,nk) = FlatVector<>(nk-j-1, &lfact[hfirstinrow[i1+j]]);
          }


        /*
          // the original version
        for (int j = 0; j < mi; j++)
          {
            for (int k = 0; k < j; k++)
              {
                double q = -a(k,k)*a(k,j);
                for (int l = j; l < nk; l++)
                  a(j,l) += q * a(k,l);
              }
            a(j,j) = 1.0 / a(j,j);
          }
        */


        /*
          // first factor A, then calc B
        for (int j = 0; j < mi; j++)
          {
            for (int k = 0; k < j; k++)
              {
                double q = -a(k,k)*a(k,j);
                for (int l = j; l < mi; l++)
                  a(j,l) += q * a(k,l);
              }
            a(j,j) = 1.0 / a(j,j);
          }

        *testout << "ldl = " << endl << a.Cols(0,mi);


        for (int l = mi; l < nk; l++)
          for (int k = 0; k < mi; k++)
            {
              double q = -a(k,k)*a(k,l);

              for (int j = k+1; j < mi; j++)
                a(j,l) += q*a(k,j);
            }

        *testout << "b-block = " << a.Cols(mi, nk);
        */


        Matrix<> a1 = a.Cols(0, mi);
        Vector<> da1(mi);

        integer na1 = a1.Width();
        integer lda = a1.Width();
        integer info;
        char uplo = 'L';
        char ch_diag = 'N';
        
        dpotrf_ (&uplo, &na1, &a1(0,0), &lda, &info);    // A = L^t L

        for (int i = 0; i < mi; i++) 
          da1(i) = a1(i,i);

        a.Cols(0,mi) = a1;

	char trans = 'N';  // 'T' 'N'
	Matrix<> b1t = Trans(a.Cols(mi,nk));
	int nrhs = nk-mi;
	int ldb = mi;
        dtrtrs_ (&uplo, &trans, &ch_diag, &na1, &nrhs, &a1(0,0), &lda, &b1t(0,0), &ldb, &info);
	Matrix<> b1 = Trans(b1t);
	a.Cols(mi,nk) = b1; 
	// *testout << "b1 with lapack: " << Trans(b1t) << endl;
	
	/*
        dtrtri_ (&uplo, &ch_diag, &na1, &a1(0,0), &lda, &info);
        Matrix<> b1 = Trans(a1) * a.Cols(mi,nk) | Lapack;
        a.Cols(mi,nk) = b1;
	*/
	// *testout << "b1 with invert: " << b1 << endl;

        for (int i = 0; i < na1; i++)
          a.Row(i) *= da1(i);
        
        for (int i = 0; i < mi; i++)
          a(i,i) = 1.0/sqr(da1(i)); 

	for (int j = 0; j < mi; j++)
	  {
            diag[i1+j] = a(j,j);
            FlatVector<>(nk-j-1, &lfact[hfirstinrow[i1+j]]) = a.Row(j).Range(j+1,nk);
          }


        timerb.Stop();
        timerc.Start();


	// merge rows
	int firsti_ri = hfirstinrow_ri[i1] + last_same-i1-1;
	int firsti = hfirstinrow[i1] + last_same-i1-1;
	int lasti = hfirstinrow[i1+1]-1;
	mi = lasti-firsti+1;

        timercla.Start();

        // Matrix<> btb = Trans(b1)*b1 | Lapack;
        Matrix<> btb(b1.Width());
        CalcAtA (b1, btb);
        
        timercla.Stop();
        // *testout << "b^t b = " << endl << btb << endl;


        timerc1.Start();


        if (mi < 100)
          {
            for (int j = 0; j < mi; j++)
              {
                FlatVector<> sum = btb.Row(j);
            
                // merge together
                int firstj = hfirstinrow[hrowindex2[firsti_ri+j]];
                int firstj_ri = hfirstinrow_ri[hrowindex2[firsti_ri+j]];
                
                for (int k = j+1; k < mi; k++)
                  {
                    int kk = hrowindex2[firsti_ri+k];
                    while (hrowindex2[firstj_ri] != kk)
                      {
                        firstj++;
                        firstj_ri++;
                      }
                    
                    lfact[firstj] -= sum[k];
                    firstj++;
                    firstj_ri++;
                  }
              }
          }

        else

          // for (int j = 0; j < mi; j++)
          ParallelFor 
            (Range(mi), [&] (int j)
             {
              FlatVector<> sum = btb.Row(j);
              
              // merge together
              int firstj = hfirstinrow[hrowindex2[firsti_ri+j]];
              int firstj_ri = hfirstinrow_ri[hrowindex2[firsti_ri+j]];
              
              for (int k = j+1; k < mi; k++)
                {
                  int kk = hrowindex2[firsti_ri+k];
                  while (hrowindex2[firstj_ri] != kk)
                    {
                      firstj++;
                      firstj_ri++;
                    }
                  
                  lfact[firstj] -= sum[k];
                  firstj++;
                    firstj_ri++;
                }
             });

        timerc1.Stop();


        timerc2.Start();
	for (int i2 = i1; i2 < last_same; i2++)
	  {
	    int first = hfirstinrow[i2] + last_same-i2-1;
	    int last = hfirstinrow[i2+1];
	    int j_ri = hfirstinrow_ri[i2] + last_same-i2-1;

	    for (int j = first; j < last; j++, j_ri++)
	      {
		double q = diag[i2] * lfact[j];
		diag[rowindex2[j_ri]] -= Trans (lfact[j]) * q;
	      }
	  }
        timerc2.Stop();
	i1 = last_same;
        timerc.Stop();
      }

    for (int i = 0, j = 0; i < n; i++)
      {
	double ai = diag[i];
	int last = hfirstinrow[i+1];

	for ( ; j < last; j++)
          lfact[j] *= ai;
      }


    if (n > 2000)
      cout << IM(4) << endl;

    // task_manager -> StartWorkers();
  }

#endif











  





  


  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {
    y = 0.0;
    MultAdd (TSCAL_VEC(1.0), x, y);
    return;
  }
  
  
  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  MultAdd (TSCAL_VEC s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseCholesky::MultAdd");
    RegionTimer reg (timer);
    timer.AddFlops (2.0*lfact.Size());

    int n = Height();
    
    const FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();

    
    Vector<TVX> hy(n);
    for (int i = 0; i < n; i++)
      hy(order[i]) = fx(i);

    TVX hv;

    // TVX * hhy = hy.Addr(0);
    FlatVector<TVX> hhy = hy;



    const TM * hlfact = &lfact[0];
    const TM * hdiag = &diag[0];
    // const int * hrowindex = &rowindex[0];
    const int * hrowindex2 = &rowindex2[0];
    const int * hfirstinrow = &firstinrow[0];
    const int * hfirstinrow_ri = &firstinrow_ri[0];
    
    int j = 0;
    for (int i = 0; i < n; i++)
      {
	TVX val = hy(i);
	int last = hfirstinrow[i+1];
	int j_ri = hfirstinrow_ri[i];
	while (j < last)
	  {
	    // hhy[hrowindex[j]] -= Trans (hlfact[j]) * val;
	    hhy[hrowindex2[j_ri]] -= Trans (hlfact[j]) * val;
	    j++;
	    j_ri++;
	  }
      }
  
    for (int i = 0; i < n; i++)
      {
	hv = hdiag[i] * hhy[i];
	hhy[i] = hv;
      }

    for (int i = n-1; i >= 0; i--)
      {
	int minj = hfirstinrow[i];
	int maxj = hfirstinrow[i+1];
	int j_ri = hfirstinrow_ri[i];

	TVX sum;
	sum = 0.0;
	
	for (j = minj; j < maxj; j++, j_ri++)
	  sum += lfact[j] * hy(rowindex2[j_ri]);
	
	hy(i) -= sum;
      }





    /*
      const TM * hlfact = &lfact[0];
      const TM * hdiag = &diag[0];
      const int * hrowindex = &rowindex[0];
      const int * hfirstinrow = &firstinrow[0];
    
      int j = 0;
      for (int i = 0; i < n; i++)
      {
      TVX val = hy(i);
      int last = hfirstinrow[i+1];
      while (j < last)
      {
      hhy[hrowindex[j]] -= Trans (hlfact[j]) * val;
      j++;
      }
      }
  
      for (int i = 0; i < n; i++)
      {
      hv = hdiag[i] * hhy[i];
      hhy[i] = hv;
      }

      for (int i = n-1; i >= 0; i--)
      {
      const int minj = hfirstinrow[i];
      const int maxj = hfirstinrow[i+1];

      TVX sum;
      sum = 0.0;
	
      for (j = minj; j < maxj; j++)
      sum += lfact[j] * hy(rowindex[j]);
	
      hy(i) -= sum;
      }
    */





    if (inner)
      {
	for (int i = 0; i < n; i++)
	  if (inner->Test(i))
	    fy(i) += s * hy(order[i]);
      }
    else if (cluster)
      {
	for (int i = 0; i < n; i++)
	  if ((*cluster)[i])
	    fy(i) += s * hy(order[i]);
      }
    else
      {
	for (int i = 0; i < n; i++)
	  fy(i) += s * hy(order[i]);
      }

  }
  





  template <>
  void SparseCholesky<double, double, double> :: 
  SolveBlock (int bnr, FlatVector<> hy) const
  {
    auto range = BlockDofs (bnr);


    // triangular solve
    for (auto i : range)
      {
        int size = range.end()-i-1;
        FlatVector<> vlfact(size, &lfact[firstinrow[i]]);
        hy.Range(i+1, range.end()) -= hy(i) * vlfact;
      }

    auto extdofs = BlockExtDofs (bnr);

    VectorMem<100> temp(extdofs.Size());
    temp = 0;

    for (auto i : range)
      {
        int first = firstinrow[i] + range.end()-i-1;

        FlatVector<> ext_lfact (extdofs.Size(), &lfact[first]);
        temp += hy(i) * ext_lfact;
      }

    for (int j : Range(extdofs))
#pragma omp atomic
      hy(extdofs[j]) -= temp(j);
  }

  template <>
  void SparseCholesky<double, double, double> :: 
  SolveBlockT (int bnr, FlatVector<> hy) const
  {
    // const double * hlfact = &lfact[0];

    // const int * hrowindex2 = &rowindex2[0];
    const int * hfirstinrow = &firstinrow[0];
    const int * hfirstinrow_ri = &firstinrow_ri[0];


    for (int i = blocks[bnr+1]-1; i >= blocks[bnr]; i--)
      {
	int minj = hfirstinrow[i];
	int maxj = hfirstinrow[i+1];
	int j_ri = hfirstinrow_ri[i];

	double sum = 0;

	for (int j = minj; j < maxj; j++, j_ri++)
	  sum += lfact[j] * hy(rowindex2[j_ri]);
	
	hy(i) -= sum;
      }

  }




  template <typename T>
  class MyQueue
  {
    Array<T> data;
    Array<atomic<int>> ok;
    atomic<int> rcnt, wcnt;
  public:
    MyQueue (int size)
      : data(size), ok(size) , rcnt(0), wcnt(0)
    {
      for (auto & d : ok) 
        d.store (0, memory_order_relaxed);
    }

    void Push (T in)
    {
      int mypos = wcnt++;
      data[mypos] = in;
      ok[mypos] = 1;
      /*
#pragma omp critical(output)
      {
        cout << "push " << in << " to " << mypos << endl;
      }
      */
    }

    bool Pop (T & out)
    {
      while (1)
        {
          if (rcnt >= data.Size()) return false;
          
          int oldval = 1;
          if (ok[rcnt].compare_exchange_weak (oldval, 0))
            {
              int mypos = rcnt;
              rcnt++;
              out = data[mypos];
              /*
#pragma omp critical(output)
              {
                cout << "pop " << out << " from " << mypos << ",   wcnt = " << wcnt << endl;
              }
              */
              return true;
            }
        }
    }
  };



  
  static TQueue queue;


  template <typename TFUNC>
  void RunParallelDependency (const Table<int> & dag, TFUNC func)
  {
    static Timer t("parallelDependency");
    static Timer t2("parallelDependency - dag");

    RegionTimer reg(t);


    Array<atomic<int>> cnt_dep(dag.Size());

    for (auto & d : cnt_dep) 
      d.store (0, memory_order_relaxed);


    ParallelFor (Range(dag),
                 [&] (int i)
                 {
                   for (int j : dag[i])
                     cnt_dep[j]++;
                 });
    

    Array<int> ready(dag.Size());
    ready.SetSize0();
    int num_final = 0;

    for (int j : Range(cnt_dep))
      {
        if (cnt_dep[j] == 0) ready.Append(j);
        if (dag[j].Size() == 0) num_final++;
      }


    /*
    while (ready.Size())
      {
        int size = ready.Size();
        int nr = ready[size-1];
        ready.SetSize(size-1);

        func(nr);

        for (int j : dag[nr];
          {
            cnt_dep[j]--;
            if (cnt_dep[j] == 0)
              ready.Append(j);
          }
      }
    */



    atomic<int> cnt_final(0);
    SharedLoop sl(Range(ready));

    Array< Vec<3> > timings(task_manager -> GetNumThreads());
    double starttime = omp_get_wtime();

    t2.Start();
    task_manager -> CreateJob 
      ([&] (const TaskInfo & ti)
       {
         // timings[ti.task_nr](0) = omp_get_wtime()-starttime;

        TPToken ptoken(queue); 
        TCToken ctoken(queue); 
        
        for (int i : sl)
          queue.enqueue (ptoken, ready[i]);

        // timings[ti.task_nr](1) = omp_get_wtime()-starttime;

        while (1)
           {
             if (cnt_final >= num_final) break;

             int nr;
             if(!queue.try_dequeue_from_producer(ptoken, nr)) 
               if(!queue.try_dequeue(ctoken, nr))  
                 continue; 
             
             if (dag[nr].Size() == 0)
               cnt_final++;

             func(nr);

             for (int j : dag[nr])
               {
                 if (--cnt_dep[j] == 0)
                   queue.enqueue (ptoken, j);
               }
           }

        // timings[ti.task_nr](2) = omp_get_wtime()-starttime;
       });
    t2.Stop();

    /*
    ofstream out("depend.out");
    for (int i = 0; i < timings.Size(); i++)
      out << 1e6 * timings[i] << endl;
    */

    /*
    MyQueue<int> queue(dag.Size());

    for (int i : ready)
      queue.Push(i);
    
    task_manager -> CreateJob 
      ([&] (const TaskInfo & ti)
       {
         while (1)
           {
             int nr;
             if (!queue.Pop(nr)) break;
             
             func(nr);

             for (int j : dag[nr])
               {
                 if (--cnt_dep[j] == 0)
                   queue.Push (j);
               }
           }
        ptqend[ti.task_nr] = omp_get_wtime();         
       }, 22);
    */
  }

















  template <>
  void SparseCholesky<double, double, double> :: 
  MultAdd (TSCAL_VEC s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseCholesky<d,d,d>::MultAdd");
    static Timer timer0("SparseCholesky<d,d,d>::MultAdd queue only");
    static Timer timer1("SparseCholesky<d,d,d>::MultAdd fac1");
    static Timer timer2("SparseCholesky<d,d,d>::MultAdd fac2");
    RegionTimer reg (timer);
    timer.AddFlops (2.0*lfact.Size());

    int n = Height();
    
    const FlatVector<> fx = x.FV<double> ();
    FlatVector<> fy = y.FV<double> ();

    Vector<> hy(n);
    for (int i = 0; i < n; i++)
      hy(order[i]) = fx(i);


    const double * hdiag = &diag[0];


    /*

    class ProfileData
    {
    public:
      double tstart, tend;
      int size, extsize;
    };

    Array<ProfileData> prof(blocks.Size()-1);
    
    double tstart = omp_get_wtime();

    timer1.Start();

    RunParallelDependency (block_dependency, 
                           [&] (int nr) 
                           {
                             int s = blocks[nr+1]-blocks[nr];
                             if (s >= 100)
                               prof[nr].tstart = omp_get_wtime();

                             SolveBlock(nr, hy); 

                             if (s >= 100)
                               prof[nr].tend = omp_get_wtime();
                             prof[nr].size = blocks[nr+1]-blocks[nr];
                             int row = blocks[nr];
                             prof[nr].extsize = firstinrow[row+1]-firstinrow[row] - prof[nr].size+1;
                           });

    timer1.Stop();

    ofstream out ("cholesky.prof");
    for (int i = 0; i < prof.Size(); i++)
      if (prof[i].size >= 100)
        out << i << "  " << prof[i].size << ", extsize = " << prof[i].extsize << ",  ts = " 
            << 1e6*(prof[i].tstart-tstart) 
            << ", te = " << 1e6*(prof[i].tend-tstart) 
            << ", delta = " << 1e6*(prof[i].tend-prof[i].tstart) << endl;
    
    */



    /*
    RunParallelDependency (block_dependency, 
                           [&] (int nr) 
                           {
                             SolveBlock(nr, hy); 
                           });
    */



    class ProfileData
    {
    public:
      double tstart, tend;
      int size, extsize;
    };

    Array<ProfileData> prof(microtasks.Size());


    timer0.Start();

    RunParallelDependency (micro_dependency, 
                           [&] (int nr) 
                           {
                             ;
                           });
    timer0.Stop();


    timer1.Start();

    double tstart = omp_get_wtime();
    RunParallelDependency (micro_dependency, 
                           [&] (int nr) 
                           {
                             MicroTask task = microtasks[nr];
                             int blocknr = task.blocknr;
                             auto range = BlockDofs (blocknr);
                            
                             /*
                             if (range.Size() > 5)
                               prof[nr].tstart = omp_get_wtime();
                             prof[nr].size = range.Size();
                             */
 
                             if (task.solveL)
                               {
                                 // SolveBlock(blocknr, hy); 
                                 
                                 for (auto i : range)
                                   {
                                     int size = range.end()-i-1;
                                     FlatVector<> vlfact(size, &lfact[firstinrow[i]]);
                                     hy.Range(i+1, range.end()) -= hy(i) * vlfact;
                                   }

                               }
                             
                             else 

                               {
                                 auto all_extdofs = BlockExtDofs (blocknr);
                                 auto myr = Range(all_extdofs).Split (task.bblock, task.nbblocks);
                                 auto extdofs = all_extdofs.Range(myr);

                                 VectorMem<520> temp(extdofs.Size());
                                 temp = 0;
                                 
                                 for (auto i : range)
                                   {
                                     int first = firstinrow[i] + range.end()-i-1;
                                     
                                     FlatVector<> ext_lfact (all_extdofs.Size(), &lfact[first]);
                                     temp += hy(i) * ext_lfact.Range(myr);
                                   }
                                 
                                 for (int j : Range(extdofs))
#pragma omp atomic
                                   hy(extdofs[j]) -= temp(j);
                               }

                             /*
                             if (range.Size() > 5)
                               prof[nr].tend = omp_get_wtime();
                             */
                           });

    double tend = omp_get_wtime();
    timer1.Stop();

    /*
    ofstream out ("cholesky.prof");
    for (int i = 0; i < prof.Size(); i++)
      if (prof[i].size > 5)
        out << i << "  " << prof[i].size << ", SolveL = " << microtasks[i].solveL 
            << ",  ts = " << 1e6*(prof[i].tstart-tstart) 
            << ", te = " << 1e6*(prof[i].tend-tstart) 
            << ", delta = " << 1e6*(prof[i].tend-prof[i].tstart) << endl;

    for (int i = 0; i < 48; i++)
      {
        out << "process " << i << " qend at " << 1e6*(ptqend[i]-tstart) << endl;
      }

    out << "queue end: " << 1e6*(tqend-tstart) << endl;
    out << "total time: " << 1e6*(tend-tstart) << endl;
    */






    /*
    for (int i = 0; i < blocks.Size()-1; i++)
      SolveBlock (i, hy);
    */

    for (int i = 0; i < n; i++)
      hy[i] *= hdiag[i];

    timer2.Start();


    /*
    for (int i = blocks.Size()-1; i >= 0; i--)
      SolveBlockT (i, hy);
    */

    RunParallelDependency (micro_dependency_trans, 
                           [&] (int nr) 
                           {
                             MicroTask task = microtasks[nr];
                             int blocknr = task.blocknr;
                             auto range = BlockDofs (blocknr);
                            
 
                             if (task.solveL)
                               {
                                 /*
                                 for (auto i : range)
                                   {
                                     int size = range.end()-i-1;
                                     FlatVector<> vlfact(size, &lfact[firstinrow[i]]);
                                     hy.Range(i+1, range.end()) -= hy(i) * vlfact;
                                   }
                                 */
                                 for (int i = range.end()-1; i >= range.begin(); i--)
                                   {
                                     int size = range.end()-i-1;
                                     FlatVector<> vlfact(size, &lfact[firstinrow[i]]);
                                     hy(i) -= InnerProduct (vlfact, hy.Range(i+1, range.end()));
                                   }

                               }
                             
                             else 

                               {
                                 auto all_extdofs = BlockExtDofs (blocknr);
                                 auto myr = Range(all_extdofs).Split (task.bblock, task.nbblocks);
                                 auto extdofs = all_extdofs.Range(myr);
                                 
                                 VectorMem<520> temp(extdofs.Size());
                                 for (int j : Range(extdofs))
                                   temp(j) = hy(extdofs[j]);

                                 for (auto i : range)
                                   {
                                     int first = firstinrow[i] + range.end()-i-1;
                                     FlatVector<> ext_lfact (all_extdofs.Size(), &lfact[first]);
                                     
                                     double val = InnerProduct (ext_lfact.Range(myr), temp);
#pragma omp atomic
                                     hy(i) -= val;
                                   }
                                 
                                 
                               }

                           });


    timer2.Stop();




    if (inner)
      {
	for (int i = 0; i < n; i++)
	  if (inner->Test(i))
	    fy(i) += s * hy(order[i]);
      }
    else if (cluster)
      {
	for (int i = 0; i < n; i++)
	  if ((*cluster)[i])
	    fy(i) += s * hy(order[i]);
      }
    else
      {
	for (int i = 0; i < n; i++)
	  fy(i) += s * hy(order[i]);
      }

  }
  






  /*
  template <>
  void SparseCholesky<double, double, double> :: 
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseCholesky::MultAdd - double");
    static Timer timerL("SparseCholesky::MultAdd - L");
    static Timer timerLt("SparseCholesky::MultAdd - Lt");
    RegionTimer reg (timer);
    timer.AddFlops (2.0*lfact.Size());

    int n = Height();
    
    const FlatVector<> fx = x.FV<double> ();
    FlatVector<> fy = y.FV<double> ();
    
    Vector<> hy(n);
    for (int i = 0; i < n; i++)
      hy(order[i]) = fx(i);

    FlatVector<> hhy = hy;

    const double * hlfact = &lfact[0];
    const double * hdiag = &diag[0];
    // const int * hrowindex = &rowindex[0];
    const int * hrowindex2 = &rowindex2[0];
    const int * hfirstinrow = &firstinrow[0];
    const int * hfirstinrow_ri = &firstinrow_ri[0];

    timerL.Start();

    enum { BS = 32 };

    Vector<> tmp1(n);

    for (int i = 0; i < n; i++)
      {
        if ( (i+BS <= n) && (blocknrs[i] == blocknrs[i+BS-1]) )
          // if (false)
          {
            // solve with trig factor
            for (int i2 = 0; i2 < BS; i2++)
              {
                double val = hy(i+i2);
                int first = hfirstinrow[i+i2];
                int j_ri = hfirstinrow_ri[i+i2];

                for (int j = first; j < first+BS-i2-1; j++, j_ri++)
                  hhy[hrowindex2[j_ri]] -= Trans (hlfact[j]) * val;
              }

            
            int nk =  hfirstinrow[i+1] - (hfirstinrow[i]+BS-1);
            FlatVector<> tmp = tmp1.Range(0,nk);
            tmp = 0.0;

            // if (nk < 100)
            if (true)
              {
                for (int i2 = 0; i2 < BS; i2+=4)
                  {
                    double val0 = hy(i+i2);
                    double val1 = hy(i+i2+1);
                    double val2 = hy(i+i2+2);
                    double val3 = hy(i+i2+3);
                    
                    FlatVector<> v0(nk, (double*)&hlfact[hfirstinrow[i+i2]]+BS-i2-1);
                    FlatVector<> v1(nk, (double*)&hlfact[hfirstinrow[i+i2+1]]+BS-i2-2);
                    FlatVector<> v2(nk, (double*)&hlfact[hfirstinrow[i+i2+2]]+BS-i2-3);
                    FlatVector<> v3(nk, (double*)&hlfact[hfirstinrow[i+i2+3]]+BS-i2-4);
                    
                    tmp += val0 * v0 + val1*v1 + val2*v2 + val3*v3;
                  }
              }
            else
              {
#pragma omp parallel
                {
                  for (int i2 = 0; i2 < BS; i2+=4)
                    {
                      double val0 = hy(i+i2);
                      double val1 = hy(i+i2+1);
                      double val2 = hy(i+i2+2);
                      double val3 = hy(i+i2+3);
                      
                      FlatVector<> v0(nk, (double*)&hlfact[hfirstinrow[i+i2]]+BS-i2-1);
                      FlatVector<> v1(nk, (double*)&hlfact[hfirstinrow[i+i2+1]]+BS-i2-2);
                      FlatVector<> v2(nk, (double*)&hlfact[hfirstinrow[i+i2+2]]+BS-i2-3);
                      FlatVector<> v3(nk, (double*)&hlfact[hfirstinrow[i+i2+3]]+BS-i2-4);
                 
#pragma omp for
                      for (int j = 0; j < nk; j++)
                        tmp(j) += val0 * v0(j) + val1*v1(j) + val2*v2(j) + val3*v3(j);
                    }
                }
              }

            int j_ri = hfirstinrow_ri[i]+BS-1;
            for (int j = 0; j < nk; j++, j_ri++)
              hhy[hrowindex2[j_ri]] -= tmp(j);
            


            i += BS-1;
          }
        else
          {
            double val = hy(i);
            int first = hfirstinrow[i];
            int last = hfirstinrow[i+1];
            int j_ri = hfirstinrow_ri[i];
            for (int j = first; j < last; j++, j_ri++)
              hhy[hrowindex2[j_ri]] -= Trans (hlfact[j]) * val;
          }
      }
    timerL.Stop();  

    for (int i = 0; i < n; i++)
      hhy[i] *= hdiag[i];

    timerLt.Start();
    for (int i = n-1; i >= 0; i--)
      {
	int minj = hfirstinrow[i];
	int maxj = hfirstinrow[i+1];
	int j_ri = hfirstinrow_ri[i];

	TVX sum;
	sum = 0.0;
	
	for (int j = minj; j < maxj; j++, j_ri++)
	  sum += lfact[j] * hy(rowindex2[j_ri]);
	
	hy(i) -= sum;
      }

    timerLt.Stop();

    if (inner)
      {
	for (int i = 0; i < n; i++)
	  if (inner->Test(i))
	    fy(i) += s * hy(order[i]);
      }
    else if (cluster)
      {
	for (int i = 0; i < n; i++)
	  if ((*cluster)[i])
	    fy(i) += s * hy(order[i]);
      }
    else
      {
	for (int i = 0; i < n; i++)
	  fy(i) += s * hy(order[i]);
      }

  }
  */  






  SparseFactorization ::     
  SparseFactorization (const BaseSparseMatrix & amatrix,
		       const BitArray * ainner,
		       const Array<int> * acluster)
    : matrix(amatrix), inner(ainner), cluster(acluster)
  { 
    smooth_is_projection = true;
    if (cluster)
      {
	int first_cluster = 0;
	for (int i = 0; i < cluster->Size(); i++)
	  if ( (*cluster)[i] != 0)
	    {
	      first_cluster = (*cluster)[i];
	      break;
	    }
	
	for (int i = 0; i < cluster->Size(); i++)
	  if ( (*cluster)[i] != 0 && (*cluster)[i] != first_cluster)
	    {
	      smooth_is_projection = false;
	      break;
	    }
      }
  }
  
  
  void SparseFactorization  :: 
  Smooth (BaseVector & u, const BaseVector & /* f */, BaseVector & y) const
  {
    auto hvec1 = u.CreateVector();
    auto hvec2 = u.CreateVector();

    hvec1 = y;
    matrix.MultAdd1 (-1, u, hvec1, inner, cluster);

    hvec2 = (*this) * hvec1;
    u += hvec2;
    
    matrix.MultAdd2 (-1, hvec2, y, inner, cluster);
  }
  


















  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: Set (int i, int j, const TM & val)
  {
    // *testout << "sparse cholesky, set (" << i << ", " << j << ") = " << val << endl;
    if (i == j)
      {
	diag[i] = val;
	return;
      }
  
    TM hval;
  
    if (i > j) 
      {
	swap (i, j);
	hval = Trans (val);
      }
    else
      hval = val;
  
    //    first = 0;
    //    if (i > 0) first = lastinrow[i-1]+1;
    // last = lastinrow[i];
    int first = firstinrow[i];
    int first_ri = firstinrow_ri[i];
    int last = firstinrow[i+1];
  
    while (first < last)
      {
	if (rowindex2[first_ri] == j)
	  {
	    lfact[first] = hval;
	    return;
	  }
	first++;
	first_ri++;
      }
    cerr << "Position " << i << ", " << j << " not found" << endl;
  }


  template <class TM, class TV_ROW, class TV_COL>
  const TM & SparseCholesky<TM, TV_ROW, TV_COL> :: Get (int i, int j) const
  {
    if (i == j)
      {
	return diag[i];
      }

    if (i > j) 
      {
	swap (i, j);
	cerr << "SparseCholesky::Get: access to upper side not available" << endl;
      }

    int first = firstinrow[i];
    int first_ri = firstinrow_ri[i];
    int last = firstinrow[i+1];
  
    while (first < last)
      {
	if (rowindex2[first_ri] == j)
	  {
	    return lfact[first];
	  }
	first++;
	first_ri++;
      }
    cerr << "Position " << i << ", " << j << " not found" << endl;
    return *new TM;
  }

  template <class TM, class TV_ROW, class TV_COL>
  ostream & SparseCholesky<TM, TV_ROW, TV_COL> :: Print (ostream & ost) const
  {
    int n = Height();

    for (int i = 0; i < n; i++)
      {
	ost << i << ": " << order[i] << " diag = "
	    << diag[i] << endl;
      }
    ost << endl;
  
    int j = 1;
    for (int i = 1; i <= n; i++)
      {
	int j_ri = firstinrow_ri[i-1];
	ost << i << ": ";
	for ( ; j < firstinrow[i]; j++, j_ri++)
	  {
	    ost << rowindex2[j_ri] << "("
		<< lfact[j] 
		<< ")  ";
	  }
	ost << endl;
      }

    return ost;
  }


  template <class TM, class TV_ROW, class TV_COL>
  SparseCholesky<TM, TV_ROW, TV_COL> :: ~SparseCholesky()
  {
    delete mdo;
  }





  template class SparseCholesky<double>;
  template class SparseCholesky<Complex>;
  template class SparseCholesky<double, Complex, Complex>;

#if MAX_SYS_DIM >= 1
  template class SparseCholesky<Mat<1,1,double> >;
  template class SparseCholesky<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SparseCholesky<Mat<2,2,double> >;
  template class SparseCholesky<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class SparseCholesky<Mat<3,3,double> >;
  template class SparseCholesky<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class SparseCholesky<Mat<4,4,double> >;
  template class SparseCholesky<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SparseCholesky<Mat<5,5,double> >;
  template class SparseCholesky<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SparseCholesky<Mat<6,6,double> >;
  template class SparseCholesky<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SparseCholesky<Mat<7,7,double> >;
  template class SparseCholesky<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SparseCholesky<Mat<8,8,double> >;
  template class SparseCholesky<Mat<8,8,Complex> >;
#endif

#ifdef CACHEBLOCKSIZE
  template class SparseCholesky<double, Vec<CACHEBLOCKSIZE>, Vec<CACHEBLOCKSIZE> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseCholesky<double, Vec<2,double>, Vec<2,double> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseCholesky<double, Vec<3,double>, Vec<3,double> >;
  template class SparseCholesky<double, Vec<4,double>, Vec<4,double> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseCholesky<double, Vec<5,double>, Vec<5,double> >;
  template class SparseCholesky<double, Vec<6,double>, Vec<6,double> >;
  template class SparseCholesky<double, Vec<7,double>, Vec<7,double> >;
  template class SparseCholesky<double, Vec<8,double>, Vec<8,double> >;
  template class SparseCholesky<double, Vec<9,double>, Vec<9,double> >;
  template class SparseCholesky<double, Vec<10,double>, Vec<10,double> >;
  template class SparseCholesky<double, Vec<11,double>, Vec<11,double> >;
  template class SparseCholesky<double, Vec<12,double>, Vec<12,double> >;
  template class SparseCholesky<double, Vec<13,double>, Vec<13,double> >;
  template class SparseCholesky<double, Vec<14,double>, Vec<14,double> >;
  template class SparseCholesky<double, Vec<15,double>, Vec<15,double> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseCholesky<double, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseCholesky<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class SparseCholesky<double, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseCholesky<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class SparseCholesky<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class SparseCholesky<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class SparseCholesky<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class SparseCholesky<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class SparseCholesky<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class SparseCholesky<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class SparseCholesky<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class SparseCholesky<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class SparseCholesky<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class SparseCholesky<double, Vec<15,Complex>, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class SparseCholesky<Complex, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class SparseCholesky<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class SparseCholesky<Complex, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class SparseCholesky<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class SparseCholesky<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class SparseCholesky<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class SparseCholesky<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class SparseCholesky<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class SparseCholesky<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class SparseCholesky<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class SparseCholesky<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class SparseCholesky<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class SparseCholesky<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class SparseCholesky<Complex, Vec<15,Complex>, Vec<15,Complex> >;
#endif



}













