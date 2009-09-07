/* *************************************************************************/
/* File:   sparseldl.cc                                                    */
/* Author: Joachim Schoeberl                                               */
/* Date:   18. Nov. 2001                                                   */
/* *************************************************************************/

#include <la.hpp>
//#ifndef MACOS
//#include <emmintrin.h>
//#endif


namespace ngla
{
  using namespace ngla;
  using namespace ngstd;

  template <class TM, class TV_ROW, class TV_COL>
  SparseCholesky<TM, TV_ROW, TV_COL> :: 
  SparseCholesky (const SparseMatrix<TM, TV_ROW, TV_COL> & a, 
		  BitArray * ainner,
		  Array<int> * acluster,
		  bool allow_refactor)
    : mat(a)
  { 
    // (*testout) << "matrix = " << a << endl;
    // (*testout) << "diag a = ";
    // for ( int i=0; i<a.Height(); i++ ) (*testout) << i << ", " << a(i,i) << endl;

    int n = a.Height();
    height = n;

    inner = ainner;
    cluster = acluster;
    
    int printstat = 0;
    
    if (printstat)
      cout << "Minimal degree ordering: N = " << n << endl;
    
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
	  FlatArray<const int> row = a.GetRowIndices(i);
	  for (int j = 0; j < row.Size(); j++)
	    {
	      int col = row[j];
	      if (col <= i)
		if ( ((*cluster)[i] == (*cluster)[col] && (*cluster)[i]) ||
		     i == col )
		  mdo->AddEdge (i, col);
	    }
	}

    if (printstat)
      cout << "start ordering" << endl;
    
    mdo->Order();
    
    endtime = clock();
    if (printstat)
      cout << "ordering time = "
	   << double (endtime - starttime) / CLOCKS_PER_SEC 
	   << " secs" << endl;
    
    starttime = endtime;
    
    if (printstat)
      cout << "," << flush;
    Allocate (mdo->order,  mdo->vertices, &mdo->blocknr[0]);

    delete mdo;
    mdo = 0;

    diag.Alloc (n);
    diag.SetName ("sparse inverse, diag");

    lfact.Alloc (nze);
    lfact.SetName ("sparse inverse, lfact");

    for (int i = 0; i < nze; i++)
      lfact[i] = 0.0;

    
    endtime = clock();
    if (printstat)
      (cout) << "allocation time = "
	     << double (endtime - starttime) / CLOCKS_PER_SEC << " secs" << endl;
    
    starttime = endtime;


    if (!inner && !cluster)
      for (int i = 0; i < n; i++)
	for (int j = 0; j < a.GetRowIndices(i).Size(); j++)
	  {
	    int col = a.GetRowIndices(i)[j];
	    if (col <= i)
	      SetOrig (i, col, a.GetRow(i)[j]);
	  }
    
    else if (inner)
      for (int i = 0; i < n; i++)
	for (int j = 0; j < a.GetRowIndices(i).Size(); j++)
	  {
	    int col = a.GetRowIndices(i)[j];
	    if (col <= i)
	      if ( (inner->Test(i) && inner->Test(col)) || i==col)
		SetOrig (i, col, a.GetRow(i)[j]);
	  }
    else
      for (int i = 0; i < n; i++)
	{
	  FlatArray<const int> row = a.GetRowIndices(i);
	  for (int j = 0; j < row.Size(); j++)
	    {
	      int col = row[j];
	      if (col <= i)
		if ( ( (*cluster)[i] == (*cluster)[col] && (*cluster)[i]) ||
		     i == col )
                  SetOrig (i, col, a.GetRow(i)[j]);
	    }
	}
    
    if (printstat)
      cout << "do factor " << flush;
    
    Factor(); 

    if (printstat)
      cout << "done" << endl;
    
    endtime = clock();

    if (printstat)
      (cout) << " factoring time = " << double(endtime - starttime) / CLOCKS_PER_SEC << " sec" << endl;
  }
  

  
  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  Allocate (const Array<int> & aorder, 
	    // const Array<CliqueEl*> & cliques,
	    const Array<MDOVertex> & vertices,
	    const int * blocknr)
  {
    int n = aorder.Size();
    
    order.Alloc (n);
    order.SetName ("sparse inverse, reorder table");

    blocknrs.Alloc (n);
    blocknrs.SetName ("sparse inverse, block nrs");
    
    DynamicMem<int> helpi(n);
    
    // order: now inverse map 
    for (int i = 0; i < n; i++)
      order[aorder[i]] = i;

    for (int i = 0; i < n; i++)
      blocknrs[i] = blocknr[i];

    int cnt = 0;
    int cnt_master = 0;

    for (int i = 0; i < n; i++)
      {
	cnt += vertices[aorder[blocknr[i]]].nconnected - (i-blocknr[i]);
	if (blocknrs[i] == i)
	  cnt_master += vertices[aorder[i]].nconnected;
      }

    nze = cnt;

    if (n > 2000)
      cout << " " << cnt*sizeof(TM)+cnt_master*sizeof(int) << " Bytes " << flush;
    //cout <<"(cnt="<<cnt<<", sizeof(TM)="<<sizeof(TM)<< ", cnt_master=" << cnt_master << ", sizeof(int)=" << sizeof(int) <<") " << flush;


    /* 
     *testout << " Sparse Cholesky mem needed " << double(cnt*sizeof(TM)+cnt_master*sizeof(int))*1e-6 << " MBytes " << endl; 
     */  

    firstinrow.Alloc(n+1);
    firstinrow.SetName ("sparse inverse, table 1a");

    firstinrow_ri.Alloc(n+1);
    firstinrow_ri.SetName ("sparse inverse, table 1a");

    rowindex2.Alloc (cnt_master);
    rowindex2.SetName ("sparse inverse, tabl 2a");
    
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
	    MergeSort (ncon, &rowindex2[firstinrow_ri[i]], &helpi[0]);

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
  }
  


  template <class TM>
  void SetIdentity( TM &identity )
  {
    for (int i=0; i<identity.Width(); i++ ) identity(i,i) = 1;
  }
  void SetIdentity( double &identity ) { identity = 1; }
  void SetIdentity( Complex &identity ) { identity = 1; }


  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  FactorNew (const SparseMatrix<TM, TV_ROW, TV_COL> & a)
  {
    if ( height != a.Height() )
      {
	cout << "SparseCholesky::FactorNew called with matrix of different size." << endl;
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
	      if ( col <= i ) SetOrig (i, col, a.GetRow(i)[j]);
	    }
	  else if (i == col)
	    SetOrig (i, i, id);
	}
    
    Factor(); 
  }




  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: Factor () 
  {
    static int factor_timer = NgProfiler::CreateTimer ("SparseCholesky::Factor");
    NgProfiler::RegionTimer reg (factor_timer);

    clock_t starttime, starttime1;
    double time1 = 0, time2 = 0;

    
    int n = Height();
    if (n > 2000)
      cout << " factor " << flush;

    // to avoid aliasing:
    int * hfirstinrow = firstinrow;
    int * hfirstinrow_ri = firstinrow_ri;
    int * hrowindex2 = rowindex2;
    TM * hlfact = &lfact[0];
    
    enum { BS = 4 };

    Array<TM> sum(BS*maxrow);

    double flops1 = 0;
    double flops2 = 0;
    starttime1 = clock();

    for (int i1 = 0; i1 < n;  )
      {
	int last_same = i1;
	while (last_same < n && blocknrs[last_same] == blocknrs[i1])
	  last_same++;


	starttime = clock();
	
	// same rows
	int mi = last_same - i1;
        int miBS = (mi / BS) * BS;

        // miBS = 0;

        // #pragma omp parallel  (was commented out)
        {
          // Array<TM> sum(BS*maxrow);
          
          for (int jj = 0; jj < miBS; jj+=4)
            {
              for (int j2 = jj; j2 < jj+4; j2++)
                if (n > 2000 && (i1+j2) % 1000 == 999)
                  {
                    if ((i1+j2) % 10000 == 9999)
                      cout << "+" << flush;
                    else
                      cout << "." << flush;
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

		  // #pragma omp parallel for
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
        }



	for (int jj = miBS ; jj < mi; jj++)
	  {
	    if (n > 2000 && (i1+jj) % 1000 == 999)
	      {
		if ((i1+jj) % 10000 == 9999)
		  cout << "+" << flush;
		else
		  cout << "." << flush;
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




	time1 += double(clock() - starttime) / CLOCKS_PER_SEC;


	starttime = clock();
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
          Array<TM> sum(BS*maxrow);
          
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

	time2 += double(clock() - starttime) / CLOCKS_PER_SEC;
	i1 = last_same;
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
      cout << endl;

    /*
    time = double(clock() - starttime1) / CLOCKS_PER_SEC;
    cout << "factorization, time = " << time << endl;

    cout << "flops1 = " << flops1 << endl;
    cout << "flops2 = " << flops2 << endl;
    cout << "time1 = " << time1 << endl;
    cout << "time2 = " << time2 << endl;
    cout << "MFLOP1 = " << flops1 / time1 * 1e-6 << endl;
    cout << "MFLOP2 = " << flops2 / time2 * 1e-6 << endl;
    cout << "MFLOP = " << (flops1+flops2) / time * 1e-6 << endl;
    */
  }




  /*
    template <>
    void SparseCholesky<double> :: Factor () 
    {
    static int factor_timer = NgProfiler::CreateTimer ("SparseCholesky::Factor");
    NgProfiler::RegionTimer reg (factor_timer);

    clock_t starttime, starttime1;
    double time, time1 = 0, time2 = 0, time3 = 0;


    double ai, aiinv, q, q1, q2, q3, q4, q5, q6, hm;
    
    int n = Height();
    if (n > 2000)
    cout << " factor " << flush;

    // to avoid aliasing:
    int * hfirstinrow = firstinrow;
    int * hfirstinrow_ri = firstinrow_ri;
    int * hrowindex2 = rowindex2;
    double * hlfact = &lfact[0];
    
    enum { BS1 = 6 };
    enum { BS = 6 };
    // Array<double> sum(BS1*maxrow)
    ;
    __m128d * dsum = (__m128d*)_mm_malloc (BS1*maxrow*sizeof(double), 16);
    FlatArray<double> sum(BS1*maxrow, (double*)dsum);

    double flops1 = 0;
    double flops2 = 0;
    starttime1 = clock();

    for (int i1 = 0; i1 < n;  )
    {
    int last_same = i1;
    while (last_same < n && blocknrs[last_same] == blocknrs[i1])
    last_same++;

    starttime = clock();
	
    // same rows
    int mi = last_same - i1;

    int jj = 0;
    for (jj = 0; jj < mi-BS1+1; jj+=BS1)
    {
    for (int j2 = jj; j2 < jj+BS1; j2++)
    if (n > 2000 && (i1+j2) % 1000 == 999)
    {
    if ((i1+j2) % 10000 == 9999)
    cout << "+" << flush;
    else
    cout << "." << flush;
    }


    int nk = hfirstinrow[i1+jj+1]-hfirstinrow[i1+jj];

    __m128d zero = _mm_set1_pd(0.0);
    for (int k = 0; k < nk*3; k++)
    dsum[k] = zero;

    for (int i2 = 0; i2 < jj; i2++)
    {
    int firsti = hfirstinrow[i1+i2] + jj-i2;
    double * hli = hlfact+firsti;
		
    __m128d dval = _mm_set1_pd (-diag[i1+i2]);
    __m128d hli12 = _mm_loadu_pd (&hli[-1]);
    __m128d hli34 = _mm_loadu_pd (&hli[1]);
    __m128d hli56 = _mm_loadu_pd (&hli[3]);

    __m128d q12 = _mm_mul_pd (hli12, dval);
    __m128d q34 = _mm_mul_pd (hli34, dval);
    __m128d q56 = _mm_mul_pd (hli56, dval);

    _mm_storeu_pd (&diag[i1+jj]  ,
    _mm_add_pd (_mm_mul_pd (hli12, q12), _mm_loadu_pd (&diag[i1+jj])));
    _mm_storeu_pd (&diag[i1+jj+2],
    _mm_add_pd (_mm_mul_pd (hli34, q34), _mm_loadu_pd (&diag[i1+jj+2])));
    _mm_storeu_pd (&diag[i1+jj+4],
    _mm_add_pd (_mm_mul_pd (hli56, q56), _mm_loadu_pd (&diag[i1+jj+4])));
		

    __m128d hd;
    dsum[0] = _mm_add_sd (_mm_mul_sd (_mm_set_sd (hli[0]), q12),  dsum[0]);
    dsum[3] = _mm_add_pd (_mm_mul_pd (_mm_set1_pd (hli[1]), q12), dsum[3]);
    dsum[6] = _mm_add_pd (_mm_mul_pd (_mm_set1_pd (hli[2]), q12), dsum[6]);
    dsum[7] = _mm_add_sd (_mm_mul_sd (_mm_set_sd (hli[2]), q34),  dsum[7]);

    dsum[9] = _mm_add_pd (_mm_mul_pd (_mm_set1_pd (hli[3]), q12), dsum[9]);
    dsum[10] = _mm_add_pd (_mm_mul_pd (_mm_set1_pd (hli[3]), q34), dsum[10]);

    dsum[12] = _mm_add_pd (_mm_mul_pd (_mm_set1_pd (hli[4]), q12), dsum[12]);
    dsum[13] = _mm_add_pd (_mm_mul_pd (_mm_set1_pd (hli[4]), q34), dsum[13]);
    dsum[14] = _mm_add_sd (_mm_mul_sd (_mm_set_sd (hli[4]), q56),  dsum[14]);

    for (int k = 5; k < nk; k++)
    {
    __m128d hv = _mm_set1_pd (hli[k]);
    dsum[3*k]   = _mm_add_pd (dsum[3*k],   _mm_mul_pd (q12, hv));
    dsum[3*k+1] = _mm_add_pd (dsum[3*k+1], _mm_mul_pd (q34, hv));
    dsum[3*k+2] = _mm_add_pd (dsum[3*k+2], _mm_mul_pd (q56, hv));
    }
    }

    for (int l = 0; l < BS1; l++)
    {
    int firstj = hfirstinrow[i1+jj+l]-l;
    for (int k = l; k < nk; k++)
    hlfact[firstj+k] += sum[k*BS1+l];
    }


    for (int j2 = jj; j2 < jj+BS1; j2++)
    {
    int firstj = hfirstinrow[i1+j2];
    for (int i2 = jj; i2 < j2; i2++)
    {
    int firsti = hfirstinrow[i1+i2] + j2-i2;
		    
    q = - diag[i1+i2] * hlfact[firsti-1];
    diag[i1+j2] += hlfact[firsti-1] * q;
		  
    for (int k = 0; k < nk-(j2-jj); k++)
    hlfact[firstj+k] += q * hlfact[firsti+k];
    }

    flops1 += (nk+1)*j2;
    CalcInverse (diag[i1+j2], aiinv);
    diag[i1+j2] = aiinv;
    }
    }


    for ( ; jj < mi; jj++)
    {
    if (n > 2000 && (i1+jj) % 1000 == 999)
    {
    if ((i1+jj) % 10000 == 9999)
    cout << "+" << flush;
    else
    cout << "." << flush;
    }


    int firstj = hfirstinrow[i1+jj];
    int nk = hfirstinrow[i1+jj+1]-firstj;

    for (int i2 = 0; i2 < jj; i2++)
    {
    int firsti = hfirstinrow[i1+i2] + jj-i2;
		  
    q = - diag[i1+i2] * hlfact[firsti-1];
    q = Trans (q);  
		  
    for (int k = 0; k < nk; k++)
    hlfact[firstj+k] += q * hlfact[firsti+k];

    diag[i1+jj] += Trans (hlfact[firsti-1]) * q;
    }

    flops1 += (nk+1)*jj;

    CalcInverse (diag[i1+jj], aiinv);
    diag[i1+jj] = aiinv;
    }



    time1 += double(clock() - starttime) / CLOCKS_PER_SEC;
    starttime = clock();
    // merge rows

    int firsti_ri = hfirstinrow_ri[i1] + last_same-i1-1;
    int firsti = hfirstinrow[i1] + last_same-i1-1;
    int lasti = hfirstinrow[i1+1]-1;
    mi = lasti-firsti+1;

    int j;
    // loop unrolling for cache
    for (j = 0; j < mi-BS+1; j+=BS)
    {
    for (int k = BS*(j+1); k < BS*mi; k++)
    sum[k] = TSCAL_MAT(0.0);

    for (int i2 = i1; i2 < last_same; i2++)
    {
    int first = hfirstinrow[i2] + last_same-i2-1;

    __m128d hv = _mm_set1_pd (diag[i2]);
    __m128d q12 = _mm_mul_pd (hv, _mm_loadu_pd (&hlfact[first+j]));
    __m128d q34 = _mm_mul_pd (hv, _mm_loadu_pd (&hlfact[first+j+2]));
    __m128d q56 = _mm_mul_pd (hv, _mm_loadu_pd (&hlfact[first+j+4]));

    dsum[3*(j+1)]   = _mm_add_sd (_mm_mul_sd (_mm_set_sd (hlfact[first+j+1]), q12),  
    dsum[3*(j+1)]);
    dsum[3*(j+3)+1] = _mm_add_sd (_mm_mul_sd (_mm_set_sd (hlfact[first+j+3]), q34),  
    dsum[3*(j+3)+1]);
    dsum[3*(j+5)+2] = _mm_add_sd (_mm_mul_sd (_mm_set_sd (hlfact[first+j+5]), q56),  
    dsum[3*(j+5)+2]);

    for (int k = j+2; k < j+BS; k++)
    dsum[3*k]   = _mm_add_pd (dsum[3*k],   
    _mm_mul_pd (q12,_mm_set1_pd (hlfact[first+k])));
    for (int k = j+4; k < j+BS; k++)
    dsum[3*k+1] = _mm_add_pd (dsum[3*k+1], 
    _mm_mul_pd (q34,_mm_set1_pd (hlfact[first+k])));
    // for (int k = j+6; k < j+BS; k++)
    //   dsum[3*k+1] = _mm_add_pd (dsum[3*k+1], _mm_mul_pd (q34,_mm_set1_pd (hlfact[first+k])));

    for (int k = j+BS; k < mi; k++)
    {
    __m128d hv = _mm_set1_pd (hlfact[first+k]);
    dsum[3*k]   = _mm_add_pd (dsum[3*k],   _mm_mul_pd (q12, hv));
    dsum[3*k+1] = _mm_add_pd (dsum[3*k+1], _mm_mul_pd (q34, hv));
    dsum[3*k+2] = _mm_add_pd (dsum[3*k+2], _mm_mul_pd (q56, hv));
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

    time2 += double(clock() - starttime) / CLOCKS_PER_SEC;
    starttime = clock();

	
    for ( ; j < mi; j++)
    {
    for (int k = j+1; k < mi; k++)
    sum[k] = TSCAL_MAT(0.0);

    for (int i2 = i1; i2 < last_same; i2++)
    {
    int first = hfirstinrow[i2] + last_same-i2-1;

    q = Trans (diag[i2] * hlfact[first+j]);
    for (int k = j+1; k < mi; k++)
    sum[k] += q * hlfact[first+k];
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
    q = diag[i2] * lfact[j];
    diag[rowindex2[j_ri]] -= Trans (lfact[j]) * q;
    }
    }

    time3 += double(clock() - starttime) / CLOCKS_PER_SEC;
    i1 = last_same;
    }

    for (int i = 0, j = 0; i < n; i++)
    {
    ai = diag[i];

    int last = hfirstinrow[i+1];
    for ( ; j < last; j++)
    {
    hm = ai * lfact[j];
    lfact[j] = hm;
    }	
    }


    _mm_free (dsum);

    if (n > 2000)
    cout << endl;

    time = double(clock() - starttime1) / CLOCKS_PER_SEC;
    cout << "factorization, time = " << time << endl;

    cout << "flops1 = " << flops1 << endl;
    cout << "flops2 = " << flops2 << endl;
    cout << "time1 = " << time1 << endl;
    cout << "time2 = " << time2 << endl;
    cout << "time3 = " << time3 << endl;
    cout << "MFLOP1 = " << flops1 / time1 * 1e-6 << endl;
    cout << "MFLOP2 = " << flops2 / time2 * 1e-6 << endl;
    cout << "MFLOP = " << (flops1+flops2) / time * 1e-6 << endl;
    }
  */








  





  


  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {
    static int timer = NgProfiler::CreateTimer ("SparseCholesky::Mult");
    NgProfiler::RegionTimer reg (timer);

    int n = Height();
    
    const FlatVector<TVX> fx = 
      dynamic_cast<const T_BaseVector<TVX> &> (x).FV();
    FlatVector<TVX> fy = 
      dynamic_cast<T_BaseVector<TVX> &> (y).FV();
    
    Vector<TVX> hy(n);
    for (int i = 0; i < n; i++)
      hy(order[i]) = fx(i);

    TVX hv;
    TVX * hhy = &hy(0);

    const TM * hlfact = &lfact[0];
    const TM * hdiag = &diag[0];

    const int * hrowindex2 = &rowindex2[0];
    const int * hfirstinrow = &firstinrow[0];
    const int * hfirstinrow_ri = &firstinrow_ri[0];


    for (int i = 0; i < n; i++)
      {
	TVX val = hy(i);
	int first = hfirstinrow[i];
	int j_ri = hfirstinrow_ri[i];
	int last = hfirstinrow[i+1];

	for (int j = first; j < last; j++, j_ri++)
	  hhy[hrowindex2[j_ri]] -= Trans (hlfact[j]) * val;
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
	
	for (int j = minj; j < maxj; j++, j_ri++)
	  sum += lfact[j] * hy(hrowindex2[j_ri]);
	
	hy(i) -= sum;
      }

    for (int i = 0; i < n; i++)
      fy(i) = hy(order[i]);

    if (inner)
      {
	for (int i = 0; i < n; i++)
	  if (!inner->Test(i))
	    fy(i) = 0;
      }

    if (cluster)
      {
	for (int i = 0; i < n; i++)
	  if (!(*cluster)[i])
	    fy(i) = 0;
      }
  }
  
  
  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  MultAdd (TSCAL_VEC s, const BaseVector & x, BaseVector & y) const
  {
    int n = Height();
    

    const FlatVector<TVX> fx = 
      dynamic_cast<const T_BaseVector<TVX> &> (x).FV();
    FlatVector<TVX> fy = 
      dynamic_cast<T_BaseVector<TVX> &> (y).FV();
    
    Vector<TVX> hy(n);
    for (int i = 0; i < n; i++)
      hy(order[i]) = fx(i);

    TVX hv;

    TVX * hhy = &hy(0);




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
  






  
  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  Smooth (BaseVector & u, const BaseVector & f, BaseVector & y) const
  {
    int n = Height();
    
    // const FlatVector<TVX> ff = 
    // dynamic_cast<const T_BaseVector<TVX> &> (f).FV();
    FlatVector<TVX> fu = 
      dynamic_cast<T_BaseVector<TVX> &> (u).FV();
    FlatVector<TVX> fy = 
      dynamic_cast<T_BaseVector<TVX> &> (y).FV();
    

    const SparseMatrixSymmetric<TM, TV_ROW> & smat = 
      dynamic_cast<const SparseMatrixSymmetric<TM,TV_ROW> &> (mat);
      

    Vector<TVX> hvec(n);
    typedef typename mat_traits<TVX>::TSCAL TSCAL;
    hvec = TSCAL(0);
    if (inner)
      {
	for (int i = 0; i < n; i++)
	  if (inner->Test(i))
	    hvec(order[i]) = fy(i) - smat.RowTimesVectorNoDiag (i, fu);
      }
    else if (cluster)
      {
	for (int i = 0; i < n; i++)
	  {
	    if ((*cluster)[i])
	      hvec(order[i]) = fy(i) - smat.RowTimesVectorNoDiag (i, fu);
	  }
      }
    else
      for (int i = 0; i < n; i++)
	hvec(order[i]) = fy(i) - smat.RowTimesVectorNoDiag (i, fu);
      

    TVX hv;
    TVX * hhvec = &hvec(0);

    const TM * hlfact = &lfact[0];
    const TM * hdiag = &diag[0];
    const int * hrowindex2 = &rowindex2[0];
    const int * hfirstinrow = &firstinrow[0];
    const int * hfirstinrow_ri = &firstinrow_ri[0];
    
    int j = 0;
    for (int i = 0; i < n; i++)
      {
	TVX val = hvec(i);
	int last = hfirstinrow[i+1];
	int j_ri = hfirstinrow_ri[i];
	while (j < last)
	  {
	    hhvec[hrowindex2[j_ri]] -= Trans (hlfact[j]) * val;
	    j++;
	    j_ri++;
	  }
      }
  
    for (int i = 0; i < n; i++)
      {
	hv = hdiag[i] * hhvec[i];
	hhvec[i] = hv;
      }

    for (int i = n-1; i >= 0; i--)
      {
	int minj = hfirstinrow[i];
	int maxj = hfirstinrow[i+1];
	int j_ri = hfirstinrow_ri[i];

	TVX sum;
	sum = 0.0;
	
	for (j = minj; j < maxj; j++, j_ri++)
	  sum += lfact[j] * hvec(rowindex2[j_ri]);
	
	hvec(i) -= sum;
      }




    if (inner)
      {
	for (int i = 0; i < n; i++)
	  if (inner->Test(i))
	    fu(i) += hvec(order[i]);
      }
    else if (cluster)
      {
	for (int i = 0; i < n; i++)
	  if ((*cluster)[i])
	    {
	      fu(i) += hvec(order[i]);
	      smat.AddRowTransToVector (i, -hvec(order[i]), fy);
	    }
      }
    else
      {
	for (int i = 0; i < n; i++)
	  fu(i) += hvec(order[i]);
      }
  }
  































  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: Set (int i, int j, const TM & val)
  {
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













