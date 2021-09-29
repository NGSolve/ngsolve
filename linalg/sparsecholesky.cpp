// #define DEBUG

/* *************************************************************************/
/* File:   sparseldl.cc                                                    */
/* Author: Joachim Schoeberl                                               */
/* Date:   18. Nov. 2001                                                   */
/* *************************************************************************/

#include <core/register_archive.hpp>
#include <la.hpp>

#include <core/concurrentqueue.h>
#include <core/taskmanager.hpp>


typedef moodycamel::ConcurrentQueue<int> TQueue; 
typedef moodycamel::ProducerToken TPToken; 
typedef moodycamel::ConsumerToken TCToken; 


namespace ngla
{
  
  static TQueue queue;


  template <typename TFUNC>
  void RunParallelDependency (const Table<int> & dag,
                              const Table<int> & trans_dag, // transposed dag
                              TFUNC func)
  {
    Array<atomic<int>> cnt_dep(dag.Size());
    for (auto i : Range(cnt_dep))
      cnt_dep[i].store (trans_dag[i].Size(), memory_order_relaxed);
    
    Array<int> ready(dag.Size());
    ready.SetSize0();
    int num_final = 0;

    for (int j : Range(cnt_dep))
      {
        if (cnt_dep[j] == 0) ready.Append(j);
        if (dag[j].Size() == 0) num_final++;
      }


    if (!task_manager)
      {
        while (ready.Size())
          {
            int size = ready.Size();
            int nr = ready[size-1];
            ready.SetSize(size-1);
            
            func(nr);
            
            for (int j : dag[nr])
              {
                cnt_dep[j]--;
                if (cnt_dep[j] == 0)
                  ready.Append(j);
              }
          }
        return;
      }

    atomic<int> cnt_final(0);
    SharedLoop sl(Range(ready));

    task_manager -> CreateJob 
      ([&] (const TaskInfo & ti)
       {
        TPToken ptoken(queue); 
        TCToken ctoken(queue); 
        
        for (int i : sl)
          queue.enqueue (ptoken, ready[i]);

        while (1)
           {
             if (cnt_final >= num_final) break;

             while (TaskManager::ProcessTask()); // do the nested tasks
             
             int nr;
             if(!queue.try_dequeue_from_producer(ptoken, nr)) 
               if(!queue.try_dequeue(ctoken, nr))  
                 continue; 
             
             func(nr);

             if (dag[nr].Size() == 0)
               cnt_final++;

             for (int j : dag[nr])
               {
                 if (--cnt_dep[j] == 0)
                   queue.enqueue (ptoken, j);
               }
           }
       });
  }






  template <class TM>
  void SetIdentity( TM &identity )
  {
    for (int i=0; i<identity.Width(); i++ ) identity(i,i) = 1;
  }
  void SetIdentity( double &identity ) { identity = 1; }
  void SetIdentity( Complex &identity ) { identity = 1; }



  void SparseFactorization::DoArchive(Archive& ar)
  {
    ar & inner & smooth_is_projection;
    if(ar.Output())
      ar << const_pointer_cast<Array<int>>(cluster);
    else
      {
        shared_ptr<Array<int>> cl;
        ar & cl;
        cluster = cl;
      }
  }

  template <class TM>
  SparseCholeskyTM<TM> :: 
  SparseCholeskyTM (shared_ptr<const SparseMatrixTM<TM>> a,
                    shared_ptr<BitArray> ainner,
                    shared_ptr<const Array<int>> acluster,
                    bool allow_refactor)
    : SparseFactorization (a, ainner, acluster)
  { 
    static Timer t("SparseCholesky - total");
    static Timer ta("SparseCholesky - allocate");
    RegionTimer reg(t);
    GetMemoryTracer().SetName("SparseCholesky");
    GetMemoryTracer().Track(order, "order",
                            inv_order, "inv_order",
                            lfact, "lfact",
                            firstinrow, "firstinrow",
                            diag, "diag",
                            rowindex2, "rowindex2",
                            firstinrow_ri, "firstinrow_ri",
                            blocknrs, "blocknrs",
                            blocks, "blocks",
                            microtasks, "microtasks",
                            block_dependency, "block_dependency",
                            micro_dependency, "micro_dependency",
                            micro_dependency_trans, "mirco_dependency_trans");

    // (*testout) << "matrix = " << a << endl;
    // (*testout) << "diag a = ";
    // for ( int i=0; i<a->Height(); i++ ) (*testout) << i << ", " << a(i,i) << endl;

    int n = a->Height();
    height = n;

    int printstat = 0;
    
    if (printstat)
      cout << IM(4) << "Minimal degree ordering: N = " << n << endl;
    
    clock_t starttime, endtime;
    starttime = clock();
    
    mdo = new MinimumDegreeOrdering (n);
    GetMemoryTracer().Track(*mdo, "MinimumDegreeOrdering");

    if (inner)
      ParallelFor (n, [&] (size_t i)
                   {
                     if (!inner->Test(i))
                       mdo->SetUnusedVertex(i);
                   });
    if (cluster)
      for (int i = 0; i < n; i++)
        if (!(*cluster)[i])
          mdo->SetUnusedVertex(i);
    

    
    if (!inner && !cluster)
      for (int i = 0; i < n; i++)
	for (int j = 0; j < a->GetRowIndices(i).Size(); j++)
	  {
	    int col = a->GetRowIndices(i)[j];
	    if (col <= i)
	      mdo->AddEdge (i, col);
	  }

    else if (inner)
      {
        for (int i = 0; i < n; i++)
          if (inner->Test(i))
            for (auto col : a->GetRowIndices(i))
              if (col <= i)
                if (inner->Test(col)) //  || i==col)
                  mdo->AddEdge (i, col);
            /*
            for (int j = 0; j < a->GetRowIndices(i).Size(); j++)
              {
                int col = a->GetRowIndices(i)[j];
                if (col <= i)
                if (inner->Test(col)) //  || i==col)
                mdo->AddEdge (i, col);
                }
            */
      }

    else 
      for (int i = 0; i < n; i++)
	{
	  FlatArray<int> row = a->GetRowIndices(i);
	  for (int j = 0; j < row.Size(); j++)
	    {
	      int col = row[j];
	      if (col <= i)
		if ( ( ((*cluster)[i] == (*cluster)[col]) && (*cluster)[i]) )
                  // || i == col )
		  mdo->AddEdge (i, col);
	    }
	}
    
    /*
    for (int i = 0; i < n; i++)
      if (a->GetPositionTest (i,i) == numeric_limits<size_t>::max())
	{
	  mdo->AddEdge (i, i);
	  *testout << "add unsused position " << i << endl;
	}
    */

    if (printstat)
      cout << IM(4) << "start ordering" << endl;
    
    // mdo -> PrintCliques ();
    mdo->Order();
    nused = mdo->nused;
    endtime = clock();
    if (printstat)
      cout << IM(4) << "ordering time = "
	   << double (endtime - starttime) / CLOCKS_PER_SEC 
	   << " secs" << endl;
    
    starttime = endtime;
    
    if (printstat)
      cout << IM(4) << "," << flush;
    ta.Start();
    Allocate (mdo->order,  mdo->vertices, mdo->blocknr.Data());
    ta.Stop();

    delete mdo;
    mdo = 0;

    diag.SetSize(nused);
    // lfact.SetSize (nze);
    lfact = NumaInterleavedArray<TM> (nze);

    // lfact = TM(0.0);     // first touch
    ParallelForRange (nze, [&] (IntRange r)
                      {
                        lfact.Range(r) = TM(0.0);
                      });
    
    endtime = clock();
    if (printstat)
      (cout) << "allocation time = "
	     << double (endtime - starttime) / CLOCKS_PER_SEC << " secs" << endl;
    
    starttime = endtime;
    FactorNew(*a);
    /*
#ifdef LAPACK
    if (a.IsSPD())
      FactorSPD();
    else
#endif
      Factor(); 
    */

    /*
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
    */

    if (printstat)
      cout << IM(4) << "done" << endl;
    
    endtime = clock();

    if (printstat)
      (cout) << " factoring time = " << double(endtime - starttime) / CLOCKS_PER_SEC << " sec" << endl;
  }
  

  
  template <class TM>
  void SparseCholeskyTM<TM> :: 
  Allocate (const Array<int> & aorder, 
	    // const Array<CliqueEl*> & cliques,
	    const Array<MDOVertex> & vertices,
	    const int * in_blocknr)
  {
    int n = aorder.Size();

    order.SetSize (n);
    blocknrs.SetSize (nused);
    
    // order: now inverse map
    ParallelForRange (order.Size(), [&] (IntRange r)
                      {
                        order.Range(r) = -1;
                      });
    for (int i = 0; i < nused; i++)
      order[aorder[i]] = i;

    inv_order.SetSize(nused);
    inv_order = aorder;
    
    for (int i = 0; i < nused; i++)
      blocknrs[i] = in_blocknr[i];

    long int cnt = 0;
    long int cnt_master = 0;

    for (int i = 0; i < nused; i++)
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

    firstinrow.SetSize(nused+1);
    firstinrow_ri.SetSize(nused+1);
    rowindex2.SetSize (cnt_master);


    cnt = 0;
    cnt_master = 0;
    maxrow = 0;

    for (int i = 0; i < nused; i++)
      {
	firstinrow[i] = cnt;
	int ii = aorder[i];
	int ncon = vertices[ii].nconnected;

	if (blocknrs[i] == i)
	  {
	    firstinrow_ri[i] = cnt_master;

	    for (int j = 0; j < ncon; j++)
	      rowindex2[firstinrow_ri[i]+j] = order[vertices[ii].connected[j]];

	    // QuickSort (FlatArray<int> (ncon, &rowindex2[firstinrow_ri[i]]));
            QuickSort (rowindex2.Part(firstinrow_ri[i], ncon));   // fixes pathologic case with ncon = 0

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

    firstinrow[nused] = cnt;
    firstinrow_ri[nused] = cnt_master;
    
    
    for (int i = 1; i < blocknrs.Size(); i++)
      if (blocknrs[i] < blocknrs[i-1])
        throw Exception ("blocknrs are unordered !!");
    /*
    for (int i = 1; i < blocknrs.Size(); i++)
      {
        if (blocknrs[i] < blocknrs[i-1]) blocknrs[i] = blocknrs[i-1];
        // limit blocksize (to 256) for better granularity in solve-phase
        if (blocknrs[i] <= i-256) blocknrs[i] = i;
      }

    // cout << "finding block-dependeny ... " << endl;
    for (int i = 0; i < nused; i++)
      if(blocknrs[i] == i) blocks.Append(i);
    blocks.Append(nused);
    */
    
    blocks.Append(0);
    for (int i = 1; i < nused; i++)
      if (blocknrs[i] == i)  //  || i >= blocks.Last()+256) // don't subdivide, for this we have the micro-blocks
        blocks.Append (i);
    if (nused > 0)
      blocks.Append(nused);


    
    // find block dependency
    Array<int> block_of_dof(nused);
    for (int i = 0; i < blocks.Size()-1; i++)
      block_of_dof[Range(blocks[i], blocks[i+1])] = i;

    DynamicTable<int> dep(blocks.Size()-1);
    for (int i = 0; i < nused; i++)
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

    // genare micro-tasks:
    Array<int> first_microtask;
    for (int i = 0; i < blocks.Size()-1; i++)
      {
        // auto extdofs = BlockExtDofs (i);
        first_microtask.Append (microtasks.Size());

        // int nb = extdofs.Size() / 256 + 1;
        // int nb = (extdofs.Size()+255) / 256;
	int nb = 0;
	if(BlockDofs(i).Size()) {
	  auto extdofs = BlockExtDofs (i);
	  nb = (extdofs.Size()+255) / 256;
	}

        if (nb == 1)
          // if (false)
          {
            MicroTask mt;
            mt.blocknr = i;
            mt.type = MicroTask::LB_BLOCK;
            mt.bblock = 0;
            mt.nbblocks = 1;
            microtasks.Append (mt);
          }
        else
          {
            MicroTask mt;
            mt.blocknr = i;
            mt.type = MicroTask::L_BLOCK;
            mt.bblock = 0;
            mt.nbblocks = 0;
            microtasks.Append (mt);
            
            for (int j = 0; j < nb; j++)
              {
                MicroTask mt;
                mt.blocknr = i;
                mt.type = MicroTask::B_BLOCK;
                mt.bblock = j;
                mt.nbblocks = nb;
                microtasks.Append (mt);
              }
          }
      }
    first_microtask.Append (microtasks.Size());

    {
      TableCreator<int> creator(microtasks.Size());
      TableCreator<int> creator_trans(microtasks.Size());
      
      for ( ; !creator.Done(); creator++, creator_trans++)
        {
          for (int i = 0; i < first_microtask.Size()-1; i++)
            {
              if (first_microtask[i+1] == first_microtask[i]+1)
                { // just one LB block
                  int b = first_microtask[i];
                  for (int o : block_dependency[i])
                    {
                      creator.Add (b, first_microtask[o]);
                      creator_trans.Add (first_microtask[o], b);
                    }
                }
              else
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

      micro_dependency = creator.MoveTable();
      micro_dependency_trans = creator_trans.MoveTable();
    }
  }
  
  template<typename TM>
  void SparseCholeskyTM<TM>::DoArchive(Archive& ar)
  {
    SparseFactorization::DoArchive(ar);
    ar & height & nused & nze & order & inv_order & lfact
      & firstinrow & diag & rowindex2 & firstinrow_ri &
      blocknrs & blocks & block_dependency & microtasks
      & micro_dependency & micro_dependency_trans & mdo
      & maxrow;
  }

  template <class TM>
  void SparseCholeskyTM<TM> :: 
  FactorNew (const SparseMatrix<TM> & a)
  {
    static Timer tf("SparseCholesky - fill factor");
    tf.Start();
    if ( height != a.Height() )
      {
	cout << IM(4) << "SparseCholesky::FactorNew called with matrix of different size." << endl;
	return;
      }
    lfact = TM(0.0);

    if (!inner && !cluster)
      ParallelFor 
        (Range(height), [&](auto i)
         {
           auto rowind = a.GetRowIndices(i);
           auto rowvals = a.GetRowValues(i);
           
           for (auto j : Range(rowind.Size()))
             if (rowind[j] <= i)
               SetOrig (i, rowind[j], rowvals[j]);
         });
        else if (inner)
      ParallelFor 
        (Range(height), [&](auto i)
         {
           if (inner->Test(i))
             for (auto j : Range(a.GetRowIndices(i)))
               {
                 auto col = a.GetRowIndices(i)[j];
                 if (col <= i)
                   {
                     if (inner->Test(col))
                       SetOrig (i, col, a.GetRowValues(i)[j]);
                   }
               }
         }, TasksPerThread(5));
    else
      for (auto i : Range(height))
	{
	  auto row = a.GetRowIndices(i);
	  for (auto j : Range(row.Size()))
	    {
	      auto col = row[j];
	      if (col <= i)
		if (((*cluster)[i] == (*cluster)[col] && (*cluster)[i]))
                  SetOrig (i, col, a.GetRowValues(i)[j]);
	    }
	}
    tf.Stop();
    FactorSPD(); 
  }
 



  template <class TM>
  void SparseCholeskyTM<TM> :: Factor () 
  {
    static Timer factor_timer("SparseCholesky::Factor");

    static Timer timerb("SparseCholesky::Factor - B", NoTracing, NoTiming);
    static Timer timerc("SparseCholesky::Factor - C", NoTracing, NoTiming);

    RegionTimer reg (factor_timer);

    
    int n = nused; // Height();
    if (n > 2000){
      cout << IM(4) << " factor " << flush;
    }

    // to avoid aliasing:
    size_t * hfirstinrow = firstinrow.Addr(0);
    size_t * hfirstinrow_ri = firstinrow_ri.Addr(0);
    int * hrowindex2 = rowindex2.Addr(0);
    TM * hlfact = lfact.Addr(0);
    
    // enum { BS = 4 };
    constexpr int BS=4;

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
            for (int j2 = 0; j2 < 4; j2++)
              if (n > 2000 && (i1+jj+j2) % 1000 == 999)
                {
                  if ((i1+jj+j2) % 10000 == 9999)
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
            
            
            for (int hj2 = 0; hj2 < 4; hj2++)
              {
                int j2 = jj+hj2;
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

    if (n > 2000){
      cout << IM(4) << endl;
    }
  }








  
  /*
  template <class TM>
  void SparseCholeskyTM<TM> :: FactorSPD () 
  {
    throw Exception ("FactorSPD called for non-double matrix");
  }
  */

  // template <>
  template <class TM>
  void SparseCholeskyTM<TM> :: FactorSPD ()
  {
    Factor();
  }

  template <>
  void SparseCholeskyTM<double> :: FactorSPD ()
  {
    FactorSPD1(5.3);
  }

  template <>
  void SparseCholeskyTM<Complex> :: FactorSPD ()
  {
    FactorSPD1(5.2);
  }
  
  template <class TM> template<typename T>
  void SparseCholeskyTM<TM> :: FactorSPD1 (T dummy) 
  {
    if (!task_manager)
      {
        RunWithTaskManager ([&] ()
                            {
                              FactorSPD1(dummy);
                            });
        return;
      }

    static Timer factor_timer("SparseCholesky::Factor SPD");
    static Timer factor_dense1("SparseCholesky::Factor SPD - setup dense cholesky");
    static Timer factor_dense("SparseCholesky::Factor SPD - dense cholesky");

    // static Timer timerb("SparseCholesky::Factor - B", NoTracing);
    // static Timer timerc("SparseCholesky::Factor - merge in rows");
    // static Timer timercla("SparseCholesky::Factor - merge(lapack)", NoTracing);
    // static Timer timerc1("SparseCholesky::Factor - merge1", NoTracing);
    // static Timer timerc2("SparseCholesky::Factor - merge2", NoTracing);

    RegionTimer reg (factor_timer);
    
    size_t n = nused; // Height();
    if (n > 2000){
      cout << IM(4) << " factor SPD " << flush;
    }

    // to avoid aliasing:
    size_t * hfirstinrow = firstinrow.Addr(0);
    size_t * hfirstinrow_ri = firstinrow_ri.Addr(0);
    int * hrowindex2 = rowindex2.Addr(0);
    TM * hlfact = lfact.Addr(0);
   
    // #define CHOLESKY_ORIGINAL
    // #define CHOLESKY_SIMPLE
    // #define CHOLESKY_PARALLEL
#define CHOLESKY_PARALLEL_ATOMIC

    
#ifdef CHOLESKY_ORIGINAL
     int percent = 0;
     Array<TM> tmpmem;
    for (size_t blocknr : Range(blocks.Size()-1))
      {
        IntRange block = BlockDofs(blocknr);
        
        size_t i1 = block.First(); // blocks[blocknr];
        size_t last_same = block.Next(); // blocks[blocknr+1];
        
        if(i1*100/n > percent)
          {
            percent = i1*100/n;
            Ng_SetThreadPercentage(percent);
          }

        if (n > 2000)
          for (auto j : block)
            {
              if (j % 1000 == 999)
                {
                  if ((j) % 10000 == 9999)
                    cout << IM(4) << "+" << flush;
                  else
                    cout << IM(4) << "." << flush;
                }
            }
        
        // timerb.Start();

	// rows in same block ...
	size_t mi = block.Size();  // last_same - i1;
        size_t nk = hfirstinrow[i1+1] - hfirstinrow[i1] + 1;
        
        // factor_dense1.Start();
        // Matrix<TM,ColMajor> tmp(nk, nk);
        tmpmem.SetSize(nk*nk);
        FlatMatrix<TM,ColMajor> tmp(nk, nk, tmpmem.Addr(0));

        bool big = nk > 1000;
        if (big)
          {
            ParallelForRange (nk, [&](IntRange r)
                              {
                                tmp.Cols(r) = TM(0.0);
                              });
          }
        else
          {
            tmp = TM(0.0);
          }

	for (size_t j = 0; j < mi; j++)
	  {
            tmp(j,j) = diag[i1+j];
            tmp.Col(j).Range(j+1,nk) = FlatVector<TM>(nk-j-1, hlfact+hfirstinrow[i1+j]);
          }

        // factor_dense1.Stop();
        
        auto A11 = tmp.Rows(0,mi).Cols(0,mi);
        auto B   = tmp.Rows(mi,nk).Cols(0,mi);
        auto A22 = tmp.Rows(mi,nk).Cols(mi,nk);
        // factor_dense.Start();
        CalcLDL (A11);
        if (mi < nk)
          {
            CalcLDL_SolveL (A11,B);
            CalcLDL_A2 (A11.Diag(),B,A22);
          }
        // factor_dense.Stop();          

        
        auto write_back_row = [&](size_t j)
          {
            diag[i1+j] = A11(j,j);
            FlatVector<TM>(nk-j-1, hlfact+hfirstinrow[i1+j]) = tmp.Col(j).Range(j+1,nk);
          };

        if (mi < 10)
          for (size_t j = 0; j < mi; j++)
            write_back_row(j);
        else
          ParallelFor (mi, write_back_row);

        // timerb.Stop();
        //timerc.Start();


	// merge rows
	size_t firsti_ri = hfirstinrow_ri[i1] + last_same-i1-1;
	size_t firsti = hfirstinrow[i1] + last_same-i1-1;
	size_t lasti = hfirstinrow[i1+1]-1;
	mi = lasti-firsti+1;

        // timerc1.Start();

        auto merge_row = [&] (size_t j)
          {
            auto sum = A22.Col(j);
            
            // merge together
            size_t firstj = hfirstinrow[hrowindex2[firsti_ri+j]];
            size_t firstj_ri = hfirstinrow_ri[hrowindex2[firsti_ri+j]];
            
            for (size_t k = j+1; k < mi; k++)
              {
                size_t kk = hrowindex2[firsti_ri+k];
                while (hrowindex2[firstj_ri] != kk)
                  {
                    firstj++;
                    firstj_ri++;
                  }
                
                lfact[firstj] += sum[k];
                firstj++;
                firstj_ri++;
              }
          };

        if (mi < 100)
          for (size_t j = 0; j < mi; j++)
            merge_row(j);
        else
          ParallelFor (Range(mi), merge_row);
          

        // timerc1.Stop();

        // timerc2.Start();
	for (size_t i2 = i1; i2 < last_same; i2++)
	  {
	    size_t first = hfirstinrow[i2] + last_same-i2-1;
	    size_t last = hfirstinrow[i2+1];
	    size_t j_ri = hfirstinrow_ri[i2] + last_same-i2-1;

	    for (auto j = first; j < last; j++, j_ri++)
	      {
		TM q = diag[i2] * lfact[j];
		diag[rowindex2[j_ri]] -= Trans (lfact[j]) * q;
	      }
	  }
        // timerc2.Stop();
	// i1 = last_same;
        // timerc.Stop();
      }
#endif


#ifdef CHOLESKY_SIMPLE
    for (size_t blocknr : Range(blocks.Size()-1))
      {
        IntRange block = BlockDofs(blocknr);
        
        size_t i1 = block.First(); 
        size_t last_same = block.Next();

	// rows in same block ...
	size_t mi = block.Size();  // last_same - i1;
        size_t nk = hfirstinrow[i1+1] - hfirstinrow[i1] + 1;

        ArrayMem<TM,10000> tmpmem(nk*nk);
        FlatMatrix<TM,ColMajor> tmp(nk, nk, tmpmem.Addr(0));

        tmp = TM(0.0);

	for (size_t j = 0; j < mi; j++)
	  {
            tmp(j,j) = diag[i1+j];
            tmp.Col(j).Range(j+1,nk) = FlatVector<TM>(nk-j-1, hlfact+hfirstinrow[i1+j]);
          }

        auto A11 = tmp.Rows(0,mi).Cols(0,mi);
        auto B   = tmp.Rows(mi,nk).Cols(0,mi);
        auto A22 = tmp.Rows(mi,nk).Cols(mi,nk);

        CalcLDL (A11);
        if (mi < nk)
          {
            CalcLDL_SolveL (A11,B);
            CalcLDL_A2 (A11.Diag(),B,A22);
          }
        
        for (size_t j = 0; j < mi; j++)
          {
            diag[i1+j] = A11(j,j);
            FlatVector<TM>(nk-j-1, hlfact+hfirstinrow[i1+j]) = tmp.Col(j).Range(j+1,nk);
          };

	// merge rows
	size_t firsti_ri = hfirstinrow_ri[i1] + last_same-i1-1;
	size_t firsti = hfirstinrow[i1] + last_same-i1-1;
	size_t lasti = hfirstinrow[i1+1]-1;
	mi = lasti-firsti+1;

        for (size_t j = 0; j < mi; j++)
          {
            auto sum = A22.Col(j);
            
            // merge together
            size_t firstj = hfirstinrow[hrowindex2[firsti_ri+j]];
            size_t firstj_ri = hfirstinrow_ri[hrowindex2[firsti_ri+j]];
            
            for (size_t k = j+1; k < mi; k++)
              {
                size_t kk = hrowindex2[firsti_ri+k];
                while (hrowindex2[firstj_ri] != kk)
                  {
                    firstj++;
                    firstj_ri++;
                  }
                
                lfact[firstj] += sum[k];
                firstj++;
                firstj_ri++;
              }
          };
          
	for (size_t i2 = i1; i2 < last_same; i2++)
	  {
	    size_t first = hfirstinrow[i2] + last_same-i2-1;
	    size_t last = hfirstinrow[i2+1];
	    size_t j_ri = hfirstinrow_ri[i2] + last_same-i2-1;

	    for (auto j = first; j < last; j++, j_ri++)
	      {
		TM q = diag[i2] * lfact[j];
		diag[rowindex2[j_ri]] -= Trans (lfact[j]) * q;
	      }
	  }
      }
#endif





#ifdef CHOLESKY_PARALLEL
    
    // find new dependency graph ...
    
    // first, find the transposed graph
    
    TableCreator<int> creator_trans(block_dependency.Size());
    for ( ; !creator_trans.Done(); creator_trans++)
      for (int i : Range(block_dependency))
        for (int j : block_dependency[i])
          creator_trans.Add(j, i);
    auto block_dep_trans = creator_trans.MoveTable();

    // sort to avoid cycles
    for (auto entry : block_dep_trans)
      QuickSort (entry);
    
    TableCreator<int> creator_transitive(block_dependency.Size());
    for ( ; !creator_transitive.Done(); creator_transitive++)
      for (int i : Range(block_dep_trans))
        if (block_dep_trans[i].Size())
          {
            for (int j = 0; j < block_dep_trans[i].Size()-1; j++)
              creator_transitive.Add(block_dep_trans[i][j], block_dep_trans[i][j+1]);
            creator_transitive.Add(block_dep_trans[i].Last(), i);
          }
    auto dep_transitive = creator_transitive.MoveTable();
    
    TableCreator<int> creator_trans2(block_dependency.Size());
    for ( ; !creator_trans2.Done(); creator_trans2++)
      for (int i : Range(dep_transitive))
        for (int j : dep_transitive[i])
          creator_trans2.Add(j, i);
    auto dep_transitive_trans = creator_trans2.MoveTable();

    static Timer tdep("paralleldep");
    static Timer tdep1("paralleldep1");
    static Timer tdep2("paralleldep2");
    
    RunParallelDependency
      (dep_transitive, dep_transitive_trans, [&] (int blocknr)
       {
         // for (size_t blocknr : Range(blocks.Size()-1))
        IntRange block = BlockDofs(blocknr);
        RegionTracer reg(TaskManager::GetThreadId(), tdep, block.Size());

        size_t i1 = block.First(); 
        size_t last_same = block.Next();

	// rows in same block ...
	size_t mi = block.Size();  // last_same - i1;
        size_t nk = hfirstinrow[i1+1] - hfirstinrow[i1] + 1;

        ArrayMem<TM,1000> tmpmem(nk*nk);
        FlatMatrix<TM,ColMajor> tmp(nk, nk, tmpmem.Addr(0));

        tmp = TM(0.0);

	for (size_t j = 0; j < mi; j++)
	  {
            tmp(j,j) = diag[i1+j];
            tmp.Col(j).Range(j+1,nk) = FlatVector<TM>(nk-j-1, hlfact+hfirstinrow[i1+j]);
          }

        auto A11 = tmp.Rows(0,mi).Cols(0,mi);
        auto B   = tmp.Rows(mi,nk).Cols(0,mi);
        auto A22 = tmp.Rows(mi,nk).Cols(mi,nk);

        {
        RegionTracer reg1(TaskManager::GetThreadId(), tdep1, block.Size());
        CalcLDL (A11);
        if (mi < nk)
          {
            CalcLDL_SolveL (A11,B);
            CalcLDL_A2 (A11.Diag(),B,A22);
          }
        }
        
        for (size_t j = 0; j < mi; j++)
          {
            diag[i1+j] = A11(j,j);
            FlatVector<TM>(nk-j-1, hlfact+hfirstinrow[i1+j]) = tmp.Col(j).Range(j+1,nk);
          };

	// merge rows
	size_t firsti_ri = hfirstinrow_ri[i1] + last_same-i1-1;
	size_t firsti = hfirstinrow[i1] + last_same-i1-1;
	size_t lasti = hfirstinrow[i1+1]-1;
	mi = lasti-firsti+1;

        {
        RegionTracer reg2(TaskManager::GetThreadId(), tdep2, block.Size());
        for (size_t j = 0; j < mi; j++)
          {
            auto sum = A22.Col(j);
            
            // merge together
            size_t firstj = hfirstinrow[hrowindex2[firsti_ri+j]];
            size_t firstj_ri = hfirstinrow_ri[hrowindex2[firsti_ri+j]];
            
            for (size_t k = j+1; k < mi; k++)
              {
                size_t kk = hrowindex2[firsti_ri+k];
                while (hrowindex2[firstj_ri] != kk)
                  {
                    firstj++;
                    firstj_ri++;
                  }
                
                lfact[firstj] += sum[k];
                firstj++;
                firstj_ri++;
              }
          };
        } 
	for (size_t i2 = i1; i2 < last_same; i2++)
	  {
	    size_t first = hfirstinrow[i2] + last_same-i2-1;
	    size_t last = hfirstinrow[i2+1];
	    size_t j_ri = hfirstinrow_ri[i2] + last_same-i2-1;

	    for (auto j = first; j < last; j++, j_ri++)
	      {
		TM q = diag[i2] * lfact[j];
		diag[rowindex2[j_ri]] -= Trans (lfact[j]) * q;
	      }
	  }
       });
#endif





#ifdef CHOLESKY_PARALLEL_ATOMIC
    
    
    // first, find the transposed graph
    
    TableCreator<int> creator_trans(block_dependency.Size());
    for ( ; !creator_trans.Done(); creator_trans++)
      // for (int i : Range(block_dependency))
      ParallelFor (block_dependency.Size(), [&] (int i)
                   {
                     for (int j : block_dependency[i])
                       creator_trans.Add(j, i);
                   });
    auto block_dep_trans = creator_trans.MoveTable();

    /*
    static Timer tdep("paralleldep");
    static Timer tdep0("paralleldep0");
    static Timer tdep1("paralleldep1");
    static Timer tdep2("paralleldep2");
    static Timer tdep3("paralleldep3");
    */
    Array<MyMutex> locks(n);
    
    RunParallelDependency
      (block_dependency, block_dep_trans, [&] (int blocknr)
       {
        IntRange block = BlockDofs(blocknr);
        // RegionTracer reg(TaskManager::GetThreadId(), tdep, block.Size());

        size_t i1 = block.First(); 
        size_t last_same = block.Next();

	// rows in same block ...
	size_t mi = block.Size();  // last_same - i1;
        size_t nk = hfirstinrow[i1+1] - hfirstinrow[i1] + 1;

        ArrayMem<TM,1000> tmpmem(nk*nk);
        FlatMatrix<TM,ColMajor> tmp(nk, nk, tmpmem.Addr(0));

        {
          // RegionTracer reg0(TaskManager::GetThreadId(), tdep0, block.Size());
        tmp = TM(0.0);
        }        
	for (size_t j = 0; j < mi; j++)
	  {
            tmp(j,j) = diag[i1+j];
            tmp.Col(j).Range(j+1,nk) = FlatVector<TM>(nk-j-1, hlfact+hfirstinrow[i1+j]);
          }

        auto A11 = tmp.Rows(0,mi).Cols(0,mi);
        auto B   = tmp.Rows(mi,nk).Cols(0,mi);
        auto A22 = tmp.Rows(mi,nk).Cols(mi,nk);

        {
          // RegionTracer reg1(TaskManager::GetThreadId(), tdep1, block.Size());
          CalcLDL (A11);
          if (mi < nk)
            {
              CalcLDL_SolveL (A11,B);
              CalcLDL_A2 (A11.Diag(),B,A22);
            }
        }
        
        for (size_t j = 0; j < mi; j++)
          {
            diag[i1+j] = A11(j,j);
            FlatVector<TM>(nk-j-1, hlfact+hfirstinrow[i1+j]) = tmp.Col(j).Range(j+1,nk);
          };

	// merge rows
	size_t firsti_ri = hfirstinrow_ri[i1] + last_same-i1-1;
	size_t firsti = hfirstinrow[i1] + last_same-i1-1;
	size_t lasti = hfirstinrow[i1+1]-1;
	mi = lasti-firsti+1;
        
        {
          // RegionTracer reg2(TaskManager::GetThreadId(), tdep2, block.Size());
          // for (size_t j = 0; j < mi; j++)
          ParallelFor (mi, [=,&locks] (size_t j)
            {
              auto other_row = hrowindex2[firsti_ri+j];
              locks[other_row].lock();
              
              auto sum = A22.Col(j);
              
              // merge together
              size_t firstj = hfirstinrow[other_row];
              size_t firstj_ri = hfirstinrow_ri[other_row];
              
              for (size_t k = j+1; k < mi; k++)
                {
                  size_t kk = hrowindex2[firsti_ri+k];
                  while (hrowindex2[firstj_ri] != kk)
                    {
                      firstj++;
                      firstj_ri++;
                    }
                  
                  lfact[firstj] += sum[k];
                  firstj++;
                  firstj_ri++;
                }
              locks[other_row].unlock();
            }, mi > 50 ? TasksPerThread(1) : 1);
       }
       

        {
          // RegionTracer reg3(TaskManager::GetThreadId(), tdep3, block.Size());

          size_t num_other = hfirstinrow[i1+1] - (hfirstinrow[i1] + last_same-i1-1);
          size_t j_ri = hfirstinrow_ri[i1] + last_same-i1-1;

          auto hdiag = diag.Addr(0);
          ParallelFor (num_other, [=,&locks] (size_t j)
            {
              auto target_row = rowindex2[j_ri+j];
              locks[target_row].lock();

              for (auto i2 : block)
                {
                  size_t first = hfirstinrow[i2] + block.Next()-i2-1;
                  
                  TM q = hdiag[i2] * hlfact[first+j];
                  hdiag[target_row] -= Trans (hlfact[first+j]) * q;
                }
              
              locks[target_row].unlock();            
            }, num_other > 50 ? TasksPerThread(1) : 1);  
          
        }
        
       });
#endif


    
    
    
    /*
    size_t j = 0;
    for (size_t i = 0; i < n; i++)
      {
	TM ai = diag[i];
	size_t last = hfirstinrow[i+1];

	for ( ; j < last; j++)
          lfact[j] = lfact[j] * ai;
      }
    */
    ParallelFor (n, [&] (size_t i)
      {
        TM ai = diag[i];
        for (auto j : Range(hfirstinrow[i], hfirstinrow[i+1]))
          lfact[j] = lfact[j] * ai;
      }, TasksPerThread(5));

    if (n > 2000){
      cout << IM(4) << endl;
    }

    // task_manager -> StartWorkers();
  }





  


  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  Mult (const BaseVector & x, BaseVector & y) const
  {
    y = 0.0;
    MultAdd (TSCAL_VEC(1.0), x, y);
    return;
  }
  


  // template <>
  // void SparseCholesky<double, double, double> :: 

  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  SolveBlock (int bnr, FlatVector<TV> hy) const
  {
    cerr << "general form of solveblock not implemented" << endl;
  }

  template <>
  void SparseCholesky<double,double,double> :: 
  SolveBlock (int bnr, FlatVector<double> hy) const
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
      {
        auto val = temp(j);
        AtomicAdd (hy(extdofs[j]), -val);
      }
  }


  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  SolveBlockT (int bnr, FlatVector<TV> hy) const
  {
    cerr << "general form of solveblock not implemented" << endl;
  }


  template <>
  void SparseCholesky<double, double, double> :: 
  SolveBlockT (int bnr, FlatVector<> hy) const
  {
    const size_t * hfirstinrow = &firstinrow[0];
    const size_t * hfirstinrow_ri = &firstinrow_ri[0];


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


  /*
  // a simple lock-free queue (but we don't need it)
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
    }

    bool Pop (T & out)
    {
      while (1)
        {
          if (rcnt >= data.Size()) return false;
          
          int oldval = 1;
          if (ok[int(rcnt)].compare_exchange_weak (oldval, 0))
            {
              int mypos = rcnt;
              rcnt++;
              out = data[mypos];
              return true;
            }
        }
    }
  };
  */





  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  SolveReordered (FlatVector<TVX> hy) const
  {
    static Timer timer1("SparseCholesky<d,d,d>::MultAdd fac1");
    static Timer timer2("SparseCholesky<d,d,d>::MultAdd fac2");

    /*
    // sequential verision 
    for (int i = 0; i < blocks.Size()-1; i++)
      SolveBlock (i, hy);
    */

    /*

    // first parallel version with dependency graph
    RunParallelDependency (block_dependency, 
                           [&] (int nr) 
                           {
                             SolveBlock(nr, hy); 
                           });
    */



    /*
    // first parallel version with dependency graph and profiling

    class ProfileData
    {
    public:
      double tstart, tend;
      int size, extsize;
    };

    Array<ProfileData> prof(blocks.Size()-1);
    
    double tstart = omp_get_wtime();

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


    // parallel version with refined tasks (micri-dependency)
    /*
    timer0.Start();    
    RunParallelDependency (micro_dependency, 
                           [&] (int nr) 
                           { ; } );
    timer0.Stop();    
    */
    timer1.Start();

    RunParallelDependency (micro_dependency, micro_dependency_trans,
                           [&,hy] (int nr) 
                           {
                             auto task = microtasks[nr];
                             size_t blocknr = task.blocknr;
                             auto range = BlockDofs (blocknr);
                             if (range.Size()==0) return;
                             
                             // if (task.solveL)
                             if (task.type == MicroTask::LB_BLOCK)
                               { // first L, then B

                                 auto extdofs = BlockExtDofs (blocknr);
                                 VectorMem<520,TVX> temp(extdofs.Size());
                                 temp = 0;

                                 for (auto i : range)
                                   {
                                     TVX hyi = hy(i);
                                     
                                     size_t size = range.end()-i-1;
                                     if (size > 0)
                                       {
                                         FlatVector<TM> vlfact(size, &lfact[firstinrow[i]]);
                                         
                                         auto hyr = hy.Range(i+1, range.end());
                                         for (size_t j = 0; j < size; j++)
                                           hyr(j) -= Trans(vlfact(j)) * hyi;
                                       }
                                     if (extdofs.Size() == 0)
                                       {
                                         // cerr << "should not be here" << endl;
                                         continue;
                                       }
                                     size_t first = firstinrow[i] + range.end()-i-1;
                                     FlatVector<TM> ext_lfact (extdofs.Size(), &lfact[first]);
                                     for (size_t j = 0; j < temp.Size(); j++)
                                       temp(j) += Trans(ext_lfact(j)) * hyi;
                                   }
                                 
                                 for (size_t j : Range(extdofs))
                                   AtomicAdd (hy(extdofs[j]), -temp(j));
                               }
                             
                             else if (task.type == MicroTask::L_BLOCK)
                               {
                                 
                                 for (auto i : range)
                                   {
                                     size_t size = range.end()-i-1;
                                     if (size == 0) continue;
                                     FlatVector<TM> vlfact(size, &lfact[firstinrow[i]]);

                                     TVX hyi = hy(i);
                                     auto hyr = hy.Range(i+1, range.end());
                                     for (size_t j = 0; j < hyr.Size(); j++)
                                       hyr(j) -= Trans(vlfact(j)) * hyi;
                                   }

                               }
                             
                             else 

                               {
                                 auto all_extdofs = BlockExtDofs (blocknr);
                                 if (all_extdofs.Size() != 0)
                                   {
                                     auto myr = Range(all_extdofs).Split (task.bblock, task.nbblocks);
                                     auto extdofs = all_extdofs.Range(myr);
 
                                     VectorMem<520,TVX> temp(extdofs.Size());
                                     temp = 0;
                                     
                                     for (auto i : range)
                                       {
                                         size_t first = firstinrow[i] + range.end()-i-1;
                                         
                                         FlatVector<TM> ext_lfact (all_extdofs.Size(), &lfact[first]);
 
                                         TVX hyi = hy(i);
                                         for (size_t j = 0; j < temp.Size(); j++)
                                           temp(j) += Trans(ext_lfact(myr.begin()+j)) * hyi;
                                       }
                                     
                                     for (size_t j : Range(extdofs))
                                       AtomicAdd (hy(extdofs[j]), -temp(j));
                                   }
                               }

                           });

    timer1.Stop();


    // solve with the diagonal
    const TM * hdiag = diag.Data();
    ParallelFor (hy.Size(), [&] (int i)
                 {
                   TVX tmp = hdiag[i] * hy[i];
                   hy[i] = tmp;
                 });


    timer2.Start();

    /*
      // sequential version 
    for (int i = blocks.Size()-2; i >= 0; i--)
      SolveBlockT (i, hy);
    */

    // advanced parallel version 
    RunParallelDependency (micro_dependency_trans, micro_dependency,
                           [&,hy] (int nr) 
                           {
                             auto task = microtasks[nr];
                             int blocknr = task.blocknr;
                             auto range = BlockDofs (blocknr);
                             if (range.Size()==0) return;
                             
                             if (task.type == MicroTask::LB_BLOCK)
                               { // first B then L

                                 auto extdofs = BlockExtDofs (blocknr);
                                 
                                 VectorMem<520,TVX> temp(extdofs.Size());
                                 for (auto j : Range(extdofs))
                                   temp(j) = hy(extdofs[j]);

                                 if (extdofs.Size())
                                   for (auto i : range)
                                     {
                                       size_t first = firstinrow[i] + range.end()-i-1;
                                       FlatVector<TM> ext_lfact (extdofs.Size(), &lfact[first]);
                                       
                                       TVX val(0.0);
                                       for (auto j : Range(extdofs))
                                         val += ext_lfact(j) * temp(j);
                                       hy(i) -= val;
                                     }
                                 for (size_t i = range.end()-1; i-- > range.begin(); )
                                   {
                                     size_t size = range.end()-i-1;
                                     if (size == 0) continue;
                                     FlatVector<TM> vlfact(size, &lfact[firstinrow[i]]);
                                     auto hyr = hy.Range(i+1, range.end());

                                     TVX hyi = hy(i);
                                     for (size_t j = 0; j < vlfact.Size(); j++)
                                       hyi -= vlfact(j) * hyr(j);
                                     hy(i) = hyi;
                                   }
                                 
                               }
                             else if (task.type == MicroTask::L_BLOCK)                               
                               {
                                 // for (int i = range.end()-1; i >= range.begin(); i--)
                                 if (range.Size() > 0) // for case [0,0)
                                 for (size_t i = range.end()-1; i-- > range.begin(); )
                                   {
                                     size_t size = range.end()-i-1;
                                     if (size == 0) continue;
                                     FlatVector<TM> vlfact(size, &lfact[firstinrow[i]]);
                                     auto hyr = hy.Range(i+1, range.end());

                                     TVX hyi = hy(i);
                                     for (size_t j = 0; j < vlfact.Size(); j++)
                                       hyi -= vlfact(j) * hyr(j);
                                     hy(i) = hyi;
                                   }

                               }
                             
                             else 

                               {
                                 auto all_extdofs = BlockExtDofs (blocknr);
                                 if (all_extdofs.Size() != 0)
                                   {
                                     auto myr = Range(all_extdofs).Split (task.bblock, task.nbblocks);
                                     auto extdofs = all_extdofs.Range(myr);
                                     
                                     VectorMem<520,TVX> temp(extdofs.Size());
                                     for (auto j : Range(extdofs))
                                       temp(j) = hy(extdofs[j]);
    
                                     for (auto i : range)
                                       {
                                         size_t first = firstinrow[i] + range.end()-i-1;
                                         FlatVector<TM> ext_lfact (all_extdofs.Size(), &lfact[first]);
    
                                         TVX val(0.0);
                                         for (auto j : Range(extdofs))
                                           val += ext_lfact(myr.begin()+j) * temp(j);
                                         AtomicAdd (hy(i), -val);
                                       }
                                   }
                               }
                           });

    timer2.Stop();

    
  }

  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  MultAdd (TSCAL_VEC s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("SparseCholesky<d,d,d>::MultAdd");
    RegionTimer reg (timer);
    timer.AddFlops (2.0*lfact.Size());

    // int n = Height();
    
    const FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();

    int nused = this->nused;
    Vector<TVX> hy1(nused);
    FlatVector<TVX> hy(hy1);

    ParallelFor (Range(height), [&] (int i)
                 {
                   if (order[i] != -1)
                     hy(order[i]) = fx(i);
                 });


    SolveReordered(hy);

    if (inner)
      {
        /*
	for (int i = 0; i < n; i++)
	  if (inner->Test(i))
	    fy(i) += s * hy(order[i]);
        */
        ParallelFor (Range(height), [&] (int i)
                     {
                       if (inner->Test(i))
                         fy(i) += s * hy(order[i]);
                     });
      }
    else if (cluster)
      {
	for (int i = 0; i < height; i++)
	  if ((*cluster)[i])
	    fy(i) += s * hy(order[i]);
      }
    else
      {
	// for (int i = 0; i < n; i++)
        // fy(i) += s * hy(order[i]);
        ParallelFor (Range(height), [&] (int i)
                     {
                       if (order[i] != -1)
                         fy(i) += s * hy(order[i]);
                     });

      }

  }
  


  template <class TM, class TV_ROW, class TV_COL>
  void SparseCholesky<TM, TV_ROW, TV_COL> :: 
  Smooth (BaseVector & u, const BaseVector & f, BaseVector & y) const
  {
    static Timer t("SparseCholesky::Smooth");
    RegionTimer reg(t);

    if(dynamic_pointer_cast<const SparseMatrixSymmetric<TM,TV>>(this->matrix.lock()))
      {
        // use the original one ...
        SparseFactorization::Smooth(u,f,y);
        return;
      }

    const FlatVector<TVX> fu = u.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();
    
    Vector<TVX> hy(this->nused);
    auto spmat = dynamic_pointer_cast<const SparseMatrix<TM,TV,TV>> (this->matrix.lock());
    if(!spmat)
      throw Exception("A matrix not available any more, needed for Smooth!");
    auto & hmat = *spmat;
    
    ParallelFor (this->nused, [&] (int i)
                 {
                   hy(i) = fy(inv_order[i]) - hmat.RowTimesVector(inv_order[i], fu);
                 });
    
    SolveReordered(hy);

    ParallelFor (this->nused, [&] (int i)
                 {
                   fu(inv_order[i]) += hy(i);
                 });
  }

  



  SparseFactorization ::     
  SparseFactorization (shared_ptr<const BaseSparseMatrix> amatrix,
		       shared_ptr<BitArray> ainner,
		       shared_ptr<const Array<int>> acluster)
    : matrix(amatrix),
      inner(ainner), cluster(acluster)
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
    static Timer t("SparseFactorization::Smooth");
    RegionTimer reg(t);
    {
    auto hvec1 = u.CreateVector();
    auto hvec2 = u.CreateVector();

    hvec1 = y;
    matrix.lock()->MultAdd1 (-1, u, hvec1, inner.get(), cluster.get());

    hvec2 = (*this) * hvec1;
    u += hvec2;
    
    matrix.lock()->MultAdd2 (-1, hvec2, y, inner.get(), cluster.get());
    }
  }
  


















  template <class TM>
  void SparseCholeskyTM<TM> :: Set (int i, int j, const TM & val)
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
    size_t first = firstinrow[i];
    size_t first_ri = firstinrow_ri[i];
    size_t last = firstinrow[i+1];
    
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


  template <class TM>
  const TM & SparseCholeskyTM<TM> :: Get (int i, int j) const
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

    size_t first = firstinrow[i];
    size_t first_ri = firstinrow_ri[i];
    size_t last = firstinrow[i+1];
  
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

  template <class TM>
  ostream & SparseCholeskyTM<TM> :: Print (ostream & ost) const
  {
    int n = Height();

    for (int i = 0; i < n; i++)
      {
	ost << i << ": " << order[i] << " diag = "
	    << diag[i] << endl;
      }
    ost << endl;
  
    size_t j = 1;
    for (int i = 1; i <= n; i++)
      {
	size_t j_ri = firstinrow_ri[i-1];
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


  template <class TM>
  SparseCholeskyTM<TM> :: ~SparseCholeskyTM()
  {
    delete mdo;
  }




  static RegisterClassForArchive<SparseCholesky<double>, SparseCholeskyTM<double>> regscd;
  static RegisterClassForArchive<SparseCholesky<Complex>, SparseCholeskyTM<Complex>> regscc;


  template class SparseCholesky<double>;
  template class SparseCholesky<Complex>;
  template class SparseCholesky<double, Complex, Complex>;

  template class SparseCholeskyTM<double>;
  template class SparseCholeskyTM<Complex>;

#if MAX_SYS_DIM >= 1
  template class SparseCholesky<Mat<1,1,double> >;
  template class SparseCholesky<Mat<1,1,Complex> >;
  template class SparseCholeskyTM<Mat<1,1,double> >;
  template class SparseCholeskyTM<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class SparseCholesky<Mat<2,2,double> >;
  template class SparseCholesky<Mat<2,2,Complex> >;
  template class SparseCholeskyTM<Mat<2,2,double> >;
  template class SparseCholeskyTM<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class SparseCholesky<Mat<3,3,double> >;
  template class SparseCholesky<Mat<3,3,Complex> >;
  template class SparseCholeskyTM<Mat<3,3,double> >;
  template class SparseCholeskyTM<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class SparseCholesky<Mat<4,4,double> >;
  template class SparseCholesky<Mat<4,4,Complex> >;
  template class SparseCholeskyTM<Mat<4,4,double> >;
  template class SparseCholeskyTM<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class SparseCholesky<Mat<5,5,double> >;
  template class SparseCholesky<Mat<5,5,Complex> >;
  template class SparseCholeskyTM<Mat<5,5,double> >;
  template class SparseCholeskyTM<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class SparseCholesky<Mat<6,6,double> >;
  template class SparseCholesky<Mat<6,6,Complex> >;
  template class SparseCholeskyTM<Mat<6,6,double> >;
  template class SparseCholeskyTM<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class SparseCholesky<Mat<7,7,double> >;
  template class SparseCholesky<Mat<7,7,Complex> >;
  template class SparseCholeskyTM<Mat<7,7,double> >;
  template class SparseCholeskyTM<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class SparseCholesky<Mat<8,8,double> >;
  template class SparseCholesky<Mat<8,8,Complex> >;
  template class SparseCholeskyTM<Mat<8,8,double> >;
  template class SparseCholeskyTM<Mat<8,8,Complex> >;
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













