/**************************************************************************/
/* File:   blockjacobi.cpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   14. Aug. 2002                                                  */
/**************************************************************************/


#include <la.hpp>


namespace ngla
{
  
  static mutex buildingblockupdate_mutex;


  BaseBlockJacobiPrecond :: 
  BaseBlockJacobiPrecond (shared_ptr<Table<int>> ablocktable)
    : blocktable(ablocktable)
  {
    maxbs = 0;
    for (auto entry : *blocktable)
      if (entry.Size() > maxbs)
        maxbs = entry.Size();
    GetMemoryTracer().Track(*blocktable, "BlockTable");
    GetMemoryTracer().Track(block_coloring, "BlockColoring");
    GetMemoryTracer().Track(color_balance, "ColorBalance");
  }


  BaseBlockJacobiPrecond :: 
  ~BaseBlockJacobiPrecond ()
  { ; }


  int BaseBlockJacobiPrecond ::
  Reorder (FlatArray<int> block, const MatrixGraph & graph,
	   FlatArray<int> block_inv,
	   LocalHeap & lh)
  {
    try
      {
	// a cheap reordering algorithm:
	
	// bool print = 0;
	
	size_t n = block.Size();
	
	void * heapp = lh.GetPointer();


	FlatArray<int> reorder(n, lh), newnum(n, lh), dist(n, lh), cluster0(n,lh);
	// FlatArray<int> block_inv(graph.Size(), lh);

	for (size_t i = 0; i < n; i++)
	  if (block[i] >= 0 && block[i] < graph.Size())
	    {
	      if (block_inv[block[i]] != -1)
		{
		  cout << IM(1) << "block has double elements " << endl;
		  cout << IM(1) << block_inv[block[i]] << " and " << i << endl;
		  cout << IM(1) << block << endl;
		}
	      block_inv[block[i]] = i;
	    }
	  else
	    {
	      cerr << "block[" << i << "] out of range" << endl;
	      cerr << "block = " << block << endl;
	      (*testout) << "block[" << i << "] out of range" << endl;
	      (*testout) << "block = " << block << endl;
	    }

	/*
	// debug only:
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	if (i != j && block[i] == block[j])
	{
	cout << "block has double elements" << endl;
	cout << "i = " << i << ", j = " << j << endl;
	cout << block << endl;
	}
	*/

	// check for sperated blocks:
	cluster0 = 0;
	cluster0[0] = 1;

	bool changed;
	do
	  {
	    changed = 0;
	    for (size_t j = 0; j < n; j++)
	      {
		FlatArray row = graph.GetRowIndices(block[j]);
		for (size_t k = 0; k < row.Size(); k++)
		  {
		    int kk = block_inv[row[k]];
		    if (kk >= 0 && unsigned(kk) < n)
		      if (block[kk] == row[k])
			{
			  if (cluster0[j] != cluster0[kk]) 
			    {
			      cluster0[j] = 1;
			      cluster0[kk] = 1;
			      changed = 1;
			    }
			  //			  if (cluster0[j])  cluster0[kk] = 1;
			  //			  if (cluster0[kk]) cluster0[j] = 1;
			}
		  }
	      }
	  }
	while (changed);
	
	/*
	// cout << "cluster0, 2 = " << cluster0 << endl;
	cluster0 = 0;
	cluster0[0] = 1;
	for (i = 0; i < n; i++)
	for (j = 0; j < n; j++)
	for (k = 0; k < j; k++)
	if ((cluster0[j] || cluster0[k]) && 
	connected.Test(j,k))
	{
	cluster0[j] = 1;
	cluster0[k] = 1;
	}

	// cout << "cluster0, 1 = " << cluster0 << endl;
	*/

	size_t cnt = 0;
	for (size_t i = 0; i < n; i++)
	  if (cluster0[i])
	    {
	      newnum[cnt] = block[i];
	      cnt++;
	    }

	if (cnt < n)
	  {
	    // separated clusters
	    size_t cnt2 = cnt;
	    for (size_t i = 0; i < n; i++)
	      if (!cluster0[i])
		newnum[cnt2++] = block[i];
	    
	    for (size_t i = 0; i < n; i++)
	      block[i] = newnum[i];

	    for (size_t i = 0; i < n; i++)
	      block_inv[block[i]] = -1;

	    // block = newnum;

	    lh.CleanUp(heapp);
	    
	    int bw1 = 
	      Reorder (FlatArray<int> (cnt, &block[0]), graph, block_inv, lh);
	    int bw2 =
	      Reorder (FlatArray<int> (cnt2-cnt, &block[cnt]), graph, block_inv, lh);
	    
	    return max2(bw1, bw2);
	  }
	

	// compute distance function

	size_t pstart = 0;
	for (int step = 0; step < 3; step++)
	  {
	    /*
	      dist = n+1;
	      dist[pstart] = 0;
	
	      bool changed;

	      do
	      {
	      changed = 0;
	      for (i = 0; i < n; i++)
	      for (j = 0; j < n; j++)
	      if (i != j && connected.Test (i,j))
	      {
	      if (dist[i] > dist[j]+1)
	      {
	      dist[i] = dist[j]+1;
	      changed = 1;
	      }
	      else if (dist[j] > dist[i]+1)
	      {
	      dist[j] = dist[i]+1;
	      changed = 1;
	      }
	      }
	      }
	      while (changed);
	    */

	    bool changed;
	    dist = n+1;
	    dist[pstart] = 0;
	    do
	      {
		changed = 0;
		for (size_t i = 0; i < n; i++)
		  {
		    FlatArray row = graph.GetRowIndices(block[i]);

		    for (size_t j = 0; j < row.Size(); j++)
		      {
			int jj = block_inv[row[j]];

			if (jj >= 0 && unsigned(jj) < n)
			  if (block[jj] == row[j])		    
			    {
			      if (dist[i] > dist[jj]+1)
				{
				  dist[i] = dist[jj]+1;
				  changed = 1;
				}
			      else if (dist[jj] > dist[i]+1)
				{
				  dist[jj] = dist[i]+1;
				  changed = 1;
				}
			    }
		      }
		  }
	      }
	    while (changed);

	    
	    int maxval = 0;
	    for (size_t i = 0; i < n; i++)
	      if (dist[i] > maxval)
		{
		  maxval = dist[i];
		  pstart = i;
		}
	    if (maxval > n)
	      {
		cerr << "Blockjacobi, reorder: separated block" << endl;
		cout << IM(1) << "block: " << block << endl;

		(*testout) << "Blockjacobi, reorder: separated block" << endl;
		(*testout) << "block: " << block << endl;
	      }	      
	  }
	

	cnt = 0;
	for (size_t i = 0; i < n; i++)
	  for (size_t j = 0; j < n; j++)
	    if (dist[j] == i)
	      {
		reorder[cnt] = j;
		cnt++;
	      }
	if (cnt != n) cerr << "BlockJac, reorder: n = " << n << " != cnt = " << cnt << endl;

	int bw = 1; 
	/*
	  for (i = 0; i < n; i++)
	  for (j = 0; j < n; j++)
	  if (connected.Test(reorder[i], reorder[j]))
	  bw = max2(bw, abs(i-j)+1);
	*/
	for (size_t i = 0; i < n; i++)
	  newnum[reorder[i]] = i;
	for (size_t j = 0; j < n; j++)
	  {
	    FlatArray row = graph.GetRowIndices(block[j]);
	    for (size_t k = 0; k < row.Size(); k++)
	      {
		int kk = block_inv[row[k]];
		if (kk >= 0 && unsigned(kk) < n)
		  if (block[kk] == row[k])
		    {
		      int inv_j = newnum[j];
		      int inv_k = newnum[kk];
		      bw = max2(bw, abs(inv_j-inv_k)+1);
		    }
	      }
	  }


	for (size_t i = 0; i < n; i++)
	  newnum[i] = block[reorder[i]];
	for (size_t i = 0; i < n; i++)
	  block[i] = newnum[i];


	lh.CleanUp(heapp);

	for (size_t i = 0; i < n; i++)
	  block_inv[block[i]] = -1;

	return bw;
      }
    catch (Exception & e)
      {
	e.Append ("in Reorder\n");
	throw;
      }


  }



  ///
  template <class TM, class TV_ROW, class TV_COL>
  BlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  BlockJacobiPrecond (shared_ptr<const SparseMatrix<TM,TV_ROW,TV_COL>> amat, 
		      shared_ptr<Table<int>> ablocktable, bool cumulate_block_diags)
    : BaseBlockJacobiPrecond(ablocktable), mat(amat), 
      invdiag(ablocktable->Size())
  {
    static Timer t("BlockJacobiPrecond ctor"); RegionTimer reg(t);
    static Timer tinv("BlockJacobiPrecond ctor inv");
    static Timer tget("BlockJacobiPrecond ctor get");
    static Timer tprep("BlockJacobiPrecond ctor prep");
    static Timer tpar("BlockJacobiPrecond ctor par");
    cout << IM(3) << "BlockJacobi Preconditioner, constructor called, #blocks = " << blocktable->Size() << endl;


    // double prevtime = WallTime();

    
    // find nze element in all blocks together 
    /*
    nze = 0;
    for (auto block : *blocktable)
      for (auto row : block)
	nze += amat.GetRowIndices(row).Size();
    */
    nze = 
      ParallelReduce (blocktable->Size(),
                      [&] (size_t i)
                      {
                        size_t nze = 0;
                        for (auto row : (*blocktable)[i])
                          nze += amat->GetRowIndices(row).Size();
                        return nze;
                      },
                      [] (size_t a, size_t b) { return a+b; },
                      size_t(0));

    /*
    size_t totmem = 0;
    for (auto i : Range (*blocktable))
      totmem += sqr ((*blocktable)[i].Size());
    */
    size_t totmem = 
      ParallelReduce (blocktable->Size(),
                      [&] (size_t i) { return sqr ((*blocktable)[i].Size()); },
                      [] (size_t a, size_t b) { return a+b; },
                      size_t(0));

    bigmem.SetSize(totmem);

    cout << IM(5) << "avg entrysize:   " << blocktable->AsArray().Size()/blocktable->Size() << endl;
    cout << IM(5) << "avg entrysize^2: " << totmem/blocktable->Size() << endl;    
    
    totmem = 0;
    for (auto i : Range (*blocktable))
      {
        size_t bs = (*blocktable)[i].Size();
        new ( & invdiag[i] ) FlatMatrix<TM> (bs, bs, bigmem.Addr(totmem));
        totmem += sqr (bs);
      }



    /** Get diagonal blocks **/
    SharedLoop2 sl1(blocktable->Size());
    ParallelJob
      ([&] (const TaskInfo & ti)
       {
         NgProfiler::StartThreadTimer (tpar, TaskManager::GetThreadId());         
         for (int i : sl1)
       {
         NgProfiler::StartThreadTimer (tprep, TaskManager::GetThreadId());

        auto blocki = (*blocktable)[i];
        QuickSort (blocki);
	if (!blocki.Size()) 
	  {
            NgProfiler::StopThreadTimer (tprep, TaskManager::GetThreadId());                             
	    invdiag[i] = 0;
	    // return;
            continue;
	  }
	
        FlatMatrix<TM> & blockmat = invdiag[i];
        NgProfiler::StopThreadTimer (tprep, TaskManager::GetThreadId());                 
        NgProfiler::StartThreadTimer (tget, TaskManager::GetThreadId());
	for (size_t j = 0; j < blocki.Size(); j++)
	  for (size_t k = 0; k < blocki.Size(); k++)
	    blockmat(j,k) = (*mat)(blocki[j], blocki[k]);
        NgProfiler::StopThreadTimer (tget, TaskManager::GetThreadId());                         
        // }, TasksPerThread(10));
       }
         NgProfiler::StopThreadTimer (tpar, TaskManager::GetThreadId());                  
       } );


    /**
       Cumulate diagonal blocks across processors - works ONLY IF:
	 - Blocks are consistent across processors
	 - The blocks are already numbered consistently
	 - For each block, all it's DOFs are shared between the same procs.
	   This excludes, for example:
	      Any blocks that reach from a subdomain interface in the interior of any domain
	      Blocks that contain all DOFs on a facet shared between two procs, as well as all
	      DOFs in that facet's BBND and BBBND nodes if any of those are also shared with some
	      third block.
        !! THIS ALSO ONLY ALLOWS FOR MULT/MULTADD, NOT FOR GSSMOOTH/GSSMOOTHBACK !!
	Calling GSSmooth/GSSmoothback leads to undefined behavior
     **/
    if (cumulate_block_diags) {
      if (auto ppds = amat->GetParallelDofs()) { // without attached ParallelDofs, we cannot cumulate blocks
	const auto & pds = *ppds;
	auto all_dps = pds.GetDistantProcs();
	const auto& btab = *blocktable;

	// iterates through all blocks that are shared with someone and calls a lambda function
	auto iterate_ex_blocks = [&](auto lam) {
	  for (auto block_num : Range(btab)) {
	    auto block = btab[block_num];
	    if (block.Size()) {
	      auto dps = pds.GetDistantProcs(block[0]);
	      for (auto p : dps) {
		lam(block_num, block, p);
		// sds[all_dps.Pos(p)] += sqr(block.Size());
	      }
	    }
	  }
	};

	// find size of data to send/recv
	Array<int> sds(all_dps.Size()); sds = 0; // shared data size
	iterate_ex_blocks([&](auto block_num, auto block, auto p){
	    sds[all_dps.Pos(p)] += sqr(block.Size());
	  });

	// allocate send/recv data
	Table<TM> send_data(sds), recv_data(sds);

	// write send-data
	sds = 0;
	iterate_ex_blocks([&](auto block_num, auto block, auto p){
	    auto pos = all_dps.Pos(p);
	    auto blocks = block.Size();
	    // auto buf_block = send_data[pos].Part(sds[pos], ds);
	    auto buf_block = FlatMatrix<TM>(blocks, blocks, send_data[pos].Addr(sds[pos]));
	    auto diag_block = invdiag[block_num];
	    sds[pos] += sqr(blocks);
	    buf_block = diag_block;
	  });

	// send-recv
        NgMPI_Requests rsend, rrecv;
	auto comm = pds.GetCommunicator();
	for (auto kp : Range(all_dps)) {
	  rsend += comm.ISend(send_data[kp], all_dps[kp], NG_MPI_TAG_SOLVE);
	  rrecv += comm.IRecv(recv_data[kp], all_dps[kp], NG_MPI_TAG_SOLVE);
	}

	// wait for recvs to finish and add to diagonal blocks
	rrecv.WaitAll();
	sds = 0;
	iterate_ex_blocks([&](auto block_num, auto block, auto p){
	    auto pos = all_dps.Pos(p);
	    auto blocks = block.Size();
	    // auto buf_block = recv_data[pos].Part(sds[pos], ds);
	    auto buf_block = FlatMatrix<TM>(blocks, blocks, recv_data[pos].Addr(sds[pos]));
	    auto diag_block = invdiag[block_num];
	    sds[pos] += sqr(blocks);
	    diag_block += buf_block;
	  });

	rsend.WaitAll(); // wait for sends to finish !!
      }
    }

    /** Invert diagonal blocks **/
    SharedLoop2 sl2(blocktable->Size());
    ParallelJob
      ([&] (const TaskInfo & ti)
       {
         NgProfiler::StartThreadTimer (tpar, TaskManager::GetThreadId());         
         for (auto i : sl2) {
	     NgProfiler::StartThreadTimer (tinv, TaskManager::GetThreadId());
	     FlatMatrix<TM> & blockmat = invdiag[i];
	     CalcInverse (blockmat);
	     NgProfiler::StopThreadTimer (tinv, TaskManager::GetThreadId());        
	   }
         NgProfiler::StopThreadTimer (tpar, TaskManager::GetThreadId());                  
       } );

    cout << IM(3) << "\rBuilding block " << blocktable->Size() << "/" << blocktable->Size() << flush;
    *testout << "block coloring";

    static Timer tcol("BlockJacobi-coloring");
    tcol.Start();

    size_t nblocks = blocktable->Size();
    Array<int> coloring(nblocks);
    coloring = -1;

    int maxcolor = 0;
    int basecol = 0;
    Array<unsigned int> mask(mat->Width());
    size_t found = 0;

    do
      {
        mask = 0;
        
        for (auto i : Range(nblocks))
          {
            if (coloring[i] >= 0) continue;

            unsigned check = 0;
	    for (int d : (*blocktable)[i] )              
              check |= mask[d];
            
            if (check != UINT_MAX) // 0xFFFFFFFF)
              {
                found++;
                unsigned checkbit = 1;
                int color = basecol;
                while (check & checkbit)
                  {
                    color++;
                    checkbit *= 2;
                  }

                coloring[i] = color;
                if (color > maxcolor) maxcolor = color;
                
                for (int d : (*blocktable)[i] )
                  for(auto coupling : mat->GetRowIndices(d))
                    mask[coupling] |= checkbit;
              }
          }
        basecol += 8*sizeof(unsigned int); // 32;
      }
    while (found < nblocks);
    tcol.Stop();    

    TableCreator<int> creator(maxcolor+1);
    for ( ; !creator.Done(); creator++)
      for (size_t i = 0; i < nblocks; i++)
          creator.Add (coloring[i], i);
    block_coloring = creator.MoveTable();

    cout << IM(4) << " using " << maxcolor+1 << " colors" << endl;

    // calc balancing:

    color_balance.SetSize (block_coloring.Size());

    for (auto c : Range (block_coloring))
      {
        color_balance[c].Calc (block_coloring[c].Size(),
                               [&] (size_t bi)
                               {
                                 int costs = 0;
                                 size_t blocknr = block_coloring[c][bi];

                                 for (auto d : (*blocktable)[blocknr])
                                   costs += mat->GetRowIndices(d).Size();
                                 return costs;
                               });

      }

    GetMemoryTracer().Track(bigmem, "InvDiag");
    cout << IM(3) << "\rBlockJacobi Preconditioner built" << endl;
  }

  ///
  template <class TM, class TV_ROW, class TV_COL>
  BlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  ~BlockJacobiPrecond () 
  {
    ;
  }


  
  template <class TM, class TV_ROW, class TV_COL>
  void BlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  MultAdd (TSCAL s, const BaseVector & x, BaseVector & y) const 
  {
    static Timer timer("BlockJacobi::MultAdd");
    RegionTimer reg (timer);
    
    x.Cumulate();
    y.Cumulate();

    FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();


    for (int c : Range(block_coloring))        
      {
        ParallelForRange
          (color_balance[c],  [&] (IntRange r) 
           {
             Vector<TVX> hxmax(maxbs);
             Vector<TVX> hymax(maxbs);
             
             for (int i : block_coloring[c].Range(r))
               {
                 int bs = (*blocktable)[i].Size();
                 if (!bs) continue;
                 
                 FlatVector<TVX> hx = hxmax.Range(0,bs); 
                 FlatVector<TVX> hy = hymax.Range(0,bs); 
                 
                 for (int j = 0; j < bs; j++)
                   hx(j) = fx((*blocktable)[i][j]);
                 
                 hy = invdiag[i] * hx;
                 
                 for (int j = 0; j < bs; j++)
                   fy((*blocktable)[i][j]) += s * hy(j);
               }
           });
      }
  }


  template <class TM, class TV_ROW, class TV_COL>
  void BlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  MultTransAdd (TSCAL s, const BaseVector & x, BaseVector & y) const 
  {
    static Timer timer("BlockJacobi::MultTransAdd");
    RegionTimer reg (timer);

    x.Cumulate();
    y.Cumulate();

    FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();

    for (size_t c = 0; c < block_coloring.Size(); c++)
      {
        ParallelForRange
          (color_balance[c], [&] (IntRange r) 
           {
             Vector<TVX> hxmax(maxbs);
             Vector<TVX> hymax(maxbs);
             
             FlatArray<int> blocks = block_coloring[c];
             for (auto ii : r) 
               {
                 size_t i = blocks[ii];
                 auto block = (*blocktable)[i];
                 size_t bs = block.Size();
                 if (!bs) continue;
                 
                 FlatVector<TVX> hx = hxmax.Range(0,bs); 
                 FlatVector<TVX> hy = hymax.Range(0,bs); 
                 
                 for (size_t j = 0; j < bs; j++)
                   hx(j) = fx(block[j]);
                 
                 hy = Trans(invdiag[i]) * hx;
                 
                 for (size_t j = 0; j < bs; j++)
                   fy(block[j]) += s * hy(j);
               }
           });
      }
  }



  template <class TM, class TV_ROW, class TV_COL>
  void BlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  GSSmooth (BaseVector & x, const BaseVector & b,
	    int steps) const 
  {
    static Timer timer ("BlockJacobiPrecond::GSSmooth");
    RegionTimer reg(timer);
    timer.AddFlops (nze);
    
    FlatVector<TVX> fb = b.FV<TVX> (); 
    FlatVector<TVX> fx = x.FV<TVX> ();

#ifdef OLD
    for (int k = 0; k < steps; k++)
      for (int c : Range(block_coloring))              
        {
          /*
          ParallelForRange
            (color_balance[c], [&] (IntRange r)
             {
               VectorMem<100,TVX> hxmax(maxbs);
               VectorMem<100,TVX> hymax(maxbs);
               
               for (size_t i : block_coloring[c].Range(r))
                 {
                   FlatArray<int> block = (*blocktable)[i];
                   size_t bs = block.Size();
                   if (!bs) continue;
                   
                   FlatVector<TVX> hx = hxmax.Range(0,bs);
                   FlatVector<TVX> hy = hymax.Range(0,bs);
                   
                   for (size_t j = 0; j < bs; j++)
                     {
                       auto jj = block[j];
                       hx(j) = fb(jj) - mat.RowTimesVector (jj, fx);
                     }
                   
                   hy = (invdiag[i]) * hx;
                   fx(block) += hy;
                 }
               
             });
          */

          auto col = block_coloring[c];
          SharedLoop2 sl(col.Range());
          
          task_manager -> CreateJob
              ( [&] (const TaskInfo & ti) 
                {
                  VectorMem<100,TVX> hxmax(maxbs);
                  VectorMem<100,TVX> hymax(maxbs);

                  for (auto mynr : sl)
                    {
                      size_t i = col[mynr];
                      FlatArray<int> block = (*blocktable)[i];
                      size_t bs = block.Size();
                      if (!bs) continue;
                      
                      FlatVector<TVX> hx = hxmax.Range(0,bs);
                      FlatVector<TVX> hy = hymax.Range(0,bs);
                      
                      for (size_t j = 0; j < bs; j++)
                        {
                          auto jj = block[j];
                          hx(j) = fb(jj) - mat->RowTimesVector (jj, fx);
                        }
                      
                      hy = (invdiag[i]) * hx;
                      fx(block) += hy;
                    }
                });
            
        }
#endif

    Array<SharedLoop2> loops(block_coloring.Size());
    
    for (int k = 0; k < steps; k++)
      {
        for (int c : Range(block_coloring))
          loops[c].Reset (block_coloring[c].Range());

        task_manager -> CreateJob
          ( [&] (const TaskInfo & ti) 
            {
              VectorMem<100,TVX> hxmax(maxbs);
              VectorMem<100,TVX> hymax(maxbs);
              
              for (int c : Range(block_coloring))              
                {
                  auto col = block_coloring[c];
                  // SharedLoop2 sl(col.Range());

                  for (auto mynr : loops[c])
                    {
                      size_t i = col[mynr];
                      FlatArray<int> block = (*blocktable)[i];
                      size_t bs = block.Size();
                      if (!bs) continue;
                      
                      FlatVector<TVX> hx = hxmax.Range(0,bs);
                      FlatVector<TVX> hy = hymax.Range(0,bs);
                      
                      for (size_t j = 0; j < bs; j++)
                        {
                          auto jj = block[j];
                          hx(j) = fb(jj) - mat->RowTimesVector (jj, fx);
                        }
                      
                      hy = (invdiag[i]) * hx;
                      fx(block) += hy;
                    }
                }
            });
      }
  }
  

  

  template <class TM, class TV_ROW, class TV_COL>
  void BlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  GSSmoothBack (BaseVector & x, const BaseVector & b,
		int steps) const 
  {
    static Timer timer ("BlockJacobiPrecond::GSSmoothBack");
    RegionTimer reg(timer);
    timer.AddFlops (nze);

    const FlatVector<TVX> fb = b.FV<TVX> (); 
    FlatVector<TVX> fx = x.FV<TVX> ();

    for (int k = 0; k < steps; k++)
      for (int c = block_coloring.Size()-1; c >=0; c--) 
        {
          ParallelForRange
            (color_balance[c], [&] (IntRange r)
             {
               VectorMem<100,TVX> hxmax(maxbs);
               VectorMem<100,TVX> hymax(maxbs);
               
               for (size_t i : block_coloring[c].Range(r))
                 {
                   FlatArray<int> block = (*blocktable)[i];
                   size_t bs = block.Size();
                   if (!bs) continue;
                   
                   FlatVector<TVX> hx = hxmax.Range(0,bs); 
                   FlatVector<TVX> hy = hymax.Range(0,bs); 
                   
                   for (size_t j = 0; j < bs; j++)
                     {
                       auto jj = block[j];
                       hx(j) = fb(jj) - mat->RowTimesVector (jj, fx);
                     }
                   
                   hy = (invdiag[i]) * hx;
                   fx(block) += hy;
                 }
             });
        }
   }








  ///
  template <class TM, class TV>
  BlockJacobiPrecondSymmetric<TM,TV> ::
  BlockJacobiPrecondSymmetric (shared_ptr<const SparseMatrixSymmetric<TM,TV>> amat, 
			       shared_ptr<Table<int>> ablocktable)
    : BaseBlockJacobiPrecond(ablocktable), mat(amat)
  {
    static Timer t("BlockJacobiPrecondSymmetric ctor"); RegionTimer reg(t);    
    cout << IM(3) << "symmetric BlockJacobi Preconditioner 2, constructor called, #blocks = " << blocktable->Size() << endl;

    lowmem = false;
    // lowmem = true;
    
    size_t maxbs = 0;
    size_t n = blocktable->Size();
	
    for (size_t i = 0; i < n; i++)
      if ((*blocktable)[i].Size() > maxbs)
	maxbs = (*blocktable)[i].Size();


    blockstart.SetSize(n);
    blocksize.SetSize(n);
    blockbw.SetSize(n);


    // int alloc = n;
    int memneed[NBLOCKS];
    for (int i = 0; i < NBLOCKS; i++)
      memneed[i] = 0;

    {
      LocalHeap lh (20000 + 5*sizeof(int)*maxbs, "blockjacobi-heap"); 
      Array<int> block_inv(amat->Height());
      block_inv = -1;

      for (size_t i = 0; i < blocktable->Size(); i++)
	{
	  int bs = (*blocktable)[i].Size();
	  
	  if (!bs) continue;
	  
	  blockbw[i] = Reorder ((*blocktable)[i], *mat, block_inv, lh);
	  blocksize[i] = bs;

	  blockstart[i] = memneed[i%NBLOCKS];
	  memneed[i%NBLOCKS] += FlatBandCholeskyFactors<TM>::RequiredMem (bs, blockbw[i]);
	  lh.CleanUp();
	}
    }

    if (!lowmem)
      {
	for (int i = 0; i < NBLOCKS; i++)
	  data[i].SetSize(memneed[i]);

        clock_t prevtime = clock();
        atomic<int> cnt(0);
        
        ParallelFor ( Range(*blocktable),
                      [&] (int i)
                      {
#ifndef __MIC__
                        cnt++;
                        if (clock()-prevtime > 0.1 * CLOCKS_PER_SEC)
                          {
                            {
                              lock_guard<mutex> guard(buildingblockupdate_mutex);
                              cout << IM(3) << "\rBuilding block " << cnt << "/" << blocktable->Size() << flush;
                              prevtime = clock();
                            }
                          }
#endif // __MIC__

                        int bs = (*blocktable)[i].Size();
                        
                        if (!bs) return;
                        int bw = blockbw[i];
                        
                        try
                          {
                            FlatBandCholeskyFactors<TM> inv (bs, bw, &data[i%NBLOCKS][blockstart[i]]);
                            ComputeBlockFactor ((*blocktable)[i], bw, inv);
                          }
                        catch (Exception & e)
                          {
                            cout << IM(1)<<  "block singular !" << endl;
                            (*testout) << "block nr = " << i << endl;
                            (*testout) << "caught: " << e.What() << endl;
                            (*testout) << "entries = " << (*blocktable)[i] << endl;
                            throw;
                          }
                      });
      }
        
    cout << IM(3) << "\rBuilding block " << blocktable->Size() << "/" << blocktable->Size() << endl;
    // cout << "\rBuilt symmetric BlockJacobi Preconditioner" << endl;




    *testout << "block coloring";
    
    int nblocks = blocktable->Size();
    Array<int> coloring(nblocks);
    Array<unsigned int> mask(mat->Width());
    int current_color = 0;
    coloring = -1;
    int colored_blocks = 0;

    
    // symmetric,   rows(block[i]) cap row(block[j]) = 0
    while(colored_blocks<nblocks) 
      {
	mask = 0;
	for (int i = 0; i < nblocks; i++) 
	  {
	    if (coloring[i]>-1) continue;
	    bool is_free = true;
	    for (int d : (*blocktable)[i] )
	      for(auto coupling : mat->GetRowIndices(d)) 
		if(mask[coupling]) 
		  {
		    is_free = false;
		    break;
		  }
	      
	    if(is_free) 
	      {
		coloring[i] = current_color;
		colored_blocks++;
		for (int d : (*blocktable)[i]) 
		  for(auto coupling : mat->GetRowIndices(d)) 
		      mask[coupling] = 1;
	      }
	  }
	current_color++;
      }


    
    TableCreator<int> creator(current_color);
    for ( ; !creator.Done(); creator++)
      for (int i=0; i<nblocks; i++)
          creator.Add(coloring[i],i);
    block_coloring = creator.MoveTable();

    cout << IM(3) <<  " using " << current_color << " colors" << endl;
    
    /*
    *testout << "matrix.h = " << mat->Height() << endl;
    *testout << "smoothing blocks = " << blocktable.Size() << endl;
    
    *testout << "coloring: " << endl << block_coloring << endl;
    *testout << "blocktable: " << endl << blocktable << endl;
    *testout << "matrix: " << endl << mat << endl;
    */

    // calc balancing:

    color_balance.SetSize (block_coloring.Size());

    for (auto c : Range (block_coloring))
      {
        color_balance[c].Calc (block_coloring[c].Size(),
                               [&] (int bi)
                               {
                                 int blocknr = block_coloring[c][bi];
                                 FlatArray<int> block = (*blocktable)[blocknr];
                                 int costs = 0;
                                 for (size_t i=0; i<block.Size(); i++)
                                   costs += mat->GetRowIndices(block[i]).Size();
                                 return costs;
                               });
      }

  

    cout << IM(3) << "\rBlockJacobi Preconditioner built" << endl;
  }




  template <class TM, class TV>
  BlockJacobiPrecondSymmetric<TM,TV> ::
  ~BlockJacobiPrecondSymmetric ()
  {
    ;
  }

  template <class TM, class TV>
  void BlockJacobiPrecondSymmetric<TM,TV> :: 
  ComputeBlockFactor (FlatArray<int> block, int bw, FlatBandCholeskyFactors<TM> & inv) const
  {
    int bs = block.Size();

    ArrayMem<TM, 10000/sizeof(TM)+1> mem(bs*bw);
    FlatSymBandMatrix<TM> blockmat(bs, bw, &mem[0]);


    blockmat = TM(0);
    for (int j = 0; j < bs; j++)
      for (int k = 0; k < bs; k++)
	if (block[j] >= block[k])
	  {
	    if (abs (j-k) < bw)
	      {
		TM val = (*mat)(block[j], block[k]);
		if (j >= k)
		  blockmat(j,k) = val;
		else
		  blockmat(k,j) = Trans (val);
	      }
	  }    

    inv.Factor (blockmat);
  } 





  template <class TM, class TV>
  void BlockJacobiPrecondSymmetric<TM,TV> :: 
  MultAdd (TSCAL s, const BaseVector & x, BaseVector & y) const 
  {
    static Timer timer("BlockJacobiSymmetric::MultAdd");
    RegionTimer reg (timer);

    x.Cumulate();
    y.Cumulate();

    FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy       = y.FV<TVX> ();

    Vector<TVX> hxmax(maxbs);
    Vector<TVX> hymax(maxbs);

    for (size_t i = 0; i < blocktable->Size(); i++)
      {
	int bs = (*blocktable)[i].Size();
	if (!bs) continue;

	FlatVector<TVX> hx = hxmax.Range (0, bs); 
	FlatVector<TVX> hy = hymax.Range (0, bs); 

	for (int j = 0; j < bs; j++)
	  hx(j) = fx((*blocktable)[i][j]);
	
	InvDiag(i).Mult (hx, hy);

	for (int j = 0; j < bs; j++)
	  fy((*blocktable)[i][j]) += s * hy(j);
      }
  }


  template <class TM, class TV>
  void BlockJacobiPrecondSymmetric<TM,TV> :: 
  MultTransAdd (TSCAL s, const BaseVector & x, BaseVector & y) const 
  {
    MultAdd (s, x, y);
  }
  


  template <class TM, class TV>
  void BlockJacobiPrecondSymmetric<TM,TV> :: 
  GSSmooth (BaseVector & x, const BaseVector & b, int steps) const 
  {
    static Timer timer ("BlockJacobiPrecondSymmetric::GSSmooth (parallel)");
    RegionTimer reg(timer);

    FlatVector<TVX> fb = b.FV<TVX> ();
    FlatVector<TVX> fx = x.FV<TVX> ();

    Vector<TVX> fy(fx.Size());

    // y = b - (D L^T) x
    fy = fb;
    for (int j = 0; j < mat->Height(); j++)
      mat->AddRowTransToVector (j, -fx(j), fy);

    
    if (task_manager)
      
      for (int k = 1; k <= steps; k++)
        for (size_t c = 0; c < block_coloring.Size(); c++)
          ParallelFor (color_balance[c], [&] (int bi)
                       {
                         SmoothBlock (block_coloring[c][bi], fx, fy);
                       });
    
    else
      
      for (int k = 1; k <= steps; k++)
        for (size_t i = 0; i < blocktable->Size(); i++)
          SmoothBlock (i, fx, fy);
  }

  
  template <class TM, class TV>
  void BlockJacobiPrecondSymmetric<TM,TV> :: 
  GSSmoothPartial (BaseVector & x, const BaseVector & b,
	    BaseVector & y) const 
  {
    static Timer timer ("BlockJacobiPrecondSymmetric::GSSmooth - partial res");
    RegionTimer reg(timer);

    FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();


    if (task_manager)
      
      for (size_t c = 0; c < block_coloring.Size(); c++)
        ParallelFor (color_balance[c], [&] (int bi)
                     {
                       SmoothBlock (block_coloring[c][bi], fx, fy);
                     });
    
    else

      for (size_t i = 0; i < blocktable->Size(); i++)
        SmoothBlock (i, fx, fy);

  }
  


  ///
  template <class TM, class TV>
  void BlockJacobiPrecondSymmetric<TM,TV> :: 
  GSSmoothResiduum (BaseVector & x, const BaseVector & b,
		    BaseVector & res, int steps) const 
  {
    static Timer timer ("BlockJacobiPrecondSymmetric::GSSmooth - residuum");
    RegionTimer reg(timer);

    // x is 0 on input
    // b is res
    
    res = b;

    // res is partial residual
    for (int k = 1; k <= steps; k++)
      GSSmoothPartial (x, b, res);

    mat->MultAdd1 (-1, x, res);
    // res = b - mat * x;
  }
  
  ///
  template <class TM, class TV>
  void BlockJacobiPrecondSymmetric<TM,TV> :: 
  GSSmoothBack (BaseVector & x, const BaseVector & b,
		int steps) const 
  {
    static Timer timer ("BlockJacobiPrecondSymmetric::SmoothBack");
    RegionTimer reg(timer);

    VVector<TVX> y(x.Size());

    // y = b - (D L^T) x
    y = 1.0*b;
    mat->MultAdd2 (-1, x, y);

    for (int k = 1; k <= steps; k++)
      GSSmoothBackPartial (x, b, y);
  }


  template <class TM, class TV>
  void BlockJacobiPrecondSymmetric<TM,TV> :: 
  GSSmoothBackPartial (BaseVector & x, const BaseVector & b,
		       BaseVector & y) const 
  {
    static Timer timer ("BlockJacobiPrecondSymmetric::GSSmoothBack - partial res");
    RegionTimer reg(timer);

    FlatVector<TVX> fx = x.FV<TVX> ();
    FlatVector<TVX> fy = y.FV<TVX> ();


    if (task_manager)
      
      for (int c = block_coloring.Size()-1; c >= 0; c--)
        ParallelFor (color_balance[c], [&] (int bi)
                     {
                       SmoothBlock (block_coloring[c][bi], fx, /* fb, */ fy);
                     });
    else

      for (int i = blocktable->Size()-1; i >= 0; i--)
        SmoothBlock (i, fx, fy);

  }



  template <class TM, class TV>
  void BlockJacobiPrecondSymmetric<TM,TV> :: 
  SmoothBlock (int i, 
	       FlatVector<TVX> & x,
	       FlatVector<TVX> & y) const
  {
    FlatArray<int> row = (*blocktable)[i];

    int bs = row.Size();
    if (bs == 0) return;

    VectorMem<1000,TVX> di (bs);
    VectorMem<1000,TVX> wi (bs);

    // di = P_i (y - L x)
    for (int j = 0; j < bs; j++)
      di(j) = y(row[j]) - mat->RowTimesVectorNoDiag (row[j], x);
    if (!lowmem)
      InvDiag(i).Mult (di, wi);
    else
      {
	int bw = blockbw[i];
	int bs = (*blocktable)[i].Size();
	ArrayMem<TM, 10000/sizeof(TM)+1> mem(bs*bw);
	FlatBandCholeskyFactors<TM> inv(bs, bw, &mem[0]);

	ComputeBlockFactor ((*blocktable)[i], bw, inv);

	inv.Mult (di, wi);
      }
    // x += P_i w
    // y -= (D L^t) P_i w
    for (int j = 0; j < bs; j++)
      {
	x(row[j]) += wi(j);
	mat->AddRowTransToVector (row[j], -wi(j), y);
      }
  }








  // compiled separately, for testing only
  template class BlockJacobiPrecond<double>;
  template class BlockJacobiPrecond<Complex>;
  template class BlockJacobiPrecond<double, Complex, Complex>;

  template class BlockJacobiPrecondSymmetric<double>;
  template class BlockJacobiPrecondSymmetric<Complex>;
  template class BlockJacobiPrecondSymmetric<double, Complex>;

#if MAX_SYS_DIM >= 1
  template class BlockJacobiPrecond<Mat<1,1,double> >;
  template class BlockJacobiPrecond<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class BlockJacobiPrecond<Mat<2,2,double> >;
  template class BlockJacobiPrecond<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class BlockJacobiPrecond<Mat<3,3,double> >;
  template class BlockJacobiPrecond<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class BlockJacobiPrecond<Mat<4,4,double> >;
  template class BlockJacobiPrecond<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class BlockJacobiPrecond<Mat<5,5,double> >;
  template class BlockJacobiPrecond<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class BlockJacobiPrecond<Mat<6,6,double> >;
  template class BlockJacobiPrecond<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class BlockJacobiPrecond<Mat<7,7,double> >;
  template class BlockJacobiPrecond<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class BlockJacobiPrecond<Mat<8,8,double> >;
  template class BlockJacobiPrecond<Mat<8,8,Complex> >;
#endif
  
#if MAX_SYS_DIM >= 1
  template class BlockJacobiPrecondSymmetric<Mat<1,1,double> >;
  template class BlockJacobiPrecondSymmetric<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class BlockJacobiPrecondSymmetric<Mat<2,2,double> >;
  template class BlockJacobiPrecondSymmetric<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class BlockJacobiPrecondSymmetric<Mat<3,3,double> >;
  template class BlockJacobiPrecondSymmetric<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class BlockJacobiPrecondSymmetric<Mat<4,4,double> >;
  template class BlockJacobiPrecondSymmetric<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class BlockJacobiPrecondSymmetric<Mat<5,5,double> >;
  template class BlockJacobiPrecondSymmetric<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class BlockJacobiPrecondSymmetric<Mat<6,6,double> >;
  template class BlockJacobiPrecondSymmetric<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class BlockJacobiPrecondSymmetric<Mat<7,7,double> >;
  template class BlockJacobiPrecondSymmetric<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class BlockJacobiPrecondSymmetric<Mat<8,8,double> >;
  template class BlockJacobiPrecondSymmetric<Mat<8,8,Complex> >;
#endif

  
#ifdef CACHEBLOCKSIZE
  template class BlockJacobiPrecond<double, Vec<CACHEBLOCKSIZE>, Vec<CACHEBLOCKSIZE> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<CACHEBLOCKSIZE> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class BlockJacobiPrecond<double, Vec<2,double>, Vec<2,double> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class BlockJacobiPrecond<double, Vec<3,double>, Vec<3,double> >;
  template class BlockJacobiPrecond<double, Vec<4,double>, Vec<4,double> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class BlockJacobiPrecond<double, Vec<5,double>, Vec<5,double> >;
  template class BlockJacobiPrecond<double, Vec<6,double>, Vec<6,double> >;
  template class BlockJacobiPrecond<double, Vec<7,double>, Vec<7,double> >;
  template class BlockJacobiPrecond<double, Vec<8,double>, Vec<8,double> >;
  template class BlockJacobiPrecond<double, Vec<9,double>, Vec<9,double> >;
  template class BlockJacobiPrecond<double, Vec<10,double>, Vec<10,double> >;
  template class BlockJacobiPrecond<double, Vec<11,double>, Vec<11,double> >;
  template class BlockJacobiPrecond<double, Vec<12,double>, Vec<12,double> >;
  template class BlockJacobiPrecond<double, Vec<13,double>, Vec<13,double> >;
  template class BlockJacobiPrecond<double, Vec<14,double>, Vec<14,double> >;
  template class BlockJacobiPrecond<double, Vec<15,double>, Vec<15,double> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class BlockJacobiPrecond<Complex, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class BlockJacobiPrecond<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class BlockJacobiPrecond<Complex, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class BlockJacobiPrecond<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class BlockJacobiPrecond<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class BlockJacobiPrecond<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class BlockJacobiPrecond<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class BlockJacobiPrecond<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class BlockJacobiPrecond<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class BlockJacobiPrecond<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class BlockJacobiPrecond<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class BlockJacobiPrecond<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class BlockJacobiPrecond<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class BlockJacobiPrecond<Complex, Vec<15,Complex>, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class BlockJacobiPrecond<double, Vec<2,Complex>, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class BlockJacobiPrecond<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class BlockJacobiPrecond<double, Vec<4,Complex>, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class BlockJacobiPrecond<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class BlockJacobiPrecond<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class BlockJacobiPrecond<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class BlockJacobiPrecond<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class BlockJacobiPrecond<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class BlockJacobiPrecond<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class BlockJacobiPrecond<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class BlockJacobiPrecond<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class BlockJacobiPrecond<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class BlockJacobiPrecond<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class BlockJacobiPrecond<double, Vec<15,Complex>, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class BlockJacobiPrecondSymmetric<double, Vec<2,double> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class BlockJacobiPrecondSymmetric<double, Vec<3,double> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<4,double> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class BlockJacobiPrecondSymmetric<double, Vec<5,double> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<6,double> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<7,double> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<8,double> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<9,double> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<10,double> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<11,double> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<12,double> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<13,double> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<14,double> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<15,double> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class BlockJacobiPrecondSymmetric<Complex, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class BlockJacobiPrecondSymmetric<Complex, Vec<3,Complex> >;
  template class BlockJacobiPrecondSymmetric<Complex, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class BlockJacobiPrecondSymmetric<Complex, Vec<5,Complex> >;
  template class BlockJacobiPrecondSymmetric<Complex, Vec<6,Complex> >;
  template class BlockJacobiPrecondSymmetric<Complex, Vec<7,Complex> >;
  template class BlockJacobiPrecondSymmetric<Complex, Vec<8,Complex> >;
  template class BlockJacobiPrecondSymmetric<Complex, Vec<9,Complex> >;
  template class BlockJacobiPrecondSymmetric<Complex, Vec<10,Complex> >;
  template class BlockJacobiPrecondSymmetric<Complex, Vec<11,Complex> >;
  template class BlockJacobiPrecondSymmetric<Complex, Vec<12,Complex> >;
  template class BlockJacobiPrecondSymmetric<Complex, Vec<13,Complex> >;
  template class BlockJacobiPrecondSymmetric<Complex, Vec<14,Complex> >;
  template class BlockJacobiPrecondSymmetric<Complex, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class BlockJacobiPrecondSymmetric<double, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class BlockJacobiPrecondSymmetric<double, Vec<3,Complex> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class BlockJacobiPrecondSymmetric<double, Vec<5,Complex> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<6,Complex> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<7,Complex> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<8,Complex> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<9,Complex> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<10,Complex> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<11,Complex> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<12,Complex> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<13,Complex> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<14,Complex> >;
  template class BlockJacobiPrecondSymmetric<double, Vec<15,Complex> >;
#endif
  

}
