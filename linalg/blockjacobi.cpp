/**************************************************************************/
/* File:   blockjacobi.cpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   14. Aug. 2002                                                  */
/**************************************************************************/


#include <la.hpp>


namespace ngla
{
  using namespace ngla;

  
  BaseBlockJacobiPrecond :: 
  BaseBlockJacobiPrecond (Table<int> & ablocktable)
    : blocktable(ablocktable)
  {
    maxbs = 0;
    for (int i = 0; i < blocktable.Size(); i++)
      if (blocktable[i].Size() > maxbs)
	maxbs = blocktable[i].Size();

  }


  BaseBlockJacobiPrecond :: 
  ~BaseBlockJacobiPrecond ()
  {
    delete &blocktable;
  }


  int BaseBlockJacobiPrecond ::
  Reorder (FlatArray<int> block, const MatrixGraph & graph,
	   FlatArray<int> block_inv,
	   LocalHeap & lh)
  {
    try
      {
	// a cheap reordering algorithm:
	
	// bool print = 0;
	
	int n = block.Size();
	
	void * heapp = lh.GetPointer();


	FlatArray<int> reorder(n, lh), newnum(n, lh), dist(n, lh), cluster0(n,lh);
	// FlatArray<int> block_inv(graph.Size(), lh);
	
	// for(int i=0; i<graph.Size(); i++)
	// block_inv[i] = -1;


	for (int i = 0; i < n; i++)
	  if (block[i] >= 0 && block[i] < graph.Size())
	    {
	      if (block_inv[block[i]] != -1)
		{
		  cout << "block has double elements " << endl;
		  cout << block_inv[block[i]] << " and " << i << endl;
		  cout << block << endl;
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
	    for (int j = 0; j < n; j++)
	      {
		FlatArray<const int> row = graph.GetRowIndices(block[j]);
		for (int k = 0; k < row.Size(); k++)
		  {
		    int kk = block_inv[row[k]];
		    if (kk >= 0 && kk < n)
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

	int cnt = 0;
	for (int i = 0; i < n; i++)
	  if (cluster0[i])
	    {
	      newnum[cnt] = block[i];
	      cnt++;
	    }

	if (cnt < n)
	  {
	    // seperated clusters
	    int cnt2 = cnt;
	    for (int i = 0; i < n; i++)
	      if (!cluster0[i])
		newnum[cnt2++] = block[i];
	    
	    for (int i = 0; i < n; i++)
	      block[i] = newnum[i];

	    for (int i = 0; i < n; i++)
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

	int pstart = 0;
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
		for (int i = 0; i < n; i++)
		  {
		    FlatArray<const int> row = graph.GetRowIndices(block[i]);

		    for (int j = 0; j < row.Size(); j++)
		      {
			int jj = block_inv[row[j]];

			if (jj >= 0 && jj < n)
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
	    for (int i = 0; i < n; i++)
	      if (dist[i] > maxval)
		{
		  maxval = dist[i];
		  pstart = i;
		}
	    if (maxval > n)
	      {
		cerr << "Blockjacobi, reorder: separated block" << endl;
		cout << "block: " << block << endl;

		(*testout) << "Blockjacobi, reorder: separated block" << endl;
		(*testout) << "block: " << block << endl;
	      }	      
	  }
	

	cnt = 0;
	for (int i = 0; i < n; i++)
	  for (int j = 0; j < n; j++)
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
	for (int i = 0; i < n; i++)
	  newnum[reorder[i]] = i;
	for (int j = 0; j < n; j++)
	  {
	    FlatArray<const int> row = graph.GetRowIndices(block[j]);
	    for (int k = 0; k < row.Size(); k++)
	      {
		int kk = block_inv[row[k]];
		if (kk >= 0 && kk < n)
		  if (block[kk] == row[k])
		    {
		      int inv_j = newnum[j];
		      int inv_k = newnum[kk];
		      bw = max2(bw, abs(inv_j-inv_k)+1);
		    }
	      }
	  }


	for (int i = 0; i < n; i++)
	  newnum[i] = block[reorder[i]];
	for (int i = 0; i < n; i++)
	  block[i] = newnum[i];


	lh.CleanUp(heapp);

	for (int i = 0; i < n; i++)
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
  BlockJacobiPrecond (const SparseMatrix<TM,TV_ROW,TV_COL> & amat, 
		      Table<int> & ablocktable)
    : BaseBlockJacobiPrecond(ablocktable), mat(amat), 
      invdiag(ablocktable.Size())
  { 
    cout << "BlockJacobi Preconditioner, constructor called, #blocks = " << blocktable.Size() << endl;


    clock_t prevtime = clock();
    for (int i = 0; i < blocktable.Size(); i++)
      {
        if (clock()-prevtime > 0.1 * CLOCKS_PER_SEC)
          {
            cout << "\rBuilding block " << i << flush;
            prevtime = clock();
          }

        /*
    for (int i = 0; i < blocktable.Size(); i++)
      {
	cout << "\rBuilding block " << i ;
        */


	int bs = blocktable[i].Size();
	if (!bs) 
	  {
	    invdiag[i] = 0;
	    continue;
	  }
	
	Matrix<TM> blockmat(bs);
	invdiag[i] = new Matrix<TM> (bs);
	
        
	for (int j = 0; j < bs; j++)
	  for (int k = 0; k < bs; k++)
	    blockmat(j,k) = mat(blocktable[i][j], blocktable[i][k]);
	
// 	(*testout) << "juhu, block " << i << " has L2norm " << L2Norm(blockmat) << endl;	
// 	(*testout) << "block = " << blocktable[i] << endl;
// 	(*testout) << "blockmat = " << endl << blockmat << endl;
	
        CalcInverse (blockmat, *invdiag[i]);
	
// 	(*testout) << "inv = " << endl << *invdiag[i] << endl;
// 	(*testout) << "prod = " << endl << (blockmat * *invdiag[i]) << endl;
      }
    cout << "\rBlockJacobi Preconditioner built" << endl;
  }

  ///
  template <class TM, class TV_ROW, class TV_COL>
  BlockJacobiPrecond<TM, TV_ROW, TV_COL> ::
  ~BlockJacobiPrecond () 
  {
    for (int i = 0; i < invdiag.Size(); i++)
      delete invdiag[i];
  }



#ifdef SYMCHOLESKY

  ///
  template <class TM, class TV>
  BlockJacobiPrecondSymmetric<TM,TV> ::
  BlockJacobiPrecondSymmetric (const SparseMatrixSymmetric<TM,TV> & amat, 
			       Table<int> & ablocktable)
    : BaseBlockJacobiPrecond(ablocktable), mat(amat), 
      invdiag(ablocktable.Size())
  { 
    cout << "symmetric BlockJacobi Preconditioner, constructor called, #blocks = " << blocktable.Size() << endl;

    int sumnn = 0;
    
    for (int i = 0; i < blocktable.Size(); i++)
      {
	cout << "\rBuilding block " << i ;
	int bs = blocktable[i].Size();
	
	if (!bs) 
	  {
	    invdiag[i] = 0;
	    continue;
	  }
	
	sumnn += bs*bs;

	//	int bw = Reorder (blocktable[i], mat);
	
	Matrix<TM> blockmat(bs);
	//	invdiag[i] = new Matrix<TM> (bs);
	
	
	for (int j = 0; j < bs; j++)
	  for (int k = 0; k < bs; k++)
	    if (blocktable[i][j] >= blocktable[i][k])
	      {
		blockmat(j,k) = 
		  mat(blocktable[i][j], blocktable[i][k]);
		if (j != k)
		  blockmat(k,j) = Trans (blockmat(j,k));
	      }
	
	try
	  {
	    invdiag[i] = new CholeskyFactors<TM> (blockmat);
	  }
	catch (Exception & e)
	  {
	    cout << "block singular !" << endl;
	    (*testout) << "caught: " << e.What() << endl;
	    (*testout) << "entries = " << blocktable[i] << endl;
	    (*testout) << "mat = " << endl << blockmat << endl;
	  }
	
	// invdiag[i]->Print(*testout);
      }
    cout << "\rBuilt symmetric BlockJacobi Preconditioner" << endl;
  }

  template <class TM, class TV>
  BlockJacobiPrecondSymmetric<TM,TV> ::
  BlockJacobiPrecondSymmetric (const SparseMatrixSymmetric<TM,TV> & amat, 
			       const FlatVector<TVX> & constraint,
			       Table<int> & ablocktable)
    : BaseBlockJacobiPrecond(ablocktable), mat(amat), 
      invdiag(ablocktable.Size())
  { 
    cout << "symmetric BlockJacobi Preconditioner, constructor called, #blocks = " << blocktable.Size() << endl;

    int sumnn = 0;

    for (int i = 0; i < blocktable.Size(); i++)
      {
	cout << "\rBuilding block " << i ;
	int bs = blocktable[i].Size();

	if (!bs) 
	  {
	    invdiag[i] = 0;
	    continue;
	  }

	sumnn += bs*bs;

	//	int bw = Reorder (blocktable[i], mat);

	Matrix<TM> blockmat(bs);
	//	invdiag[i] = new Matrix<TM> (bs);

	for (int j = 0; j < bs; j++)
	  for (int k = 0; k < bs; k++)
	    if (blocktable[i][j] >= blocktable[i][k])
	      {
		blockmat(j,k) = 
		  mat(blocktable[i][j], blocktable[i][k]);
		blockmat(k,j) = Trans (blockmat(j,k));
	      }
	
	for (int j = 0; j < bs; j++)
	  for (int k = 0; k < bs; k++)
	    blockmat(j,k) -= 1e8 * 
	      constraint(blocktable[i][j]) *
	      Trans (constraint(blocktable[i][k]));

	try
	  {
	    invdiag[i] = new CholeskyFactors<TM> (blockmat);
	  }
	catch (Exception & e)
	  {
	    cout << "block singular !" << endl;
	    (*testout) << "caught: " << e.What() << endl;
	    (*testout) << "entries = " << blocktable[i] << endl;
	    (*testout) << "mat = " << endl << blockmat << endl;
	  }
	// invdiag[i]->Print(*testout);
      }
    cout << "\rBuilt symmetric BlockJacobi Preconditioner" << endl;

  }

  ///
  template <class TM,TV>
  BlockJacobiPrecondSymmetric<TM,TV> ::
  ~BlockJacobiPrecondSymmetric () 
  {
    for (int i = 0; i < invdiag.Size(); i++)
      delete invdiag[i];
  }


#else


  ///
  template <class TM, class TV>
  BlockJacobiPrecondSymmetric<TM,TV> ::
  BlockJacobiPrecondSymmetric (const SparseMatrixSymmetric<TM,TV> & amat, 
			       Table<int> & ablocktable)
    : BaseBlockJacobiPrecond(ablocktable), mat(amat)
  { 
    cout << "symmetric BlockJacobi Preconditioner 2, constructor called, #blocks = " << blocktable.Size() << endl;

    lowmem = false;
    // lowmem = true;
    
    // int i;
    // int sumnn = 0;
    int maxbs = 0;
    
    int n = blocktable.Size();
	
    for (int i = 0; i < n; i++)
      if (blocktable[i].Size() > maxbs)
	maxbs = blocktable[i].Size();


    blockstart.Alloc(n);
    blocksize.Alloc(n);
    blockbw.Alloc(n);

    blockstart.SetName ("BlockJacobi, Table 1");
    blocksize.SetName ("BlockJacobi, Table 2");
    blockbw.SetName ("BlockJacobi, Table 3");

    for (int i = 0; i < NBLOCKS; i++)
      data[i].SetName ("BlockJacobi");

    // int alloc = n;
    int memneed[NBLOCKS];
    for (int i = 0; i < NBLOCKS; i++)
      memneed[i] = 0;
    /*
    int starti[NBLOCKS];
    for (int i = 0; i < NBLOCKS; i++)
      starti[i] = memneed[i] = 0;
    */

    {
      LocalHeap lh (20000 + 5*sizeof(int)*maxbs, "blockjacobi-heap"); 
      Array<int> block_inv(amat.Height());
      block_inv = -1;

      for (int i = 0; i < blocktable.Size(); i++)
	{
	  int bs = blocktable[i].Size();
	  
	  if (!bs) continue;
	  
	  blockbw[i] = Reorder (blocktable[i], mat, block_inv, lh);
	  blocksize[i] = bs;

	  blockstart[i] = memneed[i%NBLOCKS];
	  memneed[i%NBLOCKS] += FlatBandCholeskyFactors<TM>::RequiredMem (bs, blockbw[i]);
	  lh.CleanUp();
	}
    }

    /* 
       int tot_mem =0; 
       for(int i=0;i<NBLOCKS;i++) 
       tot_mem += memneed[i]; 
       *testout << " ******* MEMORY BlockJacobi " << endl; 
       *testout << " Doubles needed for Block-Jacobi " << double(tot_mem) << endl; 
       *testout << " Memory needed for Block-Jacobi " << double(tot_mem) * sizeof(double) * 1.e-6 << " MB " <<  endl ; 
       *testout << " NZE of Amat " << double(amat.NZE()) << endl; 
       *testout << " Memory for Amat " << double(amat.NZE())*(sizeof(int)+sizeof(double)) *1.e-6 << " MB " << endl; 
       */
       
    if (!lowmem)
      {
	for (int i = 0; i < NBLOCKS; i++)
	  data[i].Alloc (memneed[i]);

        clock_t prevtime = clock();

	
	for (int i = 0; i < blocktable.Size(); i++)
	  {
            if (clock()-prevtime > 0.1 * CLOCKS_PER_SEC)
              {
                cout << "\rBuilding block " << i << flush;
                prevtime = clock();
              }

	    int bs = blocktable[i].Size();
	    
	    if (!bs) continue;
	    
	    // blockstart[i] = starti[i%NBLOCKS];
	    // if (blockstart[i] == starti[i%NBLOCKS]) cout << "g"; else cout << "b";
	    
	    int bw = blockbw[i];
	    // int need = FlatBandCholeskyFactors<TM>::RequiredMem (bs, bw);
	    
	    /*
	      if (starti + need > alloc)
	      {
	      alloc = int (1.5 * (starti+need) + 10);
	      data.ReAlloc (alloc);
	      }
	    */
	    try
	      {
		// invdiag[i] = new BandCholeskyFactors<TM> (blockmat);
		// FlatBandCholeskyFactors<TM> inv (bs, bw, &data[i%NBLOCKS][starti[i%NBLOCKS]]);
		FlatBandCholeskyFactors<TM> inv (bs, bw, &data[i%NBLOCKS][blockstart[i]]);
		// (*testout) << "factor block " << i << endl;

		ComputeBlockFactor (blocktable[i], bw, inv);

		// inv.Print (*testout);
		// inv.Factor (blockmat);
	      }
	    catch (Exception & e)
	      {
		cout << "block singular !" << endl;
		(*testout) << "block nr = " << i << endl;
		(*testout) << "caught: " << e.What() << endl;
		(*testout) << "entries = " << blocktable[i] << endl;
		/*
		  (*testout) << "mat = " << endl;
		  blockmat.Print(*testout);
		  (*testout) << "diag = " << endl;
		  pfor (int l = 0; l < bs; l++)
		  (*testout) << l << ": " << blockmat(l,l) << endl;
		*/
		throw;
	      }
	    
	    // starti[i%NBLOCKS] += need;
	  }
      }

    cout << "\rBuilt symmetric BlockJacobi Preconditioner" << endl;
  }




  template <class TM, class TV>
  BlockJacobiPrecondSymmetric<TM,TV> ::
  BlockJacobiPrecondSymmetric (const SparseMatrixSymmetric<TM,TV> & amat, 
			       const FlatVector<TVX> & constraint,
			       Table<int> & ablocktable)
    : BaseBlockJacobiPrecond(ablocktable), mat(amat)
  {
    throw Exception ("BlockJacPrecondSym with constraints not available for banded blocks, please define SYMCHOLESKY in blocjacobi.hpp");
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
		TM val = mat(block[j], block[k]);
		if (j >= k)
		  blockmat(j,k) = val;
		else
		  blockmat(k,j) = Trans (val);
	      }
	  }    

    inv.Factor (blockmat);

//       (*testout) << "block = " << block << endl
//       << "mat = " << blockmat << endl
//       << "inv = " << endl << inv << endl;



    /*
      Matrix<TM> mat2(bs);
      mat2 = TM(0);
      for (int j = 0; j < bs; j++)
      for (int k = 0; k < bs; k++)
      if (block[j] >= block[k])
      {
      if (abs (j-k) < bw)
      {
      TM val = mat(block[j], block[k]);
      mat2(j,k) = val;
      mat2(k,j) = Trans (val);
      }
      }    
    
      CholeskyFactors<TM> inv2(mat2);
      (*testout) << "mat2 = " << endl << mat2 << endl;
      (*testout) << "inv2 = " << endl;
      inv2.Print (*testout);
      (*testout) << endl;
    */
  } 


#endif // symcholesky




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
