/*********************************************************************/
/* File:   parallel_matrices.cpp                                     */
/* Author: Joachim Schoeberl                                         */
/* Date:   June 2011                                                 */
/*********************************************************************/

#ifdef PARALLEL
 
#include <la.hpp>
#include <comp.hpp>
#include <parallelngs.hpp>


namespace ngla
{
  using namespace ngparallel;
  using namespace ngcomp;

  template <typename TM>
  MasterInverse<TM> :: MasterInverse (const SparseMatrixTM<TM> & mat, 
				      const BitArray * subset, 
				      const ParallelDofs * hpardofs)

    : loc2glob(MyMPI_GetNTasks (pardofs -> GetCommunicator())),
      pardofs(hpardofs)
  {
    inv = NULL;
    
    int id = MyMPI_GetId (pardofs -> GetCommunicator());
    int ntasks = MyMPI_GetNTasks (pardofs -> GetCommunicator());


    // consistent enumeration
    
    int ndof = pardofs->GetNDof();
    
    Array<int> global_nums(ndof);
    global_nums = -1;
    int num_master_dofs = 0;
    for (int i = 0; i < ndof; i++)
      if (pardofs -> IsMasterDof (i) && (!subset || (subset && subset->Test(i))))
	global_nums[i] = num_master_dofs++;
    
    Array<int> first_master_dof(ntasks);
    MPI_Allgather (&num_master_dofs, 1, MPI_INT, 
		   &first_master_dof[0], 1, MPI_INT, 
		   pardofs -> GetCommunicator());
    
    int num_glob_dofs = 0;
    for (int i = 0; i < ntasks; i++)
      {
	int cur = first_master_dof[i];
	first_master_dof[i] = num_glob_dofs;
	num_glob_dofs += cur;
      }
    
    for (int i = 0; i < ndof; i++)
      if (global_nums[i] != -1)
	global_nums[i] += first_master_dof[id];

    ScatterDofData (global_nums, *pardofs);

    /*
    cout << "TESTING" << endl;
    for (int i = 0; i < ndof; i++)
      *testout << "dof" << i << " = " << pardofs->GetDofNodes()[i] << " free = " << subset->Test(i) << endl;
    const MeshAccess & ma = pardofs->GetMeshAccess();

    for (int dest = 0; dest < ntasks; dest++)
      {
	*testout << "shared edges with proc " << dest << endl;
	*testout << "exdofs " << endl << pardofs->GetExchangeDofs (dest) << endl;
	Array<int> procs;
	for (int i = 0; i < ma.GetNEdges(); i++)
	  {
	    ma.GetDistantProcs (Node(NT_EDGE, i), procs);
	    if (procs.Contains (dest))
	      {
		int v1, v2;
		ma.GetEdgePNums (i, v1, v2);
		v1 = ma.GetGlobalNodeNum (Node(NT_VERTEX, v1));
		v2 = ma.GetGlobalNodeNum (Node(NT_VERTEX, v2));
		*testout << "E" << i << ": " << v1 << "-" << v2 << endl;
	      }
	  }
      }

    for (int i = 0; i < ma.GetNEdges(); i++)
      {
	*testout << "local edge " << i << endl;
	int v1, v2;
	ma.GetEdgePNums (i, v1, v2);
	*testout << "loc pnts = " << v1 << "-" << v2 << ", glob pnts = " 
		 << ma.GetGlobalNodeNum (Node(NT_VERTEX, v1)) << "-" 
		 << ma.GetGlobalNodeNum (Node(NT_VERTEX, v2)) << endl;
	Array<int> procs;
	ma.GetDistantProcs (Node(NT_EDGE, i), procs);
	*testout << "dist procs = " << procs << endl;
      }

    for (int i = 0; i < ndof; i++)
      if ( (global_nums[i] == -1) && (!subset || (subset && subset->Test(i))))
	{
	  Node node = pardofs->GetDofNodes()[i];
	  cerr << "global enumeration problem, dof = " << i << ", node = " << node << endl;
	  *testout << "global enumeration problem, dof = " << i << ", node = " << node << endl;
	}
    */

    
    if (id != 0)
      {
	// const MeshAccess & ma = nodaldofs -> GetMeshAccess();

	int ndof = pardofs->GetNDof();

	Array<int> rows, cols;
	Array<TM> vals;

	for (int row = 0; row < ndof; row++)
	  if (!subset || subset->Test(row))
	    select.Append (row);
	
	for (int row = 0; row < mat.Height(); row++)
	  if (!subset || subset->Test(row))
	    {
	      FlatArray<int> rcols = mat.GetRowIndices(row);
	      FlatVector<TM> rvals = mat.GetRowValues(row);

	      for (int j = 0; j < rcols.Size(); j++)
		if (!subset || subset->Test(rcols[j]))
		  {
		    rows.Append (global_nums[row]);
		    cols.Append (global_nums[rcols[j]]);
		    vals.Append (rvals[j]);
		  }
	    }

	MyMPI_Send (rows, 0);
	MyMPI_Send (cols, 0);
	MyMPI_Send (vals, 0);
	MyMPI_Send (global_nums, 0);

#ifdef USE_MUMPS
	if (mat.GetInverseType() == MUMPS)
	  {
	    SparseMatrixSymmetric<TM> dummy_matrix (1,1);
	    inv = new MumpsInverse<TM> (dummy_matrix, 0, 0, true);
	  }
	else
#endif 
	  inv = NULL;
      }

    else
      {
	cout << "create masterinverse" << endl;

	Array<int> rows, cols;
	Array<TM> vals;
	HashTable<INT<1>, int> ht_globdofs(100000);
	// int num_globdofs = 0; 

	for (int src = 1; src < ntasks; src++)
	  {
	    Array<int> hrows, hcols;
	    Array<TM> hvals;
	    Array<int> hglobid;

	    MyMPI_Recv (hrows, src);
	    MyMPI_Recv (hcols, src);
	    MyMPI_Recv (hvals, src);
	    MyMPI_Recv (hglobid, src);

	    for (int i = 0; i < hglobid.Size(); i ++)
	      {
		if (hglobid[i] == -1) continue;
		loc2glob.Add (src, hglobid[i]);
	      }

	    for (int i = 0; i < hrows.Size(); i++)
	      {
		if (hrows[i] < 0 || hcols[i] < 0) cerr << "Illegal (row/col)" << endl;
		rows.Append (hrows[i]);
		cols.Append (hcols[i]);
		vals.Append (hvals[i]);
	      }

	    cout << "\rmaster: got data from " << src << flush;
	  }
	cout << endl;


	cout << "now build graph" << endl;

	// build matrix
	DynamicTable<int> graph(num_glob_dofs);
	cout << "n = " << num_glob_dofs << endl;
	for (int i = 0; i < rows.Size(); i++)
	  {
	    int r = rows[i], c = cols[i];
	    if (r < c) swap (r, c);
	    graph.AddUnique (r, c);
	  }

	// *testout << "graphi = " << endl << graph << endl;

	Array<int> els_per_row(num_glob_dofs);
	for (int i = 0; i < num_glob_dofs; i++)
	  els_per_row[i] = graph[i].Size();

	cout << "now build matrix" << endl;

	SparseMatrixSymmetric<TM> matrix(els_per_row);

	for (int i = 0; i < rows.Size(); i++)
	  {
	    int r = rows[i], c = cols[i];
	    if (r < c) swap (r, c);
	    matrix.CreatePosition(r, c);
	  }
	matrix.AsVector() = 0.0;

	for (int i = 0; i < rows.Size(); i++)
	  {
	    int r = rows[i], c = cols[i];
	    if (r < c) swap (r, c);
	    matrix(r,c) += vals[i];
	  }

	cout << "have matrix, now invert" << endl;

	matrix.SetInverseType (mat.GetInverseType());
	inv = matrix.InverseMatrix ();
      }

    // MPI_Barrier (ngs_comm);
  }

  template <typename TM>
  MasterInverse<TM> :: ~MasterInverse ()
  {
    delete inv;
  }

  template <typename TM>
  void MasterInverse<TM> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    typedef typename mat_traits<TM>::TV_ROW TV;

    int id = MyMPI_GetId (pardofs -> GetCommunicator());
    int ntasks = MyMPI_GetNTasks (pardofs -> GetCommunicator());

    bool is_x_cum = (dynamic_cast<const ParallelBaseVector&> (x) . Status() == CUMULATED);
    x.Distribute();
    y.Cumulate();

    if (id > 0)
      {
	FlatVector<TV> fx = x.FV<TV> ();
	FlatVector<TV> fy = y.FV<TV> ();

	Array<TV> lx (select.Size());
	for (int i = 0; i < select.Size(); i++)
	  lx[i] = fx(select[i]);
	
	MPI_Request request = MyMPI_ISend (lx, 0, MPI_TAG_SOLVE);
	MPI_Request_free (&request);

	// only for MUMPS:
	if (inv)
	  y = (*inv) * x;
	
	request = MyMPI_IRecv (lx, 0, MPI_TAG_SOLVE);
	MPI_Wait (&request, MPI_STATUS_IGNORE);

	for (int i = 0; i < select.Size(); i++)
	  fy(select[i]) += s * lx[i];
      }

    else
      {
	VVector<TV> hx(inv->Height());
	VVector<TV> hy(inv->Height());
	hx = 0.0;

	Array<int> sizes(ntasks);
	for (int i = 0; i < ntasks; i++)
	  sizes[i] = loc2glob[i].Size();

	Table<TV> exdata(sizes);


	for (int src = 1; src < ntasks; src++)
	  {
	    FlatArray<int> selecti = loc2glob[src];

	    Array<TV> lx(selecti.Size());
	    MyMPI_Recv (lx, src);

	    for (int i = 0; i < selecti.Size(); i++)
	      hx(selecti[i]) += lx[i];
	  }

	hy = (*inv) * hx;

	Array<MPI_Request> requ;
	for (int src = 1; src < ntasks; src++)
	  {
	    FlatArray<int> selecti = loc2glob[src];
	    for (int j = 0; j < selecti.Size(); j++)
	      exdata[src][j] = hy(selecti[j]);
	    requ.Append (MyMPI_ISend (exdata[src], src, MPI_TAG_SOLVE));
	  }
	MyMPI_WaitAll (requ);
      }

    if (is_x_cum)
      dynamic_cast<const ParallelBaseVector&> (x) . Cumulate(); // AllReduce(&hoprocs);

  }



  ParallelMatrix :: ~ParallelMatrix ()
  {
    delete &mat;
  }

  void ParallelMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    x.Cumulate();
    y.Distribute();
    mat.MultAdd (s, x, y);
  }

  void ParallelMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    x.Cumulate();
    y.Distribute();
    mat.MultTransAdd (s, x, y);
  }

  BaseMatrix * ParallelMatrix :: CreateMatrix () const
  {
    return new ParallelMatrix (mat.CreateMatrix(), &pardofs);
  }

  BaseVector * ParallelMatrix :: CreateVector () const
  {
    cerr << "ParallelMatrix::CreateVector not implemented" << endl;
    return NULL;
  }

  ostream & ParallelMatrix :: Print (ostream & ost) const
  {
    ost << mat;
    return ost;
  }

  int ParallelMatrix :: VHeight() const
  {
    return mat.VHeight();
  }

  int ParallelMatrix :: VWidth() const
  {
    return mat.VWidth();
  }


  BaseMatrix * ParallelMatrix::InverseMatrix (const BitArray * subset) const
  {
    const SparseMatrixTM<double> * dmat = dynamic_cast<const SparseMatrixTM<double>*> (&mat);
    const SparseMatrixTM<Complex> * cmat = dynamic_cast<const SparseMatrixTM<Complex>*> (&mat);

#ifdef USE_MUMPS
    if (mat.GetInverseType() == MUMPS)
      {
	if (dmat) return new ParallelMumpsInverse<double> (*dmat, subset, NULL, &pardofs);
	if (cmat) return new ParallelMumpsInverse<Complex> (*cmat, subset, NULL, &pardofs);
      }
    else 
#endif
      {
	if (dmat) return new MasterInverse<double> (*dmat, subset, &pardofs);
	if (cmat) return new MasterInverse<Complex> (*cmat, subset, &pardofs);
      }
 
    cerr << "ParallelMatrix::Inverse(BitArray) not avail, typeid(mat) = " << typeid(mat).name() << endl;
    return NULL;
  }

  BaseMatrix * ParallelMatrix::InverseMatrix (const Array<int> * clusters) const
  {
    cerr << "ParallelMatrix::Inverse(ARRAY) not avail" << endl;
    return NULL;
  }

  INVERSETYPE ParallelMatrix::SetInverseType (INVERSETYPE ainversetype) const
  {
    return mat.SetInverseType (ainversetype);
  }

  INVERSETYPE ParallelMatrix::SetInverseType (string ainversetype) const
  {
    return mat.SetInverseType (ainversetype);
  }
  
  INVERSETYPE ParallelMatrix::GetInverseType () const
  {
    return mat.GetInverseType ();
  }









  template class MasterInverse<double>;
  template class MasterInverse<Complex>;

#if MAX_SYS_DIM >= 1
  template class MasterInverse<Mat<1,1,double> >;
  template class MasterInverse<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class MasterInverse<Mat<2,2,double> >;
  template class MasterInverse<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class MasterInverse<Mat<3,3,double> >;
  template class MasterInverse<Mat<3,3,Complex> >;
#endif



}

#endif
