/*********************************************************************/
/* File:   parallel_matrices.cpp                                     */
/* Author: Joachim Schoeberl                                         */
/* Date:   June 2011                                                 */
/*********************************************************************/

#ifdef PARALLEL
 
#include <la.hpp>
#include "../linalg/mumpsinverse.hpp"
#include <parallelngs.hpp>


namespace ngla
{
  
  template <typename TM> AutoVector MasterInverse<TM> :: CreateVector () const
  { return make_shared<ParallelVVector<double>> (paralleldofs->GetNDofLocal(), paralleldofs); }
  
  template <typename TM>
  MasterInverse<TM> :: MasterInverse (const SparseMatrixTM<TM> & mat, 
				      shared_ptr<BitArray> subset, 
				      shared_ptr<ParallelDofs> hpardofs)
    
    : BaseMatrix(hpardofs), loc2glob(MyMPI_GetNTasks (hpardofs -> GetCommunicator()))
  {
    inv = nullptr;
    
    MPI_Comm comm = paralleldofs->GetCommunicator();
    int id = MyMPI_GetId (comm);
    int ntasks = MyMPI_GetNTasks(comm);

    // consistent enumeration
    
    int ndof = paralleldofs->GetNDofLocal();

    Array<int> global_nums(ndof);
    global_nums = -1;
    int num_master_dofs = 0;
    for (int i = 0; i < ndof; i++)
      if (paralleldofs -> IsMasterDof (i) && (!subset || (subset && subset->Test(i))))
	global_nums[i] = num_master_dofs++;
    

    Array<int> first_master_dof(ntasks);
    MPI_Allgather (&num_master_dofs, 1, MPI_INT, 
		   &first_master_dof[0], 1, MPI_INT, 
		   comm);

    
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

    paralleldofs -> ScatterDofData (global_nums);


    /*
    cout << "TESTING" << endl;
    for (int i = 0; i < ndof; i++)
      *testout << "dof" << i << " = " << pardofs->GetDofNodes()[i] << " free = " << subset->Test(i) << endl;
    *testout << "global_nums = " << endl << global_nums << endl;


    const MeshAccess & ma = pardofs->GetMeshAccess();

    for (int dest = 0; dest < ntasks; dest++)
      {
	*testout << "shared dofs with proc " << dest << endl;
	*testout << "exdofs " << endl << pardofs->GetExchangeDofs (dest) << endl;
	Array<int> procs;
	for (int i = 0; i < ma.GetNFaces(); i++)
	  {
	    ma.GetDistantProcs (Node(NT_FACE, i), procs);
	    if (procs.Contains (dest))
	      {
		Array<int> vnums;
		ma.GetFacePNums (i, vnums);
		int v1 = ma.GetGlobalNodeNum (Node(NT_VERTEX, vnums[0]));
		int v2 = ma.GetGlobalNodeNum (Node(NT_VERTEX, vnums[1]));
		int v3 = ma.GetGlobalNodeNum (Node(NT_VERTEX, vnums[2]));
		*testout << "F" << i << ": " << v1 << "-" << v2 << "-" << v3 << endl;
	      }
	  }
      }
    */


    /*
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

	int ndof = paralleldofs->GetNDofLocal();

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
		    *testout << "send (" << row << "," << rcols[j] 
			     << "), global = (" <<  global_nums[row] << "," << global_nums[rcols[j]] << ")" << endl;
		    rows.Append (global_nums[row]);
		    cols.Append (global_nums[rcols[j]]);
		    vals.Append (rvals[j]);
		  }
	    }

	MyMPI_Send (rows, 0, MPI_TAG_SOLVE, comm);
	MyMPI_Send (cols, 0, MPI_TAG_SOLVE, comm);
	MyMPI_Send (vals, 0, MPI_TAG_SOLVE, comm);
	MyMPI_Send (global_nums, 0, MPI_TAG_SOLVE, comm);

#ifdef USE_MUMPS
	if (mat.GetInverseType() == MUMPS)
	  {
	    SparseMatrixSymmetric<TM> dummy_matrix (1,1);
	    inv = make_shared<MumpsInverse<TM>> (dummy_matrix, nullptr, nullptr, true);
	  }
	else
#endif 
	  inv = nullptr;
      }

    else
      {
	cout << "create masterinverse" << endl;

	bool symmetric = (dynamic_cast<const SparseMatrixSymmetric<TM>*>(&mat) != NULL);
	cout << "symmetric? " << symmetric << endl;

	Array<int> rows, cols;
	Array<TM> vals;
	HashTable<INT<1>, int> ht_globdofs(100000);
	// int num_globdofs = 0; 

	for (int src = 1; src < ntasks; src++)
	  {
	    Array<int> hrows, hcols;
	    Array<TM> hvals;
	    Array<int> hglobid;

	    MyMPI_Recv (hrows, src, MPI_TAG_SOLVE, comm);
	    MyMPI_Recv (hcols, src, MPI_TAG_SOLVE, comm);
	    MyMPI_Recv (hvals, src, MPI_TAG_SOLVE, comm);
	    MyMPI_Recv (hglobid, src, MPI_TAG_SOLVE, comm);

	    *testout << "got from P" << src << ":" << endl
		     << "rows " << endl << hrows << endl
		     << "cols " << endl << hcols << endl;

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
	    if (symmetric && (r < c)) swap (r, c);
	    graph.AddUnique (r, c);
	  }

	// *testout << "graphi = " << endl << graph << endl;

	Array<int> els_per_row(num_glob_dofs);
	for (int i = 0; i < num_glob_dofs; i++)
	  els_per_row[i] = graph[i].Size();

	cout << "now build matrix" << endl;

	// SparseMatrixSymmetric<TM> matrix(els_per_row);
        auto matrix = symmetric ? make_shared<SparseMatrixSymmetric<TM>> (els_per_row)
	  : make_shared<SparseMatrix<TM>> (els_per_row);

	for (int i = 0; i < rows.Size(); i++)
	  {
	    int r = rows[i], c = cols[i];
	    if (symmetric && (r < c) ) swap (r, c);
	    matrix->CreatePosition(r, c);
	  }
	matrix->AsVector() = 0.0;

	for (int i = 0; i < rows.Size(); i++)
	  {
	    int r = rows[i], c = cols[i];
	    if (symmetric && (r < c)) swap (r, c);
	    (*matrix)(r,c) += vals[i];
	  }

	cout << "have matrix, now invert" << endl;

	matrix->SetInverseType (mat.GetInverseType());
	inv = matrix->InverseMatrix ();
      }

    // MPI_Barrier (ngs_comm);
  }

  template <typename TM>
  MasterInverse<TM> :: ~MasterInverse ()
  {
    ; // delete inv;
  }

  template <typename TM>
  void MasterInverse<TM> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    typedef typename mat_traits<TM>::TV_ROW TV;
    
    MPI_Comm comm = paralleldofs->GetCommunicator();
    int id = MyMPI_GetId(comm);
    int ntasks = MyMPI_GetNTasks(comm);

    bool is_x_cum = (dynamic_cast_ParallelBaseVector(x) . Status() == CUMULATED);
    x.Distribute();
    y.Cumulate();

    if (id > 0)
      {
	FlatVector<TV> fx = x.FV<TV> ();
	FlatVector<TV> fy = y.FV<TV> ();

	Array<TV> lx (select.Size());
	for (int i = 0; i < select.Size(); i++)
	  lx[i] = fx(select[i]);
	
	MPI_Request request = MyMPI_ISend (lx, 0, MPI_TAG_SOLVE, comm);
	MPI_Request_free (&request);

	// only for MUMPS:
	if (inv)
	  y = (*inv) * x;
	
	request = MyMPI_IRecv (lx, 0, MPI_TAG_SOLVE, comm);
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
	    MyMPI_Recv (lx, src, MPI_TAG_SOLVE, comm);

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
	    requ.Append (MyMPI_ISend (exdata[src], src, MPI_TAG_SOLVE, comm));
	  }
	MyMPI_WaitAll (requ);
      }

    if (is_x_cum)
      dynamic_cast_ParallelBaseVector(x) . Cumulate(); // AllReduce(&hoprocs);

  }

  ParallelMatrix :: ParallelMatrix (shared_ptr<BaseMatrix> amat, shared_ptr<ParallelDofs> apardofs)
    : BaseMatrix(apardofs), mat(amat)
  { 
    mat->SetParallelDofs (apardofs);
#ifdef USE_MUMPS
    mat->SetInverseType(MUMPS);
#else
    mat->SetInverseType(MASTERINVERSE);
#endif
  }


  ParallelMatrix :: ~ParallelMatrix ()
  {
    ; // delete &mat;
  }

  void ParallelMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    x.Cumulate();
    y.Distribute();
    mat->MultAdd (s, x, y);
  }

  void ParallelMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    x.Cumulate();
    y.Distribute();
    mat->MultTransAdd (s, x, y);
  }

  shared_ptr<BaseMatrix> ParallelMatrix :: CreateMatrix () const
  {
    return make_shared<ParallelMatrix> (mat->CreateMatrix(), paralleldofs);
  }

  AutoVector ParallelMatrix :: CreateVector () const
  {
    if (dynamic_cast<const SparseMatrix<double>*> (mat.get()))
      return make_shared<ParallelVVector<double>> (mat->Height(), paralleldofs);

    cerr << "ParallelMatrix::CreateVector not implemented for matrix type " 
	 << typeid(mat).name()
	 << endl;
    return shared_ptr<BaseVector>();
  }

  ostream & ParallelMatrix :: Print (ostream & ost) const
  {
    ost << *mat;
    return ost;
  }

  int ParallelMatrix :: VHeight() const
  {
    return mat->VHeight();
  }

  int ParallelMatrix :: VWidth() const
  {
    return mat->VWidth();
  }


  shared_ptr<BaseMatrix> ParallelMatrix::InverseMatrix (shared_ptr<BitArray> subset) const
  {
    shared_ptr<BaseMatrix> inv;
    inv = InverseMatrixTM<double> (subset);   if (inv) return inv;
    inv = InverseMatrixTM<Complex> (subset);  if (inv) return inv;
    inv = InverseMatrixTM<Mat<2> > (subset);   if (inv) return inv;
    inv = InverseMatrixTM<Mat<3> > (subset);   if (inv) return inv;

    inv = InverseMatrixTM<Mat<2,2,Complex> > (subset);   if (inv) return inv;
    inv = InverseMatrixTM<Mat<3,3,Complex> > (subset);   if (inv) return inv;

    throw Exception ("ParallelMatrix::Inverse(BitArray) not available, typeid(mat) = " 
		     + ToString (typeid(mat).name()));
  }


  
  template <typename TM>
  shared_ptr<BaseMatrix> ParallelMatrix::InverseMatrixTM (shared_ptr<BitArray> subset) const
  {
    const SparseMatrixTM<TM> * dmat = dynamic_cast<const SparseMatrixTM<TM>*> (mat.get());
    if (!dmat) return NULL;

#ifdef USE_MUMPS
    bool symmetric = dynamic_cast<const SparseMatrixSymmetric<TM>*> (mat.get()) != NULL;
    if (mat->GetInverseType() == MUMPS)
      return make_shared<ParallelMumpsInverse<TM>> (*dmat, subset, nullptr, paralleldofs, symmetric);
    else 
#endif
      return make_shared<MasterInverse<TM>> (*dmat, subset, paralleldofs);
  }



  shared_ptr<BaseMatrix> ParallelMatrix::InverseMatrix (shared_ptr<const Array<int>> clusters) const
  {
    throw Exception ("ParallelMatrix::Inverse(cluster) not available");
  }

  INVERSETYPE ParallelMatrix::SetInverseType (INVERSETYPE ainversetype) const
  {
    return mat->SetInverseType (ainversetype);
  }

  INVERSETYPE ParallelMatrix::SetInverseType (string ainversetype) const
  {
    return mat->SetInverseType (ainversetype);
  }
  
  INVERSETYPE ParallelMatrix::GetInverseType () const
  {
    return mat->GetInverseType ();
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






  FETI_Jump_Matrix :: FETI_Jump_Matrix (shared_ptr<ParallelDofs> apardofs, shared_ptr<ParallelDofs> au_paralleldofs)
    : BaseMatrix(apardofs), u_paralleldofs(au_paralleldofs)
  {
    
    size_t njs = 0;
    for (auto p:paralleldofs->GetDistantProcs())
      njs += paralleldofs->GetExchangeDofs(p).Size();

    Array<size_t> ones(njs);
    ones = 1;
    Table<int> dps(ones);
    njs = 0;
    for (auto p:paralleldofs->GetDistantProcs()) {
      for (auto d:paralleldofs->GetExchangeDofs(p)) {
	dps[njs++][0] = p;
      }
    }    

    this->jump_paralleldofs = make_shared<ParallelDofs>(paralleldofs->GetCommunicator(), move(dps));

    return;
  }


  void FETI_Jump_Matrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    y.Distribute();
    size_t count = 0;
    /*
      for (auto p:paralleldofs->GetDistantProcs()) {
      auto exdofs = paralleldofs->GetExchangeDofs(p);
      if (p<MyMPI_GetId(paralleldofs->GetCommunicator())) {
      for (auto k:Range(exdofs.Size())) {
      y.FVDouble()[count++] -= s*x.FVDouble()[exdofs[k]];
      }
      }
      else {
      for (auto k:Range(exdofs.Size())) {
      y.FVDouble()[count++] += s*x.FVDouble()[exdofs[k]];
      }
      }
      }
    */
    auto me = MyMPI_GetId(paralleldofs->GetCommunicator());
    auto fx = x.FVDouble();
    auto fy = y.FVDouble();
    for (auto p : paralleldofs->GetDistantProcs())
      {
        double signed_s = (p < me) ? -s : s;
        for (auto exdof : paralleldofs->GetExchangeDofs(p))
          fy[count++] += signed_s * fx[exdof];
      }
  }

  void FETI_Jump_Matrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    x.Cumulate();
    size_t count = 0;
    for (auto p:paralleldofs->GetDistantProcs()) {
      auto exdofs = paralleldofs->GetExchangeDofs(p);
      if (p<MyMPI_GetId(paralleldofs->GetCommunicator())) {
	for (auto k:Range(exdofs.Size())) {
	  y.FVDouble()[exdofs[k]] -= s*x.FVDouble()[count++];
	}
      }
      else {
	for (auto k:Range(exdofs.Size())) {
	  y.FVDouble()[exdofs[k]] += s*x.FVDouble()[count++];
	}
      }
    }
  }

  
  AutoVector FETI_Jump_Matrix :: CreateRowVector () const
  {
    // throw Exception("Called FETI_Jump_Matrix :: CreateRowVector, this is not well defined");
    if (u_paralleldofs==nullptr) {
      return make_shared<VVector<double>>(VHeight());
    }
    return make_shared<ParallelVVector<double>> (u_paralleldofs->GetNDofLocal(),
						 u_paralleldofs);
  }
  
  AutoVector FETI_Jump_Matrix :: CreateColVector () const
  {
    return make_shared<ParallelVVector<double>> (jump_paralleldofs->GetNDofLocal(),
						 jump_paralleldofs);
  }
  
  
}

#endif
