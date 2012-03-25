/*********************************************************************/
/* File:   parallel_matrices.cpp                                     */
/* Author: Joachim Schoeberl                                         */
/* Date:   June 2011                                                 */
/*********************************************************************/

// #undef USE_MUMPS

#ifdef PARALLEL
 
#include <la.hpp>
#include <comp.hpp>
#include <parallelngs.hpp>


namespace ngla
{
  using namespace ngparallel;


  template <typename TM>
  MasterInverse<TM> :: MasterInverse (const SparseMatrixTM<TM> & mat, 
				      const BitArray * subset, const ParallelDofs * pardofs)
    : loc2glob(ntasks)
  {
    inv = NULL;

    if (id != 0)
      {
	const MeshAccess & ma = pardofs -> GetMeshAccess();

	int ndof = pardofs->GetNDof();

	Array<int> rows, cols, globid(3*ndof);
	Array<TM> vals;

	for (int row = 0; row < ndof; row++)
	  if (!subset || subset->Test(row))
	    select.Append (row);

	
	Array<int> compress(ndof);
	compress = -1;
	for (int i = 0; i < select.Size(); i++)
	  compress[select[i]] = i;


	globid = -1;

	/*
	Array<int> dnums;
	for (NODE_TYPE nt = NT_VERTEX; nt <= NT_CELL; nt++)
	  for (int i = 0; i < ma.GetNNodes (nt); i++)
	    {
	      fes.GetNodeDofNrs (nt, i, dnums);
	      // int distnum = parallelma->GetGlobalNodeNum (nt, i);
	      int distnum = ma.GetGlobalNodeNum (Node (nt, i));

	      for (int j = 0; j < dnums.Size(); j++)
		if (dnums[j] != -1 && (!subset || subset->Test(dnums[j])))
		  {
		    int dn = dnums[j];
		    globid[3*dn+0] = int(nt);
		    globid[3*dn+1] = distnum;
		    globid[3*dn+2] = j;
		  }
	    }
	*/

	const Array<Node> & dofnodes = pardofs -> GetDofNodes();

	Array<int> nnodes(4);
	for (int j = 0; j < 4; j++) nnodes[j] = ma.GetNNodes(NODE_TYPE(j));
	Table<int> tdofnr(nnodes);
	for (int j = 0; j < 4; j++) tdofnr[j] = 0;
	Array<int> dofnr(ndof);

	for (int i = 0; i < ndof; i++)
	  dofnr[i] = tdofnr[dofnodes[i].GetType()][dofnodes[i].GetNr()]++;

	for (int i = 0; i < ndof; i++)
	  if (!subset || subset->Test(i))
	    {
	      globid[3*i+0] = dofnodes[i].GetType();
	      globid[3*i+1] = ma.GetGlobalNodeNum (dofnodes[i]);
	      globid[3*i+2] = dofnr[i];
	    }

	for (int row = 0; row < mat.Height(); row++)
	  if (!subset || subset->Test(row))
	    {
	      FlatArray<int> rcols = mat.GetRowIndices(row);
	      FlatVector<TM> rvals = mat.GetRowValues(row);

	      for (int j = 0; j < rcols.Size(); j++)
		if (!subset || subset->Test(rcols[j]))
		  {
		    rows.Append (row);
		    cols.Append (rcols[j]);
		    vals.Append (rvals[j]);
		  }
	    }

	MyMPI_Send (rows, 0);
	MyMPI_Send (cols, 0);
	MyMPI_Send (vals, 0);
	MyMPI_Send (globid, 0);

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
	HashTable<INT<3>, int> ht_globdofs(100000);
	int num_globdofs = 0;

	for (int src = 1; src < ntasks; src++)
	  {
	    Array<int> hrows, hcols;
	    Array<TM> hvals;
	    Array<int> hglobid;

	    MyMPI_Recv (hrows, src);
	    MyMPI_Recv (hcols, src);
	    MyMPI_Recv (hvals, src);
	    MyMPI_Recv (hglobid, src);

	    /*
	    *testout << "got from proc " << src << ":" << endl
		     << "rows = " << endl << hrows << endl
		     << "cols = " << endl << hcols << endl
		     << "vals = " << endl << hvals << endl
		     << "globid = " << endl << hglobid << endl;
	    */

	    for (int i = 0; i < hrows.Size(); i++)
	      if (hglobid[3*hrows[i]] < 0)
		cout << "globid missing (rows) !!!! " << endl;
	    for (int i = 0; i < hrows.Size(); i++)
	      if (hglobid[3*hcols[i]] < 0)
		cout << "globid missing (cols) !!!! " << endl;


	    Array<int> full_loc2glob(hglobid.Size()/3);
	    full_loc2glob = -1;
	    for (int i = 0; i < hglobid.Size(); i += 3)
	      {
		if (hglobid[i] == -1) continue;
	      
		INT<3> nentry;
		nentry[0] = hglobid[i];
		nentry[1] = hglobid[i+1];
		nentry[2] = hglobid[i+2];

		int found;

		if (ht_globdofs.Used (nentry))
		  found = ht_globdofs.Get(nentry);
		else
		  {
		    found = num_globdofs;
		    num_globdofs++;
		    ht_globdofs.Set(nentry, found);
		  }
	      
		loc2glob.Add (src, found);
		full_loc2glob[i/3] = found;
	      }
	    // *testout << "full_loc2glob = " << endl << full_loc2glob << endl;
	 
	    for (int i = 0; i < hrows.Size(); i++)
	      {
		rows.Append (full_loc2glob[hrows[i]]);
		cols.Append (full_loc2glob[hcols[i]]);
		vals.Append (hvals[i]);
	      }

	    cout << "\rmaster: got data from " << src << flush;
	  }
	cout << endl;
	/*
	 *testout << "rows = " << endl << rows << endl;
	 *testout << "cols = " << endl << cols << endl;
	 *testout << "vals = " << endl << vals << endl;
	 */
	cout << "now build graph" << endl;

	// build matrix
	DynamicTable<int> graph(num_globdofs);
	cout << "n = " << num_globdofs << endl;
	for (int i = 0; i < rows.Size(); i++)
	  {
	    int r = rows[i], c = cols[i];
	    if (r < c) swap (r, c);
	    if (r < 0 || c < 0)
	      cout << "r,c = " << r << ", " << c << endl;
	    graph.AddUnique (r, c);
	  }

	// *testout << "graphi = " << endl << graph << endl;

	Array<int> els_per_row(num_globdofs);
	for (int i = 0; i < num_globdofs; i++)
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

    MPI_Barrier (ngs_comm);
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

	// *testout << "before *inv, hx = " << endl << hx << endl;

	hy = (*inv) * hx;

	// *testout << "before *inv, hy = " << endl << hy << endl;


	Array<MPI_Request> requ;
	for (int src = 1; src < ntasks; src++)
	  {
	    FlatArray<int> selecti = loc2glob[src];
	    for (int j = 0; j < selecti.Size(); j++)
	      exdata[src][j] = hy(selecti[j]);
	    requ.Append (MyMPI_ISend (exdata[src], src, MPI_TAG_SOLVE));
	  }
	MPI_Waitall (requ.Size(), &requ[0], MPI_STATUS_IGNORE);
      }

    if (is_x_cum)
      dynamic_cast<const ParallelBaseVector&> (x) . Cumulate(); // AllReduce(&hoprocs);

  }



  ParallelMatrix :: ~ParallelMatrix ()
  {
    ;
  }

  void ParallelMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    x.Cumulate();
    y.Distribute();
    if (id > 0)
      mat.MultAdd (s, x, y);
  }



  void ParallelMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    x.Cumulate();
    y.Distribute();
    if (id > 0)
      mat.MultTransAdd (s, x, y);
  }


  BaseMatrix * ParallelMatrix :: CreateMatrix () const
  {
    return new ParallelMatrix (*mat.CreateMatrix(), pardofs);
  }

  BaseVector * ParallelMatrix :: CreateVector () const
  {
    cerr << "ParallelMatrix::CreateVector not implemented" << endl;
    return NULL;
  }


  ostream & ParallelMatrix :: Print (ostream & ost) const
  {
    if (id > 0)
      ost << mat;
    return ost;
  }

  int ParallelMatrix :: VHeight() const
  {
    return (id == 0) ? 0 : mat.VHeight();
  }
  int ParallelMatrix :: VWidth() const
  {
    return (id == 0) ? 0 : mat.VWidth();
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

  INVERSETYPE ParallelMatrix::SetInverseType ( INVERSETYPE ainversetype ) const
  {
    // if (id != 0) 
    return mat.SetInverseType (ainversetype);
    // return SPARSECHOLESKY;
  }

  INVERSETYPE ParallelMatrix::SetInverseType ( string ainversetype ) const
  {
    // if (id != 0) 
    return mat.SetInverseType (ainversetype);
    // return SPARSECHOLESKY;
  }
  
  INVERSETYPE ParallelMatrix::GetInverseType () const
  {
    // if (id != 0) 
    return mat.GetInverseType ();
    // return SPARSECHOLESKY;
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
