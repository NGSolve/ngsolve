/*********************************************************************/
/* File:   elementbyelement.cpp                                      */
/* Author: Joachim Schoeberl                                         */
/* Date:   June 2010                                                 */
/*********************************************************************/

/* 
   Element by element matrix
*/

#include <la.hpp>
namespace ngla
{
  using namespace ngla;

  
  template <class SCAL> class ElementByElementMatrix;
  template <class SCAL>
  ElementByElementMatrix<SCAL> :: ElementByElementMatrix (int h, int ane, bool isymmetric) 
  {
    clone.SetSize(ane);
    clone.Clear();
    symmetric=isymmetric;
    height = h; 
    ne = ane; 
    elmats.SetSize(ne);
    rowdnums.SetSize(ne);
    coldnums.SetSize(ne);
    disjointrows = false;
    disjointcols = false;
    for (int i = 0; i < ne; i++)
      {
        elmats[i].AssignMemory (0, 0, NULL);
	new (&rowdnums[i]) FlatArray<int> (0, NULL);
	new (&coldnums[i]) FlatArray<int> (0, NULL);
	/*
        rowdnums[i] = FlatArray<int> (0, NULL);
	coldnums[i] = FlatArray<int> (0, NULL);
	*/
      }
  }
  
  template <class SCAL>
  ElementByElementMatrix<SCAL> :: ElementByElementMatrix (int h, int ane, bool isymmetric, bool adisjointrows, bool adisjointcols) 
  {
    clone.SetSize(ane);
    clone.Clear();
    symmetric=isymmetric;
    height = h; 
    ne = ane; 
    elmats.SetSize(ne);
    rowdnums.SetSize(ne);
    coldnums.SetSize(ne);
    disjointrows = adisjointrows;
    disjointcols = adisjointcols;
    for (int i = 0; i < ne; i++)
      {
        elmats[i].AssignMemory (0, 0, NULL);
	new (&rowdnums[i]) FlatArray<int> (0, NULL);
	new (&coldnums[i]) FlatArray<int> (0, NULL);
	/*
        rowdnums[i] = FlatArray<int> (0, NULL);
	coldnums[i] = FlatArray<int> (0, NULL);
	*/
      }
  }


  template <class SCAL>
  ElementByElementMatrix<SCAL> :: ~ElementByElementMatrix ()
  {
    for (int i = 0; i < ne; i++)
      if (!clone.Test(i))
	{
	  delete [] &(elmats[i](0,0));
	  delete [] &(rowdnums[i])[0];
	  delete [] &(coldnums[i])[0];
	}
  }
  
  template <class SCAL>
  void ElementByElementMatrix<SCAL> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("EBE-matrix::MultAdd");
    RegionTimer reg (timer);

    int maxs = 0;
    for (int i = 0; i < coldnums.Size(); i++)
      maxs = max2 (maxs, coldnums[i].Size());

    if (disjointrows)
      {
#pragma omp parallel
	{    
	  ArrayMem<SCAL, 100> mem1(maxs); 

	  FlatVector<SCAL> vx = x.FV<SCAL> (); 
	  FlatVector<SCAL> vy = y.FV<SCAL> (); 

#pragma omp for
	  for (int i = 0; i < rowdnums.Size(); i++) //sum over all elements
	    {
	      FlatArray<int> rdi = rowdnums[i];
	      FlatArray<int> cdi = coldnums[i];
	      
	      if (!rdi.Size() || !cdi.Size()) continue;
	      
	      FlatVector<SCAL> hv1(cdi.Size(), &mem1[0]);

	      hv1 = vx(cdi);
	      vy(rdi) += s * elmats[i] * hv1;

	      timer.AddFlops (cdi.Size()*rdi.Size());
	    }
	}//end of parallel    
	
      }
    else
      {    
	ArrayMem<SCAL, 100> mem1(maxs), mem2(maxs);
	
	FlatVector<SCAL> vx = dynamic_cast<const S_BaseVector<SCAL> & >(x).FVScal();
	FlatVector<SCAL> vy = dynamic_cast<S_BaseVector<SCAL> & >(y).FVScal();
	
	for (int i = 0; i < rowdnums.Size(); i++) //sum over all elements
	  {
	    FlatArray<int> rdi = rowdnums[i];
	    FlatArray<int> cdi = coldnums[i];
	    
	    if (!rdi.Size() || !cdi.Size()) continue;
	    
	    FlatVector<SCAL> hv(cdi.Size(), &mem1[0]);
	    
	    hv = vx(cdi);
	    vy(rdi) += s * elmats[i] * hv;
	    
	    timer.AddFlops (cdi.Size()*rdi.Size());
	  }
      }
  }

  
  template <class SCAL>
  void ElementByElementMatrix<SCAL> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
//     cout << " ElementByElementMatrix<SCAL> :: MultTansAdd here " << endl << flush;
    static Timer timer("EBE-matrix::MultTransAdd");
    RegionTimer reg (timer);
    int maxs = 0;
    for (int i = 0; i < rowdnums.Size(); i++)
      maxs = max2 (maxs, rowdnums[i].Size());
    
    if (disjointcols)
      {
#pragma omp parallel
	{
	  ArrayMem<SCAL, 100> mem1(maxs);
	  
	  FlatVector<SCAL> vx = dynamic_cast<const S_BaseVector<SCAL> & >(x).FVScal();
	  FlatVector<SCAL> vy = dynamic_cast<S_BaseVector<SCAL> & >(y).FVScal();

    #pragma omp  for    
	  for (int i = 0; i < coldnums.Size(); i++) //sum over all elements
	    {
	      FlatArray<int> rdi (rowdnums[i]);
	      FlatArray<int> cdi (coldnums[i]);
	      if (!rdi.Size() || !cdi.Size()) continue;
	      
	      FlatVector<SCAL> hv1(rdi.Size(), &mem1[0]);

	      hv1 = vx(rdi);
	      vy(cdi) += s * Trans(elmats[i]) * hv1;

	      timer.AddFlops (cdi.Size()*rdi.Size());
	    }
	}//end of parallel
      }
    else
      {
	ArrayMem<SCAL, 100> mem1(maxs);
	
	FlatVector<SCAL> vx = x.FV<SCAL> (); // dynamic_cast<const S_BaseVector<SCAL> & >(x).FVScal();
	FlatVector<SCAL> vy = y.FV<SCAL> (); // dynamic_cast<S_BaseVector<SCAL> & >(y).FVScal();
	
	for (int i = 0; i < coldnums.Size(); i++) //sum over all elements
	  {
	    FlatArray<int> rdi (rowdnums[i]);
	    FlatArray<int> cdi (coldnums[i]);
	    
	    if (!rdi.Size() || !cdi.Size()) continue;
	    
	    FlatVector<SCAL> hv1(rdi.Size(), &mem1[0]);
	    
	    hv1 = vx(rdi);
	    vy(cdi) += s * Trans(elmats[i]) * hv1;

	    timer.AddFlops (cdi.Size()*rdi.Size());
	  }
      }
  }


  template <class SCAL>
  BaseMatrix *  ElementByElementMatrix<SCAL> :: InverseMatrix ( BitArray * subset ) const
  {
    cout << "wird das tatsaechlich verwendet ???" << endl;
    throw Exception ("not available any longer!");
    return NULL;
/*    ElementByElementMatrix<SCAL> * invmat = new ElementByElementMatrix<SCAL> (height, ne);

    int maxs = 0;
    for (int i = 0; i < dnums.Size(); i++)
      maxs = max2 (maxs, dnums[i].Size());

    LocalHeap lh (maxs*maxs*sizeof(SCAL)+100);

    for ( int i = 0; i < dnums.Size(); i++ )
      {
        int nd = dnums[i] . Size();
        FlatMatrix<SCAL> mat(nd, nd, lh);
        Array<int> dnumsarray(nd);
        for ( int j = 0; j < nd; j++ )
          dnumsarray[j] = dnums[i][j];
        mat = elmats[i];

        LapackInverse(mat);

        invmat -> AddElementMatrix(i, dnumsarray, dnumsarray, mat);
        lh.CleanUp();
      }
  
    return invmat;*/
  }

  template <class SCAL>
  void ElementByElementMatrix<SCAL> :: AddElementMatrix (int elnr,
                                                         const FlatArray<int> & rowdnums_in,
                                                         const FlatArray<int> & coldnums_in,
                                                         const FlatMatrix<SCAL> & elmat)
  {
    ArrayMem<int,50> usedrows;
    for (int i = 0; i < rowdnums_in.Size(); i++)
      if (rowdnums_in[i] >= 0) usedrows.Append(i);
    int sr = usedrows.Size();

    ArrayMem<int,50> usedcols;
    for (int i = 0; i < coldnums_in.Size(); i++)
      if (coldnums_in[i] >= 0) usedcols.Append(i);
    int sc = usedcols.Size();

    FlatMatrix<SCAL> mat (sr,sc, new SCAL[sr*sc]);
    // mat = elmat.Rows(usedrows).Cols(usedcols);

    for (int i = 0; i < sr; i++)
      for (int j = 0; j < sc; j++)
        mat(i,j) = elmat(usedrows[i], usedcols[j]);

    FlatArray<int> dnr(sr, new int[sr]);
    for (int i = 0; i < sr; i++)
      dnr[i] = rowdnums_in[usedrows[i]];
    
    FlatArray<int> dnc(sc, new int[sc]);
    for (int j = 0; j < sc; j++)
      dnc[j] = coldnums_in[usedcols[j]];

    if (elnr < elmats.Size())
      {
        // rowdnums[elnr] = dnr;
	// coldnums[elnr] = dnc;

	new (&rowdnums[elnr]) FlatArray<int> (dnr);
	new (&coldnums[elnr]) FlatArray<int> (dnc);

        elmats[elnr].AssignMemory (sr, sc, &mat(0,0));
      }
  }

  template <class SCAL>
  void ElementByElementMatrix<SCAL> :: AddCloneElementMatrix (int elnr,
                                                         const FlatArray<int> & rowdnums_in,
                                                         const FlatArray<int> & coldnums_in,
                                                         int refelnr)
  {
    ArrayMem<int,50> usedrows;
    for (int i = 0; i < rowdnums_in.Size(); i++)
      if (rowdnums_in[i] >= 0) usedrows.Append(i);
    int sr = usedrows.Size();

    ArrayMem<int,50> usedcols;
    for (int i = 0; i < coldnums_in.Size(); i++)
      if (coldnums_in[i] >= 0) usedcols.Append(i);
    int sc = usedcols.Size();

//     FlatMatrix<SCAL> mat (sr,sc, new SCAL[sr*sc]);
//     for (int i = 0; i < sr; i++)
//       for (int j = 0; j < sc; j++)
//         mat(i,j) = elmat(usedrows[i], usedcols[j]);

    FlatArray<int> dnr(sr, new int[sr]);
    for (int i = 0; i < sr; i++)
      dnr[i] = rowdnums_in[usedrows[i]];
    
    FlatArray<int> dnc(sc, new int[sc]);
    for (int j = 0; j < sc; j++)
      dnc[j] = coldnums_in[usedcols[j]];

    if (elnr < elmats.Size())
      {
        // rowdnums[elnr] = dnr;
	// coldnums[elnr] = dnc;

	new (&rowdnums[elnr]) FlatArray<int> (dnr);
	new (&coldnums[elnr]) FlatArray<int> (dnc);


        elmats[elnr].AssignMemory (sr, sc, &elmats[refelnr](0,0));
	clone.Set(elnr);
      }
  }

  template <class SCAL>
  BaseBlockJacobiPrecond * ElementByElementMatrix<SCAL> :: 
  CreateBlockJacobiPrecond (Table<int> & blocks,
                            const BaseVector * constraint, int * paralleloptions) const
  { 
    return 0;//new helmholtz_exp_cpp::BlockJacobiPrecond_ElByEl<SCAL> (*this, blocks );
  }


  template <class SCAL>
  ostream & ElementByElementMatrix<SCAL> :: Print (ostream & ost) const
  {
      ost << "Element-by-Element Matrix:" << endl;
      ost << "num blocks = " << elmats.Size();
      for (int i = 0; i < elmats.Size(); i++)
	{
	  ost << "block " << i << endl;
	  ost << "rows = " << rowdnums[i] << endl;
	  ost << "cols = " << coldnums[i] << endl;
	  ost << "matrix = " << elmats[i] << endl;
	}
      return ost;
    }



  
  template class ElementByElementMatrix<double>;
  template class ElementByElementMatrix<Complex>;  
  
}
