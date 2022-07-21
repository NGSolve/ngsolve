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
    InitMemoryTracing();
    clone.SetSize(ane);
    clone.Clear();
    symmetric=isymmetric;
    height = h; 
    width = h; 
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


  template <class SCAL> class ElementByElementMatrix;
  template <class SCAL>
  ElementByElementMatrix<SCAL> :: ElementByElementMatrix (int h, int w, int ane, bool isymmetric) 
  {
    InitMemoryTracing();
    clone.SetSize(ane);
    clone.Clear();
    symmetric=isymmetric;
    height = h; 
    width = w; 
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
  ElementByElementMatrix<SCAL> :: ElementByElementMatrix (int h, int w, int ane, bool isymmetric, bool adisjointrows, bool adisjointcols) 
  {
    InitMemoryTracing();
    clone.SetSize(ane);
    clone.Clear();
    symmetric=isymmetric;
    height = h; 
    width = w; 
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
  ElementByElementMatrix<SCAL> ::
  ElementByElementMatrix (size_t h, size_t w, 
                          FlatArray<int> nrowi, FlatArray<int> ncoli,
                          bool isymmetric, bool adisjointrows, bool adisjointcols)
    : height(h), width(w), ne(nrowi.Size()), symmetric(isymmetric), disjointrows(adisjointrows), disjointcols(adisjointcols)
  {
    InitMemoryTracing();
    clone.SetSize(ne);
    clone.Clear();
    
    size_t totmem_row = 0;
    size_t totmem_col = 0;
    size_t totmem_values = 0;
    for (size_t i = 0; i < nrowi.Size(); i++)
      {
        totmem_row += nrowi[i];
        totmem_col += ncoli[i];
        totmem_values += nrowi[i]*ncoli[i];
      }

    allrow.SetSize(totmem_row+1);  // to mark as array-version
    allcol.SetSize(totmem_col+1);
    allvalues.SetSize(totmem_values+1);
    allrow = -1;
    allcol = -1;
    /*
    for (int i = 0; i < totmem_row; i++)
      allrow.get()[i] = -1;
    for (int i = 0; i < totmem_col; i++)
      allcol.get()[i] = -1;
    for (int i = 0; i < totmem_values; i++)
      allvalues.get()[i] = -1;
    */
    
    rowdnums.SetSize(ne);
    coldnums.SetSize(ne);
    elmats.SetSize(ne);
    totmem_row = 0;
    totmem_col = 0;
    totmem_values = 0;
    for (size_t i = 0; i < nrowi.Size(); i++)
      {
        new (&rowdnums[i]) FlatArray<int> (nrowi[i], allrow.Addr(totmem_row));
        new (&coldnums[i]) FlatArray<int> (ncoli[i], allcol.Addr(totmem_col));
        elmats[i].AssignMemory (nrowi[i], ncoli[i], allvalues.Addr(totmem_values));
        totmem_row += nrowi[i];
        totmem_col += ncoli[i];
        totmem_values += nrowi[i]*ncoli[i];
      }
  }
  


  
  template <class SCAL>
  ElementByElementMatrix<SCAL> :: ~ElementByElementMatrix ()
  {
    if (allvalues.Size())
      return;  // all memory in unique_ptrs 
    
    for (int i = 0; i < ne; i++)
      if (!clone.Test(i))
	{
	  delete [] &(elmats[i](0,0));
	  if (rowdnums[i].Size() > 0)
	    delete [] &(rowdnums[i])[0];
	  if (coldnums[i].Size() > 0)
	    delete [] &(coldnums[i])[0];
	}
  }
  
  template <class SCAL>
  void ElementByElementMatrix<SCAL> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("EBE-matrix::MultAdd");
    RegionTimer reg (timer);

    size_t maxs = 0;
    for (size_t i = 0; i < coldnums.Size(); i++)
      maxs = max2 (maxs, coldnums[i].Size());

    if (disjointrows)
      {
        ParallelForRange( IntRange(rowdnums.Size()), [&] ( IntRange r )
	{    
	  ArrayMem<SCAL, 100> mem1(maxs); 

	  FlatVector<SCAL> vx = x.FV<SCAL> (); 
	  FlatVector<SCAL> vy = y.FV<SCAL> (); 

	  for (int i : r)  //sum over all elements
	    {
	      FlatArray<int> rdi = rowdnums[i];
	      FlatArray<int> cdi = coldnums[i];
	      
	      if (!rdi.Size() || !cdi.Size()) continue;
	      if (rdi[0] == -1 || cdi[0] == -1) continue;  // reserved but not used
              // throw ("Illegal index in ebe mult");
	      FlatVector<SCAL> hv1(cdi.Size(), &mem1[0]);

	      hv1 = vx(cdi);
	      vy(rdi) += s * elmats[i] * hv1;

	      timer.AddFlops (cdi.Size()*rdi.Size());
	    }
	}); //end of parallel    
	
      }
    else
      {    
	ArrayMem<SCAL, 100> mem1(maxs), mem2(maxs);
	
	FlatVector<SCAL> vx = x.FV<SCAL>(); 
	FlatVector<SCAL> vy = y.FV<SCAL>(); 
	
	for (int i = 0; i < rowdnums.Size(); i++) //sum over all elements
	  {
	    FlatArray<int> rdi = rowdnums[i];
	    FlatArray<int> cdi = coldnums[i];
	    
	    if (!rdi.Size() || !cdi.Size()) continue;
            if (rdi[0] == -1 || cdi[0] == -1) continue;  // reserved but not used            
	    
	    FlatVector<SCAL> hv(cdi.Size(), &mem1[0]);
	    
	    hv = vx(cdi);
	    vy(rdi) += s * elmats[i] * hv;
	    
	    timer.AddFlops (cdi.Size()*rdi.Size());
	  }
      }
  }




  template <>
  void ElementByElementMatrix<double> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("EBE-matrix::MultAdd");
    RegionTimer reg (timer);

    size_t maxs = 0;
    for (size_t i = 0; i < coldnums.Size(); i++)
      maxs = max2 (maxs, coldnums[i].Size());


    if (task_manager)
      {
        FlatVector<> vx = x.FV<double> (); 
        FlatVector<> vy = y.FV<double> (); 
        
        task_manager -> CreateJob 
          ( [&] (const TaskInfo & ti)
            {
              ArrayMem<double, 100> mem1(max_row_size);
              ArrayMem<double, 100> mem2(max_col_size);
              
              auto myr = Range(coldnums).Split (ti.task_nr, ti.ntasks);
              
              for (int i : myr)
                {
                  FlatArray<int> rdi (rowdnums[i]);
                  FlatArray<int> cdi (coldnums[i]);
                  
                  if (!rdi.Size() || !cdi.Size()) continue;
                  if (rdi[0] == -1 || cdi[0] == -1) continue;  // reserved but not used
                  
                  FlatVector<> hv1(rdi.Size(), &mem1[0]);
                  FlatVector<> hv2(cdi.Size(), &mem2[0]);
                  
                  hv2 = vx(cdi);
                  hv1 = elmats[i] * hv2;

                  for (int j = 0; j < rdi.Size(); j++)
                  {
                      double &res = vy(rdi[j]);
                      double t = s*hv1(j);
                      AtomicAdd(res, t);
                  }
                  
                  timer.AddFlops (cdi.Size()*rdi.Size());
                }
            });
      }
    
    else
      {    
	ArrayMem<double, 100> mem1(maxs), mem2(maxs);
	
	FlatVector<double> vx = x.FV<double>(); 
	FlatVector<double> vy = y.FV<double>(); 
	
	for (int i = 0; i < rowdnums.Size(); i++) //sum over all elements
	  {
	    FlatArray<int> rdi = rowdnums[i];
	    FlatArray<int> cdi = coldnums[i];
	    
	    if (!rdi.Size() || !cdi.Size()) continue;
            if (rdi[0] == -1 || cdi[0] == -1) continue;  // reserved but not used
            
	    FlatVector<double> hv(cdi.Size(), &mem1[0]);
	    
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
    size_t maxs = 0;
    for (size_t i = 0; i < rowdnums.Size(); i++)
      maxs = max2 (maxs, rowdnums[i].Size());
    
    if (disjointcols)
      {
        ParallelForRange( IntRange(coldnums.Size()), [&] ( IntRange r )
	{
	  ArrayMem<SCAL, 100> mem1(maxs);
	  
	  FlatVector<SCAL> vx = dynamic_cast<const S_BaseVector<SCAL> & >(x).FVScal();
	  FlatVector<SCAL> vy = dynamic_cast<S_BaseVector<SCAL> & >(y).FVScal();

	  for (int i : r) //sum over all elements
	    {
	      FlatArray<int> rdi (rowdnums[i]);
	      FlatArray<int> cdi (coldnums[i]);
	      if (!rdi.Size() || !cdi.Size()) continue;
              if (rdi[0] == -1 || cdi[0] == -1) continue;  // reserved but not used
              
	      FlatVector<SCAL> hv1(rdi.Size(), &mem1[0]);

	      hv1 = vx(rdi);
	      vy(cdi) += s * Trans(elmats[i]) * hv1;

	      timer.AddFlops (cdi.Size()*rdi.Size());
	    }
	});//end of parallel
      }
    else
      {
        ArrayMem<SCAL, 100> mem1(maxs);
        
        FlatVector<SCAL> vx = x.FV<SCAL> (); 
        FlatVector<SCAL> vy = y.FV<SCAL> (); 
        
        for (int i = 0; i < coldnums.Size(); i++) //sum over all elements
          {
            FlatArray<int> rdi (rowdnums[i]);
            FlatArray<int> cdi (coldnums[i]);
            
            if (!rdi.Size() || !cdi.Size()) continue;
            if (rdi[0] == -1 || cdi[0] == -1) continue;  // reserved but not used
            
            FlatVector<SCAL> hv1(rdi.Size(), &mem1[0]);
            
            hv1 = vx(rdi);
            vy(cdi) += s * Trans(elmats[i]) * hv1;
            
            timer.AddFlops (cdi.Size()*rdi.Size());
          }
      }
  }





  template <>
  void ElementByElementMatrix<double> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer timer("EBE-matrix<double>::MultTransAdd");
    RegionTimer reg (timer);
    size_t maxs = 0;
    for (size_t i = 0; i < rowdnums.Size(); i++)
      maxs = max2 (maxs, rowdnums[i].Size());
    
    if (false) // disjointcols)
      {
        ParallelForRange( IntRange(coldnums.Size()), [&] ( IntRange r )
	{
	  ArrayMem<double, 100> mem1(maxs);
	  
	  FlatVector<> vx = dynamic_cast<const S_BaseVector<double> & >(x).FVScal();
	  FlatVector<> vy = dynamic_cast<S_BaseVector<double> & >(y).FVScal();

	  for (int i : r) //sum over all elements
	    {
	      FlatArray<int> rdi (rowdnums[i]);
	      FlatArray<int> cdi (coldnums[i]);
	      if (!rdi.Size() || !cdi.Size()) continue;
              if (rdi[0] == -1 || cdi[0] == -1) continue;  // reserved but not used
            
	      FlatVector<> hv1(rdi.Size(), &mem1[0]);

	      hv1 = vx(rdi);
	      vy(cdi) += s * Trans(elmats[i]) * hv1;

	      timer.AddFlops (cdi.Size()*rdi.Size());
	    }
	});//end of parallel
      }
    else
      {
        if (task_manager)
          {
            FlatVector<> vx = x.FV<double> (); 
            FlatVector<> vy = y.FV<double> (); 

            task_manager -> CreateJob 
              ( [&] (const TaskInfo & ti)
                {
                  ArrayMem<double, 100> mem1(max_row_size);
                  ArrayMem<double, 100> mem2(max_col_size);
                  
                  auto myr = Range(coldnums).Split (ti.task_nr, ti.ntasks);
                  
                  for (int i : myr)
                    {
                      FlatArray<int> rdi (rowdnums[i]);
                      FlatArray<int> cdi (coldnums[i]);
                
                      if (!rdi.Size() || !cdi.Size()) continue;
                      if (rdi[0] == -1 || cdi[0] == -1) continue;  // reserved but not used
                      
                      FlatVector<> hv1(rdi.Size(), &mem1[0]);
                      FlatVector<> hv2(cdi.Size(), &mem2[0]);
                      
                      hv1 = vx(rdi);
                      hv2 = Trans(elmats[i]) * hv1;

                      for (int j = 0; j < cdi.Size(); j++)
                        {
                          double &res = vy(cdi[j]);
                          double t = s * hv2(j);
                          AtomicAdd(res, t);
                        }
                      
                      timer.AddFlops (cdi.Size()*rdi.Size());
                    }
                });
          }
        
        else
          {
            ArrayMem<double, 100> mem1(maxs);
            
            FlatVector<> vx = x.FV<double> (); 
            FlatVector<> vy = y.FV<double> (); 
            
            for (int i = 0; i < coldnums.Size(); i++) //sum over all elements
              {
                FlatArray<int> rdi (rowdnums[i]);
                FlatArray<int> cdi (coldnums[i]);
                
                if (!rdi.Size() || !cdi.Size()) continue;
                if (rdi[0] == -1 || cdi[0] == -1) continue;  // reserved but not used
                
                FlatVector<> hv1(rdi.Size(), &mem1[0]);
                
                hv1 = vx(rdi);
                vy(cdi) += s * Trans(elmats[i]) * hv1;
                
                timer.AddFlops (cdi.Size()*rdi.Size());
              }
          }
      }
  }










  template <class SCAL>
  shared_ptr<BaseMatrix> ElementByElementMatrix<SCAL> :: InverseMatrix ( BitArray * subset ) const
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
                                                         FlatArray<int> rowdnums_in,
                                                         FlatArray<int> coldnums_in,
                                                         BareSliceMatrix<SCAL> elmat)
  {
    if (elnr > elmats.Size())
      throw Exception ("EBEMatrix::AddElementMatrix, illegal elnr");
    
    
    ArrayMem<int,50> usedrows;
    for (int i = 0; i < rowdnums_in.Size(); i++)
      if (rowdnums_in[i] >= 0) usedrows.Append(i);
    int sr = usedrows.Size();

    ArrayMem<int,50> usedcols;
    for (int i = 0; i < coldnums_in.Size(); i++)
      if (coldnums_in[i] >= 0) usedcols.Append(i);
    int sc = usedcols.Size();

    if (allvalues.Size())
      {
        FlatMatrix<SCAL> mat(elmats[elnr]);
        FlatArray<int> dnr(rowdnums[elnr]);
        FlatArray<int> dnc(coldnums[elnr]);

        if (dnr.Size() != sr || dnc.Size() != sc || mat.Height() != sr || mat.Width() != sc)
          {
            throw Exception ("ebe, dnr or dnc has illegal size: \n"
                             "dnr.size = "+ToString(dnr.Size()) + " sr = " +ToString(sr) + "\n"
                             "dnc.size = "+ToString(dnc.Size()) + " sc = " +ToString(sc));
          }

        for (int i = 0; i < sr; i++)
          for (int j = 0; j < sc; j++)
            mat(i,j) = elmat(usedrows[i], usedcols[j]);

        for (int i = 0; i < sr; i++)
          dnr[i] = rowdnums_in[usedrows[i]];

        for (int j = 0; j < sc; j++)
          dnc[j] = coldnums_in[usedcols[j]];
      }
    else
      {
        FlatMatrix<SCAL> mat (sr,sc, new SCAL[sr*sc]);
        
        for (int i = 0; i < sr; i++)
          for (int j = 0; j < sc; j++)
            mat(i,j) = elmat(usedrows[i], usedcols[j]);
        
        FlatArray<int> dnr(sr, new int[sr]);
        for (int i = 0; i < sr; i++)
          dnr[i] = rowdnums_in[usedrows[i]];
        
        FlatArray<int> dnc(sc, new int[sc]);
        for (int j = 0; j < sc; j++)
          dnc[j] = coldnums_in[usedcols[j]];
        
        new (&rowdnums[elnr]) FlatArray<int> (dnr);
        new (&coldnums[elnr]) FlatArray<int> (dnc);
        elmats[elnr].AssignMemory (sr, sc, &mat(0,0));
      }

    
    max_row_size = max2(max_row_size, sr);
    max_col_size = max2(max_col_size, sc);
  }

  template <class SCAL>
  void ElementByElementMatrix<SCAL> :: AddCloneElementMatrix (int elnr,
                                                         const FlatArray<int> & rowdnums_in,
                                                         const FlatArray<int> & coldnums_in,
                                                         int refelnr)
  {
    if (allvalues.Size())
      throw Exception ("AddClone + allvalues not ready");

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
	clone.SetBitAtomic(elnr);
      }
    else
      throw Exception ("EBEMatrix::AddCloneElementMatrix, illegal elnr");
    
  }

  template<class SCAL>
  void ElementByElementMatrix<SCAL>::InitMemoryTracing() const
  {
    GetMemoryTracer().Track(allrow, "allrow", allcol, "allcol",
                            allvalues, "allvalues", elmats, "elmats",
                            rowdnums, "rowdnums", coldnums, "coldnums",
                            clone, "clone");
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

  
  ConstantElementByElementMatrix ::
  ConstantElementByElementMatrix (size_t ah, size_t aw, Matrix<> amatrix,
                                  Table<int> acol_dnums, Table<int> arow_dnums)
  : h(ah), w(aw), matrix(amatrix), col_dnums(move(acol_dnums)), row_dnums(move(arow_dnums))
  {
    disjoint_cols = true;
    disjoint_rows = true;

    BitArray used_col(h);
    used_col.Clear();
    for (auto col : col_dnums)
      for (auto d : col)
        {
          if (used_col.Test(d)) disjoint_cols = false;
          used_col.SetBit(d);
        }

    BitArray used_row(w);
    used_row.Clear();
    for (auto row : row_dnums)
      for (auto d : row)
        {
          if (used_row.Test(d)) disjoint_rows = false;
          used_row.SetBit(d);
        }
    
    // cout << "disjoint_rows = " << disjoint_rows << ", disjoint_cols = " << disjoint_cols << endl;



    if (!disjoint_rows)
      {
        Array<MyMutex> locks(w);
        size_t nblocks = row_dnums.Size();
        Array<int> col(nblocks);
        col = -1;

        int maxcolor = 0;
        int basecol = 0;
        Array<unsigned int> mask(w);

        atomic<int> found(0);
        size_t cnt = row_dnums.Size();

        while (found < cnt)
          {
            ParallelForRange
              (mask.Size(),
               [&] (IntRange myrange) { mask[myrange] = 0; });

            ParallelForRange
              (nblocks, [&] (IntRange myrange)
               {
                 Array<size_t> dofs;
                 size_t myfound = 0;
                 
                 for (size_t nr : myrange)
                   {
                     if (col[nr] >= 0) continue;
                     
                     unsigned check = 0;
                     dofs = row_dnums[nr];
                     
                     QuickSort (dofs);   // sort to avoid dead-locks
                     
                     for (auto d : dofs) 
                       locks[d].lock();
                     
                     for (auto d : dofs) 
                       check |= mask[d];
                     
                     if (check != UINT_MAX) // 0xFFFFFFFF)
                       {
                         myfound++;
                         unsigned checkbit = 1;
                         int color = basecol;
                         while (check & checkbit)
                           {
                             color++;
                             checkbit *= 2;
                           }
                         
                         col[nr] = color;
                         if (color > maxcolor) maxcolor = color;
                         
                         for (auto d : dofs) 
                           mask[d] |= checkbit;
                       }
                     
                     for (auto d : dofs) 
                       locks[d].unlock();
                   }
                 found += myfound;
               });
            
            basecol += 8*sizeof(unsigned int); // 32;
          }

        Array<int> cntcol(maxcolor+1);
        cntcol = 0;
        
        for (auto nr : Range(nblocks))
          cntcol[col[nr]]++;
        // cout << "cntcol = " << cntcol << endl;
        row_coloring = Table<int> (cntcol);

	cntcol = 0;
        for (auto nr : Range(nblocks))        
          row_coloring[col[nr]][cntcol[col[nr]]++] = nr;
      }



    if (!disjoint_cols)
      {
        Array<MyMutex> locks(h);
        size_t nblocks = row_dnums.Size();
        Array<int> col(nblocks);
        col = -1;

        int maxcolor = 0;
        int basecol = 0;
        Array<unsigned int> mask(h);

        atomic<int> found(0);
        size_t cnt = row_dnums.Size();

        while (found < cnt)
          {
            ParallelForRange
              (mask.Size(),
               [&] (IntRange myrange) { mask[myrange] = 0; });

            ParallelForRange
              (nblocks, [&] (IntRange myrange)
               {
                 Array<size_t> dofs;
                 size_t myfound = 0;
                 
                 for (size_t nr : myrange)
                   {
                     if (col[nr] >= 0) continue;
                     
                     unsigned check = 0;
                     dofs = col_dnums[nr];
                     
                     QuickSort (dofs);   // sort to avoid dead-locks
                     
                     for (auto d : dofs) 
                       locks[d].lock();
                     
                     for (auto d : dofs) 
                       check |= mask[d];
                     
                     if (check != UINT_MAX) // 0xFFFFFFFF)
                       {
                         myfound++;
                         unsigned checkbit = 1;
                         int color = basecol;
                         while (check & checkbit)
                           {
                             color++;
                             checkbit *= 2;
                           }
                         
                         col[nr] = color;
                         if (color > maxcolor) maxcolor = color;
                         
                         for (auto d : dofs) 
                           mask[d] |= checkbit;
                       }
                     
                     for (auto d : dofs) 
                       locks[d].unlock();
                   }
                 found += myfound;
               });
            
            basecol += 8*sizeof(unsigned int); // 32;
          }

        Array<int> cntcol(maxcolor+1);
        cntcol = 0;
        
        for (auto nr : Range(nblocks))
          cntcol[col[nr]]++;
        col_coloring = Table<int> (cntcol);

	cntcol = 0;
        for (auto nr : Range(nblocks))        
          col_coloring[col[nr]][cntcol[col[nr]]++] = nr;
      }

    
  }

  AutoVector ConstantElementByElementMatrix :: CreateRowVector () const
  {
    return make_unique<VVector<>> (w);
  }
  
  AutoVector ConstantElementByElementMatrix :: CreateColVector () const 
  {
    return make_unique<VVector<>> (h);
  }

  
  void ConstantElementByElementMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("ConstantEBE mult add");
    static Timer tcol("ConstantEBE mult coloring");
    static Timer tpmult("ConstantEBE mult parallel mult");

    auto fx = x.FV<double>();
    auto fy = y.FV<double>();

    if (!disjoint_cols)
      {
        /*
        RegionTimer reg(ts);
        Vector<> hx(matrix.Width());
        Vector<> hy(matrix.Height());
        for (size_t i = 0; i < row_dnums.Size(); i++)
          {
            hx = fx(row_dnums[i]);
            hy = matrix * hx;
            fy(col_dnums[i]) += s * hy;
          }
        */
        RegionTimer reg(tcol);

        for (auto col : col_coloring)
          ParallelForRange
            (col.Size(), [&] (IntRange r)
             {
               constexpr size_t BS = 128;
               Matrix<> hx(BS, matrix.Width());
               Matrix<> hy(BS, matrix.Height());
               
               for (size_t bi = r.First(); bi < r.Next(); bi+= BS)
                 {
                   size_t li = min2(bi+BS, r.Next());
                   size_t num = li-bi;
                   
                   for (size_t i = 0; i < num; i++)
                     hx.Row(i) = fx(row_dnums[col[bi+i]]);
                   
                   {
                     // RegionTracer rt(TaskManager::GetThreadId(), tpmult);
                     hy.Rows(0, num) = hx.Rows(0, num) * Trans(matrix);
                   }
                   
                   for (size_t i = 0; i < num; i++)
                     fy(col_dnums[col[bi+i]]) += s * hy.Row(i);
                 }
             });
      }
    else
      {
        RegionTimer reg(t);
        ParallelForRange
          (row_dnums.Size(), [&] (IntRange r)
           {
             /*
             Vector<> hx(matrix.Width());
             Vector<> hy(matrix.Height());
             for (auto i : r)
               {
                 hx = fx(row_dnums[i]);
                 hy = matrix * hx;
                 fy(col_dnums[i]) += s * hy;
               }
             */
             constexpr size_t BS = 128;
             Matrix<> hx(BS, matrix.Width());
             Matrix<> hy(BS, matrix.Height());

             for (size_t bi = r.First(); bi < r.Next(); bi+= BS)
               {
                 size_t li = min2(bi+BS, r.Next());
                 size_t num = li-bi; 
                 for (size_t i = 0; i < num; i++)
                   hx.Row(i) = fx(row_dnums[bi+i]);
                 {
                   // NgProfiler::AddThreadFlops(tpmult, TaskManager::GetThreadId(), num*matrix.Height()*matrix.Width());
                   // RegionTimer reg(tpmult);
                   // RegionTracer rt(TaskManager::GetThreadId(), tpmult);
                   hy.Rows(0, num) = hx.Rows(0, num) * Trans(matrix);
                 }
                 for (size_t i = 0; i < num; i++)
                   fy(col_dnums[bi+i]) += s * hy.Row(i);
               }
           });
      }
  }
  
  void ConstantElementByElementMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("ConstantEBE mult trans add");
    static Timer tcol("ConstantEBE mult trans coloring");
    static Timer tpmult("ConstantEBE mult trans mult");    

    auto fx = x.FV<double>();
    auto fy = y.FV<double>();
    
    if (!disjoint_rows)
      { // use coloring
        RegionTimer reg(tcol);

        for (auto col : row_coloring)
          ParallelForRange
            (col.Size(), [&] (IntRange r)
             {
               constexpr size_t BS = 128;
               Matrix<> hx(BS, matrix.Height());
               Matrix<> hy(BS, matrix.Width());
               
               for (size_t bi = r.First(); bi < r.Next(); bi+= BS)
                 {
                   size_t li = min2(bi+BS, r.Next());
                   size_t num = li-bi;
                   
                   for (size_t i = 0; i < num; i++)
                     hx.Row(i) = fx(col_dnums[col[bi+i]]);
                   
                   {
                     // NgProfiler::AddThreadFlops(tpmult, TaskManager::GetThreadId(), num*matrix.Height()*matrix.Width());
                     // RegionTimer reg(tpmult);
                     // RegionTracer rt(TaskManager::GetThreadId(), tpmult);
                     hy.Rows(0, num) = hx.Rows(0, num) * matrix;
                   }
                   
                   for (size_t i = 0; i < num; i++)
                     fy(row_dnums[col[bi+i]]) += s * hy.Row(i);
                 }
             });
      }
    else
      {
        RegionTimer reg(t);
        ParallelForRange
          (row_dnums.Size(), [&] (IntRange r)
           {
             constexpr size_t BS = 128;
             Matrix<> hx(BS, matrix.Height());
             Matrix<> hy(BS, matrix.Width());

             for (size_t bi = r.First(); bi < r.Next(); bi+= BS)
               {
                 size_t li = min2(bi+BS, r.Next());
                 size_t num = li-bi;
                 
                 for (size_t i = 0; i < num; i++)
                   hx.Row(i) = fx(col_dnums[bi+i]);
                 
                 {
                   // NgProfiler::AddThreadFlops(tpmult, TaskManager::GetThreadId(), num*matrix.Height()*matrix.Width());
                   // RegionTimer reg(tpmult);
                   // RegionTracer rt(TaskManager::GetThreadId(), tpmult);
                   
                   hy.Rows(0, num) = hx.Rows(0, num) * matrix;
                 }
                 
                 for (size_t i = 0; i < num; i++)
                   fy(row_dnums[bi+i]) += s * hy.Row(i);
               }
           });
      }
  }

  void StructuredElementByElementMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    auto hx = x.FV<double>().AsMatrix(num, matrix.Width());
    auto hy = y.FV<double>().AsMatrix(num, matrix.Height());
    Matrix tmp = s * Trans(matrix);
    ParallelForRange (hx.Height(), [&](IntRange myrange) {
        hy.Rows(myrange) += hx.Rows(myrange) * tmp; });
  }
  
  void StructuredElementByElementMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    auto hx = x.FV<double>().AsMatrix(num, matrix.Height());
    auto hy = y.FV<double>().AsMatrix(num, matrix.Width());
    Matrix tmp = s * matrix;
    ParallelForRange (hx.Height(), [&](IntRange myrange) {
        hy.Rows(myrange) += hx.Rows(myrange) * tmp; });
  }
}
