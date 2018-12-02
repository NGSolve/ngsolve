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
                      AsAtomic(res) += t;
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
                          AsAtomic(res) += t;
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
	clone.Set(elnr);
      }
    else
      throw Exception ("EBEMatrix::AddCloneElementMatrix, illegal elnr");
    
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

  
  void ConstantElementByElementMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    auto fx = x.FV<double>();
    auto fy = y.FV<double>();
    Vector<> hx(matrix.Width());
    Vector<> hy(matrix.Height());
    for (size_t i = 0; i < row_dnums.Size(); i++)
      {
        hx = fx(row_dnums[i]);
        hy = matrix * hx;
        fy(col_dnums[i]) += s * hy;
      }
  }
  
  void ConstantElementByElementMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    auto fx = x.FV<double>();
    auto fy = y.FV<double>();
    Vector<> hx(matrix.Height());
    Vector<> hy(matrix.Width());
    for (size_t i = 0; i < row_dnums.Size(); i++)
      {
        hx = fx(col_dnums[i]);
        hy = Trans(matrix) * hx;
        fy(row_dnums[i]) += s * hy;
      }
  }


  
}
