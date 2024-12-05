/*********************************************************************/
/* File:   sparsematrix_dyn.cpp                                      */
/* Author: Joachim Schoeberl                                         */
/* Date:   July 2019                                                 */
/*********************************************************************/


#include <sparsematrix_dyn.hpp>

#include "pardisoinverse.hpp"
#include "umfpackinverse.hpp"
#include "superluinverse.hpp"
#include "mumpsinverse.hpp"


namespace ngla
{

  template <typename TSCAL>
  void SparseMatrixDynamic<TSCAL> :: Mult (const BaseVector & x, BaseVector & y) const 
  {
    y = 0.0;
    MultAdd (1, x, y);
    /*    
    auto fx = x.FV<TSCAL>();
    auto fy = y.FV<TSCAL>();
    auto matvecfunc = dispatch_addmatvec[bw];
    fy = TSCAL(0.0);
    for (size_t i = 0; i < Height(); i++)
      {
        auto rowind = GetRowIndices(i);
        TSCAL * pmat = &data[bs*firsti[i]];
        size_t my_bw = bw;
        size_t my_bh = bh;
        size_t my_bs = bs;

        FlatVector<TSCAL> yi(my_bh, &fy(i*my_bh));
        for (auto j : rowind)
          {
            FlatVector<TSCAL> xi(my_bw, &fx(j*bw));
            FlatMatrix<TSCAL> mi(my_bh, my_bw, pmat);
            // yi += mi * xi;
            // yi = mi * xi;
            (*matvecfunc) (1, mi, xi, yi);
            pmat += my_bs;
          }
      }
    */
  }
  
  template <typename TSCAL>
  void SparseMatrixDynamic<TSCAL> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const 
  {
    // for (size_t i = 0; i < Height(); i++)
    ParallelForRange
      (Height(), [&] (IntRange r)
       {
         auto fx = x.FV<TSCAL>();
         auto fy = y.FV<TSCAL>();
         auto matvecfunc = dispatch_addmatvec[bw];

         size_t my_bw = bw;
         size_t my_bh = bh;
         size_t my_bs = bs;
         double my_s = s;

         TSCAL * pmat = &data[bs*firsti[r.First()]];
         TSCAL * py = &fy(r.First()*my_bh);
         for (auto i : r)
           {
             auto rowind = GetRowIndices(i);
             FlatVector<TSCAL> yi(my_bh, py); 
             for (auto j : rowind)
               {
                 FlatVector<TSCAL> xi(my_bw, &fx(j*my_bw));
                 FlatMatrix<TSCAL> mi(my_bh, my_bw, pmat);
                 // yi += s * mi * xi;
                 (*matvecfunc) (my_s, mi, xi, yi);            
                 pmat += my_bs;
               }
             py += my_bh;
           }
       });
  }
  

  template class SparseMatrixDynamic<double>;

  template <typename TSCAL>
  SparseMatrixVariableBlocks<TSCAL> ::
  SparseMatrixVariableBlocks (const SparseMatrixTM<TSCAL> & mat)
    : height(mat.Height()), width(mat.Width())
  {
    size_t row = 0;
    nblocks = 0;
    firsti_colnr.Append(0);
    firsti_data.Append(0);
    cum_block_size.Append(0);
    while (row < mat.Height())
      {
        auto rowind = mat.GetRowIndices(row);
        // colnr.Append (rowind);
        for (auto ri : rowind) colnr.Append(ri);
        firsti_colnr.Append(colnr.Size());
        for (auto val : mat.GetRowValues(row))
          data.Append (val);
        nblocks++;
        row++;
        while (row < mat.Height() && 
               colnr.Range(firsti_colnr[nblocks-1], firsti_colnr[nblocks]) ==
               mat.GetRowIndices(row))
          {
            for (auto val : mat.GetRowValues(row))
              data.Append (val);
            row++;
          }

        firsti_data.Append(data.Size());
        cum_block_size.Append(row);

        size_t bw = rowind.Size();
        size_t bh = cum_block_size[nblocks]-cum_block_size[nblocks-1];
        TSCAL * pdata = &data[firsti_data[nblocks-1]];
        Matrix<TSCAL> tmp = FlatMatrix<TSCAL> (bh, bw, pdata);
        FlatMatrix<TSCAL> (bw, bh, pdata) = Trans(tmp);
      }

    int maxbs = 0;
    for (int i = 0; i < nblocks; i++)
      maxbs = cum_block_size[i+1]-cum_block_size[i];
    Array<int> nbs(maxbs+1);
    nbs = 0;
    for (int i = 0; i < nblocks; i++)
      nbs[cum_block_size[i+1]-cum_block_size[i]]++;
    for (int i = 0; i < nbs.Size(); i++)
      if (nbs[i])
        cout << "bs=" << i << ": " << nbs[i] << endl;
  }


  template <typename TSCAL>
  void SparseMatrixVariableBlocks<TSCAL> ::
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    auto fx = x.FV<TSCAL>();
    auto fy = y.FV<TSCAL>();
    // for (size_t i = 0; i < nblocks; i++)
    ParallelForRange
      (nblocks, [&] (IntRange myrange)
       {
         for (size_t i : myrange)
           {
             TSCAL * pdata = &data[firsti_data[i]];
             int firsti = firsti_colnr[i];
             int nexti = firsti_colnr[i+1];
             //int * pcol = &colnr[firsti];
             TSCAL * py = &fy(cum_block_size[i]);
             int bs = cum_block_size[i+1]-cum_block_size[i];
             FlatMatrix<TSCAL> mat(nexti-firsti, bs, pdata);
             auto ind = colnr.Range(firsti, nexti);
             MultAddMatTransVecIndirect (s, mat, fx, FlatVector<TSCAL>(bs, py), ind);
           }
       }, TasksPerThread(4));
        /*
        for (size_t j = 0; j < nexti-firsti; j++, pcol++)
          {
            TSCAL hx = fx[*pcol];
            for (size_t k = 0; k < bs; k++, pdata++)
              py[k] += s * *pdata * hx; 
          }
        */
  }

  template <typename TSCAL>  
  AutoVector SparseMatrixVariableBlocks<TSCAL> :: CreateRowVector () const
  {
    cout << "CreateRowVector, w = " << width << endl;
    return CreateBaseVector(width, false, 1);    
  }

  template <typename TSCAL>  
  AutoVector SparseMatrixVariableBlocks<TSCAL> :: CreateColVector () const
  {
    return CreateBaseVector(height, false, 1);        
  }

  template class SparseMatrixVariableBlocks<double>;  

}
