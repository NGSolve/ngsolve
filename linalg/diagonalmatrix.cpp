/* ************************************************************************/
/* File:   diagonalmatrix.cpp                                             */
/* Author: Joachim Schoeberl                                              */
/* Date:   14 Mar. 02                                                     */
/* ************************************************************************/


#include "diagonalmatrix.hpp"
#include "sparsematrix.hpp"

namespace ngla
{
  void Projector :: Mult (const BaseVector & x, BaseVector & y) const
  {
    y = x;
    Project (y);
  }
  
  void Projector :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    Mult (x, y);
  }


  void Projector :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("Projector::MultAdd"); RegionTimer reg(t);
    /*
    FlatSysVector<> sx = x.SV<double>();
    FlatSysVector<> sy = y.SV<double>();

    ParallelForRange
      (bits->Size(),
       [this, sx, sy, s] (IntRange myrange)
       {
         if (keep_values)
           {
             for (size_t i : myrange) //  Range(*bits))
               if ((*bits)[i])
                 sy(i) += s * sx(i);
           }
         else
           {
             for (size_t i : myrange) // Range(*bits))
               if (!(*bits)[i])
                 sy(i) += s * sx(i);
           }
       });
    */

    auto multadd = [this,s] (auto sx, auto sy)
      {
        ParallelForRange
        (bits->Size(),
         [sx, sy, s, this] (IntRange myrange)
            {
              if (keep_values)
                {
                  for (size_t i : myrange) 
                    if ((*bits)[i])
                      sy(i) += s * sx(i);
                }
              else
                {
                  for (size_t i : myrange) 
                    if (!(*bits)[i])
                      sy(i) += s * sx(i);
                }
            });
      };
    
    if (x.EntrySize() == 1)
      multadd (x.FV<double>(), y.FV<double>());
    else
      multadd (x.SV<double>(), y.SV<double>());
  }
  
  void Projector :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    return MultAdd (s, x, y);
  }

  void Projector :: Project (BaseVector & x) const
  {
    static Timer t("Projector::Project"); RegionTimer reg(t);
    /*
    FlatSysVector<> sx = x.SV<double>();
    
    ParallelForRange
      (bits->Size(),
       [this, sx] (IntRange myrange)
       {
         if (keep_values)
           {
             for (size_t i : myrange) // Range(*bits))
               if (!(*bits)[i])
                 sx(i) = 0.0;
           }
         else
           {
             for (size_t i : myrange) // Range(*bits))
               if ((*bits)[i])
                 sx(i) = 0.0;
           }
       });
    */

    auto project = [this] (BitArray & bits, auto sx)
      {
        ParallelForRange
        (bits.Size(),
         [&bits, sx, this] (IntRange myrange)
            {
              if (keep_values)
                {
                  for (size_t i : myrange) // Range(*bits))
                    if (!bits[i])
                      sx(i) = 0.0;
                }
              else
                {
                  for (size_t i : myrange) // Range(*bits))
                    if (bits[i])
                      sx(i) = 0.0;
                }
            });
      };

    if (x.EntrySize() == 1)
      project (*bits, x.FV<double>());
    else
      project (*bits, x.SV<double>());
  }
  

  void Projector :: SetValues (BaseVector & x, double value) const
  {
    auto setval = [this, value] (BitArray & bits, auto sx)
      {
        ParallelForRange
        (bits.Size(),
         [&bits, sx, this,value] (IntRange myrange)
            {
              if (keep_values)
                {
                  for (size_t i : myrange) 
                    if (bits[i])
                      sx(i) = value;
                }
              else
                {
                  for (size_t i : myrange) 
                    if (!bits[i])
                      sx(i) = value;
                }
            });
      };

    if (x.EntrySize() == 1)
      setval (*bits, x.FV<double>());
    else
      setval (*bits, x.SV<double>());
  }


  shared_ptr<BaseSparseMatrix> Projector :: CreateSparseMatrix() const
  {
    Array<int> indi(Height()), indj(Width());
    Array<double> vals(Height());
    for (int i : Range(Height()))
      {
        indi[i] = i;
        indj[i] = i;
      }
    auto mask = Mask();
    if (KeepValues())
      {
        vals = false;
        for (int i : Range(Height()))
          if ( (*mask)[i] ) vals[i] = true;
      }
    else
      {
        vals = false;
        for (int i : Range(Height()))
          if ( !(*mask)[i] ) vals[i] = true;
      }
    return SparseMatrix<double>::CreateFromCOO (indi, indj, vals, Height(), Height());           
  }


  template <typename TM>  
  DiagonalMatrix<TM> ::  DiagonalMatrix(const VVector<TM> & diag_)
    : diag(make_shared<VVector<TM>>(diag_))
  { } 
  
  template <typename TM>  
  DiagonalMatrix<TM> :: DiagonalMatrix(shared_ptr<VVector<TM>> diag_)
    : diag(diag_)
  { } 

  template <typename TM>  
  DiagonalMatrix<TM> :: ~DiagonalMatrix()
  { }
    

  template <typename TM>  
  shared_ptr<BaseSparseMatrix> DiagonalMatrix<TM> :: CreateSparseMatrix() const
  {
    Array<int> indi(Height()), indj(Width());
    Array<TM> vals(Height());
    for (int i : Range(Height()))
      {
        indi[i] = i;
        indj[i] = i;
        vals[i] = (*diag)(i);
      }
    return SparseMatrix<TM>::CreateFromCOO (indi, indj, vals, Height(), Height());           
  }

  


  template <typename TM>
  void DiagonalMatrix<TM> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DiagonalMatrix::MultAdd"); RegionTimer reg(t);    
    // if (mat_traits<TM>::WIDTH == x.EntrySize())
    if (ngbla::Width<TM>() == x.EntrySize())
      {
        typedef typename mat_traits<TM>::TV_ROW TV_ROW;
        typedef typename mat_traits<TM>::TV_COL TV_COL;
        
        auto sx = x.FV<TV_ROW>();
        auto sy = y.FV<TV_COL>();
        auto sd = this->diag->FV();
        ParallelForRange
          (Range(*diag), [sx,sy,sd,s] (IntRange myrange)
           {
             for (size_t i : myrange)
               sy(i) += s * sd(i)*sx(i);
               // sy(i) += s * this->diag(i)*sx(i);
           });
      }
    else
      {
        auto sx = x.SV<TSCAL>();
        auto sy = y.SV<TSCAL>();
        for (size_t i : Range(*diag))
          sy(i) += s * (*diag)(i)*sx(i);
      }
  }

  template <typename TM>  
  void DiagonalMatrix<TM> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    return MultAdd (s, x, y);
  }

  template <typename TM>  
  AutoVector DiagonalMatrix<TM> :: CreateRowVector () const 
  {
    // return CreateBaseVector(diag->Size(), mat_traits<TM>::IS_COMPLEX, mat_traits<TM>::WIDTH);
    return CreateBaseVector(diag->Size(), ngbla::IsComplex<TM>(), ngbla::Width<TM>());
  }

  template <typename TM>    
  AutoVector DiagonalMatrix<TM> :: CreateColVector () const 
  {
    return CreateBaseVector(diag->Size(), ngbla::IsComplex<TM>(), ngbla::Height<TM>());
  }

  template <typename TM>    
  shared_ptr<BaseMatrix> DiagonalMatrix<TM> ::
  InverseMatrix (shared_ptr<BitArray> subset) const
  {
    VVector<TM> v2(diag->Size());
    if (subset)
      {
        for (size_t i = 0; i < diag->Size(); i++)
          if (subset->Test(i))
            {
              v2(i) = (*diag)(i);
              CalcInverse(v2(i));
            }
          else
            v2(i) = TM(0.0);
      }
    else
      {
        for (size_t i = 0; i < diag->Size(); i++)
          {
            v2(i) = (*diag)(i);
            CalcInverse(v2(i));
          }
      }
    return make_shared<DiagonalMatrix<TM>> (v2);
  }

  template <typename TM>    
  ostream & DiagonalMatrix<TM> :: Print (ostream & ost) const 
  {
    return ost << diag;
  }


  template class DiagonalMatrix<double>;
  template class DiagonalMatrix<Complex>;


#if MAX_SYS_DIM >= 1
  template class DiagonalMatrix<Mat<1,1,double> >;
  template class DiagonalMatrix<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class DiagonalMatrix<Mat<2,2,double> >;
  template class DiagonalMatrix<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class DiagonalMatrix<Mat<3,3,double> >;
  template class DiagonalMatrix<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class DiagonalMatrix<Mat<4,4,double> >;
  template class DiagonalMatrix<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class DiagonalMatrix<Mat<5,5,double> >;
  template class DiagonalMatrix<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class DiagonalMatrix<Mat<6,6,double> >;
  template class DiagonalMatrix<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class DiagonalMatrix<Mat<7,7,double> >;
  template class DiagonalMatrix<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class DiagonalMatrix<Mat<8,8,double> >;
  template class DiagonalMatrix<Mat<8,8,Complex> >;
#endif
  

  template <typename T>  
  BlockDiagonalMatrix<T> ::
  BlockDiagonalMatrix(Tensor<3,T> _blockdiag)
    : blockdiag(std::move(_blockdiag))
  {
    blocks = blockdiag.GetSize();
    dimy = blockdiag.GetSubTensor().GetSize();
    dimx = blockdiag.GetSubTensor().GetSubTensor().GetSize();

    
    // check for blockdiag non-zero pattern:
    Matrix nonzero(dimy, dimx);
    for (int i = 0; i < dimy; i++)
      for (int j = 0; j < dimx; j++)
        nonzero(i,j) = L2Norm(blockdiag(STAR,i,j));

    // cout << "Blockdiag, shape = " << dimy << " x " << dimx << endl;
    // cout << "nonzero = " << endl << nonzero << endl;
    
    for (int sub = dimx; sub >= 2; sub--) // so many sub-blocks ?
      if ( (dimx % sub == 0) && (dimy % sub == 0) )
        {
          Matrix nz2 = nonzero;
          int subdimx = dimx/sub;
          int subdimy = dimy/sub;
          for (int i = 0; i < sub; i++)
            nz2.Rows(i*subdimy, (i+1)*subdimy).Cols(i*subdimx, (i+1)*subdimx) = 0.0;
          
          if (L2Norm(nz2) == 0)
            {
              cout << IM(3) << "can reduce subblocks by factor " << sub << endl;
              // cout << "nonzero = " << endl << nonzero << endl;
              
              Tensor<3,T> newblockdiag(blocks*sub, subdimy, subdimx);
              for (int i = 0; i < blocks; i++)
                for (int j = 0; j < sub; j++)
                  newblockdiag(i*sub+j,STAR,STAR) =
                    blockdiag(i,STAR,STAR).Rows(j*subdimy, (j+1)*subdimy).Cols(j*subdimx, (j+1)*subdimx);

              blocks *= sub;
              dimy = subdimy;
              dimx = subdimx;
              // blockdiag = std::move(newblockdiag);
              blockdiag.Assign (newblockdiag);
              newblockdiag.Data() = nullptr;
              return;
            }
        }

    cout << IM(5) << "Blockdiag, final shape = " << dimy << " x " << dimx << endl;    
  }
  
  template <typename T>  
  ostream & BlockDiagonalMatrix<T> :: Print (ostream & ost) const
  {
    ost << "BlockDiagmatrix, blocks = " << blocks << " of dim " << dimy << " x " << dimx << endl;
    return ost;
  }

  template <typename T>    
  AutoVector BlockDiagonalMatrix<T> :: CreateRowVector () const
  {
    return make_unique<VVector<T>>(VWidth());
  }

  template <typename T>      
  AutoVector BlockDiagonalMatrix<T> :: CreateColVector () const
  {
    return make_unique<VVector<T>>(VHeight());    
  }

  template <typename T>
  void BlockDiagonalMatrix<T> :: Mult (const BaseVector & x, BaseVector & y) const
  {
    // cout << "BlockDiagonalMult, dims = " << dimx << " x " << dimy << endl;
    static Timer t("BlockDiagonalMatrix::Mult"); RegionTimer r(t);


    if (dimx == 1 && dimy == 1)
      {
        auto fx = x.FV<T>();
        auto fy = y.FV<T>();
        ParallelFor
          (blocks, [&] (size_t i)
           {
             fy(i) = blockdiag(i,0,0) * fx(i);
           });
        return;
      }

    if (dimx == 2 && dimy == 2)
      {
        auto fx = x.FV<Vec<2,T>>();
        auto fy = y.FV<Vec<2,T>>();
        ParallelFor
          (blocks, [&] (size_t i)
           {
             fy(i) = Mat<2,2,T>(blockdiag(i,STAR,STAR)) * fx(i);
           });
        return;
      }
    
    auto fx = x.FV<T>();
    auto fy = y.FV<T>();
    ParallelFor
      (blocks, [&] (size_t i)
       {
         fy.Range(i*dimy, (i+1)*dimy) = blockdiag(i,STAR,STAR) * fx.Range(i*dimx, (i+1)*dimx);
       });
  }

  
  template <typename T>
  void BlockDiagonalMatrix<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("BlockDiagonalMatrix::MultAdd"); RegionTimer r(t);
    auto fx = x.FV<T>();
    auto fy = y.FV<T>();
    // for (size_t i = 0; i < blocks; i++)
    ParallelFor
      (blocks, [&] (size_t i)
       {
         fy.Range(i*dimy, (i+1)*dimy) += s* blockdiag(i,STAR,STAR) * fx.Range(i*dimx, (i+1)*dimx);
       });
  }

  template <typename T>  
  void BlockDiagonalMatrix<T> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("BlockDiagonalMatrix::MultTransAdd"); RegionTimer r(t);
    
    auto fx = x.FV<T>();
    auto fy = y.FV<T>();
    // for (size_t i = 0; i < blocks; i++)
    ParallelFor
      (blocks, [&] (size_t i)
       {
         fy.Range(i*dimx, (i+1)*dimx) += s* Trans(blockdiag(i,STAR,STAR)) * fx.Range(i*dimy, (i+1)*dimy);
       });
  }

  template <typename T>    
  shared_ptr<BaseMatrix> BlockDiagonalMatrix<T> :: InverseMatrix (shared_ptr<BitArray> subset) const
  {
    cout << "blockdiagmatrix, Inverse not implemented" << endl;
    return nullptr;
  }
  

  template class BlockDiagonalMatrix<double>;
  template class BlockDiagonalMatrix<Complex>;




  BlockDiagonalMatrixSoA ::
  BlockDiagonalMatrixSoA(Tensor<3> _blockdiag)
    : blockdiag(std::move(_blockdiag))
  {
    dimy = blockdiag.GetSize();
    dimx = blockdiag.GetSubTensor().GetSize();
    blocks = blockdiag.GetSubTensor().GetSubTensor().GetSize();

    // check for blockdiag non-zero pattern:
    nonzero.SetSize(dimy, dimx);
    for (int i = 0; i < dimy; i++)
      for (int j = 0; j < dimx; j++)
        nonzero(i,j) = L2Norm(blockdiag(i,j,STAR)) > 0;

    TableCreator<int> creator(dimy);
    for ( ; !creator.Done(); creator++)
      for (int i = 0; i < dimy; i++)
        for (int j = 0; j < dimx; j++)
          if (nonzero(i,j))
            creator.Add (i,j);
    sparse = creator.MoveTable();

    TableCreator<int> creatorT(dimx);
    for ( ; !creatorT.Done(); creatorT++)
      for (int i = 0; i < dimy; i++)
        for (int j = 0; j < dimx; j++)
          if (nonzero(i,j))
            creatorT.Add (j, i);
    sparseT = creatorT.MoveTable();

    // cout << "blockdiagSoA, nonzero = " << sparse.AsArray().Size() << " of " << dimx*dimy << endl;
  }

  ostream & BlockDiagonalMatrixSoA :: Print (ostream & ost) const
  {
    ost << "BlockDiagmatrixSoA, blocks = " << blocks << " of dim " << dimy << " x " << dimx << endl;
    return ost;
  }

  BaseMatrix::OperatorInfo BlockDiagonalMatrixSoA :: GetOperatorInfo () const
  {
    OperatorInfo info;
    info.name = string("BlockDiagonalMatrixSoA (bs = ") + ToString(dimy) + "x"
      + ToString (dimx) + ")";
    info.height = Height();
    info.width = Width();
    return info;
  }
    
  AutoVector BlockDiagonalMatrixSoA :: CreateRowVector () const
  {
    return make_unique<VVector<double>>(VWidth());
  }
  
  AutoVector BlockDiagonalMatrixSoA :: CreateColVector () const
  {
    return make_unique<VVector<double>>(VHeight());    
  }



  void BlockDiagonalMatrixSoA :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("BlockDiagonalMatrixSoA::Mult"); RegionTimer r(t);

    auto mx = x.FV<double>().AsMatrix(dimy, blocks);
    auto my = y.FV<double>().AsMatrix(dimx, blocks);

    int nthreads = max (blocks/1024, TasksPerThread(2));
    ParallelForRange
      (blocks, [blockdiag=FlatTensor(blockdiag),mx,my,this] (IntRange r)
       {
         for (size_t j = 0; j < sparse.Size(); j++)
           {
             switch (sparse[j].Size())
               {
               case 1:
                 {
                   my.Row(j).Range(r) = pw_mult(blockdiag(j, sparse[j][0],STAR).Range(r), mx.Row(sparse[j][0]).Range(r));
                   break;
                 }
               case 2:
                 {
                   my.Row(j).Range(r) =
                     pw_mult(blockdiag(j, sparse[j][0],STAR).Range(r), mx.Row(sparse[j][0]).Range(r))
                     +pw_mult(blockdiag(j, sparse[j][1],STAR).Range(r), mx.Row(sparse[j][1]).Range(r));
                   break;
                 }
               case 3:
                 {
                   my.Row(j).Range(r) =
                     pw_mult(blockdiag(j, sparse[j][0],STAR).Range(r), mx.Row(sparse[j][0]).Range(r))
                     +pw_mult(blockdiag(j, sparse[j][1],STAR).Range(r), mx.Row(sparse[j][1]).Range(r))
                     +pw_mult(blockdiag(j, sparse[j][2],STAR).Range(r), mx.Row(sparse[j][2]).Range(r));
                   break;
                 }
               default:
                 {
                   my.Row(j).Range(r) = 0.0;
                   for (size_t k : sparse[j])
                     my.Row(j).Range(r) += pw_mult(blockdiag(j,k,STAR).Range(r), mx.Row(k).Range(r));
                 }
               }
           }
       }, nthreads); // TasksPerThread(3));

    size_t numnz = 0;
    for (int i = 0; i < dimx; i++)
      for (int j = 0; j < dimy; j++)
        if (nonzero(j,i)) numnz++;
    t.AddFlops (blocks*numnz);
  }

  
  
  void BlockDiagonalMatrixSoA :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("BlockDiagonalMatrixSoA::MultAdd"); RegionTimer r(t);
    auto mx = x.FV<double>().AsMatrix(dimx, blocks);
    auto my = y.FV<double>().AsMatrix(dimy, blocks);
    // for (size_t i = 0; i < blocks; i++)
    ParallelForRange
      (blocks, [&] (IntRange r)
       {
         for (int i = 0; i < dimy; i++)
           for (int j = 0; j < dimx; j++)
             if (nonzero(i,j))
               my.Row(i).Range(r) += s* pw_mult (blockdiag(i,j,STAR).Range(r), mx.Row(j).Range(r));
       });

    size_t numnz = 0;
    for (int i = 0; i < dimx; i++)
      for (int j = 0; j < dimy; j++)
        if (nonzero(j,i)) numnz++;
    t.AddFlops (blocks*numnz);
  }


  void BlockDiagonalMatrixSoA :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("BlockDiagonalMatrixSoA::MultTrans"); RegionTimer r(t);

    auto mx = x.FV<double>().AsMatrix(dimy, blocks);
    auto my = y.FV<double>().AsMatrix(dimx, blocks);

    int nthreads = max (blocks/1024, TasksPerThread(2));
    ParallelForRange
      (blocks, [blockdiag=FlatTensor(blockdiag),mx,my,this] (IntRange r)
       {
         for (size_t j = 0; j < sparseT.Size(); j++)
           {
             switch (sparseT[j].Size())
               {
               case 0:
                 {
                   my.Row(j).Range(r) = 0.0;
                   break;
                 }
               case 1:
                 {
                   my.Row(j).Range(r) = pw_mult(blockdiag(sparseT[j][0],j,STAR).Range(r), mx.Row(sparseT[j][0]).Range(r));
                   break;
                 }
               case 2:
                 {
                   my.Row(j).Range(r) =
                     pw_mult(blockdiag(sparseT[j][0],j,STAR).Range(r), mx.Row(sparseT[j][0]).Range(r))
                     +pw_mult(blockdiag(sparseT[j][1],j,STAR).Range(r), mx.Row(sparseT[j][1]).Range(r));
                   break;
                 }
               case 3:
                 {
                   my.Row(j).Range(r) =
                     pw_mult(blockdiag(sparseT[j][0],j,STAR).Range(r), mx.Row(sparseT[j][0]).Range(r))
                     +pw_mult(blockdiag(sparseT[j][1],j,STAR).Range(r), mx.Row(sparseT[j][1]).Range(r))
                     +pw_mult(blockdiag(sparseT[j][2],j,STAR).Range(r), mx.Row(sparseT[j][2]).Range(r));
                   break;
                 }
               default:
                 {
                   my.Row(j).Range(r) =
                     pw_mult(blockdiag(sparseT[j][0],j,STAR).Range(r), mx.Row(sparseT[j][0]).Range(r))
                     +pw_mult(blockdiag(sparseT[j][1],j,STAR).Range(r), mx.Row(sparseT[j][1]).Range(r))
                     +pw_mult(blockdiag(sparseT[j][2],j,STAR).Range(r), mx.Row(sparseT[j][2]).Range(r))
                     +pw_mult(blockdiag(sparseT[j][3],j,STAR).Range(r), mx.Row(sparseT[j][3]).Range(r));
                   size_t k = 4;
                   for ( ; k+4 <= sparseT[j].Size(); k+=4) 
                     my.Row(j).Range(r) +=
                       pw_mult(blockdiag(sparseT[j][k+0],j,STAR).Range(r), mx.Row(sparseT[j][k+0]).Range(r))
                       +pw_mult(blockdiag(sparseT[j][k+1],j,STAR).Range(r), mx.Row(sparseT[j][k+1]).Range(r))
                       +pw_mult(blockdiag(sparseT[j][k+2],j,STAR).Range(r), mx.Row(sparseT[j][k+2]).Range(r))
                       +pw_mult(blockdiag(sparseT[j][k+3],j,STAR).Range(r), mx.Row(sparseT[j][k+3]).Range(r));
                   for ( ; k < sparseT[j].Size(); k++) // size_t k : sparseT[j])
                     my.Row(j).Range(r) += pw_mult(blockdiag(sparseT[j][k],j,STAR).Range(r), mx.Row(sparseT[j][k]).Range(r));
                 }
               }
           }
       }, nthreads); // TasksPerThread(3));

    size_t numnz = 0;
    for (int i = 0; i < dimx; i++)
      for (int j = 0; j < dimy; j++)
        if (nonzero(j,i)) numnz++;
    t.AddFlops (blocks*numnz);
  }


  
  void BlockDiagonalMatrixSoA :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("BlockDiagonalMatrixSoA::MultTransAdd"); RegionTimer r(t);

    auto mx = x.FV<double>().AsMatrix(dimy, blocks);
    auto my = y.FV<double>().AsMatrix(dimx, blocks);

    int nthreads = max (blocks/1024, TasksPerThread(2));
    ParallelForRange
      (blocks, [blockdiag=FlatTensor(blockdiag),mx,my,s,this] (IntRange r)
       {
         for (size_t i = 0; i < dimx; i++)
           for (size_t j = 0; j < dimy; j++)
             if (nonzero(j,i)) 
               my.Row(i).Range(r) += s* pw_mult(blockdiag(j,i,STAR).Range(r), mx.Row(j).Range(r));
       }, nthreads); // TasksPerThread(3));

    size_t numnz = 0;
    for (int i = 0; i < dimx; i++)
      for (int j = 0; j < dimy; j++)
        if (nonzero(j,i)) numnz++;
    t.AddFlops (blocks*numnz);
  }
  

}
