/*********************************************************************/
/* File:   special_matrix.cpp                                        */
/* Author: Joachim Schoeberl                                         */
/* Date:   14. Mar. 2002                                             */
/*********************************************************************/


#include <special_matrix.hpp>
#include "sparsematrix.hpp"

namespace ngla
{

  
  TransposeVector :: TransposeVector (int ah, int aw)
    : h(ah), w(aw)
  { ; } 

  ostream & TransposeVector :: Print (ostream & ost) const
  {
    ost << "TransposeVector, output shape is " << h << " x " << w << endl;
    return ost;
  }
    
  AutoVector TransposeVector :: CreateRowVector () const
  {
    return make_unique<VVector<double>> (h*w);
  }
  
  AutoVector TransposeVector :: CreateColVector () const
  {
    return make_unique<VVector<double>> (h*w);
  }
      
  void TransposeVector :: Mult (const BaseVector & x, BaseVector & y) const
  {
    y.FV<double>().AsMatrix(h,w) = Trans(x.FV<double>().AsMatrix(w,h));
  }

  void TransposeVector :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    y.FV<double>().AsMatrix(w,h) = Trans(x.FV<double>().AsMatrix(h,w));
  };

  
  
  
  void PermutationMatrix :: Mult (const BaseVector & x, BaseVector & y) const
  {
    auto fvx = x.FV<double>();
    auto fvy = y.FV<double>();
    for (size_t i = 0; i < ind.Size(); i++)
      fvy(i) = fvx(ind[i]);
  }
  
  void PermutationMatrix :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    auto fvx = x.FV<double>();
    auto fvy = y.FV<double>();
    y = 0;
    for (size_t i = 0; i < ind.Size(); i++)
      fvy(ind[i]) += fvx(i);
  }

  void PermutationMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    auto fvx = x.FV<double>();
    auto fvy = y.FV<double>();
    for (size_t i = 0; i < ind.Size(); i++)
      fvy(i) += s * fvx(ind[i]);
  }
  
  void PermutationMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    auto fvx = x.FV<double>();
    auto fvy = y.FV<double>();
    for (size_t i = 0; i < ind.Size(); i++)
      fvy(ind[i]) += s * fvx(i);
  }


  
  void Embedding :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("Embedding::Mult"); RegionTimer reg(t);
    y = 0.0;
    y.Range(range) = x;
  }
  
  void Embedding :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("Embedding::MultTrans"); RegionTimer reg(t);
    y = x.Range(range);    
  }
    
  void Embedding :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("Embedding::MultAdd"); RegionTimer reg(t);
    y.Range(range) += s*x;
  }
  
  void Embedding :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("Embedding::MultAddTrans"); RegionTimer reg(t);
    y += s*x.Range(range);
  }

  shared_ptr<BaseSparseMatrix> Embedding :: CreateSparseMatrix() const
  {
    Array<int> ai(range.Size()), aj(range.Size());

    for (int i = 0; i < range.Size(); i++)
      {
        ai[i] = range.First()+i;
        aj[i] = i;
      }

    if (IsComplex())
      {
        Array<Complex> vals(range.Size());
        vals = 1.;
        return SparseMatrix<Complex>::CreateFromCOO (ai, aj, vals, Height(), Width());
      }
    else
      {
        Array<double> vals(range.Size());
        vals = 1.;
        return SparseMatrix<double>::CreateFromCOO (ai, aj, vals, Height(), Width());
      }
  }

  
  
  BaseMatrix::OperatorInfo EmbeddedMatrix :: GetOperatorInfo () const
  {
    OperatorInfo info;
    info.name = "EmbeddedMatrix";
    info.height = Height();
    info.width = Width();
    info.childs += mat.get();
    return info;
  }

  
  void EmbeddedMatrix :: Mult (const BaseVector & x, BaseVector & y) const
  {
    if (Height() != y.Size()) throw Exception("Embedded matrix, h = "+ToString(Height())
                                              + " != y.size = " + ToString(y.Size()));
    if (range.Size() != mat->Height()) throw Exception("range mismatch");
    if (Width() != x.Size()) throw Exception("Embedded matrix, w = "+ToString(Width())
                                             + " != x.Size() = " + ToString(x.Size()));
    y = 0;
    y.Range(range) = *mat * x;
  }
  void EmbeddedMatrix :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    mat->MultTrans(x.Range(range), y);
  }

  void EmbeddedMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const 
  {
    y.Range(range) += s * (*mat) * x;
  }
  
  void EmbeddedMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const 
  {
    mat->MultTransAdd(s, x.Range(range), y);
  }


  shared_ptr<BaseSparseMatrix> EmbeddedMatrix :: CreateSparseMatrix() const
  {
    auto mat1 = mat->CreateSparseMatrix();
    mat1 -> EmbedHeight(range.First(), Height());
    return mat1;
  }


  void EmbeddingTranspose :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("EmbeddingTranspose::Mult"); RegionTimer reg(t);
    y = x.Range(range);    
  }
  
  void EmbeddingTranspose :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("EmbeddingTranspose::MultTrans"); RegionTimer reg(t);
    y = 0.0;
    y.Range(range) = x;
  }
    
  void EmbeddingTranspose :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("EmbeddingTranspose::MultAdd"); RegionTimer reg(t);
    y += s*x.Range(range);
  }
  
  void EmbeddingTranspose :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("EmbeddingTranspose::MultAddTrans"); RegionTimer reg(t);
    y.Range(range) += s*x;
  }


  shared_ptr<BaseSparseMatrix> EmbeddingTranspose :: CreateSparseMatrix() const
  {
    Array<int> ai(range.Size()), aj(range.Size());

    for (int i = 0; i < range.Size(); i++)
      {
        ai[i] = range.First()+i;
        aj[i] = i;
      }

    if (IsComplex())
      {
        Array<Complex> vals(range.Size());
        vals = 1.;
        return SparseMatrix<Complex>::CreateFromCOO (aj, ai, vals, Height(), Width());
      }
    else
      {
        Array<double> vals(range.Size());
        vals = 1.;
        return SparseMatrix<double>::CreateFromCOO (aj, ai, vals, Height(), Width());
      }
  }





  

  BaseMatrix::OperatorInfo EmbeddedTransposeMatrix :: GetOperatorInfo () const
  {
    OperatorInfo info;
    info.name = "EmbeddedTransposeMatrix";
    info.height = Height();
    info.width = Width();
    info.childs += mat.get();
    return info;
  }

  void EmbeddedTransposeMatrix :: Mult (const BaseVector & x, BaseVector & y) const
  {
    mat->Mult(x.Range(range), y);
  }
  
  void EmbeddedTransposeMatrix :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    y = 0;
    auto ry = y.Range(range);
    mat->MultTrans(x, ry);
  }

  void EmbeddedTransposeMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const 
  {
    mat->MultAdd (s, x.Range(range), y);
  }
  
  void EmbeddedTransposeMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const 
  {
    auto ry = y.Range(range);    
    mat->MultTransAdd(s, x, ry);
  }





  
  template <class TVR, class TVC>
  Real2ComplexMatrix<TVR,TVC> :: 
  Real2ComplexMatrix (shared_ptr<BaseMatrix> arealmatrix)
    : hx(0), hy(0)
  { 
    SetMatrix (arealmatrix); 
  }
  
  template <class TVR, class TVC>
  void Real2ComplexMatrix<TVR,TVC> ::
  SetMatrix (shared_ptr<BaseMatrix> arealmatrix)
  {
    realmatrix = arealmatrix;
    if (realmatrix)
      {
	hx.SetSize (realmatrix->Height());
	hy.SetSize (realmatrix->Width());
      }
  }

  template <class TVR, class TVC>
  AutoVector Real2ComplexMatrix<TVR,TVC> :: CreateRowVector() const
  {
    auto h = realmatrix->Width();
    return make_unique<VVector<TVC>> (h);
  }
    
  template <class TVR, class TVC>
  AutoVector Real2ComplexMatrix<TVR,TVC> :: CreateColVector() const
  {
    auto w = realmatrix->Width();
    return make_unique<VVector<TVC>> (w);
  }
    


  template <class TVR, class TVC>
  void  Real2ComplexMatrix<TVR,TVC> :: 
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    MultAdd (Complex(s), x, y);
  }

  template <class TVR, class TVC>
  void  Real2ComplexMatrix<TVR,TVC> :: 
  MultAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    try
      {
	const FlatVector<TVC> cx = x.FV<TVC> ();
	// dynamic_cast<const T_BaseVector<TVC>&> (x).FV();
	FlatVector<TVC> cy = y.FV<TVC> ();
	// dynamic_cast<T_BaseVector<TVC>&> (y).FV();

	VVector<TVR> & hhx = const_cast<VVector<TVR>&> (hx);
	VVector<TVR> & hhy = const_cast<VVector<TVR>&> (hy);

	int i, j;
	
	/*
	  for (i = 0; i < cx.Size(); i++)
	  hhx.FV()(i)(0) = cx(i)(0).real();
	*/
	for (i = 0; i < cx.Size(); i++)
	  for (j = 0; j < TVR::SIZE; j++)
	    hhx(i)(j) = cx(i)(j).real();
	realmatrix -> Mult (hhx, hhy);
	cy += s * hhy.FV();

	/*
	  for (i = 0; i < cx.Size(); i++)
	  hhx.FV()(i)(0) = cx(i)(0).imag();
	*/

	for (i = 0; i < cx.Size(); i++)
	  for (j = 0; j < TVR::SIZE; j++)
	    hhx(i)(j) = cx(i)(j).imag();

	realmatrix -> Mult (hhx, hhy);
	
	cy += (s*Complex(0,1)) * hhy.FV();
      }
    catch (Exception & e)
      {
	e.Append (string ("\nthrown by Complex2RealMatrix::MultAdd "));
	throw;
      }
  }


  template <>
  void  Real2ComplexMatrix<double,Complex> :: 
  MultAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    const FlatVector<Complex> cx =x.FV<Complex> ();
    // dynamic_cast<const T_BaseVector<Complex>&> (x).FV();
    FlatVector<Complex> cy =y.FV<Complex> ();
    // dynamic_cast<T_BaseVector<Complex>&> (y).FV();

    VVector<double> & hhx = const_cast<VVector<double>&> (hx);
    VVector<double> & hhy = const_cast<VVector<double>&> (hy);

    int i;

    /*
    for (i = 0; i < cx.Size(); i++)
      hhx.FV()(i)(0) = cx(i)(0).real();
    */
    for (i = 0; i < cx.Size(); i++)
      hhx(i) = cx(i).real();
    realmatrix -> Mult (hhx, hhy);
    cy += s * hhy.FV();

    /*
    for (i = 0; i < cx.Size(); i++)
      hhx.FV()(i)(0) = cx(i)(0).imag();
    */

    for (i = 0; i < cx.Size(); i++)
      hhx(i) = cx(i).imag();

    realmatrix -> Mult (hhx, hhy);

    cy += (s*Complex(0,1)) * hhy.FV();
  }





  template class Real2ComplexMatrix<double,Complex>;
  template class Real2ComplexMatrix<Vec<2,double>,Vec<2,Complex> >;
  template class Real2ComplexMatrix<Vec<3,double>,Vec<3,Complex> >;
  template class Real2ComplexMatrix<Vec<4,double>,Vec<4,Complex> >;








  ////////////////////////////////////////////////////////////////////////////////
  // added 08/19/2003

  template <class TVR>
  Sym2NonSymMatrix<TVR> :: 
  Sym2NonSymMatrix (const BaseMatrix * abasematrix)
    : hx(0), hy(0)
  { 
    SetMatrix (abasematrix); 
  }
  
  template <class TVR>
  void Sym2NonSymMatrix<TVR> ::
  SetMatrix (const BaseMatrix * abasematrix)
  {
    base = abasematrix;
    if (base)
      {
	hx.SetSize (base->Height());
	hy.SetSize (base->Width());
      }
  }    



  template <class TVR>
  void  Sym2NonSymMatrix<TVR> :: 
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    const FlatVector<TVR> cx = x.FV<TVR> ();
    // dynamic_cast<const T_BaseVector<TVR>&> (x).FV();
    FlatVector<TVR> cy = y.FV<TVR> ();
    // dynamic_cast<T_BaseVector<TVR>&> (y).FV();

    VVector<TVR> & hhx = const_cast<VVector<TVR>&> (hx);
    VVector<TVR> & hhy = const_cast<VVector<TVR>&> (hy);

    int i, j;

    // for each set of values (x_cos, x_sin),
    // replace them by (x_cos+x_sin,x_cos-x_sin)
    for (i = 0; i < cx.Size(); i++)
      for (j = 0; j < TVR::SIZE; j+=2)
	{
	  hhx(i)(j) = cx(i)(j) + cx(i)(j+1);
	  hhx(i)(j+1) = cx(i)(j) - cx(i)(j+1);
	}

    base -> Mult (hhx, hhy);
    cy -= (s/2) * hhy.FV(); // TODO: passt des minus eh?
  }



  template class Sym2NonSymMatrix<Vec<2,double> >;
  template class Sym2NonSymMatrix<Vec<4,double> >;
  template class Sym2NonSymMatrix<Vec<6,double> >;
  template class Sym2NonSymMatrix<Vec<8,double> >;







  ////////////////////////////////////////////////////////////////////////////////
  // added 09/02/2003

  template <class TVSMALL, class TVBIG>
  Small2BigNonSymMatrix<TVSMALL,TVBIG> :: 
  Small2BigNonSymMatrix (const BaseMatrix * abasematrix)
    : hx1(0), hx2(0), hy1(0), hy2(0)
  { 
    SetMatrix (abasematrix); 
  }
  
  template <class TVSMALL, class TVBIG>
  void Small2BigNonSymMatrix<TVSMALL,TVBIG> ::
  SetMatrix (const BaseMatrix * abasematrix)
  {
    base = abasematrix;
    if (base)
      {
	hx1.SetSize (base->Width());
	hx2.SetSize (base->Width());
	hy1.SetSize (base->Height());
	hy2.SetSize (base->Height());
      }
  }    



  template <class TVSMALL, class TVBIG>
  void  Small2BigNonSymMatrix<TVSMALL,TVBIG> :: 
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    const FlatVector<TVBIG> cx = x.FV<TVBIG> ();
    // dynamic_cast<const T_BaseVector<TVBIG>&> (x).FV();
    FlatVector<TVBIG> cy = y.FV<TVBIG> ();
    // dynamic_cast<T_BaseVector<TVBIG>&> (y).FV();

    VVector<TVSMALL> & hhx1 = const_cast<VVector<TVSMALL>&> (hx1);
    VVector<TVSMALL> & hhx2 = const_cast<VVector<TVSMALL>&> (hx2);
    VVector<TVSMALL> & hhy1 = const_cast<VVector<TVSMALL>&> (hy1);
    VVector<TVSMALL> & hhy2 = const_cast<VVector<TVSMALL>&> (hy2);

    int i, j;

    // for each set of values (x_cos, x_sin),
    // store (x_cos+x_sin) in hhx1,
    // (x_cos-x_sin) in hhx2
    for (i = 0; i < cx.Size(); i++)
      for (j = 0; j < TVSMALL::SIZE; j++)
	{
	  hhx1(i)(j) = cx(i)(2*j) + cx(i)(2*j+1);
	  hhx2(i)(j) = cx(i)(2*j) - cx(i)(2*j+1);
	}

    base -> Mult (hhx1, hhy1);
    base -> Mult (hhx2, hhy2);


    // TODO: passt des minus eh?
    for (i = 0; i < cx.Size(); i++)
      for (j = 0; j < TVSMALL::SIZE; j++)
	{
	  cy(i)(2*j) -= s/2 * hhy1(i)(j);
	  cy(i)(2*j+1) -= s/2 * hhy2(i)(j);
	}
  }




  template <>
  void  Small2BigNonSymMatrix<double, Vec<2,double> > :: 
  MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    const FlatVector<Vec<2,double> > cx = x.FV<Vec<2,double> > ();
    // dynamic_cast<const T_BaseVector<Vec<2,double> >&> (x).FV();
    FlatVector<Vec<2,double> > cy = y.FV< Vec<2,double> > ();
    // dynamic_cast<T_BaseVector<Vec<2,double> >&> (y).FV();

    VVector<double> & hhx1 = const_cast<VVector<double>&> (hx1);
    VVector<double> & hhx2 = const_cast<VVector<double>&> (hx2);
    VVector<double> & hhy1 = const_cast<VVector<double>&> (hy1);
    VVector<double> & hhy2 = const_cast<VVector<double>&> (hy2);

    // for each set of values (x_cos, x_sin),
    // store (x_cos+x_sin) in hhx1,
    // (x_cos-x_sin) in hhx2
    for (int i = 0; i < cx.Size(); i++)
      {
	hhx1(i) = cx(i)(0) + cx(i)(1);
	hhx2(i) = cx(i)(0) - cx(i)(1);
      }

    base -> Mult (hhx1, hhy1);
    base -> Mult (hhx2, hhy2);


    // TODO: passt des minus eh?
    for (int i = 0; i < cx.Size(); i++)
      {
	cy(i)(0) -= s/2 * hhy1(i);
	cy(i)(1) -= s/2 * hhy2(i);
      }
  }



  // C = C^T -> no extra implementation
  template <class TVSMALL, class TVBIG>
  void  Small2BigNonSymMatrix<TVSMALL,TVBIG> :: 
  MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    MultAdd(s,x,y);
  }


  template class Small2BigNonSymMatrix<double, Vec<2,double> >;
  template class Small2BigNonSymMatrix<Vec<2,double>, Vec<4,double> >;
  template class Small2BigNonSymMatrix<Vec<3,double>, Vec<6,double> >;
  template class Small2BigNonSymMatrix<Vec<4,double>, Vec<8,double> >;



  
  BlockMatrix :: BlockMatrix (const Array<Array<shared_ptr<BaseMatrix>>> & amats)
    : mats(amats)
  {
    h = mats.Size();
    w = (h > 0) ? mats[0].Size() : 0;
    // check if nr of blocks per row is consistent
    for (auto k:Range(h))
      if (mats[k].Size()!=h)
	throw Exception("Tried to construct a BlockMatrix with unequal nr. of blocks per row");
    // check if there is at least one block per row/col
    BitArray rowhas(h);
    rowhas.Clear();
    BitArray colhas(w);
    colhas.Clear();
    for (auto k:Range(h))
      for (auto j:Range(w))
	if (mats[k][j]!=nullptr) {
	  rowhas.SetBit(k);
	  colhas.SetBit(j);
	}
    if (rowhas.NumSet()!=h)
      throw Exception("BlockMatrix needs at least one block per row");
    if (colhas.NumSet()!=w)
      throw Exception("BlockMatrix needs at least one block per col");
    row_reps.SetSize(h);
    row_reps = nullptr;
    for (auto k:Range(h)) {
      size_t col = 0;
      while(row_reps[k]==nullptr) {
	if (mats[k][col++]!=nullptr) {
	  row_reps[k] = mats[k][col-1];
	}
      }
    }
    col_reps.SetSize(w);
    col_reps = nullptr;
    for (auto k:Range(w)) {
      size_t row = 0;
      while (col_reps[k]==nullptr) {
	if (mats[row++][k]!=nullptr) {
	  col_reps[k] = mats[row-1][k];
	}
      }
    }
  }
  
  void BlockMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    auto & bvx = dynamic_cast_BlockVector(x);
    auto & bvy = dynamic_cast_BlockVector(y);
    for (size_t i = 0; i < h; i++)
      for (size_t j = 0; j < w; j++)
        {
          auto & spmat = mats[i][j];
          if (spmat)
              spmat->MultAdd(s, *bvx[j], *bvy[i]);
        }
  }

  void BlockMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    auto & bvx = dynamic_cast_BlockVector(x);
    auto & bvy = dynamic_cast_BlockVector(y);
    for (size_t i = 0; i < h; i++)
      for (size_t j = 0; j < w; j++)
        {
          auto & spmat = mats[i][j];
          if (spmat)
              spmat->MultTransAdd(s, *bvx[i], *bvy[j]);
        }
  }

  AutoVector BlockMatrix :: CreateRowVector () const {
    Array<shared_ptr<BaseVector>> vecs(w);
    for(auto col:Range(w)) {
      vecs[col] = col_reps[col]->CreateRowVector();
    }
    return make_unique<BlockVector>(vecs);
  }
  
  AutoVector BlockMatrix :: CreateColVector () const {
    Array<shared_ptr<BaseVector>> vecs(h);
    for (auto row:Range(h)) {
      vecs[row] = row_reps[row]->CreateColVector();
    }
    return make_unique<BlockVector>(vecs);
  }




  
  BaseMatrixFromVector :: BaseMatrixFromVector (shared_ptr<BaseVector> avec)
    : vec(avec) { }
  
  void BaseMatrixFromVector :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    y += s * x.FV<double>()(0) * (*vec);
  }
    
  void BaseMatrixFromVector :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    y.FV<double>()(0) += s * InnerProduct(x, *vec);
  }
  
  AutoVector BaseMatrixFromVector :: CreateRowVector () const
  {
    // missing parallel: 1 dof for all
    shared_ptr<BaseVector> sp = make_shared<VVector<double>>(1);   
    return sp;
  }
  
  AutoVector BaseMatrixFromVector :: CreateColVector () const
  {
    return vec->CreateVector();
  }





  
  BaseMatrixFromMultiVector :: BaseMatrixFromMultiVector (shared_ptr<MultiVector> avec)
    : vec(avec) { }
  
  void BaseMatrixFromMultiVector :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    // y += s * x.FV<double>()(0) * (*vec);
    Vector<> tmp = x.FV<double>();
    tmp *= s;
    vec->AddTo(tmp, y);
  }
    
  void BaseMatrixFromMultiVector :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    Vector<> tmp = vec->InnerProductD(x);
    y.FV<double>() += s *  tmp;
  }
  
  AutoVector BaseMatrixFromMultiVector :: CreateRowVector () const
  {
    // missing parallel: 1 dof for all
    shared_ptr<BaseVector> sp = make_shared<VVector<double>>(vec->Size());   
    return sp;
  }
  
  AutoVector BaseMatrixFromMultiVector :: CreateColVector () const
  {
    return vec->RefVec()->CreateVector();
  }



  template <typename T>
  BaseMatrixFromMatrix<T> :: BaseMatrixFromMatrix (Matrix<T> amat)
    : mat(std::move(amat)) { }

  template <typename T>  
  void BaseMatrixFromMatrix<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    y.FV<T>() += s * mat * x.FV<T>();    
  }

  template <typename T>    
  void BaseMatrixFromMatrix<T> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    y.FV<T>() += s * Trans(mat) * x.FV<T>();
  }

  template <typename T>  
  void BaseMatrixFromMatrix<T> :: MultAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    y.FV<Complex>() += s * mat * x.FV<Complex>();    
  }

  template <typename T>    
  void BaseMatrixFromMatrix<T> :: MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    y.FV<Complex>() += s * Trans(mat) * x.FV<Complex>();
  }



  
  template <typename T>    
  AutoVector BaseMatrixFromMatrix<T> :: CreateRowVector () const
  {
    // missing parallel: 1 dof for all
    shared_ptr<BaseVector> sp = make_shared<VVector<T>>(mat.Width());   
    return sp;
  }

  template <typename T>    
  AutoVector BaseMatrixFromMatrix<T> :: CreateColVector () const
  {
    shared_ptr<BaseVector> sp = make_shared<VVector<T>>(mat.Height());   
    return sp;
  }

  template class BaseMatrixFromMatrix<double>;
  template class BaseMatrixFromMatrix<Complex>;




  
  
  string PS(PARALLEL_STATUS stat)
  {
    switch (stat)
      {
      case DISTRIBUTED: return "distributed";
      case CUMULATED: return "cumulated";
      default: return "sequential";
      }
  }
  LoggingMatrix :: LoggingMatrix (shared_ptr<BaseMatrix> amat, string alabel, string logfile,
                                  optional<NgMPI_Comm> acomm)
    : mat(amat), label(alabel), comm(acomm)
  {
    if(logfile=="stdout")
      out = make_unique<ostream>(cout.rdbuf());
    else if(logfile=="stderr")
      out = make_unique<ostream>(cerr.rdbuf());
    else
      out = make_unique<ofstream>(logfile);
  }
      
  LoggingMatrix :: ~LoggingMatrix () { } 
  
  BaseVector & LoggingMatrix :: AsVector()
  {
    *out << "matrix '" << label << "' AsVector called" << endl;
    return mat->AsVector();
  }
  
  const BaseVector & LoggingMatrix :: AsVector() const
  {
    *out << "matrix '" << label << "' AsVector called" << endl;
    return mat->AsVector();
  }
  
  void LoggingMatrix :: SetZero()
  {
    mat->SetZero();
  }

  AutoVector LoggingMatrix :: CreateRowVector () const 
  {
    auto vec = mat->CreateRowVector();
    *out << "matrix '" << label << "' CreateRowVector "
         << "size: " << vec.Size() << " " << PS(vec.GetParallelStatus()) << endl;
    return vec;
  }
  
  AutoVector LoggingMatrix :: CreateColVector () const
  {
    auto vec = mat->CreateColVector();
    *out << "matrix '" << label << "' CreateColVector "
         << "size: " << vec.Size() << " " << PS(vec.GetParallelStatus()) << endl;
    return vec;
  }

  
  void LoggingMatrix :: Mult (const BaseVector & x, BaseVector & y) const
  {
    if (comm.has_value()) comm->Barrier();
    const BaseMatrix & rmat = *mat;
    *out << "matrix '" << label << "' Mult: " << typeid(rmat).name() << " "
         << "x: " << x.Size() << " " << PS(x.GetParallelStatus()) << " "
         << "y: " << y.Size() << " " << PS(y.GetParallelStatus()) << endl;
    if (comm.has_value()) comm->Barrier();
    
    mat->Mult(x,y);
    
    if (comm.has_value()) comm->Barrier();
    *out << "matrix '" << label << "' Mult complete" << endl;
    if (comm.has_value()) comm->Barrier();        
  }
  
  void LoggingMatrix :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    mat->MultTrans (x,y);
  }
  void LoggingMatrix :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    if (comm.has_value()) comm->Barrier();    
    *out << "matrix '" << label << "' MultAdd: "
         << "x: " << x.Size() << " " << PS(x.GetParallelStatus()) << " "
         << "y: " << y.Size() << " " << PS(y.GetParallelStatus()) << endl;
    if (comm.has_value()) comm->Barrier();
    
    mat->MultAdd(s,x,y);
    
    if (comm.has_value()) comm->Barrier();    
    *out << "matrix '" << label << "' MultAdd complete" << endl;
    if (comm.has_value()) comm->Barrier();        
  }
  void LoggingMatrix :: MultAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    mat->MultAdd(s,x,y);
  }
  void LoggingMatrix :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    mat->MultTransAdd (s,x,y);
  }
  void LoggingMatrix :: MultTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    mat->MultTransAdd (s,x,y);
  }
  void LoggingMatrix :: MultConjTransAdd (Complex s, const BaseVector & x, BaseVector & y) const
  {
    mat->MultConjTransAdd (s,x,y);
  }
  void LoggingMatrix :: MultAdd (FlatVector<double> alpha, const MultiVector & x, MultiVector & y) const
  {
    mat->MultAdd (alpha, x, y);
  }

  
}
