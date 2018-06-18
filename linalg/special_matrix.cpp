/*********************************************************************/
/* File:   special_matrix.cpp                                        */
/* Author: Joachim Schoeberl                                         */
/* Date:   14. Mar. 2002                                             */
/*********************************************************************/

/* 
   bilinear-form and linear-form integrators
*/

#include <la.hpp>
namespace ngla
{



  void Projector :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    FlatSysVector<> sx = x.SV<double>();
    FlatSysVector<> sy = y.SV<double>();

    if (keep_values)
      {
        for (size_t i : Range(*bits))
          if ((*bits)[i])
            sy(i) += s * sx(i);
      }
    else
      {
        for (size_t i : Range(*bits))
          if (!(*bits)[i])
            sy(i) += s * sx(i);
      }
  }



  template <class TVR, class TVC>
  Real2ComplexMatrix<TVR,TVC> :: 
  Real2ComplexMatrix (const BaseMatrix * arealmatrix)
    : hx(0), hy(0)
  { 
    SetMatrix (arealmatrix); 
  }
  
  template <class TVR, class TVC>
  void Real2ComplexMatrix<TVR,TVC> ::
  SetMatrix (const BaseMatrix * arealmatrix)
  {
    realmatrix = arealmatrix;
    if (realmatrix)
      {
	hx.SetSize (realmatrix->Height());
	hy.SetSize (realmatrix->Width());
      }
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
	  rowhas.Set(k);
	  colhas.Set(j);
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
    return make_shared<BlockVector>(vecs);
  }
  
  AutoVector BlockMatrix :: CreateColVector () const {
    Array<shared_ptr<BaseVector>> vecs(h);
    for (auto row:Range(h)) {
      vecs[row] = row_reps[row]->CreateColVector();
    }
    return make_shared<BlockVector>(vecs);
  }


  
}
