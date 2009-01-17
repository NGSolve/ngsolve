#ifndef FILE_H1HOFETP
#define FILE_H1HOFETP

/*********************************************************************/
/* File:   h1hofetp.hpp                                              */
/* Author: Start                                                     */
/* Date:   11 June 2007                                              */
/*********************************************************************/

/**
  High order finite elements for H^1  in full Tensor Product basis
*/


#define JACOBI

template <int DIM> class TPFactors;

template <> class TPFactors<2>
{
public:
  Matrix<> fac_x, fac_y;
};

template <> class TPFactors<3>
{
public:
  Matrix<> fac_x, fac_y, fac_z;
};




template <ELEMENT_TYPE ET> class H1HighOrderTP;



template <ELEMENT_TYPE ET> class H1HOTPtrait { };

template<> class H1HOTPtrait<ET_SEGM>
{
public:
  typedef H1HighOrderSegm<> BASE;
};

template<> class H1HOTPtrait<ET_TRIG>
{
public:
  typedef H1HighOrderTrig<> BASE;
};

template<> class H1HOTPtrait<ET_QUAD>
{
public:
  typedef H1HighOrderQuad<> BASE;
};



template<> class H1HOTPtrait<ET_TET>
{
public:
  typedef H1HighOrderTet<> BASE;
};

template<> class H1HOTPtrait<ET_PRISM>
{
public:
  typedef H1HighOrderPrism<> BASE;
};

template<> class H1HOTPtrait<ET_HEX>
{
public:
  typedef H1HighOrderHex<> BASE;
};




template <ELEMENT_TYPE ET, int DIM = ET_trait<ET>::DIM> class BaseH1HighOrderTP;


template <ELEMENT_TYPE ET>
class BaseH1HighOrderTP<ET,1> : public H1HOTPtrait<ET>::BASE
{ 
protected:
public:
  typedef typename H1HOTPtrait<ET>::BASE BASE;
  enum { D = 1 };

  BaseH1HighOrderTP (int aorder)
    : H1HOTPtrait<ET>::BASE (aorder) { BASE::tp = true; }

  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;

  void CalcDShape (const IntegrationPoint & ip, 
                   FlatMatrix<> dshape) const;


  // int Sort (int i) const { return sort[i]; }
};


template <ELEMENT_TYPE ET>
class BaseH1HighOrderTP<ET,2> : public H1HOTPtrait<ET>::BASE
{ 
protected:
  FlatArray<int> map2dto1d;
  FlatArray<int> map1dto2d;
  FlatArray<int> pmap1dto2d;

  int ndof1d;

public:
  typedef typename H1HOTPtrait<ET>::BASE BASE;
  enum { D = 2 };

  BaseH1HighOrderTP (int aorder)
    : H1HOTPtrait<ET>::BASE (aorder) { BASE::tp = true; }

  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;

  virtual void CalcDShape (const IntegrationPoint & ip, 
                           FlatMatrix<> dshape) const;

  const H1HighOrderTP<ET> & Spec() const
  { return static_cast<const H1HighOrderTP<ET>&> (*this); }


  virtual void EvaluateShapeGridTrans (const IntegrationRuleTP<2> & ir,
				       const FlatVector<double> gridvalues,
				       FlatVector<double> coefs,
				       LocalHeap & lh) const;


  // int Sort (int i) const { return sort[i]; }

  int GetNDof1d () const { return ndof1d; }

  FlatArray<int> Get2dTo1dMapping () const { return map2dto1d; }
  FlatArray<int> Get1dTo2dMapping () const { return map1dto2d; }
  FlatArray<int> Get1dTo2dMappingPointer () const { return pmap1dto2d; }
};



template <ELEMENT_TYPE ET>
class BaseH1HighOrderTP<ET,3> : public H1HOTPtrait<ET>::BASE
{ 
protected:
  FlatArray<int> map3dto2d;
  FlatArray<int> map2dto1d;

  FlatArray<int> map2dto3d;
  FlatArray<int> pmap2dto3d;
  FlatArray<int> map1dto2d;
  FlatArray<int> pmap1dto2d;

  int ndof1d;
  int ndof2d;

public:
  typedef typename H1HOTPtrait<ET>::BASE BASE;
  enum { D = 3 };

  BaseH1HighOrderTP (int aorder)
    : H1HOTPtrait<ET>::BASE (aorder) { BASE::tp = true; }

  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;

  virtual void CalcDShape (const IntegrationPoint & ip, 
                           FlatMatrix<> dshape) const;


  const H1HighOrderTP<ET> & Spec() const
  { return static_cast<const H1HighOrderTP<ET>&> (*this); }

  virtual void EvaluateShapeGrid (const IntegrationRuleTP<3> & ir,
				  const FlatVector<double> coefs,
				  FlatVector<double> gridvalues,
				  LocalHeap & lh) const;
				  
  virtual void EvaluateShapeGridTrans (const IntegrationRuleTP<3> & ir,
				       const FlatVector<double> gridvalues,
				       FlatVector<double> coefs,
				       LocalHeap & lh) const;


  // int Sort (int i) const { return sort[i]; }

  int GetNDof1d () const { return ndof1d; }
  int GetNDof2d () const { return ndof2d; }

  FlatArray<int> Get3dTo2dMapping () const { return map3dto2d; }
  FlatArray<int> Get2dTo1dMapping () const { return map2dto1d; }
  FlatArray<int> Get2dTo3dMapping () const { return map2dto3d; }
  FlatArray<int> Get2dTo3dMappingPointer () const { return pmap2dto3d; }
  FlatArray<int> Get1dTo2dMapping () const { return map1dto2d; }
  FlatArray<int> Get1dTo2dMappingPointer () const { return pmap1dto2d; }
};






/* ********************************** 1D Elements ******************* */


template <>
class H1HighOrderTP<ET_SEGM> : public BaseH1HighOrderTP<ET_SEGM>
{
  FlatArray<int[2]> segm2tensor;
  int sort[2];

public:
  H1HighOrderTP (int aorder)
    : BaseH1HighOrderTP<ET_SEGM> (aorder) {;}

  virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx x, TFA & shape) const;
};




/* ********************************** 2D Elements ******************* */


template <>
class H1HighOrderTP<ET_TRIG> : public BaseH1HighOrderTP<ET_TRIG>
{
  FlatArray<int[3]> trig2tensor;
  FlatArray<int[2]> segm2tensor;
  int sort[3];

public:
  H1HighOrderTP (int aorder)
    : BaseH1HighOrderTP<ET_TRIG,2> (aorder) {;}

  virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);


  template<typename Tx, typename Ty, typename TFA>  
  void T_CalcShape (Tx x, Ty y, TFA & shape) const; 


  // shape[i] = facx[i] * facy[2t1[i]]

  template <typename Tx, typename Tres>
  void CalcXFactor (Tx x,  Tres & facx, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max (2, order);
    FlatMatrix<Tx> polsx(order+1, horder, lh);

    Tx hv = 1;
    for (int i = 0; i <= order; i++)
      {
#ifdef JACOBI
        JacobiPolynomialMult (horder-2, 2*x-1, max(0,2*i-2), 0, x*hv, polsx.Row(i).Addr(1));
#else
        LegendrePolynomialMult (horder-2, 2*x-1, x*hv, polsx.Row(i).Addr(1));
#endif
        polsx(i,0) = hv;
        hv *= 1-x;
      }

    for (int jj = 0; jj < ndof; jj++)
      facx[jj] = polsx(trig2tensor[jj][1]+trig2tensor[jj][2],
                       trig2tensor[jj][0]); 
  }

  template <typename Ty, typename Tres>
  void CalcYFactor (Ty y, Tres & facy, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max (order, 2);
    FlatArray<Ty> poly(horder, lh);
    Vec<2, Ty> polyy;
    
    LegendrePolynomialMult (horder-2, 2*y-1, y, poly.Addr(1));
    poly[0] = 1;

    polyy(0) = 1;
    polyy(1) = (1-y); 

    for (int j1d = 0; j1d < segm2tensor.Size(); j1d++)
      facy[j1d] = poly[segm2tensor[j1d][0]] * polyy(segm2tensor[j1d][1]);
  }

};






















template <>
class H1HighOrderTP<ET_QUAD> : public BaseH1HighOrderTP<ET_QUAD>
{
  FlatArray<int[2]> quad2tensor;
  FlatArray<int[1]> segm2tensor;
  FlatArray<double> factor;

public:
  H1HighOrderTP<ET_QUAD> (int aorder)
    : BaseH1HighOrderTP<ET_QUAD> (aorder) { ; }
  
  virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);

  template<typename Tx, typename Ty, typename TFA>  
  void T_CalcShape (Tx x, Ty y, TFA & shape) const; 


  template <typename Tx, typename Tres>
  void CalcXFactor (Tx x,  Tres & facx, LocalHeap & lh) const
  { 
    HeapReset hr(lh);
    FlatVector<Tx> polx(order+1, lh);

    LegendrePolynomialMult (order-2, 2*x-1, x*(1-x), polx.Addr(2));
    polx[0] = 1-x;
    polx[1] = x;

    for (int jj = 0; jj < ndof; jj++)
      facx[jj] = factor[jj] * polx(quad2tensor[jj][0]);
  }

  template <typename Ty, typename Tres>
  void CalcYFactor (Ty y, Tres & facy, LocalHeap & lh) const
  { 
    HeapReset hr(lh);

    LegendrePolynomialMult (order-2, 2*y-1, y*(1-y), facy.Range(2,order+1));

    facy[0] = 1-y;
    facy[1] = y;
  }
};




/* ********************************** 3D Elements ******************* */


template <>
class H1HighOrderTP<ET_TET> : public BaseH1HighOrderTP<ET_TET>
{
  FlatArray<int[4]> tet2tensor;
  FlatArray<int[3]> trig2tensor;
  FlatArray<int[2]> segm2tensor;
  int sort[4];

public:
  H1HighOrderTP<ET_TET> (int aorder) : BaseH1HighOrderTP<ET_TET> (aorder) { ; }

  /// builds tables
  virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);



  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void T_CalcShape (Tx x, Ty y, Tz z, TFA & shape) const; 


  // shape[i] = facx[i] * facy[3t2[i]] * facz[2t1[3t2[i]]]

  template <typename Tx, typename Tres>
  void CalcXFactor (Tx x,  Tres & facx, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max (2, order);
    FlatMatrix<Tx> polsx(order+1, horder, lh);

    Tx hv = 1;
    for (int i = 0; i <= order; i++)
      {
#ifdef JACOBI
        JacobiPolynomialMult (horder-2, 2*x-1, max(0,2*i-2), 0, x*hv, polsx.Row(i).Addr(1));
#else
        LegendrePolynomialMult (horder-2, 2*x-1, x*hv, polsx.Row(i).Addr(1));
#endif
        polsx(i,0) = hv;
        hv *= 1-x;
      }

    for (int jj = 0; jj < ndof; jj++)
      facx[jj] = polsx(tet2tensor[jj][1]+tet2tensor[jj][2]+tet2tensor[jj][3],
                       tet2tensor[jj][0]); 
  }

  template <typename Ty, typename Tres>
  void CalcYFactor (Ty y, Tres & facy, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max (2, order);
    FlatMatrix<Ty> polsy(order+1, horder, lh);

    Ty hv = 1.0;
    for (int i = 0; i <= order; i++)
      {
#ifdef JACOBI
        JacobiPolynomialMult (horder-2, 2*y-1, max(0,2*i-2), 0, y*hv, polsy.Row(i).Addr(1));
#else
        LegendrePolynomialMult (horder-2, 2*y-1, y*hv, polsy.Row(i).Addr(1));
#endif
        polsy(i,0) = hv;
        hv *= 1-y;
      }

    for (int jj = 0; jj < trig2tensor.Size(); jj++)
      facy[jj] = polsy(trig2tensor[jj][1]+trig2tensor[jj][2],
                       trig2tensor[jj][0]); 
  }

  template <typename Tz, typename Tres>
  void CalcZFactor (Tz z, Tres & facz, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max (order, 2);
    FlatArray<Tz> polz(horder, lh);
    Vec<2, Tz> polzz;
    
    LegendrePolynomialMult (horder-2, 2*z-1, z, polz.Addr(1));
    polz[0] = 1;

    polzz(0) = 1;
    polzz(1) = (1-z); 

    for (int j1d = 0; j1d < segm2tensor.Size(); j1d++)
      facz[j1d] = polz[segm2tensor[j1d][0]] * polzz(segm2tensor[j1d][1]);
  }
  
};











template <> 
class H1HighOrderTP<ET_PRISM> : public BaseH1HighOrderTP<ET_PRISM>
{
  FlatArray<int[4]> prism2tensor;
  FlatArray<int[3]> trig2tensor;   // y-z quad
  FlatArray<int[2]> segm2tensor;   // z
  FlatArray<double> factor;
  int sort[6];

public:
  H1HighOrderTP<ET_PRISM> (int aorder)
    : BaseH1HighOrderTP<ET_PRISM> (aorder) { ; }

  virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);

  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void T_CalcShape (Tx x, Ty y, Tz z, TFA & shape) const; 

  // shape[i] = facx[i] * facy[3t2[i]] * facz[2t1[3t2[i]]]

  template <typename Tz, typename Tres>
  void CalcXFactor (Tz z, Tres & facz, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatArray<Tz> polz(order+1, lh);

    LegendrePolynomialMult (order-2, 2*z-1, z*(1-z), polz.Addr(2));
    polz[0] = 1-z;
    polz[1] = z;

    for (int jj = 0; jj < ndof; jj++)
      facz[jj] = factor[jj] * polz[prism2tensor[jj][3]];
  }



  template <typename Tx, typename Tres>
  void CalcYFactor (Tx x, Tres & facx, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max(2, order);
    FlatArray<Tx> powx(order+1, lh);
    FlatMatrix<Tx> polsx(order+1, horder, lh);

    for (int i = 0; i <= order; i++)
      {
        Tx * hp = &polsx(i, 1);
#ifdef JACOBI
        JacobiPolynomial (horder-2, 2*x-1, max(0,2*i-2), 0, hp);
#else
        LegendrePolynomial (horder-2, 2*x-1, hp);
#endif
        for (int j = 1; j < horder; j++) polsx(i,j) *= x;
        polsx(i,0) = 1.0;
      }

    Tx hv = 1.0;
    for (int i = 0; i <= order; i++, hv *= (1-x))
      powx[i] = hv;
    
    int nd2d = trig2tensor.Size();
    for (int jj = 0; jj < nd2d; jj++)
      {
        int power = trig2tensor[jj][1]+trig2tensor[jj][2];
        facx[jj] = polsx(power, trig2tensor[jj][0]) * powx[power];
      }
  }


  template <typename Ty, typename Tres>
  void CalcZFactor (Ty y, Tres & facy, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max(2, order);
    FlatVector<Ty> poly(horder, lh);
    Vec<2,Ty> polyy;

    LegendrePolynomialMult (horder-2, 2*y-1, y, poly.Addr(1));
    poly(0) = 1.0;

    polyy(0) = 1;
    polyy(1) = 1-y;
    
    int nd1d = segm2tensor.Size();
    for (int jj = 0; jj < nd1d; jj++)
      facy[jj] = poly(segm2tensor[jj][0]) * polyy(segm2tensor[jj][1]);
  }
};










template <>
class H1HighOrderTP<ET_HEX> : public BaseH1HighOrderTP<ET_HEX>
{

  FlatArray<int[3]> hex2tensor;
  // FlatArray<int[2]> quad2tensor;   // y-z quad
  // FlatArray<int[1]> segm2tensor;   // z
  FlatArray<double> factor;

public:
  H1HighOrderTP<ET_HEX> (int aorder)
    : BaseH1HighOrderTP<ET_HEX> (aorder) { ; }

  virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);

  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void T_CalcShape (Tx x, Ty y, Tz z, TFA & shape) const; 


  // shape[i] = facx[i] * facy[3t2[i]] * facz[2t1[3t2[i]]]

  template <typename Tx, typename Tres>
  void CalcXFactor (Tx x, Tres & facx, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatVector<Tx> polx(order+1, lh);

    LegendrePolynomialMult (order-2, 2*x-1, x*(1-x), polx.Addr(2));
    polx[0] = 1-x;
    polx[1] = x;

    for (int jj = 0; jj < ndof; jj++)
      facx[jj] = factor[jj] * polx(hex2tensor[jj][0]);
  }

  template <typename Ty, typename Tres>
  void CalcYFactor (Ty y, Tres & facy, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatVector<Ty> poly(order+1, lh);

    LegendrePolynomialMult (order-2, 2*y-1, y*(1-y), poly.Addr(2));
    poly[0] = 1-y;
    poly[1] = y;

    int nd2d = (order+1)*(order+1);
    for (int jj = 0; jj < nd2d; jj++)
      facy[jj] = poly(jj/(order+1));
  }


  template <typename Tz, typename Tres>
  void CalcZFactor (Tz z, Tres & facz, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    LegendrePolynomialMult (order-2, 2*z-1, z*(1-z), facz.Range(2,order+1));

    facz[0] = 1-z;
    facz[1] = z;
  }


};











#endif



