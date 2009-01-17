#ifdef TPFO
#include "l2hofetpfo.hpp"
#endif



#ifndef FILE_L2HOFETP
#define FILE_L2HOFETP

/*********************************************************************/
/* File:   l2hofetp.hpp                                              */
/* Author: J Schöberl, Aachen                                        */
/* Date:   4 Oct 2008                                                */
/*********************************************************************/

/**
   High order finite elements for L2  in full Tensor Product basis
*/



class L2HighOrderTrigTP : public L2HighOrderTrig
{
protected:
  LocalHeap & lh;
  FlatArray<int[3]> trig2tensor;
  FlatArray<int[2]> segm2tensor;
  FlatArray<int> map2dto1d;

  FlatArray<int> map1dto2d;
  FlatArray<int> pmap1dto2d;

  int ndof1d;
  int sort[3];
  static ARRAY<RecPol*> recpols;

  static FlatMatrix<> fac_x, fac_y;
  static FlatMatrix<AutoDiff<1> > fac_xdx, fac_ydy;

public:
  L2HighOrderTrigTP (int aorder, LocalHeap & alh);

  virtual void SetVertexNumbers (FlatArray<int> & avnums);
  // virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);

  int GetNDof1d () const { return ndof1d; }

  FlatArray<int> Get2dTo1dMapping () const { return map2dto1d; }
  FlatArray<int> Get1dTo2dMapping () const { return map1dto2d; }
  FlatArray<int> Get1dTo2dMappingPointer () const { return pmap1dto2d; }

  
  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

  // int Sort (int i) const { return sort[i]; }

  template<typename Tx, typename Ty, typename TFA>  
  void T_CalcShape (Tx x, Ty y, TFA & shape) const; 



  virtual void EvaluateShapeGrid (const IntegrationRuleTP<2> & ir,
				  const FlatVector<double> coefs,
				  FlatVector<double> gridvalues,
				  LocalHeap & lh) const;


  virtual void EvaluateShapeGridTrans (const IntegrationRuleTP<2> & ir,
				       const FlatVector<double> gridvalues,
				       FlatVector<double> coefs,
				       LocalHeap & lh) const;

  virtual void EvaluateDShapeGrid (const IntegrationRuleTP<2> & ir,
				   const FlatVector<double> coefs,
				   FlatMatrixFixWidth<2> gridvalues,
				   LocalHeap & lh) const;
				  
  virtual void EvaluateDShapeGridTrans (const IntegrationRuleTP<2> & ir,
					const FlatMatrixFixWidth<2> gridvalues,
					FlatVector<double> coefs,
					LocalHeap & lh) const;

  virtual void CalcTPFactors (const IntegrationRuleTP<2> & ir,
                              TPFactors<2> & factors) const;

  virtual void SetTPFactors (const TPFactors<2> & factors);


  template <typename Tx, typename Tres>
  void CalcXFactor (Tx x,  Tres & facx, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<Tx> polsx(order+1, order+1, lh);
    Tx hv = 1;
    for (int i = 0; i <= order; i++)
      {
#ifdef JACOBI
        recpols[i] -> Evaluate (order-i, 2*x-1, polsx.Row(i));
        polsx.Row(i) *= hv;
        // JacobiPolynomialMult (order-i, 2*x-1, 2*i+1, 0, hv, polsx.Row(i));
#else
        // LegendrePolynomialMult (horder-2, 2*x-1, x*hv, polsx.Row(i).Addr(1));
#endif
        hv *= 1-x;
      }

    for (int jj = 0; jj < ndof; jj++)
      facx[jj] = polsx(trig2tensor[jj][1], trig2tensor[jj][0]); 
  }

  template <typename Ty, typename Tres>
  void CalcYFactor (Ty y, Tres & facy, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatArray<Ty> poly(order+1, lh);
    
    LegendrePolynomial (order, 2*y-1, poly);

    for (int j1d = 0; j1d < ndof1d; j1d++)
      facy[j1d] = poly[segm2tensor[j1d][0]];
  }

};





/* ************************** TET ********************************************** */








class L2HighOrderTetTP : public L2HighOrderTet
{
protected:
  LocalHeap & lh;
  FlatArray<int[4]> tet2tensor;
  FlatArray<int[3]> trig2tensor;
  FlatArray<int[2]> segm2tensor;
  FlatArray<int> map3dto2d;
  FlatArray<int> map2dto1d;

  FlatArray<int> map2dto3d;
  FlatArray<int> pmap2dto3d;
  FlatArray<int> map1dto2d;
  FlatArray<int> pmap1dto2d;


  int ndof2d;
  int ndof1d;
  int sort[4];
  static ARRAY<RecPol*> recpols;
  FlatMatrix<> fac_x, fac_y, fac_z;
  static FlatMatrix<AutoDiff<1> > fac_xdx, fac_ydy, fac_zdz;

public:
  L2HighOrderTetTP (int aorder, LocalHeap & alh);
  // virtual ~L2HighOrderTetTP ();

  virtual void SetVertexNumbers (FlatArray<int> & avnums);
  // virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);


  int GetNDof1d () const { return ndof1d; }
  int GetNDof2d () const { return ndof2d; }

  FlatArray<int> Get3dTo2dMapping () const { return map3dto2d; }
  FlatArray<int> Get2dTo1dMapping () const { return map2dto1d; }
  FlatArray<int> Get2dTo3dMapping () const { return map2dto3d; }
  FlatArray<int> Get2dTo3dMappingPointer () const { return pmap2dto3d; }
  FlatArray<int> Get1dTo2dMapping () const { return map1dto2d; }
  FlatArray<int> Get1dTo2dMappingPointer () const { return pmap1dto2d; }


  
  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

  // int Sort (int i) const { return sort[i]; }

  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void T_CalcShape (Tx x, Ty y, Tz z, TFA & shape) const; 


  virtual void CalcTPFactors (const IntegrationRuleTP<3> & ir,
                              TPFactors<3> & factors) const;

  virtual void SetTPFactors (const TPFactors<3> & factors);


  virtual void EvaluateShapeGrid (const IntegrationRuleTP<3> & ir,
				  const FlatVector<double> coefs,
				  FlatVector<double> gridvalues,
				  LocalHeap & lh) const;

  virtual void EvaluateShapeGridTrans (const IntegrationRuleTP<3> & ir,
				       const FlatVector<double> gridvalues,
				       FlatVector<double> coefs,
				       LocalHeap & lh) const;

  virtual void EvaluateDShapeGrid (const IntegrationRuleTP<3> & ir,
				   const FlatVector<double> coefs,
				   FlatMatrixFixWidth<3> gridvalues,
				   LocalHeap & lh) const;
				  
  virtual void EvaluateDShapeGridTrans (const IntegrationRuleTP<3> & ir,
					const FlatMatrixFixWidth<3> gridvalues,
					FlatVector<double> coefs,
					LocalHeap & lh) const;


  template <typename Tx, typename Tres>
  void CalcXFactor (Tx x,  Tres & facx, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<Tx> polsx(order+1, order+1, lh);

    Tx hv = 1;
    for (int i = 0; i <= order; i++)
      {
#ifdef JACOBI
        recpols[2*i+2] -> Evaluate (order-i, 2*x-1, polsx.Row(i));
        polsx.Row(i) *= hv;
#else
        // LegendrePolynomialMult (horder-2, 2*x-1, x*hv, polsx.Row(i).Addr(1));
#endif
        hv *= 1-x;
      }

    for (int jj = 0; jj < ndof; jj++)
      facx[jj] = polsx(tet2tensor[jj][1]+tet2tensor[jj][2], tet2tensor[jj][0]); 
  }


  template <typename Ty, typename Tres>
  void CalcYFactor (Ty y,  Tres & facy, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<Ty> polsy(order+1, order+1, lh);

    Ty hv = 1;
    for (int i = 0; i <= order; i++)
      {
#ifdef JACOBI
        recpols[2*i+1] -> Evaluate (order-i, 2*y-1, polsy.Row(i));
        polsy.Row(i) *= hv;
        // JacobiPolynomialMult (order-i, 2*y-1, 2*i+1, 0, hv, polsy.Row(i));
#else
        // LegendrePolynomialMult (horder-2, 2*y-1, x*hv, polsy.Row(i).Addr(1));
#endif
        hv *= 1-y;
      }

    for (int jj = 0; jj < ndof2d; jj++)
      facy[jj] = polsy(trig2tensor[jj][1], trig2tensor[jj][0]); 
  }




  template <typename Tz, typename Tres>
  void CalcZFactor (Tz z, Tres & facz, LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatArray<Tz> polz(order+1, lh);
    
    LegendrePolynomial (order, 2*z-1, polz);

    for (int j1d = 0; j1d < ndof1d; j1d++)
      facz[j1d] = polz[segm2tensor[j1d][0]];
  }

};






#endif
