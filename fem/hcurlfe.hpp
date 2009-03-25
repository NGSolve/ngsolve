#ifndef FILE_HCURLFE
#define FILE_HCURLFE

/*********************************************************************/
/* File:   hcurlfe.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   16. Apr. 2000                                             */
/*********************************************************************/

/*
   HCurl Finite Element Definitions
*/


template <int D> class HDivFiniteElement;


template <int D>
class DIM_CURL_TRAIT
{
public:
  enum { DIM = (D*(D-1))/2 };
};

/*
template <> class DIM_CURL_TRAIT<2>
{
public:
  enum { DIM = 1 };
};
*/

template <> class DIM_CURL_TRAIT<1>
{
public:
  enum { DIM = 1 };  // should be 0; changed to make gcc34 feel ok
};


/**
   H(Curl) finite element of dimension D
 */
template <int D>
class HCurlFiniteElement : public FiniteElement
{

public:
  enum { DIM = D };
  enum { DIM_CURL = DIM_CURL_TRAIT<D>::DIM };

protected:

  class IPData
  {
  public:
    FlatMatrixFixWidth<DIM> shape;
    FlatMatrixFixWidth<DIM_CURL> curlshape;
  };

  IPData * p_ipdata;
  DynamicMem<double> * block;

public:
  ///
  HCurlFiniteElement () { p_ipdata = 0; block = 0; }

  /// 
  HCurlFiniteElement (ELEMENT_TYPE aeltype, int andof, int aorder)
    : FiniteElement (DIM, aeltype, andof, aorder) { p_ipdata = 0; block = 0; }
  
  virtual ~HCurlFiniteElement ();

  virtual string ClassName(void) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<DIM> shape) const = 0;
  
  /// compute curl of shape, default: numerical diff
  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<DIM_CURL> curlshape) const;

  /// compute shape
  virtual void CalcMappedShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
                                FlatMatrixFixWidth<DIM> shape) const;

  /// compute curl of shape
  virtual void CalcMappedCurlShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
                                    FlatMatrixFixWidth<DIM_CURL> curlshape) const;


  ///
  const FlatMatrixFixWidth<DIM> GetShape (const IntegrationPoint & ip, 
					  LocalHeap & lh) const
  {
    if (p_ipdata && ip.IPNr() >= 0)
      {
	return p_ipdata[ip.IPNr()].shape;
      }
    else
      {
	FlatMatrixFixWidth<DIM> shape(ndof, lh);
	CalcShape (ip, shape);
	return shape;
      }      
  }

  ///
  const FlatMatrixFixWidth<DIM_CURL> GetCurlShape (const IntegrationPoint & ip, 
						   LocalHeap & lh) const
  {
    if (p_ipdata && ip.IPNr() >= 0)
      {
	return p_ipdata[ip.IPNr()].curlshape;
      }
    else
      {
	FlatMatrixFixWidth<DIM_CURL> curlshape(ndof, lh);
	CalcCurlShape (ip, curlshape);
	return curlshape;
      }
  }  
 
  template <typename TVX>
  Vec<DIM_CURL, typename TVX::TSCAL> 
  EvaluateCurlShape (const IntegrationPoint & ip, 
                     const TVX & x, LocalHeap & lh) const
  {
    return Trans (GetCurlShape(ip, lh)) * x;
  } 

  virtual Vec<DIM_CURL> 
  EvaluateCurlShape (const IntegrationPoint & ip, 
                     FlatVector<double> x, LocalHeap & lh) const
  {
    return Trans (GetCurlShape(ip, lh)) * x;
  }  

  ///
  void CalcIPData (Array<IPData> & ipdata);

protected:
  ///
  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<D> shape) const
  { ; }
  ///
  virtual void CalcShape2 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<D> shape) const
  { ; }
  
  virtual void CalcShape3 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<D> shape) const
  { ; }
  
  virtual void CalcShape4 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<D> shape) const
  { ; }
  
  ///
  void ComputeEdgeMoments (int enr, ScalarFiniteElement<1> & testfe,
			   FlatMatrix<> moments, int order, int shape = 1) const;
  ///
  void ComputeFaceMoments (int fnr, HDivFiniteElement<2> & testfe,
			   FlatMatrix<> moments, int order, int shape = 1) const;
  ///
  void ComputeVolMoments (HDivFiniteElement<3> & testfe,
			  FlatMatrix<> moments, int order, int shape = 1) const;
};




template <int D>
extern void ComputeGradientMatrix (const ScalarFiniteElement<D> & h1fe,
				   const HCurlFiniteElement<D> & hcurlfe,
				   FlatMatrix<> gradient);


/* **************************** Segm Elements *************** */


///
class FE_NedelecSegm1 : public HCurlFiniteElement<1>
{
public:

  enum { NDOF = 1 };

private:
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;

public:
  ///
  FE_NedelecSegm1();
  ///
  virtual ~FE_NedelecSegm1();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<1> shape) const;
};



///
class FE_NedelecSegm2 : public HCurlFiniteElement<1>
{
public:

  enum { NDOF = 2 };

private:
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;

public:
  ///
  FE_NedelecSegm2();
  ///
  virtual ~FE_NedelecSegm2();

  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<1> shape) const;
};



///
class FE_NedelecSegm3 : public HCurlFiniteElement<1>
{
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;

public:
  ///
  FE_NedelecSegm3();
  ///
  virtual ~FE_NedelecSegm3();

  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<1> shape) const;
};










/* *********************** Quad elements ******************* */

/// Gradients of Q1
class FE_NedelecQuad1 : public HCurlFiniteElement<2>
{
public:
  enum { NDOF = 4 };

protected:
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
public:
  ///
  FE_NedelecQuad1();
  ///
  virtual ~FE_NedelecQuad1();

  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<2> shape) const;

};



/*
template <int ORDER, int ZORDER>
class FE_TNedelecQuadTraits
{
public:
  enum { NDOF = ORDER * (ZORDER+1) + (ORDER+1) * ZORDER };
  enum { NEDGEDOF = 2 * (ORDER + ZORDER) - 4 };
};
*/

template <int ORDER, int ZORDER>
class FE_TNedelecQuad : public HCurlFiniteElement<2>
{
public:
  enum { NDOF = ORDER * (ZORDER+1) + (ORDER+1) * ZORDER };
  enum { NEDGEDOF = 2 * (ORDER + ZORDER) - 4 };

protected:
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  ///
  // static Mat<FE_TNedelecQuad<ORDER,ZORDER>::NDOF> trans;

  // static Mat<FE_TNedelecQuadTraits<ORDER,ZORDER>::NDOF> trans;
  // static Mat<FE_TNedelecQuadTraits<ORDER,ZORDER>::NEDGEDOF> trans2;

  static Matrix<> trans;
  static Matrix<> trans2;

  FE_NedelecQuad1 quad1;

public:
  enum { MAXORDER = (ORDER > ZORDER) ? ORDER : ZORDER };

  ///
  FE_TNedelecQuad();
  ///
  virtual ~FE_TNedelecQuad();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<2> shape) const;
  ///
  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<2> shape) const;
  ///
  virtual void CalcShape2 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<2> shape) const;
  ///
  void Orthogonalize();  
};




/* ******************** triangular elements *********************** */

/// Lowest order Nedelec
class FE_NedelecTrig1 : public HCurlFiniteElement<2>
{
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;

public:

  ///
  FE_NedelecTrig1();
  ///
  virtual ~FE_NedelecTrig1();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<2> shape) const;
};



/// Nedelec type 2, order 1, gradients of P2
class FE_NedelecTrig2 : public HCurlFiniteElement<2>
{
public:
  enum { NDOF = 6 };

private:
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  ///
  static Mat<NDOF> trans;

public:
  ///
  FE_NedelecTrig2();
  ///
  virtual ~FE_NedelecTrig2();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<2> shape) const;

  ///
  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<2> shape) const;

  ///
  void Orthogonalize();
};







/// Nedelec type 2, order 2, gradients of P3
class FE_NedelecTrig3 : public HCurlFiniteElement<2>
{
public:
  enum { NDOF = 12 };
  enum { NEDGEDOF = 6 };
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  ///
  static Mat<NDOF> trans;
  ///
  static Mat<NEDGEDOF> trans2;
  ///
  FE_NedelecTrig2 trig1;
public:
  ///
  FE_NedelecTrig3();
  ///
  virtual ~FE_NedelecTrig3();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<2> shape) const;

  ///
  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<2> shape) const;

  ///
  virtual void CalcShape2 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<2> shape) const;
  
  ///
  void Orthogonalize();
};











/* *********************** Tetrahedral elements ********************** */


///
class FE_NedelecTet1 : public HCurlFiniteElement<3>
{
public:
  enum { NDOF = 6 };

private:
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;

public:

  ///
  FE_NedelecTet1();
  ///
  virtual ~FE_NedelecTet1();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<3> curlshape) const;
};



///
class FE_NedelecTet2 : public HCurlFiniteElement<3>
{
public:

  enum { NDOF = 12 };

private:
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  ///
  static Mat<NDOF> trans;

public:
  ///
  FE_NedelecTet2();
  ///
  virtual ~FE_NedelecTet2();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;


  ///
  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  void Orthogonalize();
};





/// 2nd order Nedelec element of class II
class FE_NedelecTet3 : public HCurlFiniteElement<3>
{
public:
  enum { NDOF = 30 };
  enum { NEDGEDOF = 12 };
  enum { NFACEDOF = 12 };

protected:
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  ///
  static Mat<NDOF> trans;
  ///
  static Mat<NEDGEDOF> trans2;
  ///
  static Mat<NFACEDOF> trans3;

  FE_NedelecTet1 tet1;
public:

  ///
  FE_NedelecTet3();
  ///
  virtual ~FE_NedelecTet3();


  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<3> curlshape) const;
  
  ///
  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcShape2 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcShape3 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcCurlShape3 (const IntegrationPoint & ip, 
			       FlatMatrixFixWidth<3> shape) const;

  ///
  void Orthogonalize();
};
/// 2nd order Nedelec element of class II, without gradient fields
class FE_NedelecTet3NoGrad : public HCurlFiniteElement<3>
{
public:
  enum { NDOF = 18 };
  enum { NFACEDOF = 12 };

protected:
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  static Mat<NFACEDOF> trans3;

  FE_NedelecTet1 tet1;
public:
  FE_NedelecTet3NoGrad();
  virtual ~FE_NedelecTet3NoGrad();


  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<3> curlshape) const;

  virtual void CalcShape3 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  virtual void CalcCurlShape3 (const IntegrationPoint & ip, 
			       FlatMatrixFixWidth<3> shape) const;

  void Orthogonalize();
};



/* *********************** Hex elements ************************ */ 


///
class FE_NedelecHex1 : public HCurlFiniteElement<3> 
{
  /// 
  static Array<IPData> ipdata; 
  bool ipdatadestructed;
  
public: 
  ///
  FE_NedelecHex1(); 
  ///
  virtual ~FE_NedelecHex1(); 
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const; 
}; 
  


/* *********************** Prism elements ********************** */

///
class FE_NedelecPrism1 : public HCurlFiniteElement<3>
{
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;

public:

  ///
  FE_NedelecPrism1();
  ///
  virtual ~FE_NedelecPrism1();

  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;
};

/// \nabla Q (2,ZORDER)
template <int ZORDER>
class FE_TNedelecPrism2 : public HCurlFiniteElement<3>
{
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  ///
  static Matrix<> trans;
  ///
  static Matrix<> trans2;
  ///
  static Matrix<> trans3;

  FE_NedelecPrism1 prism1;

public:
  enum { NDOF = 6 * (ZORDER+1) + 6 * ZORDER };
  enum { NEDGEDOF = 6 + 3 * (ZORDER-1) };
  enum { NFACEDOF = 9 * ZORDER - 6} ;
  enum { MAXORDER = (2 > ZORDER) ? 2 : ZORDER };

  ///
  FE_TNedelecPrism2();
  ///
  virtual ~FE_TNedelecPrism2();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcShape2 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcShape3 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  void Orthogonalize();

};




/// \nabla Q (3,ZORDER)
template <int ZORDER>
class FE_TNedelecPrism3 : public HCurlFiniteElement<3>
{
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  ///
  static Matrix<> trans;
  ///
  static Matrix<> trans2;
  ///
  static Matrix<> trans_quad;
  ///
  static Matrix<> trans_trig;

  FE_NedelecPrism1 prism1;
  FE_NedelecTrig3 trig3;
  FE_Trig2 h1trig2;
  FE_Trig3Pot h1trig3;
  FE_TSegmL2<ZORDER> segm;
public:
  enum { NDOF = 12 * (ZORDER+1) + 10 * ZORDER };
  enum { NEDGEDOF = 12 + 3 * (ZORDER-1) };
  enum { NQUADFACEDOF = 3 * (5*ZORDER-3) };
  enum { NTRIGFACEDOF = 6 };
  enum { MAXORDER = (3 > ZORDER) ? 3 : ZORDER };
  enum { NINNERDOF = 3 * (ZORDER-1) + ZORDER };

  ///
  FE_TNedelecPrism3();
  ///
  virtual ~FE_TNedelecPrism3();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcShape2 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  /// quad face dofs
  virtual void CalcShape3 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  /// trig face dofs
  virtual void CalcShape4 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcInner (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void GetInternalDofs (Array<int> & idofs) const;

  ///
  void Orthogonalize();
};





/// \nabla Q (3,ZORDER)
template <int ZORDER>
class FE_TNedelecPrism3NoGrad : public HCurlFiniteElement<3>
{
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  ///
  static Matrix<> trans_quad;
  ///
  static Matrix<> trans_trig;

  FE_NedelecPrism1 prism1;
  FE_NedelecTrig3 trig3;
  FE_Trig2 h1trig2;
  FE_Trig3Pot h1trig3;
  FE_TSegmL2<ZORDER> segm;
public:
  //  enum { NDOF = 12 * (ZORDER+1) + 10 * ZORDER };
  //  enum { NEDGEDOF = 12 + 3 * (ZORDER-1) };
  // 12 z + 12 + 10 z - 12 - 3z + 3 = 19 z + 3
  enum { NDOF = 19 * ZORDER + 3 };
  enum { NQUADFACEDOF = 3 * (5*ZORDER-3) };
  enum { NTRIGFACEDOF = 6 };
  enum { MAXORDER = (3 > ZORDER) ? 3 : ZORDER };
  // enum { NINNERDOF = 3 * (ZORDER-1) + ZORDER };
  enum { NINNERDOF = 3 * (ZORDER-1) + 1 };

  ///
  FE_TNedelecPrism3NoGrad();
  ///
  virtual ~FE_TNedelecPrism3NoGrad();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcShape2 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  /// quad face dofs
  virtual void CalcShape3 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  /// trig face dofs + inner dofs
  virtual void CalcShape4 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void CalcInner (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void GetInternalDofs (Array<int> & idofs) const;

  ///
  void Orthogonalize();
};














///
class FE_NedelecPyramid1 : public HCurlFiniteElement<3>
{
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  ///
  static Matrix<> trans;
public:
  ///
  FE_NedelecPyramid1();
  ///
  virtual ~FE_NedelecPyramid1();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  void Orthogonalize();

};



///
class FE_NedelecPyramid2 : public HCurlFiniteElement<3>
{
public:
  enum { NDOF = 20 };
  enum { NEDGEDOF = 8 };

private:
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  ///
  static Matrix<> trans;
  static Matrix<> trans2;
  static Matrix<> trans3;
  
  ///
  FE_Quad1 quad1;
  FE_Quad2 quad2;
  FE_NedelecPyramid1 pyramid1;
public:
  ///
  FE_NedelecPyramid2();
  ///
  virtual ~FE_NedelecPyramid2();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  virtual void CalcShape2 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;


  ///
  void Orthogonalize();

};



///
class FE_NedelecPyramid3 : public HCurlFiniteElement<3>
{
public:
  enum { NDOF = 57 };
  enum { NEDGEDOF = 16 };
  enum { NFACEDOF = 24 };
  enum { NINNERDOF = 9 };
private:
  ///
  static Array<IPData> ipdata;
  bool ipdatadestructed;
  ///
  static Mat<NDOF> trans;
  static Mat<NEDGEDOF> trans2;
  static Mat<NFACEDOF> trans3;
  
  ///
  FE_Quad1 quad1;
  FE_Quad2 quad2;
  FE_Quad3 quad3;
  FE_NedelecPyramid1 pyramid1;
public:
  ///
  FE_NedelecPyramid3();
  ///
  virtual ~FE_NedelecPyramid3();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  virtual void CalcShape2 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  virtual void CalcShape3 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void GetInternalDofs (Array<int> & idofs) const;
  ///
  void Orthogonalize();
};





/*
///
class FE_NedelecPyramid3NoGrad : public HCurlFiniteElement<3>
{
public:
  //  enum { NDOF = 57 };
  // enum { NEDGEDOF = 16 };
  enum { NDOF = 41 };
  enum { NFACEDOF = 24 };
  enum { NINNERDOF = 9 };
private:
  ///
  static Array<IPData> ipdata;
  ///
  static Mat<NDOF> trans;
  // static Mat<NEDGEDOF> trans2;
  static Mat<NFACEDOF> trans3;
  
  ///
  FE_Quad1 quad1;
  FE_Quad2 quad2;
  FE_Quad3 quad3;
  FE_NedelecPyramid1 pyramid1;
public:
  ///
  FE_NedelecPyramid3();
  ///
  virtual ~FE_NedelecPyramid3();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  virtual void CalcShape2 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  virtual void CalcShape3 (const IntegrationPoint & ip, 
			   FlatMatrixFixWidth<3> shape) const;

  ///
  virtual void GetInternalDofs (Array<int> & idofs) const;
  ///
  void Orthogonalize();
};
*/









#ifdef OLD

/// extension to Nedelec type II
class FE_NedelecPyramid1b : public HCurlFiniteElement<3>
{
  ///
  static Array<IPData> ipdata;
  ///
  FE_NedelecPyramid1 pyramid1;
  ///
  static FlatMatrix<> trans;

public:

  ///
  FE_NedelecPyramid1b();
  ///
  virtual ~FE_NedelecPyramid1b();

  ///
  virtual int SpatialDim () const { return 3; }
  ///
  virtual int GetNDof () const { return 16; }
  ///
  virtual int Order () const { return 2; }
  ///
  virtual ELEMENT_TYPE ElementType() const { return ET_PYRAMID; }
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;
  ///
  virtual void CalcShape1 (const IntegrationPoint & ip, 
			   FlatMatrix<> shape) const;

  ///
  void Orthogonalize();
};

#endif








#endif
