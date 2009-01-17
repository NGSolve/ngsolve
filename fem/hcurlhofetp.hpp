#ifndef FILE_HCURLHOFETP
#define FILE_HCURLHOFETP

/*********************************************************************/
/* File:   h1curlfetp.hpp                                            */
/* Author: Start                                                     */
/* Date:   12 July 2007                                              */
/*********************************************************************/

/**
  High order finite elements for H(curl)  in full Tensor Product basis
*/


// #define JACOBI



/*
  All H(curl) basis functions are written in the form

  ( psi_i(x) phi_j(y) phi_k(z) )
  ( phi_i(x) psi_j(y) phi_k(z) )
  ( phi_i(x) phi_j(y) psi_k(z) )

  
  the curl is

  ( phi_i(x) phi_j'(y) psi_k(z) - phi_i(x) psi_j(y) phi_k'(z) )
  ( psi_i(x) phi_j(y) phi_k'(z) - phi_i'(x) phi_j(y) psi_k(z) )
  ( phi_i'(x) psi_j(y) phi_k(z) - psi_i(x) phi_j'(y) phi_k(z) )
  

  we store (phi_i, psi_i, phi_i') as triple in the Xfactor, and correspondingly 
  in Yfactor and Zfactor


  But first, we represent the basis funcitons in the form

  w (\nabla u v - u \nabla v)

  w is only needed for the original SZ Type-III trig functons
*/



class HCurlHighOrderTetTP : public HCurlHighOrderTet<>
{
  int sort[4];

  FlatArray<int[4]> tet2tensor;
  FlatArray<int> split;

  FlatArray<int[3]> trig2tensor;
  FlatArray<int> split2d;

  FlatArray<int[2]> segm2tensor;
  FlatArray<int> split1d;

  FlatArray<int> map3dto2d;
  FlatArray<int> map2dto1d;

public:
  HCurlHighOrderTetTP (int aorder) : HCurlHighOrderTet<> (aorder) { ; }

  /// builds tables
  virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);

  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<3> shape) const;

  int Sort (int i) const { return sort[i]; }

  FlatArray<int> Get3dTo2dMapping () const { return map3dto2d; }
  FlatArray<int> Get2dTo1dMapping () const { return map2dto1d; }

  FlatArray<int[4]> GetTet2TensorMapping () const { return tet2tensor; }
  FlatArray<int> GetSplit () const { return split; }

  FlatArray<int[3]> GetTrig2TensorMapping () const { return trig2tensor; }
  FlatArray<int> GetSplit2d () const { return split2d; }

  FlatArray<int[2]> GetSegm2TensorMapping () const { return segm2tensor; }
  FlatArray<int> GetSplit1d () const { return split1d; }

  int GetNDof2d () const { return trig2tensor.Size(); }
  int GetNDof1d () const { return segm2tensor.Size(); }
  bool IsGradient (int i) const { return split[i] == 0; }


  template <typename Tres>
  void CalcXFactor (AutoDiff<1> x,  Tres & facx, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max (2, order+1);
    FlatArray<AutoDiff<1> > powx(horder, lh);
    FlatMatrix<AutoDiff<1> > polsx(horder+1, horder, lh);

    for (int i = 0; i <= horder; i++)
      {
#ifdef JACOBIxxxxx
        JacobiPolynomial (horder-2, 2*x-1, max(0,2*i-2), 0, polsx.Row(i).Addr(1));
#else
        LegendrePolynomial (horder-2, 2*x-1, polsx.Row(i).Addr(1));
#endif
        for (int j = 1; j < horder; j++)
          polsx(i,j) *= x;
        polsx(i,0) = 1.0;
      }

    AutoDiff<1> hv = 1.0;
    for (int i = 0; i < horder; i++, hv *= (1-x))
      powx[i] = hv;
    
    for (int i = 0; i < ndof; i++)
      {
        int power = tet2tensor[i][1]+tet2tensor[i][2]+tet2tensor[i][3];
  
        AutoDiff<1> ux = 1.0, vx = 1.0, wx = 1.0;
        AutoDiff<1> fac[4] = { polsx(power, tet2tensor[i][0]),
                               powx[tet2tensor[i][1]],
                               powx[tet2tensor[i][2]],
                               powx[tet2tensor[i][3]] };

        
        
        for (int j = 0; j < split[i]-1; j++)
          wx *= fac[j];
        for (int j = max(0, split[i]-1); j < split[i]; j++)
          ux *= fac[j];
        for (int j = split[i]; j < 4; j++)
          vx *= fac[j];


        AutoDiff<1> prod = ux*vx*wx;

        facx(i)(0) = prod.Value();
        facx(i)(1) = wx.Value() * (ux.Value() * vx.DValue(0) - ux.DValue(0) * vx.Value());
        facx(i)(2) = prod.DValue(0);
      }
  }





  template <typename Tres>
  void CalcYFactor (AutoDiff<1> y, Tres & facy, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max (2, order+1);
    FlatVector<AutoDiff<1> > powy(horder+1, lh);
    FlatMatrix<AutoDiff<1> > polsy(horder+1, horder, lh);

    for (int i = 0; i <= horder; i++)
      {
#ifdef JACOBIxxxx
        JacobiPolynomial (horder-2, 2*y-1, max(0,2*i-2), 0, polsy.Row(i).Addr(1));
#else
        LegendrePolynomial (horder-2, 2*y-1, polsy.Row(i).Addr(1));
#endif
        for (int j = 1; j < horder; j++)
          polsy(i,j) *= y;
        polsy(i,0) = 1.0;
      }

    AutoDiff<1> hv = 1.0;
    for (int i = 0; i <= order; i++, hv *= (1-y) )
      powy[i] = hv;
    

    int ndof2d = trig2tensor.Size();
    for (int i = 0; i < ndof2d; i++)
      {
        int power = trig2tensor[i][1]+trig2tensor[i][2];

        AutoDiff<1> uy = 1.0, vy = 1.0, wy = 1.0;
        AutoDiff<1> fac[3] = { polsy(power, trig2tensor[i][0]),
                               powy[trig2tensor[i][1]],
                               powy[trig2tensor[i][2]] };

        for (int j = 0; j < split2d[i]-1; j++)
          wy *= fac[j];
        for (int j = max (0, split2d[i]-1); j < split2d[i]; j++)
          uy *= fac[j];
        for (int j = split2d[i]; j < 3; j++)
          vy *= fac[j];

        AutoDiff<1> prod = uy*vy*wy;

        facy(i)(0) = prod.Value();
        facy(i)(1) = wy.Value() * (uy.Value() * vy.DValue(0) - uy.DValue(0) * vy.Value());
        facy(i)(2) = prod.DValue(0);
      }
  }



  template <typename Tres>
  void CalcZFactor (AutoDiff<1> z, Tres & facz, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max (order+1, 2);
    FlatVector<AutoDiff<1> > polz(horder+1, lh);
    Vec<2, AutoDiff<1> > powz;
    
    LegendrePolynomialMult (horder-2, 2*z-1, z, polz.Addr(1));
    polz[0] = 1;

    powz(0) = 1;
    powz(1) = (1-z); 

    int ndof1d = segm2tensor.Size();

    for (int i = 0; i < ndof1d; i++)
      {
        if (split1d[i] == 0)
          {
            AutoDiff<1> uz = 1.0;
            AutoDiff<1> vz = polz(segm2tensor[i][0]) * powz(segm2tensor[i][1]);
        
            facz(i)(0) =  vz.Value();
            facz(i)(1) =  vz.DValue(0);
            facz(i)(2) =  vz.DValue(0);
          }
        else
          { // additional info: segm2tensor[i][0] = segm2tensor[i][1] = 1
            facz(i) (0) = 0;
            facz(i) (1) = -1;
            facz(i) (2) = 0;
          }
      }
  }

};















class HCurlHighOrderPrismTP : public HCurlHighOrderPrism<>
{
  int sort[6];

  FlatArray<int[4]> prism2tensor;
  FlatArray<int> split;
  // FlatArray<int> splita;
  FlatArray<bool> isgradient;

  FlatArray<int[3]> trig2tensor;
  FlatArray<int> split2d;
  // FlatArray<int> split2da;

  FlatArray<int[2]> segm2tensor;
  FlatArray<int> split1d;


  FlatArray<int> map3dto2d;
  FlatArray<int> map2dto1d;

  FlatArray<double> factorxy;  
  FlatArray<double> factorz;

public:
  HCurlHighOrderPrismTP (int aorder) : HCurlHighOrderPrism<> (aorder) { ; }

  /// builds tables
  virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);

  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<3> shape) const;

  int Sort (int i) const { return sort[i]; }

  FlatArray<int> Get3dTo2dMapping () const { return map3dto2d; }
  FlatArray<int> Get2dTo1dMapping () const { return map2dto1d; }

  FlatArray<int[4]> GetPrism2TensorMapping () const { return prism2tensor; }
  FlatArray<int> GetSplit () const { return split; }

  FlatArray<int[3]> GetTrig2TensorMapping () const { return trig2tensor; }
  FlatArray<int> GetSplit2d () const { return split2d; }

  FlatArray<int[2]> GetSegm2TensorMapping () const { return segm2tensor; }
  FlatArray<int> GetSplit1d () const { return split1d; }

  int GetNDof2d () const { return trig2tensor.Size(); }
  int GetNDof1d () const { return segm2tensor.Size(); }
  bool IsGradient (int i) const { return isgradient[i]; }







  template <typename Tres>
  void CalcXFactor (AutoDiff<1> z,  Tres & facz, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max (2, order+1);
    FlatVector<AutoDiff<1> > polz(horder+1, lh);   

    LegendrePolynomialMult (horder-2, 2*z-1, z*(1-z), polz.Addr(2));
    polz(0) = 1-z;
    polz(1) = z;
    
    for (int i = 0; i < ndof; i++)
      {
        facz(i)(0) = factorxy[i] * polz(prism2tensor[i][3]).Value();
        facz(i)(1) = factorz[i] * polz(prism2tensor[i][3]).DValue(0);
        facz(i)(2) = factorxy[i] * polz(prism2tensor[i][3]).DValue(0);
      }
  }




  template <typename Tres>
  void CalcYFactor (AutoDiff<1> y, Tres & facy, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max (2, order+1);
    FlatVector<AutoDiff<1> > powy(horder+1, lh);
    FlatMatrix<AutoDiff<1> > polsy(horder+1, horder, lh);

    for (int i = 0; i <= horder; i++)
      {
#ifdef JACOBIxxxx
        JacobiPolynomial (horder-2, 2*y-1, max(0,2*i-2), 0, polsy.Row(i).Addr(1));
#else
        LegendrePolynomial (horder-2, 2*y-1, polsy.Row(i).Addr(1));
#endif
        for (int j = 1; j < horder; j++)
          polsy(i,j) *= y;
        polsy(i,0) = 1.0;
      }

    AutoDiff<1> hv = 1.0;
    for (int i = 0; i <= order; i++, hv *= (1-y) )
      powy[i] = hv;
    

    int ndof2d = trig2tensor.Size();
    for (int i = 0; i < ndof2d; i++)
      {
        int power = trig2tensor[i][1]+trig2tensor[i][2];

        AutoDiff<1> uy = 1.0, vy = 1.0, wy = 1.0;
        AutoDiff<1> fac[3] = { polsy(power, trig2tensor[i][0]),
                               powy[trig2tensor[i][1]],
                               powy[trig2tensor[i][2]] };

        /*
        for (int j = 0; j < split2da[i]; j++)
          wy *= fac[j];
        for (int j = split2da[i]; j < split2d[i]; j++)
          uy *= fac[j];
        for (int j = split2d[i]; j < 3; j++)
          vy *= fac[j];
        */
        for (int j = 0; j < split2d[i]-1; j++)
          wy *= fac[j];
        for (int j = max (0, split2d[i]-1); j < split2d[i]; j++)
          uy *= fac[j];
        for (int j = split2d[i]; j < 3; j++)
          vy *= fac[j];


        AutoDiff<1> prod = uy*vy*wy;

        facy(i)(0) = prod.Value();
        facy(i)(1) = wy.Value() * (uy.Value() * vy.DValue(0) - uy.DValue(0) * vy.Value());
        facy(i)(2) = prod.DValue(0);
      }
  }



  template <typename Tres>
  void CalcZFactor (AutoDiff<1> z, Tres & facz, LocalHeap & lh) const
  {
    HeapReset hr(lh);

    int horder = max (order+1, 2);
    FlatVector<AutoDiff<1> > polz(horder+1, lh);
    Vec<2, AutoDiff<1> > powz;
    
    LegendrePolynomialMult (horder-2, 2*z-1, z, polz.Addr(1));
    polz[0] = 1;

    powz(0) = 1;
    powz(1) = (1-z); 

    int ndof1d = segm2tensor.Size();

    for (int i = 0; i < ndof1d; i++)
      {
        if (split1d[i] == 0)
          {
            AutoDiff<1> uz = 1.0;
            AutoDiff<1> vz = polz(segm2tensor[i][0]) * powz(segm2tensor[i][1]);
        
            facz(i)(0) =  vz.Value();
            facz(i)(1) =  vz.DValue(0);
            facz(i)(2) =  vz.DValue(0);
          }
        else
          { // additional info: segm2tensor[i][0] = segm2tensor[i][1] = 1
            facz(i) (0) = 0;
            facz(i) (1) = -1;
            facz(i) (2) = 0;
          }

        /*
        AutoDiff<1> uz = 1.0;
        AutoDiff<1> vz = 1.0;
        
        if (split1d[i] > 0)
          uz *= polz(segm2tensor[i][0]);
        else
          vz *= polz(segm2tensor[i][0]);
        
        vz *= powz(segm2tensor[i][1]);
        
        facz(i)(0) =  uz.Value() * vz.Value();
        facz(i)(1) =  uz.Value() * vz.DValue(0) - uz.DValue(0) * vz.Value();
        facz(i)(2) =  uz.Value() * vz.DValue(0) + uz.DValue(0) * vz.Value();

        if (split1d[i] == 1)
          {
            facz(i)(0) = 0;
            facz(i)(2) = 0;
          }
        */
      }
  }



};









#endif
