/*********************************************************************/
/* File:   h1hofe.cpp                                                */
/* Author: Start                                                     */
/* Date:   1. June 2007                                              */
/*********************************************************************/


/*

High order integrators in tensor product form

*/

// #define NOPROFILE
// #define CHECK_RANGE

#include <fem.hpp>
namespace ngfem
{
  using namespace ngfem;



  // ******************************** source integrator **************************************






  
  template <int D>
  class SourceEdgeIntegratorTP : public SourceEdgeIntegrator<D>
  {

  public:
    SourceEdgeIntegratorTP (CoefficientFunction * ac1,
                            CoefficientFunction * ac2,
                            CoefficientFunction * ac3)
      : SourceEdgeIntegrator<D> (ac1, ac2, ac3)
    { ; }

    static Integrator * Create (ARRAY<CoefficientFunction*> & coeffs)
    {
      return new SourceEdgeIntegratorTP (coeffs[0], coeffs[1], coeffs[2]);
    }


    virtual void
    AssembleElementVector (const FiniteElement & bfel,
			   const ElementTransformation & eltrans, 
			   FlatVector<double> & elvec,
			   LocalHeap & lh) const
    {
      const HCurlHighOrderTetTP * fel_tet = dynamic_cast<const HCurlHighOrderTetTP*> (&bfel);      
      if (fel_tet)
        {
          // SourceEdgeIntegrator<D>::AssembleElementVector(bfel, eltrans, elvec, lh);
          // *testout << "std, elvec = " << endl << elvec << endl;

          TAssembleElementVector (*fel_tet, eltrans, elvec, lh);
          return;
        }

      const HCurlHighOrderPrismTP * fel_prism = dynamic_cast<const HCurlHighOrderPrismTP*> (&bfel);      
      if (fel_prism)
        if (fel_prism->Sort(0)+3 == fel_prism->Sort(3))  // same singular vertex at bottom and top
          {
            // SourceEdgeIntegrator<D>::AssembleElementVector(bfel, eltrans, elvec, lh);
            // *testout << "std, elvec = " << endl << elvec << endl;

            TAssembleElementVector (*fel_prism, eltrans, elvec, lh);

            // *testout << "fast, elvec = " << endl << elvec << endl;
            return;
          }
      
      /*
        const H1HighOrderHexTP * fel_hex = dynamic_cast<const H1HighOrderHexTP*> (&bfel);      
        if (fel_hex)
        {
        TAssembleElementVector (*fel_hex, eltrans, elvec, lh);
        return;
        }
      */

      SourceEdgeIntegrator<D>::AssembleElementVector(bfel, eltrans, elvec, lh);
    }






    template <class FEL>
    void TAssembleElementVector (const FEL & fel,
                                 const ElementTransformation & eltrans, 
                                 FlatVector<double> & elvec,
                                 LocalHeap & lh) const
    {
      int ndof = fel.GetNDof();
      

      elvec.AssignMemory (ndof, lh);
      elvec = 0;


      HeapReset hr (lh);

      int order = fel.Order();
      int horder = max(order+1,2);


      IntegrationRuleTP<3> irtp(eltrans, 2*order, 1, lh);
          
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      const IntegrationRule & irz = irtp.GetIRZ();
      
      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nipz = irz.GetNIP();
      int nip = nipx * nipy * nipz;

      FlatArray<int> map2dto1d = fel.Get2dTo1dMapping ();
      FlatArray<int> map3dto2d = fel.Get3dTo2dMapping ();

      int ndof2d = fel.GetTrig2TensorMapping().Size();
      int ndof1d = fel.GetSegm2TensorMapping().Size();


      FlatVector<Vec<3> > dvecs(nip, lh);
      
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<3,3> sip(ip, eltrans, irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Vec<3> dvec;
          (this->dvecop).GenerateVector (fel, sip, dvec, lh);

          dvec *= fabs (sip.GetJacobiDet()) * ip.Weight();

          Mat<3> trans = irtp.GetDuffyJacobian (ii) * sip.GetJacobianInverse();
          dvecs[ii] =  trans * dvec;
        }
      

      FlatMatrix<Vec<3> > xfactor(ndof, nipx, lh);
      FlatMatrix<Vec<3> > yfactor(ndof2d, nipy, lh);
      FlatMatrix<Vec<3> > zfactor(ndof1d, nipz, lh);

      for (int ix = 0; ix < nipx; ix++)
        {
          AutoDiff<1> x (irx[ix](0), 0);
          fel.CalcXFactor (x, xfactor.Col(ix), lh);
        }

      for (int iy = 0; iy < nipy; iy++)
        {
          AutoDiff<1> y (iry[iy](0), 0);
          fel.CalcYFactor (y, yfactor.Col(iy), lh);
        }

      for (int iz = 0; iz < nipz; iz++)
        {
          AutoDiff<1> z (irz[iz](0), 0);
          fel.CalcZFactor (z, zfactor.Col(iz), lh);
        }


      FlatVector<Vec<2> > sum2d(ndof2d, lh);
      FlatVector<Vec<3> > sum1d(ndof1d, lh);

      for (int i1 = 0, ii = 0; i1 < nipx; i1++)
        {
          sum2d = 0;
        
          for (int i2 = 0; i2 < nipy; i2++)
            {
              sum1d = 0;
              for (int i3 = 0; i3 < nipz; i3++, ii++)
                {
                  for (int i = 0; i < ndof1d; i++)
                    {
                      const Vec<3> & zf = zfactor(i, i3);

                      sum1d(i)(0) += dvecs(ii)(0) * zf(0);
                      sum1d(i)(1) += dvecs(ii)(1) * zf(0);
                      sum1d(i)(2) += dvecs(ii)(2) * zf(1);
                    }
                }
              // *testout << "sum1d = " << endl << sum1d << endl;

              for (int i = 0; i < ndof2d; i++)
                {
                  const Vec<3> & hv = sum1d(map2dto1d[i]);
                  const Vec<3> & yf = yfactor(i, i2);

                  sum2d(i)(0) += yf(0) * hv(0);
                  sum2d(i)(1) += yf(1) * hv(1) + yf(0) * hv(2);
                }
            }
          
          // *testout << "sum2d = " << endl << sum2d << endl;
        
          for (int i = 0; i < ndof; i++)
            {
              const Vec<2> & hv = sum2d(map3dto2d[i]);
              const Vec<3> & xf = xfactor(i, i1);
              
              elvec(i) += xf(1) * hv(0) + xf(0) * hv(1);
            }
        }

      // *testout << "fast, elvec = " << endl << elvec << endl;
    }
  };













  
  template <int D>
  class CurlEdgeIntegratorTP : public CurlEdgeIntegrator<D>
  {

  public:
    CurlEdgeIntegratorTP (CoefficientFunction * ac1,
                          CoefficientFunction * ac2,
                          CoefficientFunction * ac3)
      : CurlEdgeIntegrator<D> (ac1, ac2, ac3)
    { ; }

    static Integrator * Create (ARRAY<CoefficientFunction*> & coeffs)
    {
      return new CurlEdgeIntegratorTP (coeffs[0], coeffs[1], coeffs[2]);
    }


    virtual void
    AssembleElementVector (const FiniteElement & bfel,
			   const ElementTransformation & eltrans, 
			   FlatVector<double> & elvec,
			   LocalHeap & lh) const
    {
      const HCurlHighOrderTetTP * fel_tet = dynamic_cast<const HCurlHighOrderTetTP*> (&bfel);      
      if (fel_tet)
        {
          TAssembleElementVector (*fel_tet, eltrans, elvec, lh);
          return;
        }

      const HCurlHighOrderPrismTP * fel_prism = dynamic_cast<const HCurlHighOrderPrismTP*> (&bfel);      
      if (fel_prism)
        if (fel_prism->Sort(0)+3 == fel_prism->Sort(3))  // same singular vertex at bottom and top
          {
            TAssembleElementVector (*fel_prism, eltrans, elvec, lh);
            return;
          }
      
      /*
        const H1HighOrderHexTP * fel_hex = dynamic_cast<const H1HighOrderHexTP*> (&bfel);      
        if (fel_hex)
        {
        TAssembleElementVector (*fel_hex, eltrans, elvec, lh);
        return;
        }
      */

      CurlEdgeIntegrator<D>::AssembleElementVector(bfel, eltrans, elvec, lh);
    }






    template <class FEL>
    void TAssembleElementVector (const FEL & fel,
                                 const ElementTransformation & eltrans, 
                                 FlatVector<double> & elvec,
                                 LocalHeap & lh) const
    {
      int ndof = fel.GetNDof();
      
      elvec.AssignMemory (ndof, lh);
      elvec = 0;


      HeapReset hr (lh);

      int order = fel.Order();
      int horder = max(order+1,2);


      IntegrationRuleTP<3> irtp(eltrans, 2*order, 1, lh);
          
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      const IntegrationRule & irz = irtp.GetIRZ();
      
      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nipz = irz.GetNIP();
      int nip = nipx * nipy * nipz;

      FlatArray<int> map2dto1d = fel.Get2dTo1dMapping ();
      FlatArray<int> map3dto2d = fel.Get3dTo2dMapping ();

      int ndof2d = fel.GetTrig2TensorMapping().Size();
      int ndof1d = fel.GetSegm2TensorMapping().Size();

      
      FlatVector<Vec<3> > dvecs(nip, lh);
      
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<3,3> sip(ip, eltrans, irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Vec<3> dvec;
          (this->dvecop).GenerateVector (fel, sip, dvec, lh);

          Mat<3> djac = irtp.GetDuffyJacobian (ii);
          Mat<3> trans = sip.GetJacobian() * Inv(djac);

          dvec *= fabs (sip.GetJacobiDet()) * ip.Weight() * Det(djac)/sip.GetJacobiDet();
          
          dvecs[ii] =  Trans(trans) * dvec;
        }
      

      FlatArray<int> split2d = fel.GetSplit2d();

      int firstsplit1 = ndof2d;
      for (int i = ndof2d-1; i >= 0; i--)
        if (split2d[i] > 0) firstsplit1 = i;

      int firstsplit2 = ndof2d;
      for (int i = ndof2d-1; i >= 0; i--)
        if (split2d[i] == 2) firstsplit2 = i;



      FlatMatrix<Vec<3> > xfactor(ndof, nipx, lh);
      FlatMatrix<Vec<3> > yfactor(ndof2d, nipy, lh);
      FlatMatrix<Vec<3> > zfactor(ndof1d, nipz, lh);

      for (int ix = 0; ix < nipx; ix++)
        {
          AutoDiff<1> x (irx[ix](0), 0);
          fel.CalcXFactor (x, xfactor.Col(ix), lh);
        }

      for (int iy = 0; iy < nipy; iy++)
        {
          AutoDiff<1> y (iry[iy](0), 0);
          fel.CalcYFactor (y, yfactor.Col(iy), lh);
        }

      for (int iz = 0; iz < nipz; iz++)
        {
          AutoDiff<1> z (irz[iz](0), 0);
          fel.CalcZFactor (z, zfactor.Col(iz), lh);
        }

      FlatMatrix<Vec<3> > hxfactor(ndof, nipx, lh);
      FlatMatrix<Vec<3> > hyfactor(ndof2d, nipy, lh);
      FlatMatrix<Vec<3> > hzfactor(ndof1d, nipz, lh);

      for (int ix = 0; ix < nipx; ix++)
        for (int j = 0; j < ndof; j++)
          {
            hxfactor(j, ix)(0) = xfactor(j, ix)(0);
            hxfactor(j, ix)(1) = xfactor(j, ix)(1) - xfactor(j, ix)(2);
            hxfactor(j, ix)(2) = 0.5 * (xfactor(j, ix)(1) + xfactor(j, ix)(2));
          }

      for (int iy = 0; iy < nipy; iy++)
        for (int j = 0; j < ndof2d; j++)
          {
            hyfactor(j, iy)(0) = yfactor(j, iy)(0);
            hyfactor(j, iy)(1) = yfactor(j, iy)(1) - yfactor(j, iy)(2);
            hyfactor(j, iy)(2) = 0.5 * (yfactor(j, iy)(1) + yfactor(j, iy)(2));
          }


      FlatVector<Vec<3> > sum2d(ndof2d, lh);
      FlatVector<Vec<3> > sum1d(ndof1d, lh);
      sum2d = 0.0;

      for (int i1 = 0, ii = 0; i1 < nipx; i1++)
        {
          sum2d = 0;
        
          for (int i2 = 0; i2 < nipy; i2++)
            {
              sum1d = 0;
              for (int i3 = 0; i3 < nipz; i3++, ii++)
                {
                  for (int i = 0; i < ndof1d; i++)
                    {
                      const Vec<3> & zfvec = zfactor(i, i3);

                      sum1d(i)(0) += dvecs(ii)(0) * zfvec(1);
                      sum1d(i)(1) += dvecs(ii)(1) * zfvec(1);
                      sum1d(i)(2) += dvecs(ii)(2) * zfvec(0);
                    }
                }

              // *testout << "sum1d = " << endl << sum1d << endl;

              for (int i = 0; i < firstsplit1; i++)
                {
                  const Vec<3> & hv = sum1d(map2dto1d[i]);
                  const Vec<3> & yf = hyfactor(i, i2);
                  sum2d(i)(1) += yf(0) * hv(1) - yf(2) * hv(2);
                }

              for (int i = firstsplit1; i < firstsplit2; i++)
                {
                  const Vec<3> & hv = sum1d(map2dto1d[i]);
                  const Vec<3> & yf = hyfactor(i, i2);

                  sum2d(i)(0) += -yf(1) * hv(0);
                  sum2d(i)(1) += yf(0) * hv(1) - yf(2) * hv(2);
                  sum2d(i)(2) += yf(1) * hv(2);
                }

              for (int i = firstsplit2; i < ndof2d; i++)
                {
                  const Vec<3> & hv = sum1d(map2dto1d[i]);
                  const Vec<3> & yf = hyfactor(i, i2);

                  sum2d(i)(0) += -yf(1) * hv(0);
                  sum2d(i)(1) += yf(0) * hv(1); 
                }
            }

          // *testout << "yfac = " << endl << yfactor << endl;
          // *testout << "sum2d = " << endl << sum2d << endl;
        
          for (int i = 0; i < ndof; i++)
            {
              int i2d = map3dto2d[i];

              const Vec<3> & hv = sum2d(i2d);
              const Vec<3> & xf = hxfactor(i, i1);
              
              switch (split2d[i2d])
                {
                case 0:
                  {
                    elvec(i) += xf(1) * hv(1);
                    break;
                  }
                case 1:
                  {
                    elvec(i) += InnerProduct (xf, hv);
                    break;
                  }
                case 2:
                  {
                    elvec(i) += xf(0) * hv(0) + xf(1) * hv(1);
                    break;
                  }
                }
            }
        }
    }
  };







  // ******************************** CurlCurl **************************************









  class TrafoCurlZ
  {
    const Vec<3> & z;
  public:

    TrafoCurlZ (const Vec<3> & hz) : z(hz) { ; }

    template <typename TIN, typename TOUT>
    void Transform (const TIN & in, TOUT & out)
    {
      out(0) = z(1) * in(0);
      out(1) = z(1) * in(1);
      out(2) = z(0) * in(2);
    }

    template <typename TIN, typename TOUT>
    void TransformAdd (const TIN & in, TOUT & out)
    {
      out(0) += z(1) * in(0);
      out(1) += z(1) * in(1);
      out(2) += z(0) * in(2);
    }


    template <typename TIN, typename TOUT>
    void TransformAddT (const TIN & in, TOUT & out)
    {
      out(0) += z(1) * in(0);
      out(1) += z(1) * in(1);
      out(2) += z(0) * in(2);
    }

  };


  class TrafoCurlY
  {
    const Vec<3> & y;
  public:
    
    TrafoCurlY (const Vec<3> & hy) : y(hy) { ; }
    
    template <typename TIN, typename TOUT>
    void Transform (const TIN & in, TOUT & out)
    {
      out(0) = -y(1) * in(0);
      out(1) = y(0) * in(1) - y(2) * in(2);
      out(2) = y(1) * in(2);
    }
    
    template <typename TIN, typename TOUT>
    void TransformAdd (const TIN & in, TOUT & out)
    {
      out(0) -= y(1) * in(0);
      out(1) += y(0) * in(1) - y(2) * in(2);
      out(2) += y(1) * in(2);
    }

    template <typename TIN, typename TOUT>
    void TransformAddT (const TIN & in, TOUT & out)
    {
      out(0) -= y(1) * in(0);
      out(1) += y(0) * in(1);
      out(2) += y(1) * in(2) - y(2) * in(1);
    }
  };


  class TrafoCurlX
  {
    const Vec<3> & x;
  public:
    
    TrafoCurlX (const Vec<3> & hx) : x(hx) { ; }
    
    template <typename TIN, typename TOUT>
    void Transform (const TIN & in, TOUT & out)
    {
      out = x(0) * in(0) + x(1) * in(1) + x(2) * in(2);
    }
    
    template <typename TIN, typename TOUT>
    void TransformAdd (const TIN & in, TOUT & out)
    {
      out += x(0) * in(0) + x(1) * in(1) + x(2) * in(2);
    }

    template <typename TIN, typename TOUT>
    void TransformAddT (const TIN & in, TOUT & out)
    {
      out(0) += x(0) * in;
      out(1) += x(1) * in;
      out(2) += x(2) * in;
    }

  };








  // type 1:

  class TrafoCurlY1
  {
    const Vec<3> & y;
  public:
    
    TrafoCurlY1 (const Vec<3> & hy) : y(hy) { ; }
    
    template <typename TIN, typename TOUT>
    void Transform (const TIN & in, TOUT & out)
    {
      out = y(0) * in(1) - y(2) * in(2);
    }
    
    template <typename TIN, typename TOUT>
    void TransformAdd (const TIN & in, TOUT & out)
    {
      out += y(0) * in(1) - y(2) * in(2);
    }

    template <typename TIN, typename TOUT>
    void TransformAddT (const TIN & in, TOUT & out)
    {
      out(1) += y(0) * in;
      out(2) -= y(2) * in;
    }

  };


  class TrafoCurlX1
  {
    const Vec<3> & x;
  public:
    
    TrafoCurlX1 (const Vec<3> & hx) : x(hx) { ; }
    
    template <typename TIN, typename TOUT>
    void Transform (const TIN & in, TOUT & out)
    {
      out = x(1) * in;
    }
    
    template <typename TIN, typename TOUT>
    void TransformAdd (const TIN & in, TOUT & out)
    {
      out += x(1) * in;
    }

    template <typename TIN, typename TOUT>
    void TransformAddT (const TIN & in, TOUT & out)
    {
      out += x(1) * in;
    }

  };













  
  template <int D>
  class CurlCurlIntegratorTP : public CurlCurlEdgeIntegrator<D>
  {

  public:
    CurlCurlIntegratorTP (CoefficientFunction * acoef)
      : CurlCurlEdgeIntegrator<D> (acoef)
    { ; }

    static Integrator * Create (ARRAY<CoefficientFunction*> & coeffs)
    {
      return new CurlCurlIntegratorTP (coeffs[0]);
    }

    virtual void
    AssembleElementMatrix (const FiniteElement & bfel,
			   const ElementTransformation & eltrans, 
			   FlatMatrix<double> & elmat,
			   LocalHeap & lh) const
    {

      const HCurlHighOrderTetTP * fel_tet = dynamic_cast<const HCurlHighOrderTetTP*> (&bfel);      
      const HCurlHighOrderPrismTP * fel_prism = dynamic_cast<const HCurlHighOrderPrismTP*> (&bfel);      
      // const HCurlHighOrderHexTP * fel_hex = dynamic_cast<const HCurlHighOrderHexTP*> (&bfel);      


      if (fel_tet)
        {
          TAssembleElementMatrix (*fel_tet, eltrans, elmat, lh);
          return;
        }

      if (fel_prism)
        if (fel_prism->Sort(0)+3 == fel_prism->Sort(3))  // same singular vertex at bottom and top
          {
            TAssembleElementMatrix (*fel_prism, eltrans, elmat, lh);
            return;
          }
      /*
        if (fel_hex)
        {
        TAssembleElementMatrix (*fel_hex, eltrans, elmat, lh);
        return;
        }
      */
      CurlCurlEdgeIntegrator<D>::AssembleElementMatrix(bfel, eltrans, elmat, lh);
    }





    // #define OLDVER
#define NEWVER

#ifdef NEWVER

    template <class FEL>
    void TAssembleElementMatrix (const FEL & fel,
                                 const ElementTransformation & eltrans, 
                                 FlatMatrix<double> & elmat,
                                 LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("Fast assembling H(curl)");

      static int timeri = NgProfiler::CreateTimer ("Fast assembling - init");
      static int timerg1 = NgProfiler::CreateTimer ("Fast assembling - geom1");
      static int timerg2 = NgProfiler::CreateTimer ("Fast assembling - geom2");


      NgProfiler::RegionTimer reg (timer);


      NgProfiler::StartTimer (timerg1);      

      int ndof = fel.GetNDof();
      elmat.AssignMemory (ndof, ndof, lh);
      elmat = 0.0;

      HeapReset hr (lh);

      
      int order = fel.Order();
      // int horder = max(order+1,2);


      FlatArray<int> map2dto1d = fel.Get2dTo1dMapping ();
      FlatArray<int> map3dto2d = fel.Get3dTo2dMapping ();

      int ndof2d = fel.GetNDof2d();
      int ndof1d = fel.GetNDof1d();

      
      FlatArray<int> split2d = fel.GetSplit2d();

          
      int intorder = 2 * order;
      if (fel.ElementType() == ET_TET) intorder -= 2;
      
      IntegrationRuleTP<3> irtp(eltrans, intorder, 1, lh);
      
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      const IntegrationRule & irz = irtp.GetIRZ();
      
      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nipz = irz.GetNIP();
      int nip = nipx * nipy * nipz;
      
      NgProfiler::StopTimer (timerg1);            
      NgProfiler::StartTimer (timerg2);      

      FlatVector<Mat<3> > tdmats(nip, lh);      // transformed D-matrix
      
      // geometry calculation
      for (int i = 0; i < nip; i++)
        {
          IntegrationPoint ip (irtp.GetXi(i), irtp.GetWeight(i));
          SpecificIntegrationPoint<3,3> sip(ip, eltrans, irtp.GetPoint(i), irtp.GetJacobian(i), lh);
          
          Mat<3> dmat;
          (this->dmatop).GenerateMatrix (fel, sip, dmat, lh);
          
          Mat<3> djac = irtp.GetDuffyJacobian (i);
          Mat<3> trans = sip.GetJacobian() * Inv(djac);
          
          dmat *= fabs (sip.GetJacobiDet()) * sip.IP().Weight() * sqr (Det(djac)/sip.GetJacobiDet());

          Mat<3> hmat = dmat * trans;
          tdmats[i] = Trans (trans) * hmat; 
        }

      NgProfiler::StartTimer (timeri);
      NgProfiler::StopTimer (timerg2);      


      FlatMatrix<Vec<3> > xfactor(ndof, nipx, lh);
      FlatMatrix<Vec<3> > yfactor(ndof2d, nipy, lh);
      FlatMatrix<Vec<3> > zfactor(ndof1d, nipz, lh);
      
      for (int i1 = 0; i1 < nipz; i1++)
        {
          AutoDiff<1> z (irz[i1](0), 0);
          fel.CalcZFactor (z, zfactor.Col(i1), lh);
        }
      
      for (int i1 = 0; i1 < nipy; i1++)
        {
          AutoDiff<1> y (iry[i1](0), 0);
          fel.CalcYFactor (y, yfactor.Col(i1), lh);
        }
      
      for (int i1 = 0; i1 < nipx; i1++)
        {
          AutoDiff<1> x (irx[i1](0), 0);
          fel.CalcXFactor (x, xfactor.Col(i1), lh);
        }

      FlatMatrix<Vec<3> > hxfactor(ndof, nipx, lh);
      FlatMatrix<Vec<3> > hyfactor(ndof2d, nipy, lh);
      // FlatMatrix<Vec<3> > hzfactor(ndof1d, nipz, lh);

      for (int j = 0; j < ndof; j++)
        for (int ix = 0; ix < nipx; ix++)
          {
            hxfactor(j, ix)(0) = xfactor(j, ix)(0);
            hxfactor(j, ix)(1) = xfactor(j, ix)(1) - xfactor(j, ix)(2);
            hxfactor(j, ix)(2) = 0.5 * (xfactor(j, ix)(1) + xfactor(j, ix)(2));
          }

      for (int j = 0; j < ndof2d; j++)
        for (int iy = 0; iy < nipy; iy++)
          {
            hyfactor(j, iy)(0) = yfactor(j, iy)(0);
            hyfactor(j, iy)(1) = yfactor(j, iy)(1) - yfactor(j, iy)(2);
            hyfactor(j, iy)(2) = 0.5 * (yfactor(j, iy)(1) + yfactor(j, iy)(2));
          }



      NgProfiler::StopTimer (timeri);      


      static int timerc = NgProfiler::CreateTimer ("Fast assembling - compute");
      NgProfiler::StartTimer (timerc);

      int nipxy = nipx *nipy;

      FlatVector<Mat<3> > allsum2d(nipx, lh);
      FlatVector<Mat<3> > allsum1d(nipxy, lh);

      FlatVector<Vec<3> > mult3d(nipx, lh);
      FlatVector<Mat<3> > mult2d(nipxy, lh);
      FlatVector<Mat<3> > mult1d(nip, lh);


      FlatVector<Vec<3> > mult2da(nipxy, lh);
      FlatVector<Vec<3> > allsum2da(nipx, lh);
      FlatVector<> mult3da(nipx, lh);

      FlatVector<> mult2daa(nipxy, lh);
      FlatVector<> allsum2daa(nipx, lh);


      int firstsplit1 = ndof2d;
      for (int i = ndof2d-1; i >= 0; i--)
        if (split2d[i] > 0) firstsplit1 = i;

      int firstsplit2 = ndof2d;
      for (int i = ndof2d-1; i >= 0; i--)
        if (split2d[i] == 2) firstsplit2 = i;


      FlatArray<int> map2dto3d(ndof2d * (order+1), lh);
      FlatArray<int> nmap2dto3d(ndof2d, lh);
      nmap2dto3d = 0;
      for (int j = 0; j < ndof; j++)
        {
          if (fel.IsGradient(j)) continue; 
          int j2d = map3dto2d[j];
          map2dto3d[(order+1)*j2d + nmap2dto3d[j2d]++] = j;
        }

      FlatArray<int> map1dto2d(ndof1d * 2*(order+1), lh);
      FlatArray<int> nmap1dto2d(ndof1d, lh);
      nmap1dto2d = 0;
      for (int j = 0; j < ndof2d; j++)
        {
          int j1d = map2dto1d[j];
          map1dto2d[2*(order+1)*j1d + nmap1dto2d[j1d]++] = j;
        }

      for (int i = 0; i < ndof2d; i++)
        {
          if (nmap2dto3d[i] > order+1) 
            {
              cout << "too many nmap2dto3d" << endl;
              exit(1);
            }
        }
      for (int i = 0; i < ndof1d; i++)
        {
          if (nmap1dto2d[i] > 2*(order+1)) 
            {
              cout << "too many nmap1dto2d" << endl;
              exit(1);
            }
        }

      /*
      NgProfiler::AddFlops (timerc, 9*nip*ndof1d);
      NgProfiler::AddFlops (timerc, 9*nip*ndof1d*(ndof1d+1)/2);

      NgProfiler::AddFlops (timerc, 9*nipxy*ndof1d*ndof2d);
      NgProfiler::AddFlops (timerc, 3*nipxy*ndof2d*(ndof2d+1));

      NgProfiler::AddFlops (timerc, 4*ndof2d*ndof*nipx);
      NgProfiler::AddFlops (timerc, ndof*(ndof+1)*nipx);
      */

      for (int kz = 0; kz < ndof1d; kz++)
        {
          for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
            for (int iz = 0; iz < nipz; iz++, ii++)
              for (int l = 0; l < 3; l++)
                TrafoCurlZ(zfactor(kz,iz)) . Transform(tdmats[ii].Col(l), mult1d(ii).Col(l));
          
          for (int jz = 0; jz <= kz; jz++)
            {
              for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
                {
                  Mat<3> sum (0.0);

                  for (int iz = 0; iz < nipz; iz++, ii++)
                    for (int l = 0; l < 3; l++)
                      TrafoCurlZ(zfactor(jz,iz)) . TransformAdd(mult1d(ii).Row(l), sum.Row(l));
                  
                  allsum1d(ixy) = sum;
                }

              for (int ky = 0; ky < nmap1dto2d[kz]; ky++)
                {
                  int kyz = map1dto2d[2*(order+1)*kz+ky];
                

                  if (kyz < firstsplit1)
                    {
                      
                      for (int ix = 0, ixy=0; ix < nipx; ix++)
                        for (int iy = 0; iy < nipy; iy++, ixy++)
                          for (int l = 0; l < 3; l++)
                            TrafoCurlY1 (hyfactor(kyz,iy)) . Transform(allsum1d(ixy).Col(l), mult2da(ixy)(l));
                      
                      
                      for (int jy = 0; jy < nmap1dto2d[jz]; jy++)
                        {
                          if (kz == jz && jy > ky) continue;
                          int jyz = map1dto2d[2*(order+1)*jz+jy];

                          if (jyz < firstsplit1)
                            {
                              for (int ix = 0; ix < nipx; ix++)
                                {
                                  double hsum2d = 0;
                              
                                  for (int iy = 0; iy < nipy; iy++)
                                    TrafoCurlY1(hyfactor(jyz,iy)) . TransformAdd(mult2da(ix*nipy+iy), hsum2d);
                              
                                  allsum2daa(ix) = hsum2d;
                                }
                          
                              int n1dxj = nmap2dto3d[jyz];
                              int n1dxk = nmap2dto3d[kyz];
                          
                              for (int kx = 0; kx < n1dxk; kx++)
                                {
                                  int kk = map2dto3d[kyz*(order+1)+kx];
                              
                                  for (int ix = 0; ix < nipx; ix++)
                                    TrafoCurlX1 (hxfactor(kk,ix)) . Transform(allsum2daa(ix), mult3da(ix));
                              
                                  for (int jx = 0; jx < n1dxj; jx++)
                                    {
                                      int jj = map2dto3d[jyz*(order+1)+jx];
                                  
                                      double sum = 0;
                                      for (int ix = 0; ix < nipx; ix++)
                                        TrafoCurlX1 (hxfactor(jj,ix)) . TransformAdd(mult3da(ix), sum);
                                  
                                      elmat(kk,jj) = sum;
                                      elmat(jj,kk) = sum;
                                    }
                                }
                            }

                          else
                            {
                              for (int ix = 0; ix < nipx; ix++)
                                {
                                  Vec<3> hsum2da (0.0);
                              
                                  for (int iy = 0; iy < nipy; iy++)
                                    TrafoCurlY(hyfactor(jyz,iy)) . TransformAdd(mult2da(ix*nipy+iy), hsum2da);
                              
                                  allsum2da(ix) = hsum2da;
                                }
                          
                              int n1dxj = nmap2dto3d[jyz];
                              int n1dxk = nmap2dto3d[kyz];
                          
                              for (int kx = 0; kx < n1dxk; kx++)
                                {
                                  int kk = map2dto3d[kyz*(order+1)+kx];
                              
                                  for (int ix = 0; ix < nipx; ix++)
                                    for (int l = 0; l < 3; l++)
                                      TrafoCurlX1 (hxfactor(kk,ix)) . Transform(allsum2da(ix)(l), mult3d(ix)(l));
                              
                              
                                  for (int jx = 0; jx < n1dxj; jx++)
                                    {
                                      int jj = map2dto3d[jyz*(order+1)+jx];
                                  
                                      double sum = 0;
                                      for (int ix = 0; ix < nipx; ix++)
                                        TrafoCurlX (hxfactor(jj,ix)) . TransformAdd(mult3d(ix), sum);
                                  
                                      elmat(kk,jj) = sum;
                                      elmat(jj,kk) = sum;
                                    }
                                }
                            }
                        }
                    }

                  else
                    {
                      for (int ix = 0, ixy=0; ix < nipx; ix++)
                        for (int iy = 0; iy < nipy; iy++, ixy++)
                          for (int l = 0; l < 3; l++)
                            TrafoCurlY (hyfactor(kyz,iy)) . Transform(allsum1d(ixy).Col(l), mult2d(ixy).Row(l));
                      
                      
                      for (int jy = 0; jy < nmap1dto2d[jz]; jy++)
                        {
                          if (kz == jz && jy > ky) continue;
                          int jyz = map1dto2d[2*(order+1)*jz+jy];



                          if (jyz < firstsplit1)
                            {
                              for (int ix = 0; ix < nipx; ix++)
                                {
                                  Vec<3> hsum2d (0.0);
                              
                                  for (int iy = 0; iy < nipy; iy++)
                                    for (int l = 0; l < 3; l++)
                                      TrafoCurlY1(hyfactor(jyz,iy)) . TransformAdd(mult2d(ix*nipy+iy).Col(l), hsum2d(l));
                              
                                  allsum2da(ix) = hsum2d;
                                }
                          
                              int n1dxj = nmap2dto3d[jyz];
                              int n1dxk = nmap2dto3d[kyz];
                          
                              for (int kx = 0; kx < n1dxk; kx++)
                                {
                                  int kk = map2dto3d[kyz*(order+1)+kx];
                              
                                  for (int ix = 0; ix < nipx; ix++)
                                    // for (int l = 0; l < 3; l++)
                                      TrafoCurlX (hxfactor(kk,ix)) . Transform(allsum2da(ix), mult3da(ix));
                              
                              
                                  for (int jx = 0; jx < n1dxj; jx++)
                                    {
                                      int jj = map2dto3d[jyz*(order+1)+jx];
                                  
                                      double sum = 0;
                                      for (int ix = 0; ix < nipx; ix++)
                                        TrafoCurlX1 (hxfactor(jj,ix)) . TransformAdd(mult3da(ix), sum);
                                  
                                      elmat(kk,jj) = sum;
                                      elmat(jj,kk) = sum;
                                    }
                                }
                            }

                          else

                            {
                              for (int ix = 0; ix < nipx; ix++)
                                {
                                  Mat<3> hsum2d (0.0);
                              
                                  for (int iy = 0; iy < nipy; iy++)
                                    for (int l = 0; l < 3; l++)
                                      TrafoCurlY(hyfactor(jyz,iy)) . TransformAdd(mult2d(ix*nipy+iy).Col(l), hsum2d.Row(l));
                              
                                  allsum2d(ix) = Trans(hsum2d);
                                }
                          
                              int n1dxj = nmap2dto3d[jyz];
                              int n1dxk = nmap2dto3d[kyz];
                          
                              for (int kx = 0; kx < n1dxk; kx++)
                                {
                                  int kk = map2dto3d[kyz*(order+1)+kx];
                              
                                  for (int ix = 0; ix < nipx; ix++)
                                    for (int l = 0; l < 3; l++)
                                      TrafoCurlX (hxfactor(kk,ix)) . Transform(allsum2d(ix).Row(l), mult3d(ix)(l));
                              
                              
                                  for (int jx = 0; jx < n1dxj; jx++)
                                    {
                                      int jj = map2dto3d[jyz*(order+1)+jx];
                                  
                                      double sum = 0;
                                      for (int ix = 0; ix < nipx; ix++)
                                        TrafoCurlX (hxfactor(jj,ix)) . TransformAdd(mult3d(ix), sum);
                                  
                                      elmat(kk,jj) = sum;
                                      elmat(jj,kk) = sum;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

      NgProfiler::StopTimer (timerc);


    }

#endif
















#ifdef OLDVER


    template <class FEL>
    void TAssembleElementMatrix (const FEL & fel,
                                 const ElementTransformation & eltrans, 
                                 FlatMatrix<double> & elmat,
                                 LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("Fast assembling - H(curl)");
      static int timeri = NgProfiler::CreateTimer ("Fast assembling - init");
      static int timerg = NgProfiler::CreateTimer ("Fast assembling - geom");
      static int timerz = NgProfiler::CreateTimer ("Fast assembling z");
      static int timery = NgProfiler::CreateTimer ("Fast assembling y");
      static int timery1 = NgProfiler::CreateTimer ("Fast assembling y1");
      static int timery2 = NgProfiler::CreateTimer ("Fast assembling y2");
      static int timery3 = NgProfiler::CreateTimer ("Fast assembling y3");
      static int timerx = NgProfiler::CreateTimer ("Fast assembling x");
      static int timerx1 = NgProfiler::CreateTimer ("Fast assembling x1");
      static int timerx2 = NgProfiler::CreateTimer ("Fast assembling x2");

      NgProfiler::RegionTimer reg (timer);

      int ndof = fel.GetNDof();
      
      elmat.AssignMemory (ndof, ndof, lh);


      // for (int cnt = 0; cnt < 100; cnt++)
        {


      elmat = 0.0;

      HeapReset hr (lh);

      

      int order = fel.Order();
      int horder = max(order+1,2);

      NgProfiler::StartTimer (timerg);

      IntegrationRuleTP<3> irtp(eltrans, 2*order-2, 1, lh);
      NgProfiler::StopTimer (timerg);

      NgProfiler::StartTimer (timeri);

          
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      const IntegrationRule & irz = irtp.GetIRZ();

      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nipz = irz.GetNIP();
      int nip = nipx * nipy * nipz;

      FlatArray<int> map2dto1d = fel.Get2dTo1dMapping ();
      FlatArray<int> map3dto2d = fel.Get3dTo2dMapping ();

      int ndof2d = fel.GetNDof2d(); 
      int ndof1d = fel.GetNDof1d(); 

      
      
      FlatArray<int> split2d = fel.GetSplit2d();
      int firstsplit1 = ndof2d;
      for (int i = ndof2d-1; i >= 0; i--)
        if (split2d[i] > 0) firstsplit1 = i;

      int firstsplit2 = ndof2d;
      for (int i = ndof2d-1; i >= 0; i--)
        if (split2d[i] == 2) firstsplit2 = i;


      firstsplit1 = 0;
      firstsplit2 = ndof2d;


      FlatMatrix<Vec<3> > xfactor(ndof, nipx, lh);
      FlatMatrix<Vec<3> > yfactor(ndof2d, nipy, lh);
      FlatMatrix<Vec<3> > zfactor(ndof1d, nipz, lh);

      for (int ix = 0; ix < nipx; ix++)
        {
          AutoDiff<1> x (irx[ix](0), 0);
          fel.CalcXFactor (x, xfactor.Col(ix), lh);
        }

      for (int iy = 0; iy < nipy; iy++)
        {
          AutoDiff<1> y (iry[iy](0), 0);
          fel.CalcYFactor (y, yfactor.Col(iy), lh);
        }

      for (int iz = 0; iz < nipz; iz++)
        {
          AutoDiff<1> z (irz[iz](0), 0);
          fel.CalcZFactor (z, zfactor.Col(iz), lh);
        }
      

       FlatMatrix<> xfactora(ndof, nipx, lh);

      FlatMatrix<Vec<3> > hxfactor(ndof, nipx, lh);
      FlatMatrix<Vec<3> > hyfactor(ndof2d, nipy, lh);
      FlatMatrix<Vec<3> > hzfactor(ndof1d, nipz, lh);


      for (int j = 0; j < ndof; j++)
        for (int k = 0; k < nipx; k++)
          xfactora(j,k) = xfactor(j,k)(1)-xfactor(j,k)(2);

      for (int j = 0; j < ndof; j++)
        for (int ix = 0; ix < nipx; ix++)
          {
            hxfactor(j, ix)(0) = xfactor(j, ix)(0);
            hxfactor(j, ix)(1) = xfactor(j, ix)(1) - xfactor(j, ix)(2);
            hxfactor(j, ix)(2) = 0.5 * (xfactor(j, ix)(1) + xfactor(j, ix)(2));
          }

      for (int j = 0; j < ndof2d; j++)
        for (int iy = 0; iy < nipy; iy++)
          {
            hyfactor(j, iy)(0) = yfactor(j, iy)(0);
            hyfactor(j, iy)(1) = yfactor(j, iy)(1) - yfactor(j, iy)(2);
            hyfactor(j, iy)(2) = 0.5 * (yfactor(j, iy)(1) + yfactor(j, iy)(2));
          }










      FlatMatrix<Mat<3> > hmat2d(ndof2d, lh);
      FlatMatrix<Vec<3> > hmat2da(ndof2d, firstsplit1, lh);
      FlatMatrix<> hmat2daa(firstsplit1, lh);

      FlatMatrix<Mat<3> > hmat1d(ndof1d, ndof1d, lh);
      FlatMatrix<Mat<3,3> > hmult2d(ndof1d, ndof2d, lh);
      FlatMatrix<Vec<3> > hmult2da(ndof1d, ndof2d, lh);
      FlatMatrix<Vec<3> >   hmult3d(ndof2d, ndof, lh);
      FlatMatrix<>   hmult3da(ndof2d, ndof, lh);
      FlatVector<Mat<3,3> > dbmat1d (ndof1d, lh);
      
      FlatVector<Mat<3> > tdmats(nip, lh);      // transformed D-matrix

      for (int i = 0; i < nip; i++)
        {
          IntegrationPoint ip (irtp.GetXi(i), irtp.GetWeight(i));
          SpecificIntegrationPoint<3,3> sip(ip, eltrans, irtp.GetPoint(i), irtp.GetJacobian(i), lh);
          
          Mat<3> dmat;
          (this->dmatop).GenerateMatrix (fel, sip, dmat, lh);
          
          Mat<3> djac = irtp.GetDuffyJacobian (i);
          Mat<3> trans = sip.GetJacobian() * Inv(djac);
          
          dmat *= fabs (sip.GetJacobiDet()) * sip.IP().Weight() * sqr (Det(djac)/sip.GetJacobiDet());

          Mat<3> hmat = dmat * trans;
          tdmats[i] = Trans (trans) * hmat; 
        }


      NgProfiler::StopTimer (timeri);

      FlatMatrix<Mat<3> > allsum1d(ndof1d*ndof1d, nipx*nipy, lh);
      FlatMatrix<Mat<3> > allsum2d(ndof2d*ndof2d, nipx, lh);
      FlatMatrix<Vec<3> > allsum2da(ndof2d*firstsplit1, nipx, lh);
      FlatMatrix<> allsum2daa(sqr(firstsplit1), nipx, lh);
      
      

      NgProfiler::StartTimer (timerz);
      
      for (int ix = 0, ii = 0; ix < nipx; ix++)
        for (int iy = 0; iy < nipy; iy++)
          {
            hmat1d = 0.0;
            
            for (int iz = 0; iz < nipz; iz++, ii++)
              {
                for (int k = 0; k < ndof1d; k++)
                  {
                    const Vec<3> & zf = zfactor(k, iz);
                    
                    for (int l = 0; l < 3; l++)
                      {
                        dbmat1d(k)(l,0) = tdmats[ii](l,0) * zf(1);   
                        dbmat1d(k)(l,1) = tdmats[ii](l,1) * zf(1);   
                        dbmat1d(k)(l,2) = tdmats[ii](l,2) * zf(0);   
                      }
                  }
                
                for (int j = 0; j < ndof1d; j++)
                  for (int k = 0; k <= j; k++)
                    {
                      const Mat<3,3> & hm = dbmat1d(k);
                      const Vec<3> & zf = zfactor(j, iz);
                      
                      for (int l = 0; l < 3; l++)
                        {
                          hmat1d(j, k)(0,l) += zf(1) * hm(0,l);
                          hmat1d(j, k)(1,l) += zf(1) * hm(1,l);
                          hmat1d(j, k)(2,l) += zf(0) * hm(2,l);
                        }
                    }
              }
            
            for (int j = 0; j < ndof1d; j++)
              {
                for (int k = 0; k <= j; k++)
                  allsum1d(j*ndof1d+k, ix*nipy+iy) = hmat1d(j,k);
                for (int k = j+1; k < ndof1d; k++)
                  allsum1d(j*ndof1d+k, ix*nipy+iy) = Trans (hmat1d(k,j));
              }
          }


      NgProfiler::AddFlops (timerz, nip * ndof1d * 9);
      NgProfiler::AddFlops (timerz, nip * ndof1d*(ndof1d+1) / 2 * 9);
      NgProfiler::StopTimer (timerz);


      *testout << "tdmats = " << endl << tdmats << endl;

      *testout << "allsum1d = " << endl << allsum1d << endl;




      NgProfiler::StartTimer (timery);


      // NgProfiler::StartTimer (timery1);

      for (int ix = 0, ii = 0; ix < nipx; ix++)
        {
          hmat2da = 0.0;
          hmat2daa = 0.0;

          // 2d-type0 requires less ops
          for (int iy = 0; iy < nipy; iy++)
            {
              for (int k2d = 0; k2d < firstsplit1; k2d++)
                {
                  int k1d = map2dto1d[k2d];
                  const Vec<3> & yf = hyfactor(k2d, iy);
                  for (int j1d = 0; j1d < ndof1d; j1d++)
                    {
                      const Mat<3> & hm = allsum1d(j1d+ndof1d*k1d, ix*nipy+iy);
                      for (int l = 0; l < 3; l++)
                        hmult2da(j1d, k2d)(l) = hm(1, l) * yf(0) - hm(2, l) * yf(2);
                    }
                }

              for (int j2d = 0; j2d < firstsplit1; j2d++)
                {
                  int j1d = map2dto1d[j2d];
                  const Vec<3> & yf = hyfactor(j2d, iy);
                  for (int k2d = 0; k2d <= j2d; k2d++)
                    {
                      const Vec<3> & hm = hmult2da(j1d, k2d);
                      hmat2daa(j2d, k2d) += yf(0) * hm(1) - yf(2) * hm(2);
                    }
                }

              for (int j2d = firstsplit1; j2d < firstsplit2; j2d++)
                {
                  const Vec<3> & yf = hyfactor(j2d, iy);
                  int j1d = map2dto1d[j2d];
                  for (int k2d = 0; k2d < firstsplit1; k2d++)
                    {
                      const Vec<3> & hm = hmult2da(j1d, k2d);
                      
                      hmat2da(j2d, k2d)(0) += -yf(1) * hm(0);
                      hmat2da(j2d, k2d)(1) += yf(0) * hm(1) - yf(2) * hm(2);
                      hmat2da(j2d, k2d)(2) += yf(1) * hm(2);
                    }
                }

              for (int j2d = firstsplit2; j2d < ndof2d; j2d++)
                {
                  const Vec<3> & yf = hyfactor(j2d, iy);
                  int j1d = map2dto1d[j2d];
                  for (int k2d = 0; k2d < firstsplit1; k2d++)
                    {
                      const Vec<3> & hm = hmult2da(j1d, k2d);
                      
                      hmat2da(j2d, k2d)(0) += -yf(1) * hm(0);
                      hmat2da(j2d, k2d)(1) += yf(0) * hm(1);
                    }
                }
            }

          for (int k2d = firstsplit1; k2d < ndof2d; k2d++)
            for (int j2d = 0; j2d < firstsplit1; j2d++)
              allsum2da(j2d+firstsplit1*k2d, ix) = hmat2da(k2d, j2d);

          for (int j2d = 0; j2d < firstsplit1; j2d++)
            {
              for (int k2d = 0; k2d <= j2d; k2d++)
                allsum2daa(j2d*firstsplit1+k2d, ix) = hmat2daa(j2d, k2d);
              for (int k2d = j2d+1; k2d < firstsplit1; k2d++)
                allsum2daa(j2d*firstsplit1+k2d, ix) = hmat2daa(k2d, j2d);
            }
        }
      // NgProfiler::StopTimer (timery1);
      

      // NgProfiler::StartTimer (timery2);
      for (int ix = 0, ii = 0; ix < nipx; ix++)
        {
          hmat2d = 0.0;

          // all 2d components
          for (int iy = 0; iy < nipy; iy++)
            {
              for (int k2d = firstsplit1; k2d < firstsplit2; k2d++)
                {
                  int k1d = map2dto1d[k2d];
                  const Vec<3> & yf = hyfactor(k2d, iy);

                  for (int j1d = 0; j1d < ndof1d; j1d++)
                    {
                      const Mat<3> & hm = allsum1d(j1d+ndof1d*k1d, ix*nipy+iy);
                      for (int l = 0; l < 3; l++)
                        {
                          hmult2d(j1d, k2d)(l,0) = hm(0, l) * (-yf(1));
                          hmult2d(j1d, k2d)(l,1) = hm(1, l) * yf(0) - hm(2, l) * yf(2);
                          hmult2d(j1d, k2d)(l,2) = hm(2, l) * yf(1);
                        }
                    }
                }

              for (int k2d = firstsplit2; k2d < ndof2d; k2d++)
                {
                  int k1d = map2dto1d[k2d];
                  const Vec<3> & yf = hyfactor(k2d, iy);

                  for (int j1d = 0; j1d < ndof1d; j1d++)
                    {
                      const Mat<3> & hm = allsum1d(j1d+ndof1d*k1d, ix*nipy+iy);
                      for (int l = 0; l < 3; l++)
                        {
                          hmult2d(j1d, k2d)(l,0) = hm(0, l) * (-yf(1));
                          hmult2d(j1d, k2d)(l,1) = hm(1, l) * yf(0);
                        }
                    }
                }

              for (int j2d = firstsplit1; j2d < firstsplit2; j2d++)
                {
                  const Vec<3> & yf = hyfactor(j2d, iy);
                  int j1d = map2dto1d[j2d];
                  for (int k2d = firstsplit1; k2d <= j2d; k2d++)
                    {
                      const Mat<3,3> & hm = hmult2d(j1d, k2d); // map2dto1d[j2d]);
                      
                      for (int l = 0; l < 3; l++)
                        {
                          hmat2d(j2d, k2d)(0, l) += -yf(1) * hm(0, l);
                          hmat2d(j2d, k2d)(1, l) += yf(0) * hm(1, l) - yf(2) * hm(2, l);
                          hmat2d(j2d, k2d)(2, l) += yf(1) * hm(2, l);
                        }
                    }
                }

              for (int j2d = firstsplit2; j2d < ndof2d; j2d++)
                {
                  const Vec<3> & yf = hyfactor(j2d, iy);
                  int j1d = map2dto1d[j2d];
                  for (int k2d = firstsplit1; k2d < firstsplit2; k2d++)
                    {
                      const Mat<3,3> & hm = hmult2d(j1d, k2d); 
                      
                      for (int l = 0; l < 3; l++)
                        {
                          hmat2d(j2d, k2d)(0, l) += -yf(1) * hm(0, l);
                          hmat2d(j2d, k2d)(1, l) += yf(0) * hm(1, l);
                        }
                    }

                  for (int k2d = firstsplit2; k2d <= j2d; k2d++)
                    {
                      const Mat<3,3> & hm = hmult2d(j1d, k2d); 
                      
                      for (int l = 0; l < 2; l++)
                        {
                          hmat2d(j2d, k2d)(0, l) += -yf(1) * hm(0, l);
                          hmat2d(j2d, k2d)(1, l) += yf(0) * hm(1, l);
                        }
                    }
                }
            }

          for (int j2d = firstsplit1; j2d < ndof2d; j2d++)
            {
              for (int k2d = firstsplit1; k2d <= j2d; k2d++)
                allsum2d(j2d*ndof2d+k2d, ix) = hmat2d(j2d, k2d);
              for (int k2d = j2d+1; k2d < ndof2d; k2d++)
                allsum2d(j2d*ndof2d+k2d, ix) = Trans (hmat2d(k2d, j2d));
            }
        }
      // NgProfiler::StopTimer (timery2);


      NgProfiler::AddFlops (timery, nipx*nipy*firstsplit1*ndof1d * 6);
      NgProfiler::AddFlops (timery, nipx*nipy*(ndof2d-firstsplit1)*ndof1d * 15);

      NgProfiler::AddFlops (timery, nipx*nipy*firstsplit1*(firstsplit1+1)/2 * 2);
      NgProfiler::AddFlops (timery, nipx*nipy*firstsplit1*ndof2d * 5);
      NgProfiler::AddFlops (timery, nipx*nipy*(ndof2d-firstsplit1)*(ndof2d-firstsplit1+1)/2 * 15);


      NgProfiler::AddFlops (timery1, nipx*nipy*firstsplit1*ndof1d * 6);
      NgProfiler::AddFlops (timery2, nipx*nipy*(ndof2d-firstsplit1)*ndof1d * 15);

      NgProfiler::AddFlops (timery1, nipx*nipy*firstsplit1*(firstsplit1+1)/2 * 2);
      NgProfiler::AddFlops (timery1, nipx*nipy*firstsplit1*ndof2d * 5);
      NgProfiler::AddFlops (timery2, nipx*nipy*(ndof2d-firstsplit1)*(ndof2d-firstsplit1+1)/2 * 15);

      NgProfiler::StopTimer (timery);



      *testout << "allsum2d = " << endl << allsum2d << endl;

      /*      
      NgProfiler::StartTimer (timerx);

      for (int ix = 0, ii = 0; ix < nipx; ix++)
        {
          for (int j2d = 0; j2d < ndof2d; j2d++)
            for (int k2d = 0; k2d < ndof2d; k2d++)
              hmat2d(k2d,j2d) = allsum2d(j2d*ndof2d+k2d, ix);

          for (int j = 0; j < ndof; j++)
            {
              if (fel.IsGradient(j)) continue; 

              const Vec<3> & xf = xfactor(j, ix);
              int j2d = map3dto2d[j];
                  
              if (j2d < firstsplit1)
                {
                  double dx = xf(1) -xf(2);

                  for (int k = 0; k < firstsplit1; k++)
                    hmult3da(k,j) = dx * hmat2d(j2d,k)(1,1);
                  for (int k = firstsplit1; k < ndof2d; k++)
                    hmult3d(k,j) = dx * hmat2d(j2d,k).Col(1);

                  NgProfiler::AddFlops (timerx, firstsplit1);
                  NgProfiler::AddFlops (timerx, (ndof2d-firstsplit1) * 3);
                }
              else
                {
                  for (int k = 0; k < firstsplit1; k++)
                    hmult3da(k, j) = InnerProduct (hmat2d(j2d, k).Row(1), xf);
                  for (int k = firstsplit1; k < ndof2d; k++)
                    hmult3d(k, j) = hmat2d(j2d, k) * xf;

                  NgProfiler::AddFlops (timerx, firstsplit1 * 3);
                  NgProfiler::AddFlops (timerx, (ndof2d-firstsplit1) * 9);
                }
            }

          for (int j = 0; j < ndof; j++)
            {
              if (fel.IsGradient(j)) continue;

              const Vec<3> & xf = xfactor(j, ix);
              int j2d = map3dto2d[j];

              if (j2d < firstsplit1)
                {
                  double dx = xf(1)-xf(2);
                  for (int k = 0; k <= j; k++)
                    elmat(j,k) += dx * hmult3da(j2d,k);
                  NgProfiler::AddFlops (timerx, (j+1));
                }
              else
                {
                  for (int k = 0; k <= j; k++)
                    elmat(j,k) += InnerProduct (xf, hmult3d(j2d, k));
                  NgProfiler::AddFlops (timerx, (j+1) * 3);
                }

            }          
        }
      */


      /*
      for (int ix = 0, ii = 0; ix < nipx; ix++)
        {
          for (int j2d = 0; j2d < ndof2d; j2d++)
            for (int k2d = 0; k2d < ndof2d; k2d++)
              hmat2d(k2d,j2d) = allsum2d(j2d*ndof2d+k2d, ix);

          for (int j = 0; j < ndof; j++)
            {
              if (fel.IsGradient(j)) continue; 

              const Vec<3> & xf = xfactor(j, ix);
              int j2d = map3dto2d[j];
                  
              for (int k = 0; k < ndof2d; k++)
                hmult3d(k, j) = hmat2d(j2d, k) * xf;
            }

          for (int j = 0; j < ndof; j++)
            {
              if (fel.IsGradient(j)) continue;

              const Vec<3> & xf = xfactor(j, ix);
              int j2d = map3dto2d[j];

              for (int k = 0; k <= j; k++)
                elmat(j,k) += InnerProduct (xf, hmult3d(j2d, k));
            }          
        }
      */


      NgProfiler::StartTimer (timerx);

      FlatMatrix<Vec<3> > help(ndof, nipx, lh);
      FlatMatrix<> helpa(ndof, nipx, lh);

      FlatArray<int> map2dto3d(ndof2d * (order+1), lh);
      map2dto3d = -1;
      for (int j = 0; j < ndof; j++)
        {
          if (fel.IsGradient(j)) continue; 
          int start = (order+1) * map3dto2d[j];
          while (map2dto3d[start] != -1) start++;
          map2dto3d[start] = j;
        }



      NgProfiler::StartTimer (timerx1);

      for (int k2d = 0; k2d < firstsplit1; k2d++)
        {
          for (int jj = 0; jj < ndof; jj++)
            {
              if (fel.IsGradient(jj)) continue; 

              int j2d = map3dto2d[jj];

              if (j2d < firstsplit1)
                {
                  int ind = k2d*firstsplit1+j2d;

                  for (int i1 = 0; i1 < nipx; i1++)
                    helpa(jj, i1) = allsum2daa(ind, i1) * xfactora(jj, i1);

                  NgProfiler::AddFlops (timerx1, nipx);
                }
              else
                {
                  int ind = k2d+firstsplit1*j2d;

                  for (int i1 = 0; i1 < nipx; i1++)
                    helpa(jj, i1) = InnerProduct (allsum2da(ind, i1), hxfactor(jj, i1));
                  NgProfiler::AddFlops (timerx1, nipx * 3);
                }
            }
              
          for (int k1d = 0; k1d <= order; k1d++)
            {
              int kk = map2dto3d[k2d*(order+1)+k1d];
              if (kk == -1) continue;

              for (int j = 0; j <= kk; j++)
                {
                  if (fel.IsGradient(j)) continue; 

                  elmat(kk*ndof+j) = InnerProduct (helpa.Row(j), xfactora.Row(kk));

                  NgProfiler::AddFlops (timerx1, nipx);
                }
            }
        }              


      NgProfiler::StopTimer (timerx1);
      NgProfiler::StartTimer (timerx2);

      for (int k2d = firstsplit1; k2d < ndof2d; k2d++)
        {
          for (int jj = 0; jj < ndof; jj++)
            {
              if (fel.IsGradient(jj)) continue; 

              int j2d = map3dto2d[jj];

              if (j2d < firstsplit1)
                {
                  int ind = k2d*firstsplit1+j2d;

                  for (int i1 = 0; i1 < nipx; i1++)
                    help(jj, i1) = xfactora(jj, i1) * allsum2da(ind, i1);
                  NgProfiler::AddFlops (timerx2, nipx * 3);
                }
              else
                {
                  int ind = k2d*ndof2d+j2d;

                  for (int i1 = 0; i1 < nipx; i1++)
                    help(jj, i1) = allsum2d(ind, i1) * hxfactor(jj, i1);
                  NgProfiler::AddFlops (timerx2, nipx * 9);
                }
            }
              
          for (int k1d = 0; k1d <= order; k1d++)
            {
              int kk = map2dto3d[k2d*(order+1)+k1d];
              if (kk == -1) continue;

              for (int j = 0; j <= kk; j++)
                {
                  if (fel.IsGradient(j)) continue; 

                  elmat(kk*ndof+j) = InnerProduct (help.Row(j), hxfactor.Row(kk));
                  
                  NgProfiler::AddFlops (timerx2, nipx * 3);
                }
            }
        }              

      NgProfiler::StopTimer (timerx2);


      for (int j = 0; j < ndof; j++)
        for (int k = 0; k < j; k++)
          elmat(k,j) = elmat(j,k);
      NgProfiler::StopTimer (timerx);

        }
    }

#endif



  };






















  // ******************************** MassEdge **************************************



  
  template <int D>
  class MassEdgeIntegratorTP : public MassEdgeIntegrator<D>
  {

  public:
    MassEdgeIntegratorTP (CoefficientFunction * acoef)
      : MassEdgeIntegrator<D> (acoef)
    { ; }

    static Integrator * Create (ARRAY<CoefficientFunction*> & coeffs)
    {
      return new MassEdgeIntegratorTP (coeffs[0]);
    }

    virtual void
    AssembleElementMatrix (const FiniteElement & bfel,
			   const ElementTransformation & eltrans, 
			   FlatMatrix<double> & elmat,
			   LocalHeap & lh) const
    {

      const HCurlHighOrderTetTP * fel_tet = dynamic_cast<const HCurlHighOrderTetTP*> (&bfel);      
      const HCurlHighOrderPrismTP * fel_prism = dynamic_cast<const HCurlHighOrderPrismTP*> (&bfel);      
      // const HCurlHighOrderHexTP * fel_hex = dynamic_cast<const HCurlHighOrderHexTP*> (&bfel);      


      if (fel_tet)
        {
          TAssembleElementMatrix (*fel_tet, eltrans, elmat, lh);
          return;
        }

      if (fel_prism)
        if (fel_prism->Sort(0)+3 == fel_prism->Sort(3))  // same singular vertex at bottom and top
          {
            TAssembleElementMatrix (*fel_prism, eltrans, elmat, lh);
            return;
          }
      
      /*
        if (fel_hex)
        {
        TAssembleElementMatrix (*fel_hex, eltrans, elmat, lh);
        return;
        }
      */
      MassEdgeIntegrator<D>::AssembleElementMatrix(bfel, eltrans, elmat, lh);
    }




    template <class FEL>
    void TAssembleElementMatrix (const FEL & fel,
                                 const ElementTransformation & eltrans, 
                                 FlatMatrix<double> & elmat,
                                 LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("Fast assembling - Massedge");
      static int timeri = NgProfiler::CreateTimer ("Fast assembling - init");
      static int timerz = NgProfiler::CreateTimer ("Fast assembling z");
      static int timery = NgProfiler::CreateTimer ("Fast assembling y");
      static int timerx = NgProfiler::CreateTimer ("Fast assembling x");

      NgProfiler::RegionTimer reg (timer);

      int ndof = fel.GetNDof();
      
      elmat.AssignMemory (ndof, ndof, lh);
      elmat = 0.0;

      HeapReset hr (lh);

      NgProfiler::StartTimer (timeri);
      

      int order = fel.Order();
      int horder = max(order+1,2);

      IntegrationRuleTP<3> irtp(eltrans, 2*order, 1, lh);

          
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      const IntegrationRule & irz = irtp.GetIRZ();

      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nipz = irz.GetNIP();
      int nip = nipx * nipy * nipz;

      FlatArray<int> map2dto1d = fel.Get2dTo1dMapping ();
      FlatArray<int> map3dto2d = fel.Get3dTo2dMapping ();

      int ndof2d = fel.GetNDof2d();
      int ndof1d = fel.GetNDof1d();

      FlatMatrix<Vec<3> > xfactor(ndof, nipx, lh);
      FlatMatrix<Vec<3> > yfactor(ndof2d, nipy, lh);
      FlatMatrix<Vec<3> > zfactor(ndof1d, nipz, lh);

      for (int ix = 0; ix < nipx; ix++)
        {
          AutoDiff<1> x (irx[ix](0), 0);
          fel.CalcXFactor (x, xfactor.Col(ix), lh);
        }

      for (int iy = 0; iy < nipy; iy++)
        {
          AutoDiff<1> y (iry[iy](0), 0);
          fel.CalcYFactor (y, yfactor.Col(iy), lh);
        }

      for (int iz = 0; iz < nipz; iz++)
        {
          AutoDiff<1> z (irz[iz](0), 0);
          fel.CalcZFactor (z, zfactor.Col(iz), lh);
        }


      FlatMatrix<Mat<2> > hmat2d(ndof2d, ndof2d, lh);
      FlatMatrix<Mat<3> > hmat1d(ndof1d, ndof1d, lh);
      FlatMatrix<Mat<3,2> > hmult2d(ndof2d, ndof1d, lh);
      FlatMatrix<Vec<2> >   hmult3d(ndof2d, ndof, lh);
      FlatVector<Mat<3,3> > dbmat1d (ndof1d, lh);


      FlatVector<Mat<3> > tdmats(nip, lh);      // transformed D-matrix
      for (int i = 0; i < nip; i++)
        {
          IntegrationPoint ip (irtp.GetXi(i), irtp.GetWeight(i));
          SpecificIntegrationPoint<3,3> sip(ip, eltrans, irtp.GetPoint(i), irtp.GetJacobian(i), lh);
          
          Mat<3> dmat;
          (this->dmatop).GenerateMatrix (fel, sip, dmat, lh);
          dmat *= fabs (sip.GetJacobiDet()) * sip.IP().Weight(); 

          Mat<3> trans = irtp.GetDuffyJacobian (i) * sip.GetJacobianInverse();
          Mat<3> hmat = dmat * Trans (trans);
          tdmats[i] = trans * hmat; 
        }


      NgProfiler::StopTimer (timeri);

      for (int ix = 0, ii = 0; ix < nipx; ix++)
        {
          hmat2d = 0.0;
          
          for (int iy = 0; iy < nipy; iy++)
            {
              NgProfiler::StartTimer (timerz);

              hmat1d = 0.0;
              
              for (int iz = 0; iz < nipz; iz++, ii++)
                {
                  for (int k = 0; k < ndof1d; k++)
                    {
                      const Vec<3> & zf2 = zfactor(k, iz);
                      for (int l = 0; l < 3; l++)
                        {
                          dbmat1d(k)(l,0) = tdmats[ii](l,0) * zf2(0);   
                          dbmat1d(k)(l,1) = tdmats[ii](l,1) * zf2(0);   
                          dbmat1d(k)(l,2) = tdmats[ii](l,2) * zf2(1);
                        }
                    }

                  for (int j = 0; j < ndof1d; j++)
                    for (int k = 0; k <= j; k++)
                      {
                        const Mat<3,3> & hm = dbmat1d(k);
                        const Vec<3> & zf2 = zfactor(j, iz);

                        for (int l = 0; l < 3; l++)
                          {
                            hmat1d(k,j)(0,l) += zf2(0) * hm(0,l);
                            hmat1d(k,j)(1,l) += zf2(0) * hm(1,l);
                            hmat1d(k,j)(2,l) += zf2(1) * hm(2,l);
                          }
                      }
                }

              NgProfiler::AddFlops (timerz, nipz * ndof1d * 9);
              NgProfiler::AddFlops (timerz, nipz * ndof1d*(ndof1d+1) / 2 * 9);

              for (int j = 0; j < ndof1d; j++)
                for (int k = 0; k < j; k++)
                  hmat1d(j, k) = Trans(hmat1d(k,j));

              NgProfiler::StopTimer (timerz);
              NgProfiler::StartTimer (timery);
              
              for (int j1d = 0; j1d < ndof1d; j1d++)
                for (int k2d = 0; k2d < ndof2d; k2d++)
                  {
                    const Mat<3,3> & hm = hmat1d(map2dto1d[k2d], j1d);
                    const Vec<3> & yf2 = yfactor(k2d, iy);

                    for (int l = 0; l < 3; l++)
                      {
                        hmult2d(k2d, j1d)(l,0) =  hm(l,0) * yf2(0);
                        hmult2d(k2d, j1d)(l,1) =  yf2(1) * hm(l,1) + yf2(0) * hm(l,2);
                      }
                  }

              for (int j2d = 0; j2d < ndof2d; j2d++)
                {
                  const Vec<3> & yf2 = yfactor(j2d, iy);

                  for (int k2d = 0; k2d <= j2d; k2d++)
                    {
                      const Mat<3,2> & hm = hmult2d(k2d, map2dto1d[j2d]);
                      
                      for (int l = 0; l < 2; l++)
                        {
                          hmat2d(k2d, j2d)(0, l) += yf2(0) * hm(0, l);
                          hmat2d(k2d, j2d)(1, l) += yf2(1) * hm(1, l) + yf2(0) * hm(2, l);
                        }
                    }
                }

              NgProfiler::StopTimer (timery);

              NgProfiler::AddFlops (timery, ndof2d*ndof1d * 9);
              NgProfiler::AddFlops (timery, ndof2d*(ndof2d+1) * 6 / 2);
            }

          NgProfiler::StartTimer (timery);
          for (int j2d = 0; j2d < ndof2d; j2d++)
            for (int k2d = 0; k2d < j2d; k2d++)
              hmat2d(j2d, k2d) = Trans(hmat2d(k2d,j2d));

          NgProfiler::StopTimer (timery);


          NgProfiler::RegionTimer regx (timerx);

          for (int j = 0; j < ndof; j++)
            {
              const Vec<3> & hxf = xfactor(j, ix);
              int j2d = map3dto2d[j];

              for (int k = 0; k < ndof2d; k++)
                for (int l = 0; l < 2; l++)
                  hmult3d(k, j)(l) = hmat2d(j2d, k)(l,0) * hxf(1) + hmat2d(j2d,k)(l,1) * hxf(0);

              NgProfiler::AddFlops (timerx, ndof2d * 4);
            }

          for (int j = 0; j < ndof; j++)
            {
              const Vec<3> & hxf = xfactor(j, ix);
              int j2d = map3dto2d[j];

              for (int k = 0; k <= j; k++)
                elmat(j,k) += hxf(1) * hmult3d(j2d, k)(0) + hxf(0) * hmult3d(j2d,k)(1);

              NgProfiler::AddFlops (timerx, (j+1) * 2);
            }
        }


      for (int j = 0; j < ndof; j++)
        for (int k = 0; k < j; k++)
          elmat(k,j) = elmat(j,k);
    }
  };

















  
  
  namespace
  {
    class Init
    { 
    public: 
      Init ();
    };        
    
    Init::Init()
    {
      GetIntegrators().AddLFIntegrator ("sourceedgetp", 3, 3,
                                        SourceEdgeIntegratorTP<3>::Create);
      GetIntegrators().AddLFIntegrator ("curledgetp", 3, 3,
                                        CurlEdgeIntegratorTP<3>::Create);
      GetIntegrators().AddBFIntegrator ("curlcurledgetp", 3, 1,
                                        CurlCurlIntegratorTP<3>::Create);
      GetIntegrators().AddBFIntegrator ("massedgetp", 3, 1,
                                        MassEdgeIntegratorTP<3>::Create);

    }
    
    Init init;
  }

}
