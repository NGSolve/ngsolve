/*********************************************************************/
/* File:   tpequations.cpp                                           */
/* Author: Joachim Schöberl                                          */
/* Date:   1. June 2007                                              */
/*********************************************************************/


/*

High order integrators in tensor product form

*/



#include <fem.hpp>
namespace ngfem
{
  using namespace ngfem;



  // ******************************** source integrator **************************************

  

  template <int D>
  class SourceIntegratorTP : public SourceIntegrator<D>
  {
  public:
    SourceIntegratorTP (CoefficientFunction * acoef)
      : SourceIntegrator<D> (acoef)
    { ; }

    static Integrator * Create (ARRAY<CoefficientFunction*> & coeffs)
    {
      return new SourceIntegratorTP (coeffs[0]);
    }

    virtual void
    AssembleElementVector (const FiniteElement & bfel,
			   const ElementTransformation & eltrans, 
			   FlatVector<double> & elvec,
			   LocalHeap & lh) const
    {
      const NodalFiniteElement<D> & fel = 
        dynamic_cast<const NodalFiniteElement<D>&> (bfel);

      int ndof = fel.GetNDof();
      
      elvec.AssignMemory (ndof, lh);

      HeapReset hr (lh);

      IntegrationRuleTP<D> irtp(eltrans, 2*fel.Order(), 1, lh);
      int nip = irtp.GetNIP(); 
      
      FlatVector<double> dvecs(nip, lh);
      
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<D,D> sip(ip, eltrans, irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Vec<1> dvec;
          (this->dvecop).GenerateVector (fel, sip, dvec, lh);
          dvecs[ii] = dvec(0) * fabs (sip.GetJacobiDet()) * ip.Weight();
        }

      fel.EvaluateShapeGridTrans (irtp, dvecs, elvec, lh);
    }
  };




  template <int D>
  class NeumannIntegratorTP : public NeumannIntegrator<D>
  {
  public:
    NeumannIntegratorTP (CoefficientFunction * acoef)
      : NeumannIntegrator<D> (acoef)
    { ; }

    static Integrator * Create (ARRAY<CoefficientFunction*> & coeffs)
    {
      return new NeumannIntegratorTP (coeffs[0]);
    }

    virtual void
    AssembleElementVector (const FiniteElement & bfel,
			   const ElementTransformation & eltrans, 
			   FlatVector<double> & elvec,
			   LocalHeap & lh) const
    {
      const NodalFiniteElement<D-1> & fel = 
        dynamic_cast<const NodalFiniteElement<D-1>&> (bfel);

      int ndof = fel.GetNDof();
      
      elvec.AssignMemory (ndof, lh);

      HeapReset hr (lh);

      IntegrationRuleTP<D-1> irtp(eltrans, 2*fel.Order(), 0, lh);
      int nip = irtp.GetNIP(); 
      
      FlatVector<double> dvecs(nip, lh);
      
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<D-1,D> sip(ip, eltrans, lh); // irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Vec<1> dvec;
          (this->dvecop).GenerateVector (fel, sip, dvec, lh);
          dvecs[ii] = dvec(0) * fabs (sip.GetJacobiDet()) * ip.Weight();
        }

      fel.EvaluateShapeGridTrans (irtp, dvecs, elvec, lh);
    }
  };













  // ******************************** Laplace **************************************





  
  template <int D>
  class LaplaceIntegratorTP : public LaplaceIntegrator<D>
  { };







  /* *********************** Laplace 2D *************************************** */





  class TrafoGradient2dY
  {
    const AutoDiff<1> & ady;
  public:

    TrafoGradient2dY (const AutoDiff<1> & hady) : ady(hady) { ; }

    template <typename TIN, typename TOUT>
    void Transform (const TIN & in, TOUT & out)
    {
      out(1) = ady.Value() * in(0);
      out(0) = ady.DValue(0) * in(1);
    }

    template <typename TIN, typename TOUT>
    void TransformAdd (const TIN & in, TOUT & out)
    {
      out(1) += ady.Value() * in(0);
      out(0) += ady.DValue(0) * in(1);
    }

    template <typename TIN, typename TOUT>
    void TransformAddT (const TIN & in, TOUT & out)
    {
      out(0) += ady.Value() * in(1);
      out(1) += ady.DValue(0) * in(0);
    }
  };


  class TrafoGradient2dX
  {
    const AutoDiff<1> & adx;
  public:
    
    TrafoGradient2dX (const AutoDiff<1> & hx) : adx(hx) { ; }
    
    template <typename TIN, typename TOUT>
    void Transform (const TIN & in, TOUT & out)
    {
      out = adx.Value() * in(0) + adx.DValue(0) * in(1);
    }
    
    template <typename TIN, typename TOUT>
    void TransformAdd (const TIN & in, TOUT & out)
    {
      out += adx.Value() * in(0) + adx.DValue(0) * in(1);
    }

    template <typename TIN, typename TOUT>
    void TransformAddT (const TIN & in, TOUT & out)
    {
      out(0) += adx.Value() * in;
      out(1) += adx.DValue(0) * in;
    }

  };





  
  template <>
  class LaplaceIntegratorTP<2> : public LaplaceIntegrator<2>
  {
    enum { D = 2 };
  public:
    LaplaceIntegratorTP (CoefficientFunction * acoef)
      : LaplaceIntegrator<D> (acoef)
    { ; }

    static Integrator * Create (ARRAY<CoefficientFunction*> & coeffs)
    {
      return new LaplaceIntegratorTP (coeffs[0]);
    }

    virtual void
    AssembleElementMatrix (const FiniteElement & bfel,
			   const ElementTransformation & eltrans, 
			   FlatMatrix<double> & elmat,
			   LocalHeap & lh) const
    {

      const H1HighOrderTP<ET_TRIG> * fel_trig = dynamic_cast<const H1HighOrderTP<ET_TRIG>*> (&bfel);      

      if (fel_trig)
        {
          TAssembleElementMatrix (*fel_trig, eltrans, elmat, lh);
          return;
        }

      const H1HighOrderTP<ET_QUAD> * fel_quad = dynamic_cast<const H1HighOrderTP<ET_QUAD>*> (&bfel);      

      if (fel_quad)
        {
          TAssembleElementMatrix (*fel_quad, eltrans, elmat, lh);
          return;
        }



      LaplaceIntegrator<D>::AssembleElementMatrix(bfel, eltrans, elmat, lh);
    }




    virtual void
    ApplyElementMatrix (const FiniteElement & bfel,
                        const ElementTransformation & eltrans, 
                        const FlatVector<double> & elx, 
                        FlatVector<double> & ely,
                        void * precomputed,
                        LocalHeap & lh) const
    {
      const H1HighOrderTP<ET_TRIG> * fel_trig = dynamic_cast<const H1HighOrderTP<ET_TRIG>*> (&bfel);      


      if (fel_trig)
        {
          TApplyElementMatrix (*fel_trig, eltrans, elx, ely, precomputed, lh);
          return;
        }

      LaplaceIntegrator<D>::ApplyElementMatrix(bfel, eltrans, elx, ely, precomputed, lh);
    }




    template <class FEL>
    void TAssembleElementMatrix (const FEL & fel,
                                 const ElementTransformation & eltrans, 
                                 FlatMatrix<double> & elmat,
                                 LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("Fast assembling - Laplace");
      static int timerg = NgProfiler::CreateTimer ("Fast assembling - geom");
      static int timeri = NgProfiler::CreateTimer ("Fast assembling - init");
      static int timerc = NgProfiler::CreateTimer ("Fast assembling - compute");

      NgProfiler::RegionTimer reg (timer);
      NgProfiler::StartTimer (timerg);      


      int ndof = fel.GetNDof();
      elmat.AssignMemory (ndof, ndof, lh);

      HeapReset hr (lh);

      int intorder = 2 * fel.Order();
      if (fel.ElementType() == ET_TRIG) intorder -= 2;

      IntegrationRuleTP<D> irtp(eltrans, intorder, 1, lh);
      
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      
      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nip = nipx * nipy;
      
      FlatVector<Mat<2> > tdmat(nip, lh);      // transformed D-matrix
      
      // geometry calculation
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<2,2> sip(ip, eltrans, irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Mat<2> dmat, hmat, trans;
          this->dmatop.GenerateMatrix (fel, sip, dmat, lh);
          dmat *= fabs (sip.GetJacobiDet()) * ip.Weight();
          
          trans = irtp.GetDuffyJacobian(ii) * sip.GetJacobianInverse();
          hmat = trans * dmat;
          tdmat[ii] = hmat * Trans (trans);
        }

      NgProfiler::StopTimer (timerg);      

      NgProfiler::StartTimer (timeri);

      int ndof1d = fel.GetNDof1d();

      FlatArray<int> pmap1dto2d = fel.Get1dTo2dMappingPointer ();
      FlatArray<int> map1dto2d  = fel.Get1dTo2dMapping ();
      FlatArray<int> map2dto1d  = fel.Get2dTo1dMapping ();

      
      FlatMatrix<AutoDiff<1> > fac_xdx(ndof,   nipx, lh);    
      FlatMatrix<AutoDiff<1> > fac_ydy(ndof1d, nipy, lh);    
      
      for (int i1 = 0; i1 < nipy; i1++)
        {
          AutoDiff<1> y (iry[i1](0), 0);
          fel.CalcYFactor (y, fac_ydy.Col(i1), lh);
        }
      
      for (int i1 = 0; i1 < nipx; i1++)
        {
          AutoDiff<1> x (irx[i1](0), 0);
          fel.CalcXFactor (x, fac_xdx.Col(i1), lh);
        }

      FlatVector<Mat<2> > allsum1d(nipx, lh);
      FlatVector<Vec<2> > mult2d(nipx, lh);
      FlatVector<Mat<2> > mult1d(nip, lh);

      NgProfiler::StopTimer (timeri);      
      NgProfiler::StartTimer (timerc);

      for (int ky = 0; ky < ndof1d; ky++)
        {
          for (int ix = 0, ii = 0; ix < nipx; ix++)
            for (int iy = 0; iy < nipy; iy++, ii++)
              for (int l = 0; l < 2; l++)
                TrafoGradient2dY(fac_ydy(ky,iy)) . Transform(tdmat[ii].Row(l), mult1d(ii).Row(l));
          
          for (int jy = 0; jy <= ky; jy++)
            {
              for (int ix = 0, ii = 0; ix < nipx; ix++)
                {
                  Mat<2> sum = 0.0;

                  for (int iy = 0; iy < nipy; iy++, ii++)
                    for (int l = 0; l < 2; l++)
                      TrafoGradient2dY(fac_ydy(jy,iy)) . TransformAdd(mult1d(ii).Col(l), sum.Col(l));
                  
                  allsum1d(ix) = sum;
                }

              for (int kx = pmap1dto2d[ky]; kx < pmap1dto2d[ky+1]; kx++)
                {
                  int kk = map1dto2d[kx];

                  for (int ix = 0; ix < nipx; ix++)
                    for (int l = 0; l < 2; l++)
                      TrafoGradient2dX (fac_xdx(kk,ix)) . Transform(allsum1d(ix).Row(l), mult2d(ix)(l));
                  
                  for (int jx = pmap1dto2d[jy]; jx < pmap1dto2d[jy+1]; jx++)
                    {
                      int jj = map1dto2d[jx];

                      double sum = 0;
                      for (int ix = 0; ix < nipx; ix++)
                        TrafoGradient2dX (fac_xdx(jj,ix)) . TransformAdd(mult2d(ix), sum);
                      
                      elmat(kk,jj) = sum;
                      elmat(jj,kk) = sum;
                    }
                }
            }
        }

      NgProfiler::StopTimer (timerc);

      NgProfiler::AddFlops (timerc, 4*nip*ndof1d);
      NgProfiler::AddFlops (timerc, 4*nip*ndof1d*(ndof1d+1)/2);

      NgProfiler::AddFlops (timerc, 4*ndof1d*ndof*nipx);
      NgProfiler::AddFlops (timerc, ndof*(ndof+1)*nipx);

    }








    template <class FEL>
    void TApplyElementMatrix (const FEL & fel,
                              const ElementTransformation & eltrans, 
                              const FlatVector<double> & elx, 
                              FlatVector<double> & ely,
                              void * precomputed,
                              LocalHeap & lh) const
    {
#ifdef NONE
      static int timer = NgProfiler::CreateTimer ("Fast apply");


      NgProfiler::RegionTimer reg (timer);

      int ndof = fel.GetNDof();

      HeapReset hr (lh);

      
      int order = fel.Order();

      FlatArray<int> map2dto1d = fel.Get2dTo1dMapping ();
      FlatArray<int> map3dto2d = fel.Get3dTo2dMapping ();

      int ndof2d = fel.GetNDof2d();
      int ndof1d = fel.GetNDof1d();
          
      int intorder = 2 * order;
      if (fel.ElementType() == ET_TET) intorder -= 2;
      
      IntegrationRuleTP<D> irtp(eltrans, intorder, lh);
      
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      const IntegrationRule & irz = irtp.GetIRZ();
      
      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nipz = irz.GetNIP();
      int nip = nipx * nipy * nipz;
      
      FlatVector<Mat<3> > tdmat(nip, lh);      // transformed D-matrix
      
      // geometry calculation
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<3,3> sip(ip, eltrans, irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Mat<3> dmat, hmat, trans;
          this->dmatop.GenerateMatrix (fel, sip, dmat, lh);
          dmat *= fabs (sip.GetJacobiDet()) * ip.Weight();
          
          trans = irtp.GetDuffyJacobian(ii) * sip.GetJacobianInverse();
          hmat = trans * dmat;
          tdmat[ii] = hmat * Trans (trans);
        }

      
      FlatMatrix<AutoDiff<1> > fac_xdx(ndof, nipx, lh);    // ndof   x nipx
      FlatMatrix<AutoDiff<1> > fac_ydy(ndof2d, nipy, lh);    // ndof2d x nipy
      FlatMatrix<AutoDiff<1> > fac_zdz(ndof1d, nipz, lh);    // ndof1d x nipz

      
      for (int i1 = 0; i1 < nipz; i1++)
        {
          AutoDiff<1> z (irz[i1](0), 0);
          fel.CalcZFactor (z, fac_zdz.Col(i1), lh);
        }
      
      for (int i1 = 0; i1 < nipy; i1++)
        {
          AutoDiff<1> y (iry[i1](0), 0);
          fel.CalcYFactor (y, fac_ydy.Col(i1), lh);
        }
      
      for (int i1 = 0; i1 < nipx; i1++)
        {
          AutoDiff<1> x (irx[i1](0), 0);
          fel.CalcXFactor (x, fac_xdx.Col(i1), lh);
        }


      static int timerc = NgProfiler::CreateTimer ("Fast apply - compute");
      NgProfiler::StartTimer (timerc);
      
      int nipxy = nipx *nipy;
      
      FlatMatrix<Mat<2> > allsum2d(nipx, lh);
      FlatVector<Mat<3> > allsum1d(nipxy, lh);

      FlatVector<Vec<2> > mult3d(nipx, lh);
      FlatVector<Mat<3,2> > mult2d(nipxy, lh);
      FlatVector<Mat<3> > mult1d(nip, lh);

      FlatMatrix<Vec<2> > hmat1(ndof2d, nipx, lh);
      FlatMatrix<Vec<3> > hmat2(ndof1d, nipx*nipy, lh);
      FlatVector<Vec<3> > hmat3(nip, lh);

      /*
        cout << "ndof = " << ndof << endl;
        int flops = 2 * (ndof*nipx*2 +  ndof2d*nipx*nipy*3 + ndof1d*nip*3 ) +  nip*9;
        cout << "flops = " << flops << ", flops/n = " << flops/ndof << endl;
      */

      // for (int cnt = 0; cnt < ndof; cnt++)
      {

        hmat1 = 0.0;
        for (int j = 0; j < ndof; j++)
          {
            int j2d = map3dto2d[j];
            for (int ix = 0; ix < nipx; ix++)
              TrafoGradientX(fac_xdx(j, ix)) . TransformAddT (elx(j), hmat1(j2d, ix));
          }
        NgProfiler::AddFlops (timerc, ndof*nipx*2);

        hmat2 = 0.0;
        for (int j = 0; j < ndof2d; j++)
          {
            int j1d = map2dto1d[j];
            for (int ix = 0; ix < nipx; ix++)
              for (int iy = 0; iy < nipy; iy++)
                TrafoGradientY(fac_ydy(j, iy)) . TransformAddT (hmat1(j, ix), hmat2(j1d, ix*nipy+iy));
          }
        NgProfiler::AddFlops (timerc, ndof2d*nipx*nipy*3);

        hmat3 = 0.0;
        for (int j = 0; j < ndof1d; j++)
          for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
            for (int iz = 0; iz < nipz; iz++, ii++)
              TrafoGradientZ(fac_zdz(j, iz)) . TransformAddT (hmat2(j, ixy), hmat3(ii));


        NgProfiler::AddFlops (timerc, ndof1d*nip*3);

        for (int ii = 0; ii < nip; ii++)
          {
            Vec<3> hv = hmat3(ii);
            hmat3(ii) = tdmat[ii] * hv;
          }
        NgProfiler::AddFlops (timerc, nip*9);

        hmat2 = 0.0;
        for (int j = 0; j < ndof1d; j++)
          for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
            for (int iz = 0; iz < nipz; iz++, ii++)
              TrafoGradientZ(fac_zdz(j, iz)) . TransformAdd (hmat3(ii), hmat2(j, ixy));

        NgProfiler::AddFlops (timerc, ndof1d*nip*3);

        hmat1 = 0.0;
        for (int j = 0; j < ndof2d; j++)
          {
            int j1d = map2dto1d[j];
            for (int ix = 0; ix < nipx; ix++)
              for (int iy = 0; iy < nipy; iy++)
                TrafoGradientY(fac_ydy(j, iy)) . TransformAdd (hmat2(j1d, ix*nipy+iy), hmat1(j, ix));
          }
        NgProfiler::AddFlops (timerc, ndof2d*nipx*nipy*3);

        ely = 0.0;
        for (int j = 0; j < ndof; j++)
          {
            int j2d = map3dto2d[j];
            for (int ix = 0; ix < nipx; ix++)
              TrafoGradientX(fac_xdx(j, ix)) . TransformAdd (hmat1(j2d, ix), ely(j));
          }
        NgProfiler::AddFlops (timerc, ndof*nipx*2);
      }
      

      NgProfiler::StopTimer (timerc);
#endif      
    }

  };






  /* *********************** Laplace 3D *************************************** */





  class TrafoGradientZ
  {
    const AutoDiff<1> & adz;
  public:

    TrafoGradientZ (const AutoDiff<1> & hadz) : adz(hadz) { ; }

    template <typename TIN, typename TOUT>
    void Transform (const TIN & in, TOUT & out)
    {
      out(0) = adz.Value() * in(0);
      out(1) = adz.Value() * in(1);
      out(2) = adz.DValue(0) * in(2);
    }

    template <typename TIN, typename TOUT>
    void TransformAdd (const TIN & in, TOUT & out)
    {
      out(0) += adz.Value() * in(0);
      out(1) += adz.Value() * in(1);
      out(2) += adz.DValue(0) * in(2);
    }


    template <typename TIN, typename TOUT>
    void TransformAddT (const TIN & in, TOUT & out)
    {
      out(0) += adz.Value() * in(0);
      out(1) += adz.Value() * in(1);
      out(2) += adz.DValue(0) * in(2);
    }

  };


  class TrafoGradientY
  {
    const AutoDiff<1> & ady;
  public:
    
    TrafoGradientY (const AutoDiff<1> & hady) : ady(hady) { ; }
    
    template <typename TIN, typename TOUT>
    void Transform (const TIN & in, TOUT & out)
    {
      out(0) = ady.Value() * in(2) + ady.DValue(0) * in(1);
      out(1) = ady.Value() * in(0);
    }
    
    template <typename TIN, typename TOUT>
    void TransformAdd (const TIN & in, TOUT & out)
    {
      out(0) += ady.Value() * in(2) + ady.DValue(0) * in(1);
      out(1) += ady.Value() * in(0);
    }

    template <typename TIN, typename TOUT>
    void TransformAddT (const TIN & in, TOUT & out)
    {
      out(0) += ady.Value() * in(1);
      out(1) += ady.DValue(0) * in(0);
      out(2) += ady.Value() * in(0);
    }

  };


  class TrafoGradientX
  {
    const AutoDiff<1> & adx;
  public:
    
    TrafoGradientX (const AutoDiff<1> & hx) : adx(hx) { ; }
    
    template <typename TIN, typename TOUT>
    void Transform (const TIN & in, TOUT & out)
    {
      out = adx.Value() * in(0) + adx.DValue(0) * in(1);
    }
    
    template <typename TIN, typename TOUT>
    void TransformAdd (const TIN & in, TOUT & out)
    {
      out += adx.Value() * in(0) + adx.DValue(0) * in(1);
    }

    template <typename TIN, typename TOUT>
    void TransformAddT (const TIN & in, TOUT & out)
    {
      out(0) += adx.Value() * in;
      out(1) += adx.DValue(0) * in;
    }

  };





  
  template <>
  class LaplaceIntegratorTP<3> : public LaplaceIntegrator<3>
  {
    enum { D = 3 };
  public:
    LaplaceIntegratorTP (CoefficientFunction * acoef)
      : LaplaceIntegrator<D> (acoef)
    { ; }

    static Integrator * Create (ARRAY<CoefficientFunction*> & coeffs)
    {
      return new LaplaceIntegratorTP (coeffs[0]);
    }

    virtual void
    AssembleElementMatrix (const FiniteElement & bfel,
			   const ElementTransformation & eltrans, 
			   FlatMatrix<double> & elmat,
			   LocalHeap & lh) const
    {

      if (bfel.IsTPElement())
        {
          const H1HighOrderTP<ET_TET> * fel_tet = dynamic_cast<const H1HighOrderTP<ET_TET>*> (&bfel);      
          const H1HighOrderTP<ET_PRISM> * fel_prism = dynamic_cast<const H1HighOrderTP<ET_PRISM>*> (&bfel);      
          const H1HighOrderTP<ET_HEX> * fel_hex = dynamic_cast<const H1HighOrderTP<ET_HEX>*> (&bfel);      
          
          if (fel_tet)
            TAssembleElementMatrix (*fel_tet, eltrans, elmat, lh);
          
          if (fel_prism)
            TAssembleElementMatrix (*fel_prism, eltrans, elmat, lh);
          
          if (fel_hex)
            TAssembleElementMatrix (*fel_hex, eltrans, elmat, lh);

          return;
        }

      LaplaceIntegrator<D>::AssembleElementMatrix(bfel, eltrans, elmat, lh);
    }




    virtual void
    ApplyElementMatrix (const FiniteElement & bfel,
                        const ElementTransformation & eltrans, 
                        const FlatVector<double> & elx, 
                        FlatVector<double> & ely,
                        void * precomputed,
                        LocalHeap & lh) const
    {
      if (bfel.IsTPElement())
        {
          const H1HighOrderTP<ET_TET> * fel_tet = dynamic_cast<const H1HighOrderTP<ET_TET>*> (&bfel);      
          const H1HighOrderTP<ET_PRISM> * fel_prism = dynamic_cast<const H1HighOrderTP<ET_PRISM>*> (&bfel);      
          const H1HighOrderTP<ET_HEX> * fel_hex = dynamic_cast<const H1HighOrderTP<ET_HEX>*> (&bfel);      
          
          
          if (fel_tet)
            TApplyElementMatrix (*fel_tet, eltrans, elx, ely, precomputed, lh);
          
          if (fel_prism)
            TApplyElementMatrix (*fel_prism, eltrans, elx, ely, precomputed, lh);
      
          if (fel_hex)
            TApplyElementMatrix (*fel_hex, eltrans, elx, ely, precomputed, lh);

          return;
        }

      LaplaceIntegrator<D>::ApplyElementMatrix(bfel, eltrans, elx, ely, precomputed, lh);
    }






    template <class FEL>
    void TAssembleElementMatrix (const FEL & fel,
                                 const ElementTransformation & eltrans, 
                                 FlatMatrix<double> & elmat,
                                 LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("Fast assembling - Laplace");
      static int timerg = NgProfiler::CreateTimer ("Fast assembling - geom");
      static int timeri = NgProfiler::CreateTimer ("Fast assembling - init");
      static int timerc = NgProfiler::CreateTimer ("Fast assembling - compute");

      NgProfiler::RegionTimer reg (timer);
      NgProfiler::StartTimer (timerg);      


      int ndof = fel.GetNDof();
      elmat.AssignMemory (ndof, ndof, lh);

      HeapReset hr (lh);

      int intorder = 2 * fel.Order();
      if (fel.ElementType() == ET_TET) intorder -= 2;

      IntegrationRuleTP<D> irtp(eltrans, intorder, 1, lh);
      
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      const IntegrationRule & irz = irtp.GetIRZ();
      
      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nipz = irz.GetNIP();
      int nip = nipx * nipy * nipz;
      
      FlatVector<Mat<3> > tdmat(nip, lh);      // transformed D-matrix
      
      // geometry calculation
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<3,3> sip(ip, eltrans, irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Mat<3> dmat, hmat, trans;
          this->dmatop.GenerateMatrix (fel, sip, dmat, lh);
          dmat *= fabs (sip.GetJacobiDet()) * ip.Weight();
          
          trans = irtp.GetDuffyJacobian(ii) * sip.GetJacobianInverse();
          hmat = trans * dmat;
          tdmat[ii] = hmat * Trans (trans);
        }

      NgProfiler::StopTimer (timerg);      

      NgProfiler::StartTimer (timeri);

      int ndof2d = fel.GetNDof2d();
      int ndof1d = fel.GetNDof1d();

      FlatArray<int> pmap1dto2d = fel.Get1dTo2dMappingPointer ();
      FlatArray<int> map1dto2d  = fel.Get1dTo2dMapping ();

      FlatArray<int> pmap2dto3d = fel.Get2dTo3dMappingPointer ();
      FlatArray<int> map2dto3d  = fel.Get2dTo3dMapping ();


      
      FlatMatrix<AutoDiff<1> > fac_xdx(ndof,   nipx, lh);    
      FlatMatrix<AutoDiff<1> > fac_ydy(ndof2d, nipy, lh);    
      FlatMatrix<AutoDiff<1> > fac_zdz(ndof1d, nipz, lh);    
      
      for (int i1 = 0; i1 < nipz; i1++)
        {
          AutoDiff<1> z (irz[i1](0), 0);
          fel.CalcZFactor (z, fac_zdz.Col(i1), lh);
        }
      
      for (int i1 = 0; i1 < nipy; i1++)
        {
          AutoDiff<1> y (iry[i1](0), 0);
          fel.CalcYFactor (y, fac_ydy.Col(i1), lh);
        }
      
      for (int i1 = 0; i1 < nipx; i1++)
        {
          AutoDiff<1> x (irx[i1](0), 0);
          fel.CalcXFactor (x, fac_xdx.Col(i1), lh);
        }


      int nipxy = nipx *nipy;
      FlatVector<Mat<2> > allsum2d(nipx, lh);
      FlatVector<Mat<3> > allsum1d(nipxy, lh);
      FlatVector<Vec<2> > mult3d(nipx, lh);
      FlatVector<Mat<2,3> > mult2d(nipxy, lh);
      FlatVector<Mat<3> > mult1d(nip, lh);

      NgProfiler::StopTimer (timeri);      



      NgProfiler::StartTimer (timerc);

      for (int kz = 0; kz < ndof1d; kz++)
        {
          for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
            for (int iz = 0; iz < nipz; iz++, ii++)
              for (int l = 0; l < 3; l++)
                TrafoGradientZ(fac_zdz(kz,iz)) . Transform(tdmat[ii].Col(l), mult1d(ii).Col(l));
          
          for (int jz = 0; jz <= kz; jz++)
            {
              for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
                {
                  Mat<3> sum = 0.0;

                  for (int iz = 0; iz < nipz; iz++, ii++)
                    for (int l = 0; l < 3; l++)
                      TrafoGradientZ(fac_zdz(jz,iz)) . TransformAdd(mult1d(ii).Row(l), sum.Row(l));
                  
                  allsum1d(ixy) = sum;
                }
              
              for (int ky = pmap1dto2d[kz]; ky < pmap1dto2d[kz+1]; ky++)
                {
                  int kyz = map1dto2d[ky];
                  for (int ix = 0, ixy=0; ix < nipx; ix++)
                    for (int iy = 0; iy < nipy; iy++, ixy++)
                      for (int l = 0; l < 3; l++)
                        TrafoGradientY (fac_ydy(kyz,iy)) . Transform(allsum1d(ixy).Col(l), mult2d(ixy).Col(l));


                  for (int jy = pmap1dto2d[jz]; jy < pmap1dto2d[jz+1]; jy++)
                    {
                      if (jy > ky) continue;
                      int jyz = map1dto2d[jy];

                      for (int ix = 0; ix < nipx; ix++)
                        {
                          Mat<2> hsum2d = 0.0;
                        
                          for (int iy = 0; iy < nipy; iy++)
                            for (int l = 0; l < 2; l++)
                              TrafoGradientY(fac_ydy(jyz,iy)) . TransformAdd(mult2d(ix*nipy+iy).Row(l), hsum2d.Row(l));
                        
                          allsum2d(ix) = Trans(hsum2d);
                        }
                    
                      for (int kx = pmap2dto3d[kyz]; kx < pmap2dto3d[kyz+1]; kx++)
                        {
                          int kk = map2dto3d[kx];

                          for (int ix = 0; ix < nipx; ix++)
                            for (int l = 0; l < 2; l++)
                              TrafoGradientX (fac_xdx(kk,ix)) . Transform(allsum2d(ix).Row(l), mult3d(ix)(l));

                          for (int jx = pmap2dto3d[jyz]; jx < pmap2dto3d[jyz+1]; jx++)
                            {
                              int jj = map2dto3d[jx];

                              double sum = 0;
                              for (int ix = 0; ix < nipx; ix++)
                                TrafoGradientX (fac_xdx(jj,ix)) . TransformAdd(mult3d(ix), sum);
                            
                              elmat(kk,jj) = sum;
                              elmat(jj,kk) = sum;
                            }
                        }
                    }
                }
            }
        }

      NgProfiler::StopTimer (timerc);


      NgProfiler::AddFlops (timerc, 9*nip*ndof1d);
      NgProfiler::AddFlops (timerc, 9*nip*ndof1d*(ndof1d+1)/2);

      NgProfiler::AddFlops (timerc, 9*nipxy*ndof1d*ndof2d);
      NgProfiler::AddFlops (timerc, 3*nipxy*ndof2d*(ndof2d+1));

      NgProfiler::AddFlops (timerc, 4*ndof2d*ndof*nipx);
      NgProfiler::AddFlops (timerc, ndof*(ndof+1)*nipx);
    }








    template <class FEL>
    void TApplyElementMatrix (const FEL & fel,
                              const ElementTransformation & eltrans, 
                              const FlatVector<double> & elx, 
                              FlatVector<double> & ely,
                              void * precomputed,
                              LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("Fast apply");


      NgProfiler::RegionTimer reg (timer);

      int ndof = fel.GetNDof();

      HeapReset hr (lh);

      
      int order = fel.Order();

      FlatArray<int> map2dto1d = fel.Get2dTo1dMapping ();
      FlatArray<int> map3dto2d = fel.Get3dTo2dMapping ();

      int ndof2d = fel.GetNDof2d();
      int ndof1d = fel.GetNDof1d();
          
      int intorder = 2 * order;
      if (fel.ElementType() == ET_TET) intorder -= 2;
      
      IntegrationRuleTP<D> irtp(eltrans, intorder, 1, lh);
      
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      const IntegrationRule & irz = irtp.GetIRZ();
      
      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nipz = irz.GetNIP();
      int nip = nipx * nipy * nipz;
      
      FlatVector<Mat<3> > tdmat(nip, lh);      // transformed D-matrix
      
      // geometry calculation
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<3,3> sip(ip, eltrans, irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Mat<3> dmat, hmat, trans;
          this->dmatop.GenerateMatrix (fel, sip, dmat, lh);
          dmat *= fabs (sip.GetJacobiDet()) * ip.Weight();
          
          trans = irtp.GetDuffyJacobian(ii) * sip.GetJacobianInverse();
          hmat = trans * dmat;
          tdmat[ii] = hmat * Trans (trans);
        }

      
      FlatMatrix<AutoDiff<1> > fac_xdx(ndof, nipx, lh);    // ndof   x nipx
      FlatMatrix<AutoDiff<1> > fac_ydy(ndof2d, nipy, lh);    // ndof2d x nipy
      FlatMatrix<AutoDiff<1> > fac_zdz(ndof1d, nipz, lh);    // ndof1d x nipz

      
      for (int i1 = 0; i1 < nipz; i1++)
        {
          AutoDiff<1> z (irz[i1](0), 0);
          fel.CalcZFactor (z, fac_zdz.Col(i1), lh);
        }
      
      for (int i1 = 0; i1 < nipy; i1++)
        {
          AutoDiff<1> y (iry[i1](0), 0);
          fel.CalcYFactor (y, fac_ydy.Col(i1), lh);
        }
      
      for (int i1 = 0; i1 < nipx; i1++)
        {
          AutoDiff<1> x (irx[i1](0), 0);
          fel.CalcXFactor (x, fac_xdx.Col(i1), lh);
        }


      static int timerc = NgProfiler::CreateTimer ("Fast apply - compute");
      NgProfiler::StartTimer (timerc);
      
      int nipxy = nipx *nipy;

      FlatMatrix<Vec<2> > hmat1(ndof2d, nipx, lh);
      FlatMatrix<Vec<3> > hmat2(ndof1d, nipx*nipy, lh);
      FlatVector<Vec<3> > hmat3(nip, lh);

      /*
        cout << "ndof = " << ndof << endl;
        int flops = 2 * (ndof*nipx*2 +  ndof2d*nipx*nipy*3 + ndof1d*nip*3 ) +  nip*9;
        cout << "flops = " << flops << ", flops/n = " << flops/ndof << endl;
      */

      // for (int cnt = 0; cnt < ndof; cnt++)
      {

        hmat1 = 0.0;
        for (int j = 0; j < ndof; j++)
          {
            int j2d = map3dto2d[j];
            for (int ix = 0; ix < nipx; ix++)
              TrafoGradientX(fac_xdx(j, ix)) . TransformAddT (elx(j), hmat1(j2d, ix));
          }
        NgProfiler::AddFlops (timerc, ndof*nipx*2);

        hmat2 = 0.0;
        for (int j = 0; j < ndof2d; j++)
          {
            int j1d = map2dto1d[j];
            for (int ix = 0; ix < nipx; ix++)
              for (int iy = 0; iy < nipy; iy++)
                TrafoGradientY(fac_ydy(j, iy)) . TransformAddT (hmat1(j, ix), hmat2(j1d, ix*nipy+iy));
          }
        NgProfiler::AddFlops (timerc, ndof2d*nipx*nipy*3);

        hmat3 = 0.0;
        for (int j = 0; j < ndof1d; j++)
          for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
            for (int iz = 0; iz < nipz; iz++, ii++)
              TrafoGradientZ(fac_zdz(j, iz)) . TransformAddT (hmat2(j, ixy), hmat3(ii));


        NgProfiler::AddFlops (timerc, ndof1d*nip*3);

        for (int ii = 0; ii < nip; ii++)
          {
            Vec<3> hv = hmat3(ii);
            hmat3(ii) = tdmat[ii] * hv;
          }
        NgProfiler::AddFlops (timerc, nip*9);

        hmat2 = 0.0;
        for (int j = 0; j < ndof1d; j++)
          for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
            for (int iz = 0; iz < nipz; iz++, ii++)
              TrafoGradientZ(fac_zdz(j, iz)) . TransformAdd (hmat3(ii), hmat2(j, ixy));

        NgProfiler::AddFlops (timerc, ndof1d*nip*3);

        hmat1 = 0.0;
        for (int j = 0; j < ndof2d; j++)
          {
            int j1d = map2dto1d[j];
            for (int ix = 0; ix < nipx; ix++)
              for (int iy = 0; iy < nipy; iy++)
                TrafoGradientY(fac_ydy(j, iy)) . TransformAdd (hmat2(j1d, ix*nipy+iy), hmat1(j, ix));
          }
        NgProfiler::AddFlops (timerc, ndof2d*nipx*nipy*3);

        ely = 0.0;
        for (int j = 0; j < ndof; j++)
          {
            int j2d = map3dto2d[j];
            for (int ix = 0; ix < nipx; ix++)
              TrafoGradientX(fac_xdx(j, ix)) . TransformAdd (hmat1(j2d, ix), ely(j));
          }
        NgProfiler::AddFlops (timerc, ndof*nipx*2);
      }
      

      NgProfiler::StopTimer (timerc);
      
    }

  };
















  // ******************************** Mass **************************************

  
  template <int D>
  class MassIntegratorTP : public MassIntegrator<D>
  { };







  template <>
  class MassIntegratorTP<2> : public MassIntegrator<2>
  {
    enum { D = 2 };
  public:
    MassIntegratorTP (CoefficientFunction * acoef)
      : MassIntegrator<D> (acoef)
    { ; }

    static Integrator * Create (ARRAY<CoefficientFunction*> & coeffs)
    {
      return new MassIntegratorTP (coeffs[0]);
    }

    virtual void
    AssembleElementMatrix (const FiniteElement & bfel,
			   const ElementTransformation & eltrans, 
			   FlatMatrix<double> & elmat,
			   LocalHeap & lh) const
    {
      const H1HighOrderTP<ET_TRIG> * fel_trig = dynamic_cast<const H1HighOrderTP<ET_TRIG>*> (&bfel);      

      const L2HighOrderTrigTP * fel_trigl2 = dynamic_cast<const L2HighOrderTrigTP*> (&bfel);      

      if (fel_trig)
        {
          TAssembleElementMatrix (*fel_trig, eltrans, elmat, lh);
          return;
        }

      if (fel_trigl2)
        {
          TAssembleElementMatrix (*fel_trigl2, eltrans, elmat, lh);
          return;
        }

      MassIntegrator<D>::AssembleElementMatrix(bfel, eltrans, elmat, lh);
    }


    /*
    virtual void
    ApplyElementMatrix (const FiniteElement & bfel,
                        const ElementTransformation & eltrans, 
                        const FlatVector<double> & elx, 
                        FlatVector<double> & ely,
                        void * precomputed,
                        LocalHeap & lh) const
    {
      
      const H1HighOrderTP<ET_TET> * fel_tet = dynamic_cast<const H1HighOrderTP<ET_TET>*> (&bfel);      
      const H1HighOrderTP<ET_PRISM> * fel_prism = dynamic_cast<const H1HighOrderTP<ET_PRISM>*> (&bfel);      
      const H1HighOrderTP<ET_HEX> * fel_hex = dynamic_cast<const H1HighOrderTP<ET_HEX>*> (&bfel);      


      if (fel_tet)
        {
          TApplyElementMatrix (*fel_tet, eltrans, elx, ely, precomputed, lh);
          return;
        }

      
      if (fel_prism)
        if (fel_prism->Sort(0)+3 == fel_prism->Sort(3))  // same singular vertex at bottom and top
          {
            TApplyElementMatrix (*fel_prism, eltrans, elx, ely, precomputed, lh);
            return;
          }
      
      if (fel_hex)
        {
          TApplyElementMatrix (*fel_hex, eltrans, elx, ely, precomputed, lh);
          return;
        }

      MassIntegrator<D>::ApplyElementMatrix(bfel, eltrans, elx, ely, precomputed, lh);
    }
    */



    template <class FEL>
    void TAssembleElementMatrix (const FEL & fel,
                                 const ElementTransformation & eltrans, 
                                 FlatMatrix<double> & elmat,
                                 LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("Fast assembling - Mass");
      static int timerg = NgProfiler::CreateTimer ("Fast assembling - geom");
      static int timeri = NgProfiler::CreateTimer ("Fast assembling - init");
      static int timerc = NgProfiler::CreateTimer ("Fast assembling - compute");

      NgProfiler::RegionTimer reg (timer);
      NgProfiler::StartTimer (timerg);      


      int ndof = fel.GetNDof();
      elmat.AssignMemory (ndof, ndof, lh);

      HeapReset hr (lh);

      int intorder = 2 * fel.Order();
      IntegrationRuleTP<D> irtp(eltrans, intorder, 1, lh);
      
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      
      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nip = nipx * nipy;
      
      FlatVector<Mat<1> > tdmat(nip, lh);      // transformed D-matrix
      
      // geometry calculation
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<D,D> sip(ip, eltrans, irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Mat<1> dmat;
          this->dmatop.GenerateMatrix (fel, sip, dmat, lh);
          dmat *= fabs (sip.GetJacobiDet()) * ip.Weight();
          tdmat[ii] = dmat;
        }

      NgProfiler::StopTimer (timerg);      
      NgProfiler::StartTimer (timeri);

      int ndof1d = fel.GetNDof1d();

      FlatArray<int> pmap1dto2d = fel.Get1dTo2dMappingPointer ();
      FlatArray<int> map1dto2d  = fel.Get1dTo2dMapping ();

      FlatMatrix<> fac_x(ndof,   nipx, lh); 
      FlatMatrix<> fac_y(ndof1d, nipy, lh); 
      
      for (int i1 = 0; i1 < nipy; i1++)
        fel.CalcYFactor (iry[i1](0), fac_y.Col(i1), lh);

      for (int i1 = 0; i1 < nipx; i1++)
        fel.CalcXFactor (irx[i1](0), fac_x.Col(i1), lh);

      FlatVector<double> allsum1d(nipx, lh);
      FlatVector<double> mult2d(nipx, lh);
      FlatVector<double> mult1d(nip, lh);

      NgProfiler::StopTimer (timeri);      



      NgProfiler::StartTimer (timerc);

      for (int ky = 0; ky < ndof1d; ky++)
        {
          for (int ix = 0, ii = 0; ix < nipx; ix++)
            for (int iy = 0; iy < nipy; iy++, ii++)
              mult1d(ii) = fac_y(ky,iy) * tdmat[ii](0,0);
          
          for (int jy = 0; jy <= ky; jy++)
            {
              for (int ix = 0, ii = 0; ix < nipx; ix++)
                {
                  double sum = 0.0;

                  for (int iy = 0; iy < nipy; iy++, ii++)
                    sum += fac_y(jy,iy) * mult1d(ii);
                  
                  allsum1d(ix) = sum;
                }
              
              for (int kx = pmap1dto2d[ky]; kx < pmap1dto2d[ky+1]; kx++)
                {
                  int kk = map1dto2d[kx];
                  
                  for (int ix = 0; ix < nipx; ix++)
                    mult2d(ix) = fac_x(kk,ix) * allsum1d(ix);
                          
                  for (int jx = pmap1dto2d[jy]; jx < pmap1dto2d[jy+1]; jx++)
                    {
                      int jj = map1dto2d[jx];
                              
                      double sum = 0;
                      for (int ix = 0; ix < nipx; ix++)
                        sum += fac_x(jj,ix) * mult2d(ix);
                              
                      elmat(kk,jj) = sum;
                      elmat(jj,kk) = sum;
                    }
                }
            }
        }

      NgProfiler::StopTimer (timerc);
    }




    /*
    template <class FEL>
    void TApplyElementMatrix (const FEL & fel,
                              const ElementTransformation & eltrans, 
                              const FlatVector<double> & elx, 
                              FlatVector<double> & ely,
                              void * precomputed,
                              LocalHeap & lh) const
    {
      int ndof = fel.GetNDof();

      HeapReset hr (lh);

      
      int order = fel.Order();

      FlatArray<int> map2dto1d = fel.Get2dTo1dMapping ();
      FlatArray<int> map3dto2d = fel.Get3dTo2dMapping ();

      int ndof2d = fel.GetNDof2d();
      int ndof1d = fel.GetNDof1d();
          
      int intorder = 2 * order;
      
      IntegrationRuleTP<D> irtp(eltrans, intorder, lh);
      
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      const IntegrationRule & irz = irtp.GetIRZ();
      
      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nipz = irz.GetNIP();
      int nip = nipx * nipy * nipz;
      
      FlatVector<Mat<1> > tdmat(nip, lh);      // transformed D-matrix
      

      // geometry calculation
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<3,3> sip(ip, eltrans, irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Mat<1> dmat;
          this->dmatop.GenerateMatrix (fel, sip, dmat, lh);
          dmat *= fabs (sip.GetJacobiDet()) * ip.Weight();
          tdmat[ii] = dmat;
        }


      FlatMatrix<> fac_x(ndof,   nipx, lh); 
      FlatMatrix<> fac_y(ndof2d, nipy, lh); 
      FlatMatrix<> fac_z(ndof1d, nipz, lh); 
      
      for (int i1 = 0; i1 < nipz; i1++)
        fel.CalcZFactor (irz[i1](0), fac_z.Col(i1), lh);

      for (int i1 = 0; i1 < nipy; i1++)
        fel.CalcYFactor (iry[i1](0), fac_y.Col(i1), lh);

      for (int i1 = 0; i1 < nipx; i1++)
        fel.CalcXFactor (irx[i1](0), fac_x.Col(i1), lh);

      int nipxy = nipx *nipy;
      
      FlatMatrix<double> hmat1(ndof2d, nipx, lh);
      FlatMatrix<double> hmat2(ndof1d, nipx*nipy, lh);
      FlatVector<double> hmat3(nip, lh);

      hmat1 = 0.0;
      for (int j = 0; j < ndof; j++)
        {
          int j2d = map3dto2d[j];
          for (int ix = 0; ix < nipx; ix++)
            hmat1(j2d,ix) += fac_x(j, ix) * elx(j);
        }

      hmat2 = 0.0;
      for (int j = 0; j < ndof2d; j++)
        {
          int j1d = map2dto1d[j];
          for (int ix = 0; ix < nipx; ix++)
            for (int iy = 0; iy < nipy; iy++)
              hmat2(j1d, ix*nipy+iy) += fac_y(j, iy) * hmat1(j, ix);
        }
      
      hmat3 = 0.0;
      for (int j = 0; j < ndof1d; j++)
        for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
          for (int iz = 0; iz < nipz; iz++, ii++)
            hmat3(ii) += fac_z(j, iz) * hmat2(j, ixy);
      
      for (int ii = 0; ii < nip; ii++)
        hmat3(ii) *= tdmat[ii](0,0);
      
      hmat2 = 0.0;
      for (int j = 0; j < ndof1d; j++)
        for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
          for (int iz = 0; iz < nipz; iz++, ii++)
            hmat2(j,ixy) += fac_z(j, iz) * hmat3(ii);
      
      hmat1 = 0.0;
      for (int j = 0; j < ndof2d; j++)
        {
          int j1d = map2dto1d[j];
          for (int ix = 0; ix < nipx; ix++)
            for (int iy = 0; iy < nipy; iy++)
              hmat1(j, ix) += fac_y(j, iy) * hmat2(j1d, ix*nipy+iy);
        }
      
      ely = 0.0;
      for (int j = 0; j < ndof; j++)
        {
          int j2d = map3dto2d[j];
          for (int ix = 0; ix < nipx; ix++)
            ely(j) += fac_x(j, ix) * hmat1(j2d, ix);
        }
    }
    */
  };
















  template <>
  class MassIntegratorTP<3> : public MassIntegrator<3>
  {
    enum { D = 3 };
  public:
    MassIntegratorTP (CoefficientFunction * acoef)
      : MassIntegrator<D> (acoef)
    { ; }

    static Integrator * Create (ARRAY<CoefficientFunction*> & coeffs)
    {
      return new MassIntegratorTP (coeffs[0]);
    }

    virtual void
    AssembleElementMatrix (const FiniteElement & bfel,
			   const ElementTransformation & eltrans, 
			   FlatMatrix<double> & elmat,
			   LocalHeap & lh) const
    {
      if (bfel.IsTPElement())
        {
          const H1HighOrderTP<ET_TET> * fel_tet = dynamic_cast<const H1HighOrderTP<ET_TET>*> (&bfel);      
          const H1HighOrderTP<ET_PRISM> * fel_prism = dynamic_cast<const H1HighOrderTP<ET_PRISM>*> (&bfel);      
          const H1HighOrderTP<ET_HEX> * fel_hex = dynamic_cast<const H1HighOrderTP<ET_HEX>*> (&bfel);      

          const L2HighOrderTetTP * fel_tetl2 = dynamic_cast<const L2HighOrderTetTP*> (&bfel);                

          if (fel_tet)
            TAssembleElementMatrix (*fel_tet, eltrans, elmat, lh);
          
          if (fel_prism)
            TAssembleElementMatrix (*fel_prism, eltrans, elmat, lh);

          if (fel_hex)
            TAssembleElementMatrix (*fel_hex, eltrans, elmat, lh);

          if (fel_tetl2)
            TAssembleElementMatrix (*fel_tetl2, eltrans, elmat, lh);

          return;
        }

      MassIntegrator<D>::AssembleElementMatrix(bfel, eltrans, elmat, lh);
    }



    virtual void
    ApplyElementMatrix (const FiniteElement & bfel,
                        const ElementTransformation & eltrans, 
                        const FlatVector<double> & elx, 
                        FlatVector<double> & ely,
                        void * precomputed,
                        LocalHeap & lh) const
    {
      if (bfel.IsTPElement())
        {
          const H1HighOrderTP<ET_TET> * fel_tet = dynamic_cast<const H1HighOrderTP<ET_TET>*> (&bfel);      
          const H1HighOrderTP<ET_PRISM> * fel_prism = dynamic_cast<const H1HighOrderTP<ET_PRISM>*> (&bfel);      
          const H1HighOrderTP<ET_HEX> * fel_hex = dynamic_cast<const H1HighOrderTP<ET_HEX>*> (&bfel);      
          
          if (fel_tet)
            TApplyElementMatrix (*fel_tet, eltrans, elx, ely, precomputed, lh);

          if (fel_prism)
            TApplyElementMatrix (*fel_prism, eltrans, elx, ely, precomputed, lh);
      
          if (fel_hex)
            TApplyElementMatrix (*fel_hex, eltrans, elx, ely, precomputed, lh);

          return;
        }

      MassIntegrator<D>::ApplyElementMatrix(bfel, eltrans, elx, ely, precomputed, lh);
    }



    template <class FEL>
    void TAssembleElementMatrix (const FEL & fel,
                                 const ElementTransformation & eltrans, 
                                 FlatMatrix<double> & elmat,
                                 LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("Fast assembling - Mass");
      static int timerg = NgProfiler::CreateTimer ("Fast assembling - geom");
      static int timeri = NgProfiler::CreateTimer ("Fast assembling - init");
      static int timerc = NgProfiler::CreateTimer ("Fast assembling - compute");

      NgProfiler::RegionTimer reg (timer);
      NgProfiler::StartTimer (timerg);      


      int ndof = fel.GetNDof();
      elmat.AssignMemory (ndof, ndof, lh);

      HeapReset hr (lh);

      int intorder = 2 * fel.Order();
      IntegrationRuleTP<D> irtp(eltrans, intorder, 1, lh);
      
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      const IntegrationRule & irz = irtp.GetIRZ();
      
      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nipz = irz.GetNIP();
      int nip = nipx * nipy * nipz;
      
      FlatVector<Mat<1> > tdmat(nip, lh);      // transformed D-matrix
      
      // geometry calculation
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<3,3> sip(ip, eltrans, irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Mat<1> dmat;
          this->dmatop.GenerateMatrix (fel, sip, dmat, lh);
          dmat *= fabs (sip.GetJacobiDet()) * ip.Weight();
          tdmat[ii] = dmat;
        }

      NgProfiler::StopTimer (timerg);      
      NgProfiler::StartTimer (timeri);

      int ndof2d = fel.GetNDof2d();
      int ndof1d = fel.GetNDof1d();

      FlatArray<int> pmap1dto2d = fel.Get1dTo2dMappingPointer ();
      FlatArray<int> map1dto2d  = fel.Get1dTo2dMapping ();

      FlatArray<int> pmap2dto3d = fel.Get2dTo3dMappingPointer ();
      FlatArray<int> map2dto3d  = fel.Get2dTo3dMapping ();

      FlatMatrix<> fac_x(ndof,   nipx, lh); 
      FlatMatrix<> fac_y(ndof2d, nipy, lh); 
      FlatMatrix<> fac_z(ndof1d, nipz, lh); 
      
      for (int i1 = 0; i1 < nipz; i1++)
        fel.CalcZFactor (irz[i1](0), fac_z.Col(i1), lh);

      for (int i1 = 0; i1 < nipy; i1++)
        fel.CalcYFactor (iry[i1](0), fac_y.Col(i1), lh);

      for (int i1 = 0; i1 < nipx; i1++)
        fel.CalcXFactor (irx[i1](0), fac_x.Col(i1), lh);

      int nipxy = nipx *nipy;
      FlatVector<double> allsum2d(nipx, lh);
      FlatVector<double> allsum1d(nipxy, lh);
      FlatVector<double> mult3d(nipx, lh);
      FlatVector<double> mult2d(nipxy, lh);
      FlatVector<double> mult1d(nip, lh);

      NgProfiler::StopTimer (timeri);      



      NgProfiler::StartTimer (timerc);

      for (int kz = 0; kz < ndof1d; kz++)
        {
          for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
            for (int iz = 0; iz < nipz; iz++, ii++)
              mult1d(ii) = fac_z(kz,iz) * tdmat[ii](0,0);
          
          for (int jz = 0; jz <= kz; jz++)
            {
              for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
                {
                  double sum = 0.0;

                  for (int iz = 0; iz < nipz; iz++, ii++)
                    sum += fac_z(jz,iz) * mult1d(ii);
                  
                  allsum1d(ixy) = sum;
                }
              
              for (int ky = pmap1dto2d[kz]; ky < pmap1dto2d[kz+1]; ky++)
                {
                  int kyz = map1dto2d[ky];
                  for (int ix = 0, ixy=0; ix < nipx; ix++)
                    for (int iy = 0; iy < nipy; iy++, ixy++)
                      mult2d(ixy) = fac_y(kyz,iy) * allsum1d(ixy);
              
                  for (int jy = pmap1dto2d[jz]; jy < pmap1dto2d[jz+1]; jy++)
                    {
                      if (jy > ky) continue;
                      int jyz = map1dto2d[jy];

                      for (int ix = 0; ix < nipx; ix++)
                        {
                          double hsum = 0.0;
                          
                          for (int iy = 0; iy < nipy; iy++)
                            hsum += fac_y(jyz,iy) * mult2d(ix*nipy+iy);
                          
                          allsum2d(ix) = hsum;
                        }
                      
                      for (int kx = pmap2dto3d[kyz]; kx < pmap2dto3d[kyz+1]; kx++)
                        {
                          int kk = map2dto3d[kx];

                          for (int ix = 0; ix < nipx; ix++)
                            mult3d(ix) = fac_x(kk,ix) * allsum2d(ix);
                          
                          for (int jx = pmap2dto3d[jyz]; jx < pmap2dto3d[jyz+1]; jx++)
                            {
                              int jj = map2dto3d[jx];
                              
                              double sum = 0;
                              for (int ix = 0; ix < nipx; ix++)
                                sum += fac_x(jj,ix) * mult3d(ix);
                              
                              elmat(kk,jj) = sum;
                              elmat(jj,kk) = sum;
                            }
                        }
                    }
                }
            }
        }

      NgProfiler::StopTimer (timerc);

      NgProfiler::AddFlops (timerc, 9*nip*ndof1d);
      NgProfiler::AddFlops (timerc, 9*nip*ndof1d*(ndof1d+1)/2);

      NgProfiler::AddFlops (timerc, 9*nipxy*ndof1d*ndof2d);
      NgProfiler::AddFlops (timerc, 3*nipxy*ndof2d*(ndof2d+1));

      NgProfiler::AddFlops (timerc, 4*ndof2d*ndof*nipx);
      NgProfiler::AddFlops (timerc, ndof*(ndof+1)*nipx);
    }





    template <class FEL>
    void TApplyElementMatrix (const FEL & fel,
                              const ElementTransformation & eltrans, 
                              const FlatVector<double> & elx, 
                              FlatVector<double> & ely,
                              void * precomputed,
                              LocalHeap & lh) const
    {
      int ndof = fel.GetNDof();

      HeapReset hr (lh);

      
      int order = fel.Order();

      FlatArray<int> map2dto1d = fel.Get2dTo1dMapping ();
      FlatArray<int> map3dto2d = fel.Get3dTo2dMapping ();

      int ndof2d = fel.GetNDof2d();
      int ndof1d = fel.GetNDof1d();
          
      int intorder = 2 * order;
      
      IntegrationRuleTP<D> irtp(eltrans, intorder, 1, lh);
      
      const IntegrationRule & irx = irtp.GetIRX();
      const IntegrationRule & iry = irtp.GetIRY();
      const IntegrationRule & irz = irtp.GetIRZ();
      
      int nipx = irx.GetNIP();
      int nipy = iry.GetNIP();
      int nipz = irz.GetNIP();
      int nip = nipx * nipy * nipz;
      
      FlatVector<Mat<1> > tdmat(nip, lh);      // transformed D-matrix
      

      // geometry calculation
      for (int ii = 0; ii < nip; ii++)
        {
          IntegrationPoint ip(irtp.GetXi(ii), irtp.GetWeight(ii));
          SpecificIntegrationPoint<3,3> sip(ip, eltrans, irtp.GetPoint(ii), irtp.GetJacobian(ii), lh);
          
          Mat<1> dmat;
          this->dmatop.GenerateMatrix (fel, sip, dmat, lh);
          dmat *= fabs (sip.GetJacobiDet()) * ip.Weight();
          tdmat[ii] = dmat;
        }


      FlatMatrix<> fac_x(ndof,   nipx, lh); 
      FlatMatrix<> fac_y(ndof2d, nipy, lh); 
      FlatMatrix<> fac_z(ndof1d, nipz, lh); 
      
      for (int i1 = 0; i1 < nipz; i1++)
        fel.CalcZFactor (irz[i1](0), fac_z.Col(i1), lh);

      for (int i1 = 0; i1 < nipy; i1++)
        fel.CalcYFactor (iry[i1](0), fac_y.Col(i1), lh);

      for (int i1 = 0; i1 < nipx; i1++)
        fel.CalcXFactor (irx[i1](0), fac_x.Col(i1), lh);

      int nipxy = nipx *nipy;
      
      FlatMatrix<double> hmat1(ndof2d, nipx, lh);
      FlatMatrix<double> hmat2(ndof1d, nipx*nipy, lh);
      FlatVector<double> hmat3(nip, lh);

      hmat1 = 0.0;
      for (int j = 0; j < ndof; j++)
        {
          int j2d = map3dto2d[j];
          for (int ix = 0; ix < nipx; ix++)
            hmat1(j2d,ix) += fac_x(j, ix) * elx(j);
        }

      hmat2 = 0.0;
      for (int j = 0; j < ndof2d; j++)
        {
          int j1d = map2dto1d[j];
          for (int ix = 0; ix < nipx; ix++)
            for (int iy = 0; iy < nipy; iy++)
              hmat2(j1d, ix*nipy+iy) += fac_y(j, iy) * hmat1(j, ix);
        }
      
      hmat3 = 0.0;
      for (int j = 0; j < ndof1d; j++)
        for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
          for (int iz = 0; iz < nipz; iz++, ii++)
            hmat3(ii) += fac_z(j, iz) * hmat2(j, ixy);
      
      for (int ii = 0; ii < nip; ii++)
        hmat3(ii) *= tdmat[ii](0,0);
      
      hmat2 = 0.0;
      for (int j = 0; j < ndof1d; j++)
        for (int ixy = 0, ii = 0; ixy < nipxy; ixy++)
          for (int iz = 0; iz < nipz; iz++, ii++)
            hmat2(j,ixy) += fac_z(j, iz) * hmat3(ii);
      
      hmat1 = 0.0;
      for (int j = 0; j < ndof2d; j++)
        {
          int j1d = map2dto1d[j];
          for (int ix = 0; ix < nipx; ix++)
            for (int iy = 0; iy < nipy; iy++)
              hmat1(j, ix) += fac_y(j, iy) * hmat2(j1d, ix*nipy+iy);
        }
      
      ely = 0.0;
      for (int j = 0; j < ndof; j++)
        {
          int j2d = map3dto2d[j];
          for (int ix = 0; ix < nipx; ix++)
            ely(j) += fac_x(j, ix) * hmat1(j2d, ix);
        }
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
      GetIntegrators().AddLFIntegrator ("sourcetp", 2, 1,
                                        SourceIntegratorTP<2>::Create);
      GetIntegrators().AddLFIntegrator ("sourcetp", 3, 1,
                                        SourceIntegratorTP<3>::Create);

      GetIntegrators().AddLFIntegrator ("neumanntp", 2, 1,
					NeumannIntegratorTP<2>::Create);
      GetIntegrators().AddLFIntegrator ("neumanntp", 3, 1,
					NeumannIntegratorTP<3>::Create);

      GetIntegrators().AddBFIntegrator ("laplacetp", 2, 1,
                                        LaplaceIntegratorTP<2>::Create);
      GetIntegrators().AddBFIntegrator ("laplacetp", 3, 1,
                                        LaplaceIntegratorTP<3>::Create);

      GetIntegrators().AddBFIntegrator ("masstp", 2, 1,
                                        MassIntegratorTP<2>::Create);
      GetIntegrators().AddBFIntegrator ("masstp", 3, 1,
                                        MassIntegratorTP<3>::Create);

    }
    
    Init init;
  }

}
