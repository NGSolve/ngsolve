#ifdef TPFO
#include "l2hofetpfo.cpp"
#else


/*********************************************************************/
/* File:   l2hofetp.cpp                                              */
/* Author: J. Schöberl, Aachen                                       */
/* Date:   4. Oct 2008                                               */
/*********************************************************************/

#define xxNOPROFILE
#define noDEBUG

#include <fem.hpp>
namespace ngfem
{
  using namespace ngfem;





  ARRAY<RecPol*> L2HighOrderTrigTP :: recpols;

  FlatMatrix<> L2HighOrderTrigTP :: fac_x (0, 0, 0);
  FlatMatrix<> L2HighOrderTrigTP :: fac_y (0, 0, 0);
  FlatMatrix<AutoDiff<1> > L2HighOrderTrigTP :: fac_xdx (0, 0, 0);
  FlatMatrix<AutoDiff<1> > L2HighOrderTrigTP :: fac_ydy (0, 0, 0);


  L2HighOrderTrigTP :: L2HighOrderTrigTP (int aorder, LocalHeap & alh)
    : L2HighOrderTrig (aorder), lh(alh)
  {
    if (recpols.Size() < order+1)
      {
        for (int i = 0; i < recpols.Size(); i++)
          delete recpols[i];

        recpols.SetSize(order+1);
        for (int i = 0; i <= order; i++)
          recpols[i] = new RecPolJacobi(order+1, 2*i+1, 0);
      }
  }

  
  void L2HighOrderTrigTP :: SetVertexNumbers (FlatArray<int> & avnums) // , LocalHeap & lh)
  {
    L2HighOrderTrig :: SetVertexNumbers (avnums, lh);

    for (int i = 0; i < 3; i++) sort[i] = i;
 
    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
    if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
    if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
    // vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]] 



    trig2tensor = FlatArray<int[3]> (ndof, lh);

    // element dofs
    int pi = order_inner[0];
    int ii = 0;
    for (int j2 = 0; j2 <= pi; j2++)
      for (int i2 = 0; i2 <= pi-j2; i2++, ii++)
        {
          trig2tensor[ii][0] = j2;
          trig2tensor[ii][1] = i2;
          trig2tensor[ii][2] = 0;
        }


    int horder = order+1;
    int nd1d = 2 * horder;

    FlatArray<int> tensor2segm (nd1d, lh);
    tensor2segm = -1;
    
    for (int ii = 0; ii < ndof; ii++)
      tensor2segm[trig2tensor[ii][1] * 2 + trig2tensor[ii][2]] = 0;

    ndof1d = 0;
    for (int i = 0; i < nd1d; i++)
      if (tensor2segm[i] != -1)
        ndof1d++;

    segm2tensor = FlatArray<int[2]> (ndof1d, lh);

    ndof1d = 0;
    for (int i3 = 0, i1d = 0; i3 < horder; i3++)
      for (int i4 = 0; i4 < 2; i4++, i1d++)
        if (tensor2segm[i1d] != -1)
          {
            tensor2segm[i1d] = ndof1d;
            segm2tensor[ndof1d][0] = i3;
            segm2tensor[ndof1d][1] = i4;
            ndof1d++;
          }

    map2dto1d = FlatArray<int> (ndof, lh);
    for (int i = 0; i < ndof; i++)
      map2dto1d[i] = tensor2segm[trig2tensor[i][1] * 2 + trig2tensor[i][2]];

    // invert mapping:
    pmap1dto2d = FlatArray<int> (ndof1d+2, lh);
    map1dto2d = FlatArray<int> (ndof, lh);
    pmap1dto2d = 0;
    for (int j = 0; j < ndof; j++)
      pmap1dto2d[map2dto1d[j]+2]++;

    for (int j = 2; j <= ndof1d; j++)
      pmap1dto2d[j] += pmap1dto2d[j-1]; 

    for (int j = 0; j < ndof; j++)
      map1dto2d[pmap1dto2d[map2dto1d[j]+1]++] = j;
  }
  
  void L2HighOrderTrigTP :: CalcShape (const IntegrationPoint & ip, 
                                       FlatVector<> shape) const
  {
    T_CalcShape (ip(0), ip(1), shape); 
  }
  
  void L2HighOrderTrigTP :: CalcDShape (const IntegrationPoint & ip, 
                                        FlatMatrix<> dshape) const
  {
    AutoDiff<2> x(ip(0), 0);
    AutoDiff<2> y(ip(1), 1);

    ArrayMem<AutoDiff<2>,40> sds(ndof);
    T_CalcShape (x,y,sds);

    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 2; j++)
	dshape(i, j) = sds[i].DValue(j);
  }

  template<typename Tx, typename Ty, typename TFA>  
  void  L2HighOrderTrigTP :: T_CalcShape (Tx x, Ty y, TFA & sds) const
  {
    Tx lami[3] = { x, y,  1-x-y };

    int sort[] = { 0, 1, 2 };
    for (int i = 0; i < 3; i++) sort[i] = i;
 
    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
    if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
    if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
    // vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]] 

    Tx lamis[3];
    for (int i = 0; i < 3; i++)
      lamis[i] = lami[sort[i]];


#ifdef JACOBI
    int horder = order+1;

    ArrayMem<Tx, 100> memx((order+1)*horder);

    FlatMatrix<Tx> polsx(order+1, horder, &memx[0]);
    VectorMem<20, Tx> polsy(horder);
    
    for (int i = 0; i <= order; i++)
      JacobiPolynomial (horder-1, 2*lamis[0]-1, max(0,2*i+1), 0, polsx.Row(i));

    ScaledLegendrePolynomial (horder-1, lamis[1]-lamis[2], lamis[1]+lamis[2], polsy);

    for (int i = 0; i < ndof; i++)
      sds[i] =
        polsx(trig2tensor[i][1]+trig2tensor[i][2], trig2tensor[i][0]) * 
        polsy(trig2tensor[i][1]);
#else

    cout << "not implemented" << endl;

    int horder = max(2, order);
    ArrayMem<Tx, 20> polx(horder), poly(horder), polyy(2); 

    Tx * hp = &polx[1];
    LegendrePolynomialMult (horder-2, 2*lamis[0]-1, lamis[0], hp);
    polx[0] = 1.0;
    
    hp = &poly[1];
    ScaledLegendrePolynomialMult (horder-2, lamis[1]-lamis[2], lamis[1] + lamis[2], lamis[1], hp);
    poly[0] = 1.0;

    polyy[0] = 1.0;
    polyy[1] = lamis[2];

    ArrayMem<int[3], 100> trig2tensor(ndof);
    GetTrig2TensorMapping (trig2tensor);

    for (int i = 0; i < ndof; i++)
      sds[i] =
        polx[trig2tensor[i][0]] * 
        poly[trig2tensor[i][1]] * 
        polyy[trig2tensor[i][2]];
#endif
  }


  void L2HighOrderTrigTP :: CalcTPFactors (const IntegrationRuleTP<2> & irtp,
                                           TPFactors<2> & factors) const
  {
    HeapReset hr(lh);

    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
    
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();

    factors.fac_x.SetSize (ndof, nipx);
    factors.fac_y.SetSize (ndof1d, nipy);

    for (int i1 = 0; i1 < nipy; i1++)
      CalcYFactor (iry[i1](0), factors.fac_y.Col(i1), lh);
    
    for (int i1 = 0; i1 < nipx; i1++)
      CalcXFactor (irx[i1](0), factors.fac_x.Col(i1), lh);
  }

  void L2HighOrderTrigTP :: SetTPFactors (const TPFactors<2> & factors)
  {
    fac_x.AssignMemory (factors.fac_x.Height(),
                        factors.fac_x.Width(),
                        &factors.fac_x(0,0));
                          
    fac_y.AssignMemory (factors.fac_y.Height(),
                        factors.fac_y.Width(),
                        &factors.fac_y(0,0));
  }





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




  void L2HighOrderTrigTP ::
  EvaluateShapeGrid (const IntegrationRuleTP<2> & irtp,
                     const FlatVector<double> coefs,
                     FlatVector<double> gridvalues,
                     LocalHeap & lh) const
  {
    HeapReset hr (lh);
    
    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
      
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();
    int nip = nipx * nipy;

    /*
      FlatMatrix<> fac_x(ndof,   nipx, lh); 
      FlatMatrix<> fac_y(ndof1d, nipy, lh); 
    
      for (int i1 = 0; i1 < nipy; i1++)
      CalcYFactor (iry[i1](0), fac_y.Col(i1), lh);
    
      for (int i1 = 0; i1 < nipx; i1++)
      CalcXFactor (irx[i1](0), fac_x.Col(i1), lh);
    */

    if (fac_x.Height() != ndof)
      {
        FlatMatrix<> hfac_x(ndof,   nipx, new double[ndof*nipx]); 
        FlatMatrix<> hfac_y(ndof1d, nipy, new double[ndof1d*nipy]); 
        
        for (int i1 = 0; i1 < nipy; i1++)
          CalcYFactor (iry[i1](0), hfac_y.Col(i1), lh);
        
        for (int i1 = 0; i1 < nipx; i1++)
          CalcXFactor (irx[i1](0), hfac_x.Col(i1), lh);

        
        const_cast<FlatMatrix<> &> (fac_x).AssignMemory (ndof  , nipx, &hfac_x(0,0));
        const_cast<FlatMatrix<> &> (fac_y).AssignMemory (ndof1d, nipy, &hfac_y(0,0));
      }


    FlatMatrix<> mgrid(nipx, nipy, &gridvalues(0));
    FlatMatrix<> gridx(ndof1d, nipx, lh);

    gridx = 0.0;
    for (int i = 0; i < ndof; i++)
      gridx.Row(map2dto1d[i]) += coefs(i) * fac_x.Row(i);
    
    mgrid = 0.0;
    for (int ix = 0; ix < nipx; ix++)
      {
        FlatVector<> row = mgrid.Row(ix);
        for (int j = 0; j < ndof1d; j++)
          row += gridx(j,ix) * fac_y.Row(j);
      }
  }


  void L2HighOrderTrigTP ::
  EvaluateShapeGridTrans (const IntegrationRuleTP<2> & irtp,
                          const FlatVector<double> gridvalues,
                          FlatVector<double> coefs,
                          LocalHeap & lh) const
  {
    HeapReset hr (lh);
    
    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
      
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();
    int nip = nipx * nipy;

    if (fac_x.Height() != ndof || fac_x.Width() != nipx || fac_y.Width() != nipy)
      {
        FlatMatrix<> hfac_x(ndof,   nipx, new double[ndof*nipx]); 
        FlatMatrix<> hfac_y(ndof1d, nipy, new double[ndof1d*nipy]); 
        
        for (int i1 = 0; i1 < nipy; i1++)
          CalcYFactor (iry[i1](0), hfac_y.Col(i1), lh);
        
        for (int i1 = 0; i1 < nipx; i1++)
          CalcXFactor (irx[i1](0), hfac_x.Col(i1), lh);

        const_cast<FlatMatrix<> &> (fac_x).AssignMemory (ndof  , nipx, &hfac_x(0,0));
        const_cast<FlatMatrix<> &> (fac_y).AssignMemory (ndof1d, nipy, &hfac_y(0,0));
      }

    /*
      FlatMatrix<> fac_x(ndof,   nipx, lh); 
      FlatMatrix<> fac_y(ndof1d, nipy, lh); 
      
      for (int i1 = 0; i1 < nipy; i1++)
      CalcYFactor (iry[i1](0), fac_y.Col(i1), lh);

      for (int i1 = 0; i1 < nipx; i1++)
      CalcXFactor (irx[i1](0), fac_x.Col(i1), lh);
    */

    FlatMatrix<> mgrid(nipx, nipy, &gridvalues(0));
    FlatMatrix<> gridx(ndof1d, nipx, lh);

    for (int j = 0; j < ndof1d; j++)
      for (int ix = 0; ix < nipx; ix++)
        gridx(j, ix) = InnerProduct (mgrid.Row(ix), fac_y.Row(j));

    for (int i = 0; i < ndof; i++)
      coefs(i) = InnerProduct (fac_x.Row(i), gridx.Row(map2dto1d[i]));
  }



  void  L2HighOrderTrigTP ::
  EvaluateDShapeGrid (const IntegrationRuleTP<2> & irtp,
                      const FlatVector<double> coefs,
                      FlatMatrixFixWidth<2> gridvalues,
                      LocalHeap & lh) const
  {
    HeapReset hr (lh);

    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
      
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();
    int nip = nipx * nipy;

    // FlatMatrix<AutoDiff<1> > fac_xdx(ndof, nipx, lh);  
    // FlatMatrix<AutoDiff<1> > fac_ydy(ndof1d, nipy, lh);
    

    if (fac_xdx.Height() != ndof)
      {
        FlatMatrix<AutoDiff<1> > & hfac_xdx = const_cast< FlatMatrix<AutoDiff<1> > & > (fac_xdx);
        FlatMatrix<AutoDiff<1> > & hfac_ydy = const_cast< FlatMatrix<AutoDiff<1> > & > (fac_ydy);
        
        hfac_xdx.AssignMemory (ndof, nipx, new AutoDiff<1>[ndof*nipx]);
        hfac_ydy.AssignMemory (ndof1d, nipy, new AutoDiff<1>[ndof1d*nipy]);
        
        for (int i1 = 0; i1 < nipy; i1++)
          {
            AutoDiff<1> y (iry[i1](0), 0);
            CalcYFactor (y, hfac_ydy.Col(i1), lh);
          }
        
        for (int i1 = 0; i1 < nipx; i1++)
          {
            AutoDiff<1> x (irx[i1](0), 0);
            CalcXFactor (x, hfac_xdx.Col(i1), lh);
          }
      }

    FlatMatrix<Vec<2> > gridx(ndof1d, nipx, lh);
    FlatMatrix<Vec<2> > mgrid(nipx, nipy, lh);
    
    gridx = 0.0;
    for (int i = 0; i < ndof; i++)
      {
        int i2d = map2dto1d[i];
        double coefi = coefs(i);
        for (int ix = 0; ix < nipx; ix++)
          TrafoGradient2dX(fac_xdx(i, ix)) . TransformAddT (coefi, gridx(i2d, ix));
      }
    
    mgrid = 0.0;
    for (int ix = 0; ix < nipx; ix++)
      for (int j = 0; j < ndof1d; j++)
        {
          Vec<2> hv = gridx(j,ix);
          for (int iy = 0; iy < nipy; iy++)
            TrafoGradient2dY(fac_ydy(j, iy)) . TransformAddT (hv /* gridx(j, ix) */, mgrid(ix, iy));
        }
    for (int ix = 0, ii = 0; ix < nipx; ix++)
      for (int iy = 0; iy < nipy; iy++, ii++)
        gridvalues.Row(ii) = Trans (irtp.GetDuffyJacobian(ii)) * mgrid(ix, iy);
  }
  
  void  L2HighOrderTrigTP ::
  EvaluateDShapeGridTrans (const IntegrationRuleTP<2> & irtp,
                           const FlatMatrixFixWidth<2> gridvalues,
                           FlatVector<double> coefs,
                           LocalHeap & lh) const
  {
    HeapReset hr (lh);

    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
      
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();
    int nip = nipx * nipy;

    // FlatMatrix<AutoDiff<1> > fac_xdx(ndof, nipx, lh);  
    // FlatMatrix<AutoDiff<1> > fac_ydy(ndof1d, nipy, lh);

    
    if (fac_xdx.Height() != ndof)
      {
        FlatMatrix<AutoDiff<1> > & hfac_xdx = const_cast< FlatMatrix<AutoDiff<1> > & > (fac_xdx);
        FlatMatrix<AutoDiff<1> > & hfac_ydy = const_cast< FlatMatrix<AutoDiff<1> > & > (fac_ydy);
        
        hfac_xdx.AssignMemory (ndof, nipx, new AutoDiff<1>[ndof*nipx]);
        hfac_ydy.AssignMemory (ndof1d, nipy, new AutoDiff<1>[ndof1d*nipy]);
        
        for (int i1 = 0; i1 < nipy; i1++)
          {
            AutoDiff<1> y (iry[i1](0), 0);
            CalcYFactor (y, hfac_ydy.Col(i1), lh);
          }
        
        for (int i1 = 0; i1 < nipx; i1++)
          {
            AutoDiff<1> x (irx[i1](0), 0);
            CalcXFactor (x, hfac_xdx.Col(i1), lh);
          }
      }

    FlatMatrix<Vec<2> > gridx(ndof1d, nipx, lh);
    FlatMatrix<Vec<2> > mgrid(nipx, nipy, lh);
    

    
    for (int ix = 0, ii = 0; ix < nipx; ix++)
      for (int iy = 0; iy < nipy; iy++, ii++)
        mgrid(ix, iy) = irtp.GetDuffyJacobian(ii) * gridvalues.Row(ii);

    for (int ix = 0; ix < nipx; ix++)
      for (int j = 0; j < ndof1d; j++)
        {
          Vec<2> sum(0,0);
          for (int iy = 0; iy < nipy; iy++)
            TrafoGradient2dY(fac_ydy(j, iy)) . TransformAdd (mgrid(ix, iy), sum); 
          gridx(j,ix) = sum;
        }

    for (int i = 0; i < ndof; i++)
      {
        int i2d = map2dto1d[i];
        double sum = 0;
        for (int ix = 0; ix < nipx; ix++)
          TrafoGradient2dX(fac_xdx(i, ix)) . TransformAdd (gridx(i2d, ix), sum);
        coefs(i) = sum;
      }
  }










  

  ////////////////// L2 Tet



  ARRAY<RecPol*> L2HighOrderTetTP :: recpols;


  // FlatMatrix<> L2HighOrderTetTP :: fac_x (0, 0, 0);
  // FlatMatrix<> L2HighOrderTetTP :: fac_y (0, 0, 0);
  // FlatMatrix<> L2HighOrderTetTP :: fac_z (0, 0, 0);
  FlatMatrix<AutoDiff<1> > L2HighOrderTetTP :: fac_xdx (0, 0, 0);
  FlatMatrix<AutoDiff<1> > L2HighOrderTetTP :: fac_ydy (0, 0, 0);
  FlatMatrix<AutoDiff<1> > L2HighOrderTetTP :: fac_zdz (0, 0, 0);



  L2HighOrderTetTP :: L2HighOrderTetTP (int aorder, LocalHeap & alh)
    : L2HighOrderTet (aorder), lh(alh)
  {
    tp = true;

    if (recpols.Size() < 2*order+3)
      {
        for (int i = 0; i < recpols.Size(); i++)
          delete recpols[i];

        recpols.SetSize(2*order+3);
        for (int i = 0; i <= 2*order+2; i++)
          recpols[i] = new RecPolJacobi(order+1, i, 0);
      }
  }

  void L2HighOrderTetTP :: SetVertexNumbers (FlatArray<int> & avnums)
  {
    L2HighOrderTet :: SetVertexNumbers (avnums);

    for (int i = 0; i < 4; i++) sort[i] = i;

    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
    if (vnums[sort[2]] > vnums[sort[3]]) Swap (sort[2], sort[3]);
    if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
    if (vnums[sort[1]] > vnums[sort[3]]) Swap (sort[1], sort[3]);
    if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);


    tet2tensor = FlatArray<int[4]> (ndof, lh);

    // element dofs
    int pi = order_inner[0];
    int ii = 0;
    for (int j2 = 0; j2 <= pi; j2++)
      for (int i2 = 0; i2 <= pi-j2; i2++)
        for (int k2 = 0; k2 <= pi-j2-i2; k2++, ii++)
          {
            tet2tensor[ii][0] = j2;
            tet2tensor[ii][1] = i2;
            tet2tensor[ii][2] = k2;
            tet2tensor[ii][3] = 0;
          }


    int horder = order+1;

    int nd1d = 2 * horder;
    int nd2d = horder * nd1d;

    FlatArray<int> tensor2trig (nd2d, lh);
    tensor2trig = -1;

    for (int ii = 0; ii < ndof; ii++)
      {
        int j2d = tet2tensor[ii][1] * nd1d + tet2tensor[ii][2] * 2 + tet2tensor[ii][3];
        tensor2trig[j2d] = 0;
      }
    
    int ndtrig = 0;
    for (int i = 0; i < nd2d; i++)
      if (tensor2trig[i] != -1)
        ndtrig++;
    
    trig2tensor = FlatArray<int[3]> (ndtrig, lh);

    ndtrig = 0;
    for (int i2 = 0, i2d = 0; i2 < horder; i2++)
      for (int i3 = 0; i3 < horder; i3++)
        for (int i4 = 0; i4 < 2; i4++, i2d++)
          if (tensor2trig[i2d] != -1)
            {
              tensor2trig[i2d] = ndtrig;
              trig2tensor[ndtrig][0] = i2;
              trig2tensor[ndtrig][1] = i3;
              trig2tensor[ndtrig][2] = i4;
              ndtrig++;
            }

    FlatArray<int> tet2trig (ndof, lh);
    for (int i = 0; i < ndof; i++)
      tet2trig[i] = tensor2trig[tet2tensor[i][1] * nd1d + tet2tensor[i][2] * 2 + tet2tensor[i][3]];


    FlatArray<int> tensor2segm ( nd1d, lh);
    tensor2segm = -1;
    
    for (int ii = 0; ii < ndof; ii++)
      {
        int j1d = tet2tensor[ii][2] * 2 + tet2tensor[ii][3];
        tensor2segm[j1d] = 0;
      }
    
    int ndsegm = 0;
    for (int i = 0; i < nd1d; i++)
      if (tensor2segm[i] != -1)
        ndsegm++;

    segm2tensor = FlatArray<int[2]> (ndsegm, lh);

    ndsegm = 0;
    for (int i3 = 0, i1d = 0; i3 < horder; i3++)
      for (int i4 = 0; i4 < 2; i4++, i1d++)
        if (tensor2segm[i1d] != -1)
          {
            tensor2segm[i1d] = ndsegm;
            segm2tensor[tensor2segm[i1d]][0] = i3;
            segm2tensor[tensor2segm[i1d]][1] = i4;
            ndsegm++;
          }    

    FlatArray<int> trig2segm (ndtrig, lh);
    for (int i = 0; i < ndtrig; i++)
      trig2segm[i] = tensor2segm[trig2tensor[i][1] * 2 + trig2tensor[i][2]];

    
    map2dto1d = trig2segm;
    map3dto2d = tet2trig;

    ndof2d = ndtrig;
    ndof1d = ndsegm;

    /*    
    *testout << "tet2tensor = " << endl;
    for (int i = 0; i < ndof; i++)
      *testout << tet2tensor[i][0] << ", " << tet2tensor[i][1] << ", " << tet2tensor[i][2] << endl;

    *testout << "map3dto2d = " << endl << map3dto2d << endl;
    */

    FlatArray<int> cnt(ndof2d, lh);

    pmap1dto2d = FlatArray<int> (ndof1d+1, lh);
    map1dto2d = FlatArray<int> (ndof2d, lh);
    cnt = 0;
    for (int j = 0; j < ndof2d; j++)
      cnt[map2dto1d[j]]++;
    pmap1dto2d[0] = 0;
    for (int j = 0; j < ndof1d; j++)
      pmap1dto2d[j+1] = pmap1dto2d[j] + cnt[j];
    cnt = 0;
    for (int j = 0; j < ndof2d; j++)
      {
        int j1d = map2dto1d[j];
        map1dto2d[pmap1dto2d[j1d]+cnt[j1d]++] = j;
      }
    
    pmap2dto3d = FlatArray<int> (ndof2d+1, lh);
    map2dto3d = FlatArray<int> (ndof, lh);
    cnt = 0;
    for (int j = 0; j < ndof; j++)
      cnt[map3dto2d[j]]++;
    pmap2dto3d[0] = 0;
    for (int j = 0; j < ndof2d; j++)
      pmap2dto3d[j+1] = pmap2dto3d[j] + cnt[j];
    cnt = 0;
    for (int j = 0; j < ndof; j++)
      {
        int j2d = map3dto2d[j];
        map2dto3d[pmap2dto3d[j2d]+cnt[j2d]++] = j;
      }
  }





  void L2HighOrderTetTP :: CalcShape (const IntegrationPoint & ip, 
                                      FlatVector<> shape) const
  {
    T_CalcShape (ip(0), ip(1), ip(2), shape); 
  }
  
  void L2HighOrderTetTP :: CalcDShape (const IntegrationPoint & ip, 
                                       FlatMatrix<> dshape) const
  {
    AutoDiff<3> x(ip(0), 0);
    AutoDiff<3> y(ip(1), 1);
    AutoDiff<3> z(ip(2), 2);

    ArrayMem<AutoDiff<3>,40> sds(ndof);
    T_CalcShape (x, y, z, sds);

    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 3; j++)
        dshape(i, j) = sds[i].DValue(j);
  }


  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void  L2HighOrderTetTP :: T_CalcShape (Tx x, Ty y, Tz z, TFA & sds) const
  {
    Tx lami[4] = { x, y, z, 1-x-y-z };

    Tx lamis[4];
    for (int i = 0; i < 4; i++)
      lamis[i] = lami[sort[i]];

    ArrayMem<Tx, 20> memx(sqr(order+1));
    ArrayMem<Ty, 20> memy(sqr(order+1));

    FlatMatrix<Tx> polsx(order+1, &memx[0]);
    FlatMatrix<Ty> polsy(order+1, &memy[0]);
    VectorMem<10, Tz> polsz(order+1);
    
    for (int i = 0; i <= order; i++)
      JacobiPolynomial (order, 2*lamis[0]-1, 2*i+2, 0, polsx.Row(i));
    for (int i = 0; i <= order; i++)
      ScaledJacobiPolynomial (order, lamis[1]-lamis[2]-lamis[3], 1-lamis[0], 2*i+1, 0, polsy.Row(i));

    ScaledLegendrePolynomial (order, lamis[2]-lamis[3], lamis[2]+lamis[3], polsz);

    for (int i = 0; i < ndof; i++)
      sds[i] =
        polsx(tet2tensor[i][1]+tet2tensor[i][2]+tet2tensor[i][3], tet2tensor[i][0]) 
        * polsy(tet2tensor[i][2]+tet2tensor[i][3], tet2tensor[i][1])
        * polsz(tet2tensor[i][2]) ;
  }


  void L2HighOrderTetTP :: CalcTPFactors (const IntegrationRuleTP<3> & irtp,
                                          TPFactors<3> & factors) const
  {
    HeapReset hr(lh);

    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
    const IntegrationRule & irz = irtp.GetIRZ();
    
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();
    int nipz = irz.GetNIP();

    factors.fac_x.SetSize (ndof, nipx);
    factors.fac_y.SetSize (ndof2d, nipy);
    factors.fac_z.SetSize (ndof1d, nipz);

    for (int i1 = 0; i1 < nipz; i1++)
      CalcZFactor (irz[i1](0), factors.fac_z.Col(i1), lh);
    
    for (int i1 = 0; i1 < nipy; i1++)
      CalcYFactor (iry[i1](0), factors.fac_y.Col(i1), lh);
    
    for (int i1 = 0; i1 < nipx; i1++)
      CalcXFactor (irx[i1](0), factors.fac_x.Col(i1), lh);
  }

  void L2HighOrderTetTP :: SetTPFactors (const TPFactors<3> & factors)
  {
    fac_x.AssignMemory (factors.fac_x.Height(),
                        factors.fac_x.Width(),
                        &factors.fac_x(0,0));
                          
    fac_y.AssignMemory (factors.fac_y.Height(),
                        factors.fac_y.Width(),
                        &factors.fac_y(0,0));

    fac_z.AssignMemory (factors.fac_z.Height(),
                        factors.fac_z.Width(),
                        &factors.fac_z(0,0));
  }




  void L2HighOrderTetTP ::
  EvaluateShapeGrid (const IntegrationRuleTP<3> & irtp,
                     const FlatVector<double> coefs,
                     FlatVector<double> gridvalues,
                     LocalHeap & lh) const
  {
    HeapReset hr (lh);

    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
    const IntegrationRule & irz = irtp.GetIRZ();
      
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();
    int nipz = irz.GetNIP();

    int nipxy = nipx*nipy;
    int nip = nipxy * nipz;

    /*
      FlatMatrix<> fac_x(ndof,   nipx, lh); 
      FlatMatrix<> fac_y(ndof2d, nipy, lh); 
      FlatMatrix<> fac_z(ndof1d, nipz, lh); 

      for (int i1 = 0; i1 < nipz; i1++)
      CalcZFactor (irz[i1](0), fac_z.Col(i1), lh);

      for (int i1 = 0; i1 < nipy; i1++)
      CalcYFactor (iry[i1](0), fac_y.Col(i1), lh);

      for (int i1 = 0; i1 < nipx; i1++)
      CalcXFactor (irx[i1](0), fac_x.Col(i1), lh);
    */

    if (fac_x.Height() != ndof)
      {
        FlatMatrix<> hfac_x(ndof,   nipx, new double[ndof*nipx]); 
        FlatMatrix<> hfac_y(ndof2d, nipy, new double[ndof2d*nipy]); 
        FlatMatrix<> hfac_z(ndof1d, nipz, new double[ndof1d*nipz]); 
        
        for (int i1 = 0; i1 < nipz; i1++)
          CalcZFactor (irz[i1](0), hfac_z.Col(i1), lh);

        for (int i1 = 0; i1 < nipy; i1++)
          CalcYFactor (iry[i1](0), hfac_y.Col(i1), lh);
        
        for (int i1 = 0; i1 < nipx; i1++)
          CalcXFactor (irx[i1](0), hfac_x.Col(i1), lh);

        
        const_cast<FlatMatrix<> &> (fac_x).AssignMemory (ndof  , nipx, &hfac_x(0,0));
        const_cast<FlatMatrix<> &> (fac_y).AssignMemory (ndof2d, nipy, &hfac_y(0,0));
        const_cast<FlatMatrix<> &> (fac_z).AssignMemory (ndof1d, nipz, &hfac_z(0,0));
      }

    FlatMatrix<> mgrid(nipxy, nipz, &gridvalues(0));
    FlatMatrix<> gridx(ndof2d, nipx, lh);
    FlatMatrix<> gridxy(ndof1d*nipx, nipy, lh);

    gridx = 0.0;
    for (int i = 0; i < ndof; i++)
      gridx.Row(map3dto2d[i]) += coefs(i) * fac_x.Row(i);
    
    gridxy = 0.0;
    for (int i = 0; i < ndof2d; i++)
      {
        int i1d = map2dto1d[i]*nipx;
        for (int ix = 0; ix < nipx; ix++)
          gridxy.Row(i1d+ix) += gridx(i, ix) * fac_y.Row(i);
      }

    mgrid = 0.0;
    for (int ixy = 0; ixy < nipxy; ixy++)
      {
        FlatVector<> row = mgrid.Row(ixy);
        for (int j = 0; j < ndof1d; j++)
          row += gridxy(j*nipxy+ixy) * fac_z.Row(j);
      }
  }

  void L2HighOrderTetTP ::
  EvaluateShapeGridTrans (const IntegrationRuleTP<3> & irtp,
                          const FlatVector<double> gridvalues,
                          FlatVector<double> coefs,
                          LocalHeap & lh) const
  {
    static int timer = NgProfiler::CreateTimer ("tet:shapegridtrans");
    NgProfiler::StartTimer (timer);

    HeapReset hr (lh);
    
    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
    const IntegrationRule & irz = irtp.GetIRZ();
      
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();
    int nipz = irz.GetNIP();
    int nipxy = nipx*nipy;
    int nip = nipxy * nipz;

    /*
      FlatMatrix<> fac_x(ndof,   nipx, lh); 
      FlatMatrix<> fac_y(ndof2d, nipy, lh); 
      FlatMatrix<> fac_z(ndof1d, nipz, lh); 

      for (int i1 = 0; i1 < nipz; i1++)
      CalcZFactor (irz[i1](0), fac_z.Col(i1), lh);
    
      for (int i1 = 0; i1 < nipy; i1++)
      CalcYFactor (iry[i1](0), fac_y.Col(i1), lh);
    
      for (int i1 = 0; i1 < nipx; i1++)
      CalcXFactor (irx[i1](0), fac_x.Col(i1), lh);
    */

    if (fac_x.Height() != ndof || fac_x.Width() != nipx ||
        fac_y.Width() != nipy || fac_z.Width() != nipz)
      {
        FlatMatrix<> hfac_x(ndof,   nipx, new double[ndof*nipx]); 
        FlatMatrix<> hfac_y(ndof2d, nipy, new double[ndof2d*nipy]); 
        FlatMatrix<> hfac_z(ndof1d, nipz, new double[ndof1d*nipz]); 
        
        for (int i1 = 0; i1 < nipz; i1++)
          CalcZFactor (irz[i1](0), hfac_z.Col(i1), lh);

        for (int i1 = 0; i1 < nipy; i1++)
          CalcYFactor (iry[i1](0), hfac_y.Col(i1), lh);
        
        for (int i1 = 0; i1 < nipx; i1++)
          CalcXFactor (irx[i1](0), hfac_x.Col(i1), lh);

        
        const_cast<FlatMatrix<> &> (fac_x).AssignMemory (ndof  , nipx, &hfac_x(0,0));
        const_cast<FlatMatrix<> &> (fac_y).AssignMemory (ndof2d, nipy, &hfac_y(0,0));
        const_cast<FlatMatrix<> &> (fac_z).AssignMemory (ndof1d, nipz, &hfac_z(0,0));
      }


    FlatMatrix<> mgrid(nipxy, nipz, &gridvalues(0));
    FlatMatrix<> gridx(ndof2d, nipx, lh);
    FlatMatrix<> gridxy(ndof1d*nipx, nipy, lh);

    for (int j = 0; j < ndof1d; j++) 
      for (int ixy = 0; ixy < nipxy; ixy++)
        gridxy(j*nipxy+ixy) = InnerProduct (fac_z.Row(j), mgrid.Row(ixy));

    for (int i = 0; i < ndof2d; i++)
      for (int ix = 0; ix < nipx; ix++)
        gridx(i, ix) = InnerProduct (fac_y.Row(i), gridxy.Row(map2dto1d[i]*nipx+ix));

    for (int i = 0; i < ndof; i++)
      coefs(i) = InnerProduct (fac_x.Row(i), gridx.Row(map3dto2d[i]));

    NgProfiler::AddFlops (timer, ndof1d*nipx*nipy*nipz);
    NgProfiler::AddFlops (timer, ndof2d*nipx*nipy);
    NgProfiler::AddFlops (timer, ndof*nipx);

    NgProfiler::StopTimer (timer);
  }







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

  void L2HighOrderTetTP ::
  EvaluateDShapeGrid (const IntegrationRuleTP<3> & irtp,
                      const FlatVector<double> coefs,
                      FlatMatrixFixWidth<3> gridvalues,
                      LocalHeap & lh) const
  {
    HeapReset hr (lh);

    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
    const IntegrationRule & irz = irtp.GetIRZ();
      
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();
    int nipz = irz.GetNIP();
    int nip = nipx * nipy * nipz;

    FlatMatrix<AutoDiff<1> > fac_xdx(ndof,   nipx, lh);  
    FlatMatrix<AutoDiff<1> > fac_ydy(ndof2d, nipy, lh);
    FlatMatrix<AutoDiff<1> > fac_zdz(ndof1d, nipz, lh);


    for (int i1 = 0; i1 < nipz; i1++)
      {
        AutoDiff<1> z (irz[i1](0), 0);
        CalcZFactor (z, fac_zdz.Col(i1), lh);
      }
    
    for (int i1 = 0; i1 < nipy; i1++)
      {
        AutoDiff<1> y (iry[i1](0), 0);
        CalcYFactor (y, fac_ydy.Col(i1), lh);
      }
    
    for (int i1 = 0; i1 < nipx; i1++)
      {
        AutoDiff<1> x (irx[i1](0), 0);
        CalcXFactor (x, fac_xdx.Col(i1), lh);
      }
    
    FlatMatrix<Vec<3> > mgrid(nipx*nipy, nipz, lh);
    FlatMatrix<Vec<2> > gridx(ndof2d, nipx, lh);
    FlatMatrix<Vec<3> > gridxy(ndof1d*nipx, nipy, lh);
    
    gridx = 0.0;
    for (int i = 0; i < ndof; i++)
      {
        int i2d = map3dto2d[i];
        double coefi = coefs(i);
        for (int ix = 0; ix < nipx; ix++)
          TrafoGradientX(fac_xdx(i, ix)) . TransformAddT (coefi, gridx(i2d, ix));
      }
    
    gridxy = 0.0;
    for (int i = 0; i < ndof2d; i++)
      {
        int i1d = map2dto1d[i];
        for (int ix = 0; ix < nipx; ix++)
          for (int iy = 0; iy < nipy; iy++)
            TrafoGradientY(fac_ydy(i, iy)) . TransformAddT (gridx(i, ix), gridxy(i1d*nipx+ix, iy));
      }
    
    mgrid = 0.0;
    for (int ix = 0; ix < nipx; ix++)
      for (int iy = 0; iy < nipy; iy++)
        for (int j = 0; j < ndof1d; j++)
          {
            Vec<3> hv = gridxy(j*nipx+ix,iy);
            for (int iz = 0; iz < nipz; iz++)
              TrafoGradientZ(fac_zdz(j, iz)) . TransformAddT (hv, mgrid(ix*nipy+iy, iz));
          }
    
    for (int ix = 0, ii = 0; ix < nipx; ix++)
      for (int iy = 0; iy < nipy; iy++)
        for (int iz = 0; iz < nipz; iz++, ii++)
          gridvalues.Row(ii) = Trans (irtp.GetDuffyJacobian(ii)) * mgrid(ix*nipy+iy, iz);
  }
  
  




  void L2HighOrderTetTP ::
  EvaluateDShapeGridTrans (const IntegrationRuleTP<3> & irtp,
                           const FlatMatrixFixWidth<3> gridvalues,
                           FlatVector<double> coefs,
                           LocalHeap & lh) const
  {
    HeapReset hr (lh);

    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
    const IntegrationRule & irz = irtp.GetIRZ();
      
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();
    int nipz = irz.GetNIP();
    int nipxy = nipx*nipy;
    int nip = nipx * nipy * nipz;

    FlatMatrix<AutoDiff<1> > fac_xdx(ndof,   nipx, lh);  
    FlatMatrix<AutoDiff<1> > fac_ydy(ndof2d, nipy, lh);
    FlatMatrix<AutoDiff<1> > fac_zdz(ndof1d, nipz, lh);

    for (int i1 = 0; i1 < nipz; i1++)
      {
        AutoDiff<1> z (irz[i1](0), 0);
        CalcZFactor (z, fac_zdz.Col(i1), lh);
      }
    
    for (int i1 = 0; i1 < nipy; i1++)
      {
        AutoDiff<1> y (iry[i1](0), 0);
        CalcYFactor (y, fac_ydy.Col(i1), lh);
      }
    
    for (int i1 = 0; i1 < nipx; i1++)
      {
        AutoDiff<1> x (irx[i1](0), 0);
        CalcXFactor (x, fac_xdx.Col(i1), lh);
      }

    FlatMatrix<Vec<3> > mgrid(nipxy, nipz, lh);
    FlatMatrix<Vec<2> > gridx(ndof2d, nipx, lh);
    FlatMatrix<Vec<3> > gridxy(ndof1d*nipx, nipy, lh);

    for (int ii = 0; ii < nip; ii++)
      mgrid(ii) = irtp.GetDuffyJacobian(ii) * gridvalues.Row(ii);

    for (int ixy = 0; ixy < nipxy; ixy++)
      for (int j = 0; j < ndof1d; j++)
        {
          Vec<3> sum(0.0);
          for (int iz = 0; iz < nipz; iz++)
            TrafoGradientZ(fac_zdz(j, iz)) . TransformAdd (mgrid(ixy, iz), sum);
          gridxy(j*nipxy+ixy) = sum;
        }
    
    
    for (int i = 0; i < ndof2d; i++)
      {
        int i1d = map2dto1d[i];
        for (int ix = 0; ix < nipx; ix++)
          {
            Vec<2> sum(0.0);
            for (int iy = 0; iy < nipy; iy++)
              TrafoGradientY(fac_ydy(i, iy)) . TransformAdd (gridxy(i1d*nipx+ix, iy), sum);
            gridx(i, ix) = sum;
          }
      }
    
    for (int i = 0; i < ndof; i++)
      {
        int i2d = map3dto2d[i];
        double sum = 0.0;
        for (int ix = 0; ix < nipx; ix++)
          TrafoGradientX(fac_xdx(i, ix)) . TransformAdd (gridx(i2d, ix), sum);
        coefs(i) = sum;
      }
  }
}


#endif
