/*********************************************************************/
/* File:   h1hofetp.cpp                                              */
/* Author: Start                                                     */
/* Date:   1. June 2007                                              */
/*********************************************************************/

#define xxNOPROFILE
#define noDEBUG

#include <fem.hpp>
namespace ngfem
{
  using namespace ngfem;





  /// compute shape

  template <ELEMENT_TYPE ET>
  void BaseH1HighOrderTP<ET,1> ::CalcShape (const IntegrationPoint & ip, 
                                            FlatVector<> shape) const
  {
    static_cast<const H1HighOrderTP<ET>&> (*this) . T_CalcShape (ip(0), shape); 
  }
  

  template <ELEMENT_TYPE ET>
  void BaseH1HighOrderTP<ET,1> ::CalcDShape (const IntegrationPoint & ip, 
                                             FlatMatrix<> dshape) const
  {
    AutoDiff<D> x(ip(0), 0);

    ArrayMem<AutoDiff<D>,40> sds( BASE::ndof);

    static_cast<const H1HighOrderTP<ET>&> (*this) . T_CalcShape (x, sds); 

    for (int i = 0; i < BASE::ndof; i++)
      for (int j = 0; j < D; j++)
	dshape(i, j) = sds[i].DValue(j);
  }


  template <ELEMENT_TYPE ET>
  void BaseH1HighOrderTP<ET,2> ::CalcShape (const IntegrationPoint & ip, 
                                            FlatVector<> shape) const
  {
    static_cast<const H1HighOrderTP<ET>&> (*this) . T_CalcShape (ip(0), ip(1), shape); 
  }


  template <ELEMENT_TYPE ET>
  void BaseH1HighOrderTP<ET,2> ::CalcDShape (const IntegrationPoint & ip, 
                                             FlatMatrix<> dshape) const
  {
    AutoDiff<D> x(ip(0), 0);
    AutoDiff<D> y(ip(1), 1);

    ArrayMem<AutoDiff<D>,40> sds( BASE::ndof);

    static_cast<const H1HighOrderTP<ET>&> (*this) . T_CalcShape (x, y, sds); 

    for (int i = 0; i < BASE::ndof; i++)
      for (int j = 0; j < D; j++)
	dshape(i, j) = sds[i].DValue(j);
  }




  template <ELEMENT_TYPE ET>
  void BaseH1HighOrderTP<ET,2> :: 
  EvaluateShapeGridTrans (const IntegrationRuleTP<2> & irtp,
                          const FlatVector<double> gridvalues,
                          FlatVector<double> coefs,
                          LocalHeap & lh) const
  {
    int ndof = BASE::GetNDof();

    HeapReset hr (lh);
    
    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
      
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();
    int nip = nipx * nipy;

    /*
      FlatVector<> fac_x(ndof, lh); 
      FlatVector<> fac_y(ndof1d, lh); 
      FlatMatrix<> gridx(nipx, ndof1d, lh);
    
      gridx = 0;
      for (int iy = 0; iy < nipy; iy++)
      {
      Spec().CalcYFactor (iry[iy](0), fac_y, lh);
      for (int ix = 0; ix < nipx; ix++)
      gridx.Row(ix) += gridvalues(ix*nipy+iy) * fac_y;
      }

      coefs = 0;
      for (int ix = 0; ix < nipx; ix++)
      {
      Spec().CalcXFactor (irx[ix](0), fac_x, lh);
      for (int i = 0; i < ndof; i++)
      coefs(i) += fac_x(i) * gridx(ix, map2dto1d[i]);
      }
    */


    FlatMatrix<> fac_x(ndof,   nipx, lh); 
    FlatMatrix<> fac_y(ndof1d, nipy, lh); 
      
    for (int i1 = 0; i1 < nipy; i1++)
      Spec().CalcYFactor (iry[i1](0), fac_y.Col(i1), lh);

    for (int i1 = 0; i1 < nipx; i1++)
      Spec().CalcXFactor (irx[i1](0), fac_x.Col(i1), lh);

    FlatMatrix<> mgrid(nipx, nipy, &gridvalues(0));
    FlatMatrix<> gridx(ndof1d, nipx, lh);

    for (int ix = 0; ix < nipx; ix++)
      for (int j = 0; j < ndof1d; j++)
        gridx(j, ix) = InnerProduct (mgrid.Row(ix), fac_y.Row(j));

    for (int i = 0; i < ndof; i++)
      coefs(i) = InnerProduct (fac_x.Row(i), gridx.Row(map2dto1d[i]));

    
    /*
      FlatMatrix<> fac_x(ndof,   nipx, lh); 
      FlatMatrix<> fac_y(ndof1d, nipy, lh); 
      
      for (int i1 = 0; i1 < nipy; i1++)
      Spec().CalcYFactor (iry[i1](0), fac_y.Col(i1), lh);

      for (int i1 = 0; i1 < nipx; i1++)
      Spec().CalcXFactor (irx[i1](0), fac_x.Col(i1), lh);

      FlatVector<> sum1d(ndof1d, lh);
      
      coefs = 0;
      for (int i1 = 0, ii = 0; i1 < nipx; i1++)
      {
      sum1d = 0;
        
      for (int i2 = 0; i2 < nipy; i2++, ii++)
      {
      for (int i = 0; i < ndof1d; i++)
      sum1d(i) += gridvalues(ii) * fac_y(i, i2);
      }
        
      for (int i = 0; i < ndof; i++)
      coefs(i) += sum1d(map2dto1d[i]) * fac_x(i, i1);
      }
    */
  }












  template <ELEMENT_TYPE ET>
  void BaseH1HighOrderTP<ET,3> ::CalcShape (const IntegrationPoint & ip, 
                                            FlatVector<> shape) const
  {
    static_cast<const H1HighOrderTP<ET>&> (*this) . T_CalcShape (ip(0), ip(1), ip(2), shape); 
  }


  template <ELEMENT_TYPE ET>
  void BaseH1HighOrderTP<ET,3> ::CalcDShape (const IntegrationPoint & ip, 
                                             FlatMatrix<> dshape) const
  {
    AutoDiff<D> x(ip(0), 0);
    AutoDiff<D> y(ip(1), 1);
    AutoDiff<D> z(ip(2), 2);

    ArrayMem<AutoDiff<D>,40> sds( BASE::ndof);
    
    static_cast<const H1HighOrderTP<ET>&> (*this) . T_CalcShape (x, y, z, sds); 

    for (int i = 0; i < BASE::ndof; i++)
      for (int j = 0; j < D; j++)
	dshape(i, j) = sds[i].DValue(j);
  }










  template <ELEMENT_TYPE ET>
  void BaseH1HighOrderTP<ET,3> ::
  EvaluateShapeGrid (const IntegrationRuleTP<3> & ir,
                     const FlatVector<double> coefs,
                     FlatVector<double> gridvalues,
                     LocalHeap & lh) const
  {
    cout << "h1highordertet, evaluategrid not implemented" << endl;
  }


				  
  template <ELEMENT_TYPE ET>
  void BaseH1HighOrderTP<ET,3> ::
  EvaluateShapeGridTrans (const IntegrationRuleTP<3> & irtp,
                          const FlatVector<double> gridvalues,
                          FlatVector<double> coefs,
                          LocalHeap & lh) const
  {
    if (!BASE::IsTPElement())
      {
        BASE::EvaluateShapeGridTrans (irtp, gridvalues, coefs, lh);
        return;
      }

    int ndof = BASE::GetNDof();
    

    HeapReset hr (lh);
    
    const IntegrationRule & irx = irtp.GetIRX();
    const IntegrationRule & iry = irtp.GetIRY();
    const IntegrationRule & irz = irtp.GetIRZ();
      
    int nipx = irx.GetNIP();
    int nipy = iry.GetNIP();
    int nipz = irz.GetNIP();
    int nip = nipx * nipy * nipz;
      
    FlatMatrix<> fac_x(ndof,   nipx, lh); 
    FlatMatrix<> fac_y(ndof2d, nipy, lh); 
    FlatMatrix<> fac_z(ndof1d, nipz, lh); 
      
    for (int i1 = 0; i1 < nipz; i1++)
      Spec().CalcZFactor (irz[i1](0), fac_z.Col(i1), lh);

    for (int i1 = 0; i1 < nipy; i1++)
      Spec().CalcYFactor (iry[i1](0), fac_y.Col(i1), lh);

    for (int i1 = 0; i1 < nipx; i1++)
      Spec().CalcXFactor (irx[i1](0), fac_x.Col(i1), lh);

    FlatVector<> sum2d(ndof2d, lh);
    FlatVector<> sum1d(ndof1d, lh);
      
    coefs = 0;
    for (int i1 = 0, ii = 0; i1 < nipx; i1++)
      {
        sum2d = 0;
        
        for (int i2 = 0; i2 < nipy; i2++)
          {
            sum1d = 0;
            for (int i3 = 0; i3 < nipz; i3++, ii++)
              {
                for (int i = 0; i < ndof1d; i++)
                  sum1d(i) += gridvalues(ii) * fac_z(i, i3);
              }
            
            for (int i = 0; i < ndof2d; i++)
              sum2d(i) += sum1d(map2dto1d[i]) * fac_y(i, i2);
          }
        
        for (int i = 0; i < ndof; i++)
          coefs(i) += sum2d(map3dto2d[i]) * fac_x(i, i1);
      }
  }



  template class BaseH1HighOrderTP<ET_SEGM>;
  template class BaseH1HighOrderTP<ET_TRIG>;
  template class BaseH1HighOrderTP<ET_QUAD>;
  template class BaseH1HighOrderTP<ET_TET>;
  template class BaseH1HighOrderTP<ET_PRISM>;
  template class BaseH1HighOrderTP<ET_HEX>;





  // ************************** Segm ************************************ 

  void H1HighOrderTP<ET_SEGM> :: SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    H1HighOrderSegm<> :: SetVertexNumbers (avnums, lh);

    for (int i = 0; i < 2; i++) sort[i] = i;

    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);

    segm2tensor = FlatArray<int[2]> (ndof, lh);

    int isort[2];
    for (int i = 0; i < 2; i++) isort[sort[i]] = i;
    
    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 2; j++)
	segm2tensor[i][j] = 0;
    
    // vertex shapes
    for (int i = 0; i < 2; i++)
      segm2tensor[i][isort[i]] = 1;
    
    int ii = 2;
    
    //  edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_SEGM);
    if (order_inner[0] >= 2)
      {
        int es = edges[0][0], ee = edges[0][1];
        if (vnums[es] > vnums[ee]) swap (es, ee);
        int nde = order_inner[0]-1;
        for (int j = 1; j <= nde; j++, ii++)
          {
            segm2tensor[ii][isort[es]] = j; // lami_s * ''pol(lami_s-lami_e,lami_s+lami_e)'' 
            segm2tensor[ii][isort[ee]] = 1; // lami_e
          }
      }
  }

  
  template<typename Tx, typename TFA>  
  void  H1HighOrderTP<ET_SEGM> :: T_CalcShape (Tx x, TFA & sds) const
  {
    Tx lami[2] = { x, 1-x };
    
    Tx lamis[2];
    for (int i = 0; i < 2; i++)
      lamis[i] = lami[sort[i]];

    int horder = max(2, order);
    ArrayMem<Tx, 20> polx(horder), polxx(2); 

    LegendrePolynomialMult (horder-2, 2*lamis[0]-1, lamis[0], polx.Addr(1));
    polx[0] = 1.0;
    
    polxx[0] = 1.0;
    polxx[1] = lamis[1];


    for (int i = 0; i < ndof; i++)
      sds[i] =
        polx[segm2tensor[i][0]] * 
        polxx[segm2tensor[i][1]];
  }


  


  // ************************** TRIG ************************************ 


  // template <>
  void H1HighOrderTP<ET_TRIG> :: SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    H1HighOrderTrig<> :: SetVertexNumbers (avnums, lh);

    for (int i = 0; i < 3; i++) sort[i] = i;

    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
    if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);

    // vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]] < vnums[sort[3]]
    
    // build tables

    trig2tensor = FlatArray<int[3]> (ndof, lh);

    int isort[3];
    for (int i = 0; i < 3; i++) isort[sort[i]] = i;
    
    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 3; j++)
	trig2tensor[i][j] = 0;
    
    // vertex shapes
    for (int i = 0; i < 3; i++)
      trig2tensor[i][isort[i]] = 1;
    
    int ii = 3;
    
    //  edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
    for (int i = 0; i < 3; i++)
      if (order_edge[i] >= 2)
	{
	  int es = edges[i][0], ee = edges[i][1];
	  if (vnums[es] > vnums[ee]) swap (es, ee);
	  int nde = order_edge[i]-1;
	  for (int j = 1; j <= nde; j++, ii++)
	    {
	      trig2tensor[ii][isort[es]] = j; // lami_s * ''pol(lami_s-lami_e,lami_s+lami_e)'' 
	      trig2tensor[ii][isort[ee]] = 1; // lami_e
	    }
	}
    
    // face dofs
    if (order_inner[0] >= 3)
      {
	int fav[3] = { 0, 1, 2 };
	
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
	if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
	if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
          
	int pi = order_inner[0];
	for (int j2 = 0; j2 <= pi-3; j2++)
	  for (int i2 = 0; i2 <= pi-3-j2; i2++, ii++)
	    {
	      trig2tensor[ii][isort[fav[0]]] = j2+1;
	      trig2tensor[ii][isort[fav[1]]] = i2+1;
	      trig2tensor[ii][isort[fav[2]]] = 1;
	    }
      }


    int horder = max(2, order);
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


  template<typename Tx, typename Ty, typename TFA>  
  void  H1HighOrderTP<ET_TRIG> :: T_CalcShape (Tx x, Ty y, TFA & sds) const
  {
    Tx lami[3] = { x, y,  1-x-y };

    Tx lamis[3];
    for (int i = 0; i < 3; i++)
      lamis[i] = lami[sort[i]];

#ifdef JACOBI
    int horder = max(2, order);

    ArrayMem<Tx, 100> memx((order+1)*horder);
    FlatMatrix<Tx> polsx(order+1, horder, &memx[0]);

    VectorMem<100, Tx> polsy(horder); 
    Vec<2, Tx> polsyy;

    for (int i = 0; i <= order; i++)
      {
        JacobiPolynomialMult (horder-2, 2*lamis[0]-1, max(0,2*i-2), 0, 
                              lamis[0], 
                              polsx.Row(i).Addr(1));
        polsx(i,0) = 1.0;
      }

    ScaledLegendrePolynomialMult (horder-2, lamis[1]-lamis[2], lamis[1]+lamis[2], 
                                  lamis[1], 
                                  polsy.Addr(1));
    polsy(0) = 1.0;

    polsyy(0) = 1.0;
    polsyy(1) = lamis[2];

    for (int i = 0; i < ndof; i++)
      sds[i] =
        polsx(trig2tensor[i][1]+trig2tensor[i][2], trig2tensor[i][0]) * 
        polsy(trig2tensor[i][1]) *
        polsyy(trig2tensor[i][2]);

#else
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

    for (int i = 0; i < ndof; i++)
      sds[i] =
        polx[trig2tensor[i][0]] * 
        poly[trig2tensor[i][1]] * 
        polyy[trig2tensor[i][2]];
#endif
  }




  // ******************************** QUAD **************************************

  void H1HighOrderTP<ET_QUAD> :: 
  SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    H1HighOrderQuad<> :: SetVertexNumbers (avnums, lh);

    quad2tensor = FlatArray<int[2]> (ndof, lh);
    factor = FlatArray<double> (ndof, lh);
    
    static int quad_points[][2] = 
      { 
        { 0, 0 },
        { 1, 0 },
        { 1, 1 },
        { 0, 1 }
      };


    // vertex shapes

    for (int i = 0; i < 4; i++)
      {
        quad2tensor[i][0] = quad_points[i][0];
        quad2tensor[i][1] = quad_points[i][1];
        factor[i] = 1;
      }
        
    int ii = 4;

    // edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_QUAD);
    for (int i = 0; i < 4; i++)
      if (order_edge[i] >= 2)
        {
          int p = order_edge[i];
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es, ee);

          double fac = 1, mult = 1;
          int dir = 0;

          for (int k = 0; k < 2; k++)
            if (quad_points[es][k] != quad_points[ee][k])
              {
                dir = k;
                if (quad_points[es][k] < quad_points[ee][k])
                  mult = -1;
              }

          for (int l = 0; l <= p-2; l++)
            {
              for (int k = 0; k < 2; k++)
                quad2tensor[ii][k] = quad_points[ee][k];
              quad2tensor[ii][dir] = l+2;

              factor[ii] = fac;
              fac *= mult;
              ii++;
            }
        }


     
    const FACE * faces = ElementTopology::GetFaces (ET_QUAD);
    for (int i = 0; i < 1; i++)
      if (order_inner[0] >= 2 && order_inner[1] >= 2)
        {
          int px = order_inner[0];
          int py = order_inner[1];
          int p = max2(px, py);

          int fmax = 0;
          for (int j=1; j<4; j++) 
            if (vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;  
	  
          int f1 = faces[i][(fmax+3)%4];
          int f2 = faces[i][(fmax+1)%4]; 
          fmax = faces[i][fmax]; 
	  
          if(vnums[f2]>vnums[f1]) swap(f1,f2);  // fmax > f1 > f2 
	  

          int ind[3];
          double fac1 = 1, fac2 = 1, mult1 = 1, mult2 = 1;
          for (int l = 0; l <= px-2; l++)
            {
              fac2 = 1;
              for (int m = 0; m <= py-2; m++)
                {
                  for (int k = 0; k < 2; k++)
                    {
                      ind[k] = quad_points[fmax][k];
                      if (quad_points[fmax][k] != quad_points[f2][k])
                        {
                          ind[k] = m+2;
                          if (quad_points[fmax][k] > quad_points[f2][k])
                            mult2 = -1;
                        }
                      if (quad_points[fmax][k] != quad_points[f1][k])
                        {
                          ind[k] = l+2;
                          if (quad_points[fmax][k] > quad_points[f1][k])
                            mult1 = -1;
                        }
                    }

                  factor[ii] = fac1*fac2;
                  quad2tensor[ii][0] = ind[0];
                  quad2tensor[ii][1] = ind[1];
                  ii++;

                  fac2 *= mult2;
                }
              fac1 *= mult1;
            }
          
        }
    ndof1d = order+1;
    int ndof2d = ndof;

    map2dto1d = FlatArray<int> (ndof, lh);
    for (int i = 0; i < ndof; i++)
      map2dto1d[i] = quad2tensor[i][1];


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
  }




  template<typename Tx, typename Ty, typename TFA>  
  void  H1HighOrderTP<ET_QUAD> :: T_CalcShape (Tx x, Ty y, TFA & sds) const
  {
    Tx xi[2] = { 2*x-1, 2*y-1 };
    
    ArrayMem<Tx, 80> polxi( 2 * (order+1) );
    Tx * poli[2];
    for (int i = 0; i < 2; i++) poli[i] = &polxi[i*(order+1)];

    for (int i = 0; i < 2; i++)
      {
        Tx * hp = poli[i]+2;
        poli[i][0] = 0.5 * (1-xi[i]);
        poli[i][1] = 0.5 * (1+xi[i]);
        LegendrePolynomialMult (order-2, xi[i], 0.25 * (1-xi[i]*xi[i]), hp);
      }

    // ArrayMem<int[2], 100> quad2tensor(ndof);
    // ArrayMem<double, 100> fac(ndof);

    // GetQuad2TensorMapping (quad2tensor, fac);

    for (int i = 0; i < ndof; i++)
      sds[i] = factor[i] * 
        poli[0][quad2tensor[i][0]] * 
        poli[1][quad2tensor[i][1]];
  }







  // ******************************** TET **************************************

  void H1HighOrderTP<ET_TET> :: SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    H1HighOrderTet<> :: SetVertexNumbers (avnums, lh);

    for (int i = 0; i < 4; i++) sort[i] = i;

    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
    if (vnums[sort[2]] > vnums[sort[3]]) Swap (sort[2], sort[3]);
    if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
    if (vnums[sort[1]] > vnums[sort[3]]) Swap (sort[1], sort[3]);
    if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);

    // vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]] < vnums[sort[3]]


    // build tables

    tet2tensor = FlatArray<int[4]> (ndof, lh);

    int isort[4];
    for (int i = 0; i < 4; i++)
      isort[sort[i]] = i;

    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 4; j++)
        tet2tensor[i][j] = 0;

    // vertex shapes
    for (int i = 0; i < 4; i++)
      tet2tensor[i][isort[i]] = 1;
    
    int ii = 4;
    
    //  edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_TET);
    for (int i = 0; i < 6; i++)
      if (order_edge[i] >= 2)
        {
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es, ee);
          int nde = order_edge[i]-1;
          for (int j = 1; j <= nde; j++, ii++)
            {
              tet2tensor[ii][isort[es]] = j;
              tet2tensor[ii][isort[ee]] = 1;
            }
        }
    
    // face dofs
    const FACE * faces = ElementTopology::GetFaces (ET_TET);
    for (int i = 0; i < 4; i++)
      if (order_face[i][0] >= 3)
        {
          int fav[3] = { faces[i][0], faces[i][1], faces[i][2] };
          
          if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
          if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
          if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
          
          int n = order_face[i][0];
          for (int j2 = 0; j2 <= n-3; j2++)
            for (int i2 = 0; i2 <= n-3-j2; i2++, ii++)
              {
                tet2tensor[ii][isort[fav[0]]] = j2+1;
                tet2tensor[ii][isort[fav[1]]] = i2+1;
                tet2tensor[ii][isort[fav[2]]] = 1;
              }
        }
    
    // inner dofs
    if (order_inner[0] >= 4)
      {
        int n = order_inner[0];
        for (int i = 0; i <= n-4; i++)
          for (int j = 0; j <= n-4-i; j++)
            for (int k = 0; k <= n-4-i-j; k++, ii++)
              {
                tet2tensor[ii][0] = k+1;
                tet2tensor[ii][1] = j+1;
                tet2tensor[ii][2] = i+1;
                tet2tensor[ii][3] = 1;
              }
      }





    int horder = max(2, order);

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
    


  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void  H1HighOrderTP<ET_TET> :: T_CalcShape (Tx x, Ty y, Tz z, TFA & sds) const
  {
    Tx lami[4] = { x, y, z, 1-x-y-z };

    Tx lamis[4];
    for (int i = 0; i < 4; i++)
      lamis[i] = lami[sort[i]];

#ifdef JACOBI
    int horder = max(2, order);

    ArrayMem<Tx, 100> memx((order+1)*horder);
    ArrayMem<Tx, 100> memy((order+1)*horder);
    ArrayMem<Tx, 20> memz(2*horder);
    ArrayMem<Tx, 20> memzz(2);

    FlatMatrix<Tx> polsx(order+1, horder, &memx[0]);
    FlatMatrix<Tx> polsy(order+1, horder, &memy[0]);
    FlatMatrix<Tx> polsz(2, horder, &memz[0]);
    FlatMatrix<Tx> polszz(1, 2, &memzz[0]);

    for (int i = 0; i <= order; i++)
      {
        JacobiPolynomialMult (horder-2, 2*lamis[0]-1, 
                              max(0,2*i-2), 0, lamis[0], polsx.Row(i).Addr(1));
        polsx(i,0) = 1.0;
      }

    for (int i = 0; i <= order; i++)
      {
        ScaledJacobiPolynomialMult (horder-2, lamis[1]-lamis[2]-lamis[3], 1-lamis[0], 
                                    max(0,2*i-2), 0, lamis[1], polsy.Row(i).Addr(1));
        polsy(i, 0) = 1.0;
      }

    for (int i = 0; i < 2; i++)
      {
        ScaledLegendrePolynomialMult (horder-2, lamis[2]-lamis[3], lamis[2]+lamis[3], 
                                      lamis[2], polsz.Row(i).Addr(1));
        polsz(i,0) = 1.0;
      }

    polszz(0,0) = 1.0;
    polszz(0,1) = lamis[3];

    for (int i = 0; i < ndof; i++)
      sds[i] =
        polsx(tet2tensor[i][1]+tet2tensor[i][2]+tet2tensor[i][3], tet2tensor[i][0]) * 
        polsy(tet2tensor[i][2]+tet2tensor[i][3], tet2tensor[i][1]) * 
        polsz(0, tet2tensor[i][2]) *
        polszz(0, tet2tensor[i][3]);

#else
    int horder = max(2, order);

    ArrayMem<Tx, 20> polx(horder), poly(horder), polz(horder), polzz(2);

    LegendrePolynomialMult (horder-2, 2*lamis[0]-1, lamis[0], polx.Addr(1));
    polx[0] = 1.0;

    ScaledLegendrePolynomialMult (horder-2, lamis[1]-lamis[2]-lamis[3], 1-lamis[0], lamis[1], poly.Addr(1));
    poly[0] = 1.0;

    ScaledLegendrePolynomialMult (horder-2, lamis[2]-lamis[3], lamis[2]+lamis[3], lamis[2], polz.Addr(1));
    polz[0] = 1.0;

    polzz[0] = 1.0;
    polzz[1] = lamis[3];

    for (int i = 0; i < ndof; i++)
      sds[i] =
        polx[tet2tensor[i][0]] * 
        poly[tet2tensor[i][1]] * 
        polz[tet2tensor[i][2]] *
        polzz[tet2tensor[i][3]];
#endif
  }




  // ******************************** PRISM **************************************


  
  void H1HighOrderTP<ET_PRISM> :: SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    H1HighOrderPrism<> :: SetVertexNumbers (avnums, lh);

    for (int i = 0; i < 6; i++) sort[i] = i;

    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
    if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);

    if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);
    if (vnums[sort[4]] > vnums[sort[5]]) Swap (sort[4], sort[5]);
    if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);

    BASE::tp = sort[0] + 3 == sort[3];
    if (!BASE::tp) return;
    // if (sort[0]+3 != sort[3]) return;


    // build tables

    prism2tensor = FlatArray<int[4]> (ndof, lh);
    factor = FlatArray<double> (ndof, lh);

    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 4; j++)
        prism2tensor[i][j] = 0;
    for (int i = 0; i < ndof; i++)
      factor[i] = 0;


    int isort[6];
    for (int i = 0; i < 6; i++) isort[sort[i]] = i;

    // vertex shapes
    for (int i = 0; i < 6; i++)
      {
        prism2tensor[i][isort[i%3]] = 1;
        prism2tensor[i][3] = i / 3;
        factor[i] = 1;
      }    

    int ii = 6;
    
    const EDGE * edges = ElementTopology::GetEdges (ET_PRISM);
    // horizontal edge dofs
    for (int i = 0; i < 6; i++)
      if (order_edge[i] >= 2)
        {
          int p = order_edge[i];
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es%3] > vnums[ee%3]) swap (es, ee);
          int mult = ( (vnums[es] > vnums[ee]) == (vnums[es%3] > vnums[ee%3])) ? 1 : -1;
            
          int fac = 1;
          int nde = order_edge[i]-1;
          for (int j = 1; j <= nde; j++, ii++)
            {
              prism2tensor[ii][isort[es%3]] = j;
              prism2tensor[ii][isort[ee%3]] = 1;
              prism2tensor[ii][3] = es / 3;
              factor[ii] = fac;
              fac *= mult;
            }
        }
    
    // vertical edges
    for (int i = 6; i < 9; i++)
      if (order_edge[i] >= 2)
        {
          int p = order_edge[i];
          int es = edges[i][0], ee = edges[i][1];

          int mult = 1, fac = 1;
          if (vnums[es] < vnums[ee]) mult = -1;   
          
          int nde = order_edge[i]-1;
          for (int j = 1; j <= nde; j++, ii++)
            {
              prism2tensor[ii][isort[es%3]] = 1;
              prism2tensor[ii][3] = j+1;
              factor[ii] = fac;
              fac *= mult;
            }
        }

    
    const FACE * faces = ElementTopology::GetFaces (ET_PRISM);

    // trig faces
    for (int i = 0; i < 2; i++)
      if (order_face[i][0] >= 3)
        {
          int fav[3] = { faces[i][0], faces[i][1], faces[i][2] };

          if(vnums[fav[0]%3] > vnums[fav[1]%3]) swap(fav[0],fav[1]); 
          if(vnums[fav[1]%3] > vnums[fav[2]%3]) swap(fav[1],fav[2]);
          if(vnums[fav[0]%3] > vnums[fav[1]%3]) swap(fav[0],fav[1]); 	

          int mult = ( (vnums[fav[1]] > vnums[fav[2]]) == (vnums[fav[1]%3] > vnums[fav[2]%3])) ? 1 : -1;
          
          int n = order_face[i][0];

          for (int j2 = 0; j2 <= n-3; j2++)
            {
              int fac = 1;
              for (int i2 = 0; i2 <= n-3-j2; i2++, ii++)
                {
                  prism2tensor[ii][isort[fav[0]%3]] = j2+1;
                  prism2tensor[ii][isort[fav[1]%3]] = i2+1;
                  prism2tensor[ii][isort[fav[2]%3]] = 1;
                  prism2tensor[ii][3] = fav[0] / 3;                  
                  factor[ii] = fac;
                  fac *= mult;
                }
            }
        }


    // quad face dofs
    for (int i = 2; i < 5; i++)
      if (order_face[i][0] >= 2 && order_face[i][1] >= 2)
        {
          INT<2> p = order_face[i];
          int pp = int(max2(p[0], p[1]))+1;
          
          int fmax = 0;
          for(int j = 1; j < 4; j++) 
            if(vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;
          
          int fz = 3-fmax; 
          int ftrig = fmax^1; 

          fmax = faces[i][fmax];
          fz = faces[i][fz];
          ftrig = faces[i][ftrig];

          int multtrig = ( (vnums[fmax] > vnums[ftrig]) == (vnums[fmax%3] > vnums[ftrig%3])) ? 1 : -1;
          int multz = (fmax < 3) ? 1 : -1;

          int es = fmax%3, ee = ftrig%3;
          if (vnums[es] > vnums[ee]) swap (es, ee);

          // global x-direction is towards second-largest vertex (JS)
          if (vnums[ftrig] > vnums[fz])  
            {
              int factrig = 1;
              for (int k = 0; k < p[0]-1; k++)
                {
                  int facz = 1;
                  for (int j = 0; j < p[1]-1; j++, ii++)
                    {
                      prism2tensor[ii][isort[es]] = k+1;
                      prism2tensor[ii][isort[ee]] = 1;
                      prism2tensor[ii][3] = j+2;

                      factor[ii] = facz * factrig;
                      facz *= multz;
                    }
                  factrig *= multtrig;
                }

            }
          else
            {
              int facz = 1;
              for (int j = 0; j < p[0]-1; j++)
                {
                  int factrig = 1;
                  for (int k = 0; k < p[1]-1; k++, ii++)
                    {
                      prism2tensor[ii][isort[es]] = k+1;
                      prism2tensor[ii][isort[ee]] = 1;
                      prism2tensor[ii][3] = j+2;
                      factor[ii] = facz * factrig;
                      factrig *= multtrig;
                    }
                  facz *= multz;
                }
            }
        }

    // volume dofs:
    if (order_inner[0] > 2 && order_inner[2] > 1)
      {
        int pxy = order_inner[0];
        int pz = order_inner[2];

        for (int j = 0; j <= pxy-3; j++)
          for (int k = 0; k <= pxy-3-j; k++)
            for (int l = 0; l <= pz-2; l++, ii++)
              {
                prism2tensor[ii][0] = j+1;
                prism2tensor[ii][1] = k+1;
                prism2tensor[ii][2] = 1;
                prism2tensor[ii][3] = l+2;
                factor[ii] = 1;
              }
      }


    /*
      int horder = max(2, order);
      int nd1d = order+1;      // zorder
      int nd2d = (order+1) * (horder+1) * 2;  // poly, polyy, polz


      FlatArray<int> tensor2trig ( nd2d, lh);
      tensor2trig = -1;

      for (int ii = 0; ii < ndof; ii++)
      {
      int j2d = prism2tensor[ii][1] * (2*nd1d) + prism2tensor[ii][2] * nd1d + prism2tensor[ii][3];
      tensor2trig[j2d] = 0;
      }
      
      int ndof2d =0;
      for (int i = 0; i < nd2d; i++)
      if (tensor2trig[i] != -1)
      {
      tensor2trig[i] = ndof2d;
      ndof2d++;
      }

      // ArrayMem<int[3], 100> trig2tensor(ndof2d);
      trig2tensor = FlatArray<int[3]> (ndof2d, lh);
      for (int i2 = 0; i2 <= horder; i2++)
      for (int i3 = 0; i3 < 2; i3++)
      for (int i4 = 0; i4 <= order; i4++)
      {
      int i2d = i2 * 2*nd1d + i3 * nd1d + i4;
      if (tensor2trig[i2d] != -1)
      {
      trig2tensor[tensor2trig[i2d]][0] = i2;
      trig2tensor[tensor2trig[i2d]][1] = i3;
      trig2tensor[tensor2trig[i2d]][2] = i4;
      }
      }
    
      FlatArray<int> prism2trig (ndof, lh);
      for (int i = 0; i < ndof; i++)
      prism2trig[i] = tensor2trig[prism2tensor[i][1] * 2*nd1d + prism2tensor[i][2] * nd1d + prism2tensor[i][3]];
    

      FlatArray<int> tensor2segm ( nd1d, lh);
      for (int i = 0; i < nd1d; i++)
      tensor2segm[i] = i;
    
      // int ndof1d = order+1;
    



      FlatArray<int> trig2segm (ndof2d, lh);
      for (int i = 0; i < ndof2d; i++)
      trig2segm[i] = tensor2segm[trig2tensor[i][2]];
    
      map2dto1d = trig2segm;
      map3dto2d = prism2trig;
    */


    int horder = max(2, order);
    int nd1d = 2*horder;                     // zorder
    int nd2d = (order+1) * (horder+1) * 2;  // polx, poly, polyy


    FlatArray<int> tensor2trig ( nd2d, lh);
    tensor2trig = -1;

    for (int ii = 0; ii < ndof; ii++)
      {
        int j2d = prism2tensor[ii][0] * 2*(horder+1) + prism2tensor[ii][1] * 2 + prism2tensor[ii][2];
        tensor2trig[j2d] = 0;
      }
      
    ndof2d =0;
    for (int i = 0; i < nd2d; i++)
      if (tensor2trig[i] != -1)
        {
          tensor2trig[i] = ndof2d;
          ndof2d++;
        }

    // ArrayMem<int[3], 100> trig2tensor(ndof2d);
    trig2tensor = FlatArray<int[3]> (ndof2d, lh);
    for (int i1 = 0; i1 <= horder; i1++)
      for (int i2 = 0; i2 <= horder; i2++)
        for (int i3 = 0; i3 < 2; i3++)
          {
            int i2d = i1 * 2*(horder+1) + i2 * 2 + i3;
            if (tensor2trig[i2d] != -1)
              {
                trig2tensor[tensor2trig[i2d]][0] = i1;
                trig2tensor[tensor2trig[i2d]][1] = i2;
                trig2tensor[tensor2trig[i2d]][2] = i3;
              }
          }
    
    FlatArray<int> prism2trig (ndof, lh);
    for (int i = 0; i < ndof; i++)
      prism2trig[i] = tensor2trig[prism2tensor[i][0] * 2*(horder+1) + 
                                  prism2tensor[i][1] * 2 + prism2tensor[i][2]];



    FlatArray<int> tensor2segm ( nd1d, lh);
    tensor2segm = -1;
    
    for (int ii = 0; ii < ndof; ii++)
      {
        int j1d = prism2tensor[ii][1] * 2 + prism2tensor[ii][2];
        tensor2segm[j1d] = 0;
      }
    
    int ndsegm =0;
    for (int i = 0; i < nd1d; i++)
      if (tensor2segm[i] != -1)
        {
          tensor2segm[i] = ndsegm;
          ndsegm++;
        }
    
    FlatArray<int> trig2segm (ndof2d, lh);
    for (int i = 0; i < ndof2d; i++)
      trig2segm[i] = tensor2segm[trig2tensor[i][1] * 2 + trig2tensor[i][2]];


    segm2tensor = FlatArray<int[2]> (ndsegm, lh);

    for (int i3 = 0; i3 < horder; i3++)
      for (int i4 = 0; i4 < 2; i4++)
        {
          int i1d = i3 * 2 + i4;
          if (tensor2segm[i1d] != -1)
            {
              segm2tensor[tensor2segm[i1d]][0] = i3;
              segm2tensor[tensor2segm[i1d]][1] = i4;
            }
        }    

    
    map2dto1d = trig2segm;
    map3dto2d = prism2trig;




    ndof1d = ndsegm;
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
    


  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void  H1HighOrderTP<ET_PRISM> :: T_CalcShape (Tx x, Ty y, Tz z, TFA & sds) const
  {
    if (sort[0]+3 == sort[3])  // same singular vertex at bottom and top
      {
        Tx lami2d[3] = { x, y, 1-x-y };

        Tx lami2ds[3];
        for (int i = 0; i < 3; i++)
          lami2ds[i] = lami2d[sort[i]];
        
        int horder = max(2, order);
        ArrayMem<Tx, 60> polx (horder+1), 
          poly(horder+1), 
          polyy(2);

        ArrayMem<Tx, 100> memx((order+1)*horder);
        FlatMatrix<Tx> polsx(order+1, horder, &memx[0]);

        for (int i = 0; i <= order; i++)
          {
            Tx * hp = &polsx(i, 1);
#ifdef JACOBI
            JacobiPolynomial (horder-2, 2*lami2ds[0]-1, max(0,2*i-2), 0, hp);
#else
            LegendrePolynomial (horder-2, 2*lami2ds[0]-1, hp);
#endif
            for (int j = 1; j < horder; j++)
              polsx(i,j) *= lami2ds[0];
            polsx(i,0) = 1.0;
          }


        Tx * hp = &poly[1];
        ScaledLegendrePolynomialMult (horder-2, lami2ds[1]-lami2ds[2], 1-lami2ds[0], lami2ds[1], hp);
        poly[0] = 1;

        polyy[0] = 1;
        polyy[1] = lami2ds[2];


        ArrayMem<Tx,10> polz(order+1);
        polz[0] = 1-z;
        polz[1] = z;
        hp = &polz[2];
        LegendrePolynomialMult (order-2, 2*z-1, z*(1-z), hp);


        for (int i = 0; i < ndof; i++)
          {
            sds[i] = factor[i] * 
              polsx(prism2tensor[i][1]+prism2tensor[i][2], prism2tensor[i][0]) * 
              poly[prism2tensor[i][1]] * 
              polyy[prism2tensor[i][2]] *
              polz[prism2tensor[i][3]];
          }

        return;
      }



    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

    Tx lami[6] = { x, y, 1-x-y, x, y, 1-x-y };
    Tx muz[6]  = { 1-z, 1-z, 1-z, z, z, z };

    
    Tx lamisum[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        {
          if (vnums[j] > vnums[i]) lamisum[i] += lami[j];
          if (vnums[j+3] > vnums[i+3]) lamisum[i+3] += lami[j+3];
        }

    // sds = 0.0;

    // vertex shapes
    for (int i = 0; i < 6; i++)
      sds[i] = lami[i] * muz[i];

    int ii = 6;

    // horizontal edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_PRISM);
    for (int i = 0; i < 6; i++)
      if (order_edge[i] >= 2)
        {
          int p = order_edge[i];
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es, ee);
	  
          Tx * hp = &sds[ii];
          ScaledLegendrePolynomialMult (p-2, lami[es]-lamisum[es], lami[es]+lamisum[es], lami[es]*lami[ee]*muz[ee], hp);
          ii += p-1;
        }
    
    // vertical edges
    for (int i = 6; i < 9; i++)
      if (order_edge[i] >= 2)
        {
          int p = order_edge[i];
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es, ee);
	  
          Tx * hp = &sds[ii];
          LegendrePolynomialMult (p-2, muz[es]-muz[ee], muz[es]*muz[ee]*lami[ee], hp);
          ii += p-1;
        }
    

    ArrayMem<Tx,20> polx(order+1), poly(order+1), polz(order+1);

    const FACE * faces = ElementTopology::GetFaces (ET_PRISM);
    // trig face dofs

#ifdef JACOBI
    ArrayMem<Tx, 100> memx((order+1)*order);
    FlatMatrix<Tx> polsx(order+1, order, &memx[0]);

    for (int i = 0; i < 2; i++)
      if (order_face[i][0] >= 3)
        {
          int fav[3] = { faces[i][0], faces[i][1], faces[i][2] };
          if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
          if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
          if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
	  
          for (int j = 0; j <= order; j++)
            {
              Tx * hp = &polsx(j, 0);
              JacobiPolynomial (order-2, lami[fav[0]]-lamisum[fav[0]], max(0,2*j-2), 0, hp);
            }
          ScaledLegendrePolynomial (order-2, lami[fav[1]]-lamisum[fav[1]], lamisum[fav[0]], poly);

          Tx bub = lami[fav[0]] * lami[fav[1]] * lami[fav[2]] * muz[fav[2]];
          
          int n = order_face[i][0];
          for (int j2 = 0; j2 <= n-3; j2++)
            for (int i2 = 0; i2 <= n-3-j2; i2++, ii++)
              sds[ii] = bub * polsx(i2+2,j2) * poly[i2];
        }
#else
    for (int i = 0; i < 2; i++)
      if (order_face[i][0] >= 3)
        {
          int fav[3] = { faces[i][0], faces[i][1], faces[i][2] };
          if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
          if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
          if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
	  
          LegendrePolynomialMult (order-2, lami[fav[0]]-lamisum[fav[0]], lami[fav[0]], polx);
          ScaledLegendrePolynomialMult (order-2, lami[fav[1]]-lamisum[fav[1]], lamisum[fav[0]], lami[fav[1]], poly);

          Tx bub = lami[fav[2]] * muz[fav[2]];
          
          int n = order_face[i][0];
          for (int j2 = 0; j2 <= n-3; j2++)
            {
              Tx hv = bub * polx[j2];
              for (int i2 = 0; i2 <= n-3-j2; i2++, ii++)
                sds[ii] = hv * poly[i2];
            }
        }
#endif


    // quad face dofs
    for (int i = 2; i < 5; i++)
      if (order_face[i][0] >= 2 && order_face[i][1] >= 2)
        {
          INT<2> p = order_face[i];
          int pp = int(max2(p[0], p[1]))+1;
	 
          int fmax = 0;
          for(int j = 1; j < 4; j++) 
            if(vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;
	  
          int fz = 3-fmax; 
          int ftrig = fmax^1; 
	  
          Tx xi = lami[faces[i][fmax]] - lami[faces[i][ftrig]]; 
          Tx eta = 1-lami[faces[i][fmax]]-lami[faces[i][ftrig]]; 
          Tx zeta = muz[faces[i][fmax]]-muz[faces[i][fz]]; 

          Tx l1 = lami[faces[i][fmax]];
          Tx l2 = lami[faces[i][ftrig]];
          ScaledLegendrePolynomialMult (pp-2, l2-l1, l1+l2, l1*l2, polx);
          LegendrePolynomialMult (pp-2, -zeta, 0.25 * (1-zeta*zeta), polz);

          // global x-direction is towards second-largest vertex (JS)
          if (vnums[faces[i][ftrig]] > vnums[faces[i][fz]])  
            for (int k = 0; k < p[0]-1; k++)
              for (int j = 0; j < p[1]-1; j++)
                sds[ii++] = polx[k] * polz[j];
          else
            for (int j = 0; j < p[0]-1; j++)
              for (int k = 0; k < p[1]-1; k++)
                sds[ii++] = polx[k] * polz[j];
        }
    
    // volume dofs:
    if (order_inner[0] > 2 && order_inner[2] > 1)
      {
        ArrayMem<Tx,20> pol_trig((order_inner[0]-1)*(order_inner[0]-2)/2);
        int nf = T_TRIGFACESHAPES::Calc (order_inner[0], x-y, 1-x-y, pol_trig);

        T_ORTHOPOL:: Calc(order_inner[2], 2*z-1, polz);
        for (int i = 0; i < nf; i++)
          for (int k = 0; k < order_inner[2]-1; k++)
            sds[ii++] = pol_trig[i] * polz[k];
      }
  }





  // ******************************** HEXX **************************************

  void H1HighOrderTP<ET_HEX> :: SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    H1HighOrderHex<> :: SetVertexNumbers (avnums, lh);

    hex2tensor = FlatArray<int[3]> (ndof, lh);
    factor = FlatArray<double> (ndof, lh);

    static int hex_points[][3] = 
      { 
        { 0, 0, 0 },
        { 1, 0, 0 },
        { 1, 1, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
        { 1, 0, 1 },
        { 1, 1, 1 },
        { 0, 1, 1 }
      };


    // vertex shapes

    for (int i = 0; i < 8; i++)
      {
        hex2tensor[i][0] = hex_points[i][0];
        hex2tensor[i][1] = hex_points[i][1];
        hex2tensor[i][2] = hex_points[i][2];
        factor[i] = 1;
      }
        
    int ii = 8;

    // edge dofs
    const EDGE * edges = ElementTopology::GetEdges (ET_HEX);
    for (int i = 0; i < 12; i++)
      if (order_edge[i] >= 2)
        {
          int p = order_edge[i];
          int es = edges[i][0], ee = edges[i][1];
          if (vnums[es] > vnums[ee]) swap (es, ee);

          double fac = 1, mult = 1;
          int dir = 0;

          for (int k = 0; k < 3; k++)
            if (hex_points[es][k] != hex_points[ee][k])
              {
                dir = k;
                if (hex_points[es][k] < hex_points[ee][k])
                  mult = -1;
              }

          for (int l = 0; l <= p-2; l++)
            {
              for (int k = 0; k < 3; k++)
                hex2tensor[ii][k] = hex_points[ee][k];
              hex2tensor[ii][dir] = l+2;

              factor[ii] = fac;
              fac *= mult;
              ii++;
            }
        }


     
    const FACE * faces = ElementTopology::GetFaces (ET_HEX);
    for (int i = 0; i < 6; i++)
      if (order_face[i][0] >= 2 && order_face[i][1] >= 2)
        {
          int px = order_face[i][0];
          int py = order_face[i][1];
          int p = max2(px, py);

          int fmax = 0;
          for (int j=1; j<4; j++) 
            if (vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;  
	  
          int f1 = faces[i][(fmax+3)%4];
          int f2 = faces[i][(fmax+1)%4]; 
          fmax = faces[i][fmax]; 
	  
          if(vnums[f2]>vnums[f1]) swap(f1,f2);  // fmax > f1 > f2 
	  

          int ind[3];
          double fac1 = 1, fac2 = 1, mult1 = 1, mult2 = 1;
          for (int l = 0; l <= px-2; l++)
            {
              fac2 = 1;
              for (int m = 0; m <= py-2; m++)
                {
                  for (int k = 0; k < 3; k++)
                    {
                      ind[k] = hex_points[fmax][k];
                      if (hex_points[fmax][k] != hex_points[f2][k])
                        {
                          ind[k] = m+2;
                          if (hex_points[fmax][k] > hex_points[f2][k])
                            mult2 = -1;
                        }
                      if (hex_points[fmax][k] != hex_points[f1][k])
                        {
                          ind[k] = l+2;
                          if (hex_points[fmax][k] > hex_points[f1][k])
                            mult1 = -1;
                        }
                    }

                  factor[ii] = fac1*fac2;
                  hex2tensor[ii][0] = ind[0];
                  hex2tensor[ii][1] = ind[1];
                  hex2tensor[ii][2] = ind[2];
                  ii++;

                  fac2 *= mult2;
                }
              fac1 *= mult1;
            }

        }
    
    // volume dofs:
    if (order_inner[0] >= 2 && order_inner[1] >= 2 && order_inner[2] >= 2)
      {
        for (int i = 2; i <= order_inner[0]; i++)
          for (int j = 2; j <= order_inner[1]; j++)
            for (int k = 2; k <= order_inner[2]; k++)
              {
                factor[ii] = 1;
                hex2tensor[ii][0] = k;
                hex2tensor[ii][1] = j;
                hex2tensor[ii][2] = i;
                ii++;
              }
      }    


    ndof2d = (order+1)*(order+1);

    map3dto2d = FlatArray<int> (ndof, lh);
    for (int jj = 0; jj < ndof; jj++)
      map3dto2d[jj] = hex2tensor[jj][1] * (order+1) + hex2tensor[jj][2];

    map2dto1d = FlatArray<int> (ndof2d, lh);
    for (int jj = 0; jj < ndof2d; jj++)
      map2dto1d[jj] = jj % (order+1);






    ndof1d = (order+1);
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


  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void  H1HighOrderTP<ET_HEX> :: T_CalcShape (Tx x, Ty y, Tz z, TFA & sds) const
  {
    Tx xi[3] = { 2*x-1, 2*y-1, 2*z-1 };
    
    ArrayMem<Tx, 90> mem(3*(order+1));
    FlatMatrix<Tx> poli(3, order+1, &mem[0]);

    for (int i = 0; i < 3; i++)
      {
        poli(i,0) = 0.5 * (1-xi[i]);
        poli(i,1) = 0.5 * (1+xi[i]);
        LegendrePolynomialMult (order-2, xi[i], 0.25 * (1-xi[i]*xi[i]), poli.Row(i).Addr(2));
      }

    for (int i = 0; i < ndof; i++)
      sds[i] = factor[i] * 
        poli(0,hex2tensor[i][0]) * 
        poli(1,hex2tensor[i][1]) * 
        poli(2,hex2tensor[i][2]);
  }
}
