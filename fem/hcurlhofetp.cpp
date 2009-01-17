/*********************************************************************/
/* File:   h1hofe.cpp                                                */
/* Author: Start                                                     */
/* Date:   1. June 2007                                              */
/*********************************************************************/

// #define NOPROFILE

 
#include <fem.hpp>
namespace ngfem
{
  using namespace ngfem;


  void HCurlHighOrderTetTP :: SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    HCurlHighOrderTet<> :: SetVertexNumbers (avnums, lh);

    for (int i = 0; i < 4; i++) sort[i] = i;

    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
    if (vnums[sort[2]] > vnums[sort[3]]) Swap (sort[2], sort[3]);
    if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
    if (vnums[sort[1]] > vnums[sort[3]]) Swap (sort[1], sort[3]);
    if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);

    // vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]] < vnums[sort[3]]


    // build tables

    tet2tensor = FlatArray<int[4]> (ndof, lh);
    split = FlatArray<int> (ndof, lh);

    int isort[4];
    for (int i = 0; i < 4; i++)
      isort[sort[i]] = i;

    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 4; j++)
        tet2tensor[i][j] = 0;

    
    const EDGE * edges = ElementTopology::GetEdges (ET_TET);
    const FACE * faces = ElementTopology::GetFaces (ET_TET); 
    

    int ii = 6; 
    for (int i = 0; i < 6; i++)
      { 
	int p = order_edge[i]; 
	int es = edges[i][0]; 
	int ee = edges[i][1]; 

	if (vnums[es] > vnums[ee]) swap (es, ee);
	
	//Nedelec low order edge shape function 
        tet2tensor[i][isort[es]] = 1;
        tet2tensor[i][isort[ee]] = 1;

        // split[i] = max (isort[es], isort[ee]);
        split[i] = isort[es]+1;

	//HO-Edge shape functions (Gradient Fields) 	
	if(p>0 && usegrad_edge[i]) 
	  {
	    for(int j=0; j< p;j++,ii++) 	      
              {
                tet2tensor[ii][isort[es]] = 1+j;
                tet2tensor[ii][isort[ee]] = 1;
                split[ii] = 0;
              }
	  }
      }


    // face dofs
    for (int i = 0; i < 4; i++)
      if (order_face[i][0] >= 2)
        {
          int fav[3] = { faces[i][0], faces[i][1], faces[i][2] };
          
          if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
          if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
          if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]);

          int p = order_face[i][0];
          
          // gradients 
          if (usegrad_face[i])
            for (int j = 0; j <= p-2; j++)
              for (int k = 0; k <= p-2-j; k++, ii++)
                {
                  tet2tensor[ii][isort[fav[0]]] = j+1;
                  tet2tensor[ii][isort[fav[1]]] = k+1;
                  tet2tensor[ii][isort[fav[2]]] = 1;
                  split[ii] = 0;
                }

          for (int j = 0; j <= p-2; j++)
            for (int k = 0; k <= p-2-j; k++, ii++)
              {
                tet2tensor[ii][isort[fav[0]]] = j+1;
                tet2tensor[ii][isort[fav[1]]] = k+1;
                tet2tensor[ii][isort[fav[2]]] = 1;
                split[ii] = isort[fav[0]]+1;
              }

          for (int j = 0; j <= p-2; j++, ii++)
            // for (int k = 0; k <= p-2-j; k++, ii++)
              {
                
                tet2tensor[ii][isort[fav[0]]] = j+1;
                tet2tensor[ii][isort[fav[1]]] = 1;
                tet2tensor[ii][isort[fav[2]]] = 1;
                split[ii] = isort[fav[1]]+1;
              }
        }


    // inner dofs
    if (order_inner[0] >= 3)
      {
        int n = order_inner[0];

        if(usegrad_cell)
          for (int i = 0; i <= n-3; i++)
            for (int j = 0; j <= n-3-i; j++)
              for (int k = 0; k <= n-3-i-j; k++, ii++)
                {
                  tet2tensor[ii][0] = k+1;
                  tet2tensor[ii][1] = j+1;
                  tet2tensor[ii][2] = i+1;
                  tet2tensor[ii][3] = 1;
                  split[ii] = 0;
                }


        for (int i = 0; i <= n-3; i++)
          for (int j = 0; j <= n-3-i; j++)
            for (int k = 0; k <= n-3-i-j; k++)
              for (int l = 1; l <= 2; l++, ii++)
                {
                  tet2tensor[ii][0] = k+1;
                  tet2tensor[ii][1] = j+1;
                  tet2tensor[ii][2] = i+1;
                  tet2tensor[ii][3] = 1;
                  split[ii] = l;
                }

        for (int i = 0; i <= n-3; i++)
          for (int j = 0; j <= n-3-i; j++, ii++)
            {
              tet2tensor[ii][0] = i+1;
              tet2tensor[ii][1] = j+1;
              tet2tensor[ii][2] = 1;
              tet2tensor[ii][3] = 1;
              split[ii] = 3;
            }
        
      }


    int horder = max(order+1,2);
    int nd2d = horder * horder * 2 * 3;
    
    FlatArray<int> tensor2trig ( nd2d, lh);
    tensor2trig = -1;
    
    for (int ii = 0; ii < ndof; ii++)
      {
        int hsplit = max(1, split[ii])-1;

        int j2d = ((tet2tensor[ii][1] * horder + tet2tensor[ii][2]) * 2 + tet2tensor[ii][3]) + hsplit * nd2d/3;
        tensor2trig[j2d] = 0;
      }
    
    int ndof2d =0;
    for (int i = 0; i < nd2d; i++)
      if (tensor2trig[i] != -1)
        {
          tensor2trig[i] = ndof2d;
          ndof2d++;
        }

    map3dto2d = FlatArray<int> (ndof, lh);
    trig2tensor = FlatArray<int[3]> (ndof2d, lh);
    split2d = FlatArray<int> (ndof2d, lh);
    
    for (int ii = 0; ii < ndof; ii++)
      {
        int hsplit = max(1, split[ii])-1;

        int j2d = ((tet2tensor[ii][1] * horder + tet2tensor[ii][2]) * 2 + tet2tensor[ii][3]) + hsplit * nd2d/3;
        int jtrig = tensor2trig[j2d];
        
        map3dto2d[ii] = jtrig;
        trig2tensor[jtrig][0] = tet2tensor[ii][1];
        trig2tensor[jtrig][1] = tet2tensor[ii][2];
        trig2tensor[jtrig][2] = tet2tensor[ii][3];
        split2d[jtrig] = hsplit;
      }


    map2dto1d = FlatArray<int> (ndof2d, lh);
    FlatArray<int[2]> hsegm2tensor(ndof2d, lh);
    FlatArray<int> hsplit1d (ndof2d, lh);
    
    int ndof1d = 0;
    for (int i = 0; i < ndof2d; i++)
      {
        bool has = 0;
        for (int j = 0; j < ndof1d; j++)
          if (hsegm2tensor[j][0] == trig2tensor[i][1] &&
              hsegm2tensor[j][1] == trig2tensor[i][2] &&
              hsplit1d[j] == max(split2d[i]-1, 0))
            {
              map2dto1d[i] = j;
              has = 1;
              break;
            }
        if (!has)
          {
            map2dto1d[i] = ndof1d;
            hsegm2tensor[ndof1d][0] = trig2tensor[i][1];
            hsegm2tensor[ndof1d][1] = trig2tensor[i][2];
            hsplit1d[ndof1d] = max(split2d[i]-1, 0);
            ndof1d++;
          }
      }
      
    segm2tensor = hsegm2tensor.Range (0, ndof1d);
    split1d = hsplit1d.Range (0, ndof1d);

    /*
    *testout << "map3dto2d = " << endl << map3dto2d << endl;
    *testout << "map2dto1d = " << endl << map2dto1d << endl;
   
    for (int i = 0; i < ndof; i++)
      {
        *testout << "tet-shape " << i << ". ";
        for (int j = 0; j < 4; j++)
          *testout << " " << tet2tensor[i][j];
        *testout << ", split = " << split[i]  << endl;
      }

    for (int i = 0; i < ndof2d; i++)
      {
        *testout << "trig-shape " << i << ". ";
        for (int j = 0; j < 3; j++)
          *testout << " " << trig2tensor[i][j];
        *testout << ", split = " << split2d[i] << endl;
      }

    for (int i = 0; i < ndof1d; i++)
      {
        *testout << "segm-shape " << i << ". ";
        for (int j = 0; j < 2; j++)
          *testout << " " << segm2tensor[i][j];
        *testout << ", split = " << split1d[i] << endl;
      }
    */
  }
    



  void HCurlHighOrderTetTP :: CalcShape (const IntegrationPoint & ip, 
                                         FlatMatrixFixWidth<3> shape) const
  {
    AutoDiff<3> x(ip(0), 0);
    AutoDiff<3> y(ip(1), 1);
    AutoDiff<3> z(ip(2), 2);

    AutoDiff<3> lami[4] = { x, y, z, 1-x-y-z };

    AutoDiff<3> lamis[4];
    for (int i = 0; i < 4; i++)
      lamis[i] = lami[sort[i]];

    
    int horder = max(2, order+1);
    ArrayMem<AutoDiff<3>, 20> polx(horder), poly(horder), polz(horder), polzz(2);

    AutoDiff<3> * hp = &polx[1];
    LegendrePolynomialMult (horder-2, 2*lamis[0]-1, lamis[0], hp);
    polx[0] = 1.0;

    hp = &poly[1];
    ScaledLegendrePolynomialMult (horder-2, lamis[1]-lamis[2]-lamis[3], 1-lamis[0], lamis[1], hp);
    poly[0] = 1.0;

    hp = &polz[1];
    ScaledLegendrePolynomialMult (horder-2, lamis[2]-lamis[3], lamis[2]+lamis[3], lamis[2], hp);
    polz[0] = 1.0;

    polzz[0] = 1.0;
    polzz[1] = lamis[3];


    AutoDiff<3> *pols[4];
    pols[0] = &polx[0];
    pols[1] = &poly[0];
    pols[2] = &polz[0];
    pols[3] = &polzz[0];

    for (int i = 0; i < 4; i++)
      {
        pols[i][0] = 1.0;
        pols[i][1] = lamis[i];
      }

    for (int i = 0; i < ndof; i++)
      {
        AutoDiff<3> u = 1.0, v = 1.0, w = 1.0;

        for (int j = 0; j < split[i]-1; j++)
          w *= pols[j][tet2tensor[i][j]];
        for (int j = max (0, split[i]-1); j < split[i]; j++)
          u *= pols[j][tet2tensor[i][j]];
        for (int j = split[i]; j < 4; j++)
          v *= pols[j][tet2tensor[i][j]];
        
        for (int j = 0; j < 3; j++)
          shape(i,j) = w.Value() * (u.Value() * v.DValue(j) - u.DValue(j) * v.Value());
      }
  }


  inline AutoDiff<3> Cross (const AutoDiff<3> & u,
			    const AutoDiff<3> & v)
  {
    AutoDiff<3> hv;
    hv.Value() = 0.0;
    hv.DValue(0) = u.DValue(1)*v.DValue(2)-u.DValue(2)*v.DValue(1);
    hv.DValue(1) = u.DValue(2)*v.DValue(0)-u.DValue(0)*v.DValue(2);
    hv.DValue(2) = u.DValue(0)*v.DValue(1)-u.DValue(1)*v.DValue(0);
    return hv;
  }


  
  void HCurlHighOrderTetTP :: CalcCurlShape (const IntegrationPoint & ip, 
                                             FlatMatrixFixWidth<3> shape) const
  {
    // HCurlFiniteElementD<3> :: CalcCurlShape (ip, shape);
    // return;

    // *testout << "num curl shape = " << endl << shape << endl;

    AutoDiff<3> x(ip(0), 0);
    AutoDiff<3> y(ip(1), 1);
    AutoDiff<3> z(ip(2), 2);

    AutoDiff<3> lami[4] = { x, y, z, 1-x-y-z };

    AutoDiff<3> lamis[4];
    for (int i = 0; i < 4; i++)
      lamis[i] = lami[sort[i]];

    
    int horder = max(2, order+1);
    ArrayMem<AutoDiff<3>, 20> polx(horder), poly(horder), polz(horder), polzz(2);

    AutoDiff<3> * hp = &polx[1];
    LegendrePolynomialMult (horder-2, 2*lamis[0]-1, lamis[0], hp);
    polx[0] = 1.0;

    hp = &poly[1];
    ScaledLegendrePolynomialMult (horder-2, lamis[1]-lamis[2]-lamis[3], 1-lamis[0], lamis[1], hp);
    poly[0] = 1.0;

    hp = &polz[1];
    ScaledLegendrePolynomialMult (horder-2, lamis[2]-lamis[3], lamis[2]+lamis[3], lamis[2], hp);
    polz[0] = 1.0;

    polzz[0] = 1.0;
    polzz[1] = lamis[3];

    AutoDiff<3> *pols[4];
    pols[0] = &polx[0];
    pols[1] = &poly[0];
    pols[2] = &polz[0];
    pols[3] = &polzz[0];


    for (int i = 0; i < ndof; i++)
      {
        int spl = split[i];

        if (spl == 0)
          {
            for (int j = 0; j < 3; j++)
              shape(i,j) = 0;
            continue;
          }

        AutoDiff<3> u = 1.0, v = 1.0, w = 1.0;;

        for (int j = 0; j < spl-1; j++)
          w *= pols[j][tet2tensor[i][j]];
        u = pols[spl-1][tet2tensor[i][spl-1]];
        for (int j = spl; j < 4; j++)
          v *= pols[j][tet2tensor[i][j]];

        AutoDiff<3> cs = Cross (u*w, v) + Cross (u, v*w);
        for (int j = 0; j < 3; j++)
          shape(i,j) = cs.DValue(j);
      }
    // *testout << "new shape = " << endl << shape << endl;
  }












  /*    *************************************  PRISM ********************************** */




  void HCurlHighOrderPrismTP :: SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh)
  {
    HCurlHighOrderPrism<> :: SetVertexNumbers (avnums, lh);

    for (int i = 0; i < 6; i++) sort[i] = i;

    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
    if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
    if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);

    if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);
    if (vnums[sort[4]] > vnums[sort[5]]) Swap (sort[4], sort[5]);
    if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);

    if (sort[0]+3 != sort[3]) return;


    // build tables

    prism2tensor = FlatArray<int[4]> (ndof, lh);
    factorxy = FlatArray<double> (ndof, lh);
    factorz = FlatArray<double> (ndof, lh);
    split = FlatArray<int> (ndof, lh);
    // splita = FlatArray<int> (ndof, lh);
    isgradient = FlatArray<bool> (ndof, lh);

    // splita = 0;
    isgradient = false;

    int isort[6];
    for (int i = 0; i < 6; i++)
      isort[sort[i]] = i;

    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 4; j++)
        prism2tensor[i][j] = 0;

    
    const EDGE * edges = ElementTopology::GetEdges (ET_PRISM);
    const FACE * faces = ElementTopology::GetFaces (ET_PRISM); 
    

    int ii = 9; 
    // horizontal edge dofs
    for (int i = 0; i < 6; i++)
      { 
	int p = order_edge[i]; 
	int es = edges[i][0]; 
	int ee = edges[i][1]; 

	if (vnums[es%3] > vnums[ee%3]) swap (es, ee);
        int mult = ( (vnums[es] > vnums[ee]) == (vnums[es%3] > vnums[ee%3])) ? 1 : -1;
	
	//Nedelec low order edge shape function 
        prism2tensor[i][isort[es%3]] = 1;
        prism2tensor[i][isort[ee%3]] = 1;
        prism2tensor[i][3] = i / 3;

        // split[i] = isort[ee%3];
        split[i] = isort[es%3]+1;
        // splita[i] = isort[es%3];

        factorxy[i] = mult;
        factorz[i] = 0;
        
	//HO-Edge shape functions (Gradient Fields) 	
        int fac = 1;
	if(p>0 && usegrad_edge[i]) 
	  {
	    for(int j=0; j< p;j++,ii++) 	      
              {
                prism2tensor[ii][isort[es%3]] = 1+j;
                prism2tensor[ii][isort[ee%3]] = 1;
                prism2tensor[ii][3] = i / 3;
                split[ii] = 0;
                isgradient[ii] = true;
                factorxy[ii] = fac;
                factorz[ii] = fac;
                fac *= mult;
              }
	  }
      }


    // vertical edges
    for (int i = 6; i < 9; i++)
      {
        int p = order_edge[i];
        int es = edges[i][0];
        int ee = edges[i][1];
        
        int mult = 1, fac = 1;
        if (vnums[es] < vnums[ee]) mult = -1;   
        
        //Nedelec low order edge shape function 
        prism2tensor[i][isort[es]%3] = 1;
        prism2tensor[i][3] = 1;  // the constant 1 
        split[i] = 0;
        factorxy[i] = 0;
        factorz[i] = -mult;
        
        if (p>0 && usegrad_edge[i])
          {
            int nde = order_edge[i]-1;
            for (int j = 0; j < p; j++, ii++)
              {
                prism2tensor[ii][isort[es%3]] = 1;
                prism2tensor[ii][3] = j+2;
                split[ii] = 0;
                isgradient[ii] = true;
                factorxy[ii] = fac;
                factorz[ii] = fac;
                fac *= mult;
              }
          }
      }
    


    // trig face dofs
    for (int i = 0; i < 2; i++)
      if (order_face[i][0] >= 2)
        {
          int fav[3] = { faces[i][0], faces[i][1], faces[i][2] };
          
          if(vnums[fav[0]%3] > vnums[fav[1]%3]) swap(fav[0],fav[1]); 
          if(vnums[fav[1]%3] > vnums[fav[2]%3]) swap(fav[1],fav[2]);
          if(vnums[fav[0]%3] > vnums[fav[1]%3]) swap(fav[0],fav[1]);

          int mult = ( (vnums[fav[1]] > vnums[fav[2]]) == (vnums[fav[1]%3] > vnums[fav[2]%3])) ? 1 : -1;

          int p = order_face[i][0];
          
          // gradients 
          if (usegrad_face[i])
            for (int j = 0; j <= p-2; j++)
              {
                int fac = 1;
                for (int k = 0; k <= p-2-j; k++, ii++)
                  {
                    prism2tensor[ii][isort[fav[0]%3]] = j+1;
                    prism2tensor[ii][isort[fav[1]%3]] = k+1;
                    prism2tensor[ii][isort[fav[2]%3]] = 1;
                    prism2tensor[ii][3] = i;
                    split[ii] = 0;
                    isgradient[ii] = true;
                    
                    factorxy[ii] = fac;
                    factorz[ii] = fac;
                    fac *= mult;
                  }
              }

          for (int j = 0; j <= p-2; j++)
            {            
              int fac = 1;
              for (int k = 0; k <= p-2-j; k++, ii++)
                {
                  prism2tensor[ii][isort[fav[0]%3]] = j+1;
                  prism2tensor[ii][isort[fav[1]%3]] = k+1;
                  prism2tensor[ii][isort[fav[2]%3]] = 1;
                  prism2tensor[ii][3] = i;
                  split[ii] = 1; // isort[fav[1]%3];

                  factorxy[ii] = fac;
                  factorz[ii] = 0;
                  fac *= mult;
                }
            }

          for (int j = 0; j <= p-2; j++, ii++)
            {
              prism2tensor[ii][isort[fav[0]%3]] = j+1;
              prism2tensor[ii][isort[fav[1]%3]] = 1;
              prism2tensor[ii][isort[fav[2]%3]] = 1;
              prism2tensor[ii][3] = i;
              split[ii] = 2; // isort[fav[2]%3]; 
              // splita[ii] = 1;
              
              factorxy[ii] = mult;
              factorz[ii] = 0;
            }
        }






    // quad faces
    for (int i = 2; i < 5; i++)
      {
	INT<2> p = order_face[i];
	 
	int fmax = 0;
	for (int j = 1; j < 4; j++)
	  if (vnums[faces[i][j]] > vnums[faces[i][fmax]]) fmax = j;

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
            int facz, factrig;

            if (usegrad_face[i])
              {
                factrig = 1;
                for (int k = 0; k < p[0]; k++)
                  {
                    facz = 1;
                    for (int j = 0; j < p[1]; j++, ii++)
                      {
                        prism2tensor[ii][isort[es%3]] = k+1;
                        prism2tensor[ii][isort[ee%3]] = 1;
                        prism2tensor[ii][3] = j+2;
                        
                        factorxy[ii] = facz * factrig;
                        factorz[ii] = facz * factrig;
                        split[ii] = 0;
                        isgradient[ii] = true;

                        facz *= multz;
                      }
                    factrig *= multtrig;
                  }
              }

            factrig = 1;
            for (int k = 0; k < p[0]; k++)
              {
                facz = 1;
                for (int j = 0; j < p[1]; j++, ii++)
                  {
                    prism2tensor[ii][isort[es%3]] = k+1;
                    prism2tensor[ii][isort[ee%3]] = 1;
                    prism2tensor[ii][3] = j+2;

                    factorxy[ii] = facz * factrig;
                    factorz[ii] = -facz * factrig;                    
                    split[ii] = 0;                        
                    facz *= multz;
                  }
                factrig *= multtrig;
              }
            
            factrig = -multz;    // nabla(trig) * e_z
            for (int k = 0; k < p[0]; k++, ii++)
              {
                prism2tensor[ii][isort[es%3]] = k+1;
                prism2tensor[ii][isort[ee%3]] = 1;
                prism2tensor[ii][3] = 1;
                    
                split[ii] = 0;        
                factorxy[ii] = 0;
                factorz[ii] = factrig;
                factrig *= multtrig;
              }

            facz = multtrig;
            for (int j = 0; j < p[1]; j++, ii++)
              {
                prism2tensor[ii][isort[es%3]] = 1;
                prism2tensor[ii][isort[ee%3]] = 1;
                prism2tensor[ii][3] = j+2;
                
                // split[ii] = isort[ee%3];
                split[ii] = isort[es%3]+1;                
                // splita[ii] = isort[es%3];

                factorxy[ii] = facz;
                factorz[ii] = 0;
                facz *= multz;
              }

            
          }
        else
          {
            int facz, factrig;

            if (usegrad_face[i])
              {
                facz = 1;
                for (int j = 0; j < p[0]; j++)
                  {
                    factrig = 1;
                    for (int k = 0; k < p[1]; k++, ii++)
                      {
                        prism2tensor[ii][isort[es%3]] = k+1;
                        prism2tensor[ii][isort[ee%3]] = 1;
                        prism2tensor[ii][3] = j+2;

                        factorxy[ii] = facz * factrig;
                        factorz[ii] = facz * factrig;
                        split[ii] = 0;                        
                        isgradient[ii] = true;

                        factrig *= multtrig;
                      }
                    facz *= multz;
                  }
              }

            facz = 1;
            for (int j = 0; j < p[0]; j++)
              {
                int factrig = 1;
                for (int k = 0; k < p[1]; k++, ii++)
                  {
                    prism2tensor[ii][isort[es%3]] = k+1;
                    prism2tensor[ii][isort[ee%3]] = 1;
                    prism2tensor[ii][3] = j+2;

                    factorxy[ii] = -facz * factrig;
                    factorz[ii] =  facz * factrig;                    
                    split[ii] = 0;                        
                    factrig *= multtrig;
                  }
                facz *= multz;
              }


            facz = multtrig;
            for (int j = 0; j < p[1]; j++, ii++)
              {
                prism2tensor[ii][isort[es%3]] = 1;
                prism2tensor[ii][isort[ee%3]] = 1;
                prism2tensor[ii][3] = j+2;
                
                // split[ii] = isort[ee%3];
                split[ii] = isort[es%3]+1;
                // splita[ii] = isort[es%3];
                factorxy[ii] = facz; // * factrig;
                factorz[ii] = 0;
                facz *= multz;
              }

            factrig = -multz;    // nabla(trig) * const
            for (int k = 0; k < p[0]; k++, ii++)
              {
                prism2tensor[ii][isort[es%3]] = k+1;
                prism2tensor[ii][isort[ee%3]] = 1;
                prism2tensor[ii][3] = 1;
                    
                split[ii] = 0;        
                factorxy[ii] = 0;
                factorz[ii] = factrig;
                factrig *= multtrig;
              }
          }
      }
    
    if(order_inner[0] > 1&& order_inner[2]>0) 
      {
	// gradientfields
	if(usegrad_cell)
	  for (int i=0;i<=order_inner[0]-2;i++)
	    for (int j=0;j<=order_inner[0]-2-i;j++)
	      for (int k=0;k<=order_inner[2]-1;k++)
		{
                  prism2tensor[ii][0] = i+1;
                  prism2tensor[ii][1] = j+1;
                  prism2tensor[ii][2] = 1;
                  prism2tensor[ii][3] = k+2;
                  
                  split[ii] = 0;        
                  factorxy[ii] = 1;
                  factorz[ii] = 1;
                  isgradient[ii] = true;
                  ii++;
		}
	 

	// Rotations of gradientfields
	for (int i=0;i<=order_inner[0]-2;i++)
	  for (int j=0;j<=order_inner[0]-2-i;j++)
	    for (int k=0;k<=order_inner[2]-1;k++)
	      {
                prism2tensor[ii][0] = i+1;
                prism2tensor[ii][1] = j+1;
                prism2tensor[ii][2] = 1;
                prism2tensor[ii][3] = k+2;
                
                split[ii] = 0;        
                factorxy[ii] = 1;
                factorz[ii] = -1;
                ii++;



                prism2tensor[ii][0] = i+1;
                prism2tensor[ii][1] = j+1;
                prism2tensor[ii][2] = 1;
                prism2tensor[ii][3] = k+2;
                
                split[ii] = 1;        
                factorxy[ii] = 1;
                factorz[ii] = -1;
                ii++;
	      }

	// Type 3 
	// ned0(trig) * polxy2[j]*polz 
	// z.DValue(0) * polxy1[i] * polxy2[j] 
	// double ned_trig[2] = {y.Value(),-x.Value()};  
	for (int j=0;j<=order_inner[0]-2;j++) 
	  for (int k=0;k<=order_inner[2]-1;k++) 
            {
              prism2tensor[ii][0] = j+1;
              prism2tensor[ii][1] = 1;
              prism2tensor[ii][2] = 1;
              prism2tensor[ii][3] = k+2;
              
              split[ii] = 2;        
              // splita[ii] = 1;        
              factorxy[ii] = 1;
              //  factorz[ii] = -1;
              factorz[ii] = 0;
              ii++;
            }
	    
    	for (int i=0;i<=order_inner[0]-2;i++) 
	  for (int j=0;j<=order_inner[0]-2-i;j++) 
            {
              prism2tensor[ii][0] = i+1;
              prism2tensor[ii][1] = j+1;
              prism2tensor[ii][2] = 1;
              prism2tensor[ii][3] = 1;
              
              split[ii] = 0;        
              // splita[ii] = 0;        
              factorxy[ii] = 0;
              factorz[ii] = 1;
              ii++;
            }
	
      }

    /*
    for (int i = 0; i < ndof; i++)
      splita[i] = max(0, split[i]-1);
    */

    int horder = max(order+1,2);
    int nd2d = horder * horder * 2 * 6;
    
    FlatArray<int> tensor2trig ( nd2d, lh);
    tensor2trig = -1;
    
    for (int ii = 0; ii < ndof; ii++)
      {
        int hsplit = split[ii];
        // int hsplita = splita[ii];

        int j2d = ((prism2tensor[ii][0] * horder + prism2tensor[ii][1]) * 2 + prism2tensor[ii][2]) + (hsplit /* *2+hsplita */) * nd2d/6;
        tensor2trig[j2d] = 0;
      }
    
    int ndof2d =0;
    for (int i = 0; i < nd2d; i++)
      if (tensor2trig[i] != -1)
        {
          tensor2trig[i] = ndof2d;
          ndof2d++;
        }

    map3dto2d = FlatArray<int> (ndof, lh);
    trig2tensor = FlatArray<int[3]> (ndof2d, lh);
    split2d = FlatArray<int> (ndof2d, lh);
    // split2da = FlatArray<int> (ndof2d, lh);
    
    for (int ii = 0; ii < ndof; ii++)
      {
        int hsplit = split[ii];
        // int hsplita = splita[ii];
        int j2d = ((prism2tensor[ii][0] * horder + prism2tensor[ii][1]) * 2 + prism2tensor[ii][2]) + (hsplit /* *2+hsplita */) * nd2d/6;
        int jtrig = tensor2trig[j2d];
        
        map3dto2d[ii] = jtrig;
        trig2tensor[jtrig][0] = prism2tensor[ii][0];
        trig2tensor[jtrig][1] = prism2tensor[ii][1];
        trig2tensor[jtrig][2] = prism2tensor[ii][2];
        split2d[jtrig] = hsplit;
        // split2da[jtrig] = hsplita;
      }


    map2dto1d = FlatArray<int> (ndof2d, lh);
    FlatArray<int[2]> hsegm2tensor(ndof2d, lh);
    FlatArray<int> hsplit1d (ndof2d, lh);
    
    int ndof1d = 0;
    for (int i = 0; i < ndof2d; i++)
      {
        bool has = 0;
        for (int j = 0; j < ndof1d; j++)
          if (hsegm2tensor[j][0] == trig2tensor[i][1] &&
              hsegm2tensor[j][1] == trig2tensor[i][2] &&
              hsplit1d[j] == max(split2d[i]-1, 0))
            {
              map2dto1d[i] = j;
              has = 1;
              break;
            }
        if (!has)
          {
            map2dto1d[i] = ndof1d;
            hsegm2tensor[ndof1d][0] = trig2tensor[i][1];
            hsegm2tensor[ndof1d][1] = trig2tensor[i][2];
            hsplit1d[ndof1d] = max(split2d[i]-1, 0);
            ndof1d++;
          }
      }
      
    segm2tensor = hsegm2tensor.Range (0, ndof1d);
    split1d = hsplit1d.Range (0, ndof1d);

    /*
    *testout << "map3dto2d = " << endl << map3dto2d << endl;
    *testout << "map2dto1d = " << endl << map2dto1d << endl;
   
    for (int i = 0; i < ndof; i++)
      {
        *testout << "prism-shape " << i << ". ";
        for (int j = 0; j < 4; j++)
          *testout << " " << prism2tensor[i][j];
        *testout << ", split = " << split[i] << ", factorxy = " << factorxy[i] << ", factorz = " << factorz[i] << endl;
      }

    for (int i = 0; i < ndof2d; i++)
      {
        *testout << "trig-shape " << i << ". ";
        for (int j = 0; j < 3; j++)
          *testout << " " << trig2tensor[i][j];
        *testout << ", split = " << split2d[i] << endl; // ", splita = " << split2da[i] << endl;
      }

    for (int i = 0; i < ndof1d; i++)
      {
        *testout << "segm-shape " << i << ". ";
        for (int j = 0; j < 2; j++)
          *testout << " " << segm2tensor[i][j];
        *testout << ", split = " << split1d[i] << endl;
      }
    */
  }
    



  void HCurlHighOrderPrismTP :: CalcShape (const IntegrationPoint & ip, 
                                           FlatMatrixFixWidth<3> shape) const
  {
    /*
    HCurlHighOrderPrism<>::CalcShape (ip, shape);
    *testout << "ip = " << ip << endl;
    *testout << "basis, shape = " << shape << endl;
    */

    if (sort[0]+3 == sort[3])  // same singular vertex at bottom and top
      {
        AutoDiff<3> x(ip(0), 0);
        AutoDiff<3> y(ip(1), 1);
        AutoDiff<3> z(ip(2), 2);

        AutoDiff<3>  lami2d[3] = { x, y, 1-x-y };

        AutoDiff<3>  lami2ds[3];
        for (int i = 0; i < 3; i++)
          lami2ds[i] = lami2d[sort[i]];
        
        int horder = max(2, order+1);
        ArrayMem<AutoDiff<3> , 60> polx (horder+1), 
          poly(horder+1), 
          polyy(2);

        ArrayMem<AutoDiff<3> , 100> memx((horder+1)*horder);
        FlatMatrix<AutoDiff<3> > polsx(horder+1, horder, &memx[0]);

        for (int i = 0; i <= horder; i++)
          {
            AutoDiff<3>  * hp = &polsx(i, 1);
#ifdef JACOBI
            JacobiPolynomial (horder-2, 2*lami2ds[0]-1, max(0,2*i-2), 0, hp);
#else
            LegendrePolynomial (horder-2, 2*lami2ds[0]-1, hp);
#endif
            for (int j = 1; j < horder; j++)
              polsx(i,j) *= lami2ds[0];
            polsx(i,0) = 1.0;
          }

        for (int i = 0; i < horder; i++)
          polx[i] = polsx(0,i);


        AutoDiff<3>  * hp = &poly[1];
        ScaledLegendrePolynomialMult (horder-2, lami2ds[1]-lami2ds[2], 1-lami2ds[0], lami2ds[1], hp);
        poly[0] = 1;

        polyy[0] = 1;
        polyy[1] = lami2ds[2];


        
        ArrayMem<AutoDiff<3> ,10> polzmem(order+2);
        AutoDiff<3> * polz = &polzmem[0];
        polz[0] = 1-z;
        polz[1] = z;
        hp = &polz[2];
        LegendrePolynomialMult (order-1, 2*z-1, z*(1-z), hp);



        AutoDiff<3> *pols[4];
        pols[0] = &polx[0];
        pols[1] = &poly[0];
        pols[2] = &polyy[0];
        pols[3] = &polz[0];

        for (int i = 0; i < ndof; i++)
          {
            AutoDiff<3> u = 1.0, v = 1.0, uu = 1.0;

            /*
            for (int j = 0; j < splita[i]; j++)
              uu *= pols[j][prism2tensor[i][j]];
            for (int j = splita[i]; j < split[i]; j++)
              u *= pols[j][prism2tensor[i][j]];
            for (int j = split[i]; j < 3; j++)
              v *= pols[j][prism2tensor[i][j]];
            */
            for (int j = 0; j < split[i]-1; j++)
              uu *= pols[j][prism2tensor[i][j]];
            for (int j = max (0, split[i]-1); j < split[i]; j++)
              u *= pols[j][prism2tensor[i][j]];
            for (int j = split[i]; j < 3; j++)
              v *= pols[j][prism2tensor[i][j]];

            AutoDiff<3> w = polz[prism2tensor[i][3]];

            for (int j = 0; j < 2; j++)
              shape(i,j) = factorxy[i] * (u.Value() * v.DValue(j) - u.DValue(j) * v.Value()) * uu.Value() * w.Value();
            shape(i,2) = factorz[i] * uu.Value() * u.Value() * v.Value() * w.DValue(2);
          }

        return;
      }
    
    throw Exception ("H(curl) prism does not match");
  }


  
  void HCurlHighOrderPrismTP :: CalcCurlShape (const IntegrationPoint & ip, 
                                             FlatMatrixFixWidth<3> shape) const
  {
    HCurlFiniteElementD<3> :: CalcCurlShape (ip, shape);
    return;
  }







}
