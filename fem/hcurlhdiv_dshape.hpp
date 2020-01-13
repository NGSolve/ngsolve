#ifndef FILE_HCURLHDIV_DSHAPE
#define FILE_HCURLHDIV_DSHAPE


namespace ngfem
{

  /** calculates [du1/dx1 du2/dx1 (du3/dx1) du1/dx2 du2/dx2 (du3/dx2) (du1/dx3 du2/dx3 du3/dx3)] */
  template<typename FEL, int DIMSPACE, int DIM, int DIM_STRESS>
  void CalcDShapeFE(const FEL & fel, const MappedIntegrationPoint<DIM,DIMSPACE>& sip, SliceMatrix<> bmatu, LocalHeap& lh, double eps = 1e-4){
    HeapReset hr(lh);
    // bmatu = 0;
    // evaluate dshape by numerical diff
    //fel, eltrans, sip, returnval, lh
    int nd_u = fel.GetNDof();
    const IntegrationPoint& ip = sip.IP();//volume_ir[i];
    const ElementTransformation & eltrans = sip.GetTransformation();
    FlatMatrixFixWidth<DIM_STRESS> shape_ul(nd_u, lh);
    FlatMatrixFixWidth<DIM_STRESS> shape_ur(nd_u, lh);
    FlatMatrixFixWidth<DIM_STRESS> shape_ull(nd_u, lh);
    FlatMatrixFixWidth<DIM_STRESS> shape_urr(nd_u, lh);
    FlatMatrixFixWidth<DIM_STRESS> dshape_u_ref(nd_u, lh);//(shape_ur); ///saves "reserved lh-memory"

    FlatMatrixFixWidth<DIM> dshape_u_ref_comp(nd_u, lh);
    FlatMatrixFixWidth<DIMSPACE> dshape_u(nd_u, lh);//(shape_ul);///saves "reserved lh-memory"
    
    for (int j = 0; j < DIM; j++)   // d / dxj
      {
        IntegrationPoint ipl(ip);
        ipl(j) -= eps;
        IntegrationPoint ipr(ip);
        ipr(j) += eps;
        IntegrationPoint ipll(ip);
        ipll(j) -= 2*eps;
        IntegrationPoint iprr(ip);
        iprr(j) += 2*eps;
        
        MappedIntegrationPoint<DIM,DIMSPACE> sipl(ipl, eltrans);
        MappedIntegrationPoint<DIM,DIMSPACE> sipr(ipr, eltrans);
        MappedIntegrationPoint<DIM,DIMSPACE> sipll(ipll, eltrans);
        MappedIntegrationPoint<DIM,DIMSPACE> siprr(iprr, eltrans);
        
        fel.CalcMappedShape (sipl, shape_ul);
        fel.CalcMappedShape (sipr, shape_ur);
        fel.CalcMappedShape (sipll, shape_ull);
        fel.CalcMappedShape (siprr, shape_urr);
        
        dshape_u_ref = (1.0/(12.0*eps)) * (8.0*shape_ur-8.0*shape_ul-shape_urr+shape_ull);
        
        // dshape_u_ref = (1.0/(2*eps)) * (shape_ur-shape_ul);
        // dshape_u_ref = (1.0/(4*eps)) * (shape_urr-shape_ull);
        
        /*
	  for (int k = 0; k < nd_u; k++)
          for (int l = 0; l < D; l++)
          bmatu(k, j*D+l) = dshape_u_ref(k,l);
        */
        for (int l = 0; l < DIM_STRESS; l++)
          bmatu.Col(j*DIM_STRESS+l) = dshape_u_ref.Col(l);
      }
    
    for (int j = 0; j < DIM_STRESS; j++)
      {
        for (int k = 0; k < nd_u; k++)
          for (int l = 0; l < DIM; l++)
            dshape_u_ref_comp(k,l) = bmatu(k, l*DIM_STRESS+j);
        
        dshape_u = dshape_u_ref_comp * sip.GetJacobianInverse();
        
        for (int k = 0; k < nd_u; k++)
          for (int l = 0; l < DIMSPACE; l++)
            bmatu(k, l*DIM_STRESS+j) = dshape_u(k,l);
      }
  }
}

#endif 
