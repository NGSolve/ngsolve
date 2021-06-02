#ifndef FILE_HCURLHDIV_DSHAPE
#define FILE_HCURLHDIV_DSHAPE


namespace ngfem
{

  /** calculates [du1/dx1 du2/dx1 (du3/dx1) du1/dx2 du2/dx2 (du3/dx2) (du1/dx3 du2/dx3 du3/dx3)] */
  template<typename FEL, int DIMSPACE, int DIM, int DIM_STRESS>
  void CalcDShapeFE(const FEL & fel, const MappedIntegrationPoint<DIM,DIMSPACE>& mip, SliceMatrix<> bmatu, LocalHeap& lh, double eps = 1e-4)
  {
    HeapReset hr(lh);
    // bmatu = 0;
    // evaluate dshape by numerical diff
    //fel, eltrans, mip, returnval, lh
    int nd_u = fel.GetNDof();
    const IntegrationPoint& ip = mip.IP();//volume_ir[i];
    const ElementTransformation & eltrans = mip.GetTransformation();
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
        
        MappedIntegrationPoint<DIM,DIMSPACE> mipl(ipl, eltrans);
        MappedIntegrationPoint<DIM,DIMSPACE> mipr(ipr, eltrans);
        MappedIntegrationPoint<DIM,DIMSPACE> mipll(ipll, eltrans);
        MappedIntegrationPoint<DIM,DIMSPACE> miprr(iprr, eltrans);
        
        fel.CalcMappedShape (mipl, shape_ul);
        fel.CalcMappedShape (mipr, shape_ur);
        fel.CalcMappedShape (mipll, shape_ull);
        fel.CalcMappedShape (miprr, shape_urr);
        
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
        
        dshape_u = dshape_u_ref_comp * mip.GetJacobianInverse();
        
        for (int k = 0; k < nd_u; k++)
          for (int l = 0; l < DIMSPACE; l++)
            bmatu(k, l*DIM_STRESS+j) = dshape_u(k,l);
      }
  }


  template<typename FEL, int DIMSPACE, int DIM, int DIM_STRESS, class TVX, class TVY>
  void ApplyDShapeFE(const FEL & fel, const MappedIntegrationPoint<DIM,DIMSPACE>& mip, const TVX & x, TVY & y, LocalHeap& lh, double eps = 1e-4)
  {
    const IntegrationPoint& ip = mip.IP();
    const ElementTransformation & eltrans = mip.GetTransformation();
    Mat<DIM_STRESS,1> shape_ul;
    Mat<DIM_STRESS,1> shape_ur;
    Mat<DIM_STRESS,1> shape_ull;
    Mat<DIM_STRESS,1> shape_urr;
    Mat<DIM_STRESS,1> dshape_u_ref;

    Vec<DIM> dshape_u_ref_comp;
    Vec<DIMSPACE> dshape_u;
    
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
        
        MappedIntegrationPoint<DIM,DIMSPACE> mipl(ipl, eltrans);
        MappedIntegrationPoint<DIM,DIMSPACE> mipr(ipr, eltrans);
        MappedIntegrationPoint<DIM,DIMSPACE> mipll(ipll, eltrans);
        MappedIntegrationPoint<DIM,DIMSPACE> miprr(iprr, eltrans);
        
        fel.EvaluateMappedShape (mipl,  x, shape_ul);
        fel.EvaluateMappedShape (mipr,  x, shape_ur);
        fel.EvaluateMappedShape (mipll, x, shape_ull);
        fel.EvaluateMappedShape (miprr, x, shape_urr);
        
        dshape_u_ref = (1.0/(12.0*eps)) * (8.0*shape_ur-8.0*shape_ul-shape_urr+shape_ull);
        
        for (int l = 0; l < DIM_STRESS; l++)
          y(j*DIM_STRESS+l) = dshape_u_ref(l);
      }
    
    for (int j = 0; j < DIM_STRESS; j++)
      {
        for (int l = 0; l < DIM; l++)
          dshape_u_ref_comp(l) = y(l*DIM_STRESS+j);
        
        dshape_u = Trans(mip.GetJacobianInverse()) * dshape_u_ref_comp;
        
        for (int l = 0; l < DIMSPACE; l++)
          y(l*DIM_STRESS+j) = dshape_u(l);
      }
  }

  template<typename FEL, int DIMSPACE, int DIM, int DIM_STRESS, class TVX, class TVY>
  void ApplyTransDShapeFE(const FEL & fel_u, const MappedIntegrationPoint<DIM,DIMSPACE>& mip, const TVX & x, TVY & by, LocalHeap & lh, double eps = 1e-4)
  {
    typedef typename TVX::TSCAL TSCALX;

    HeapReset hr(lh);
    int nd_u = fel_u.GetNDof();
    FlatMatrixFixWidth<DIM_STRESS*DIMSPACE> bmatu(nd_u,lh);
    
    auto y = by.Range(0, nd_u);
    const IntegrationPoint& ip = mip.IP();
    const ElementTransformation & eltrans = mip.GetTransformation();
    FlatMatrixFixWidth<DIM_STRESS> shape_ul(nd_u, lh);
    FlatMatrixFixWidth<DIM_STRESS> shape_ur(nd_u, lh);
    FlatMatrixFixWidth<DIM_STRESS> shape_ull(nd_u, lh);
    FlatMatrixFixWidth<DIM_STRESS> shape_urr(nd_u, lh);
    FlatMatrixFixWidth<DIM_STRESS> dshape_u_ref(nd_u, lh);
    FlatMatrixFixWidth<DIM_STRESS> dshape_u(nd_u, lh);
    
    FlatMatrix<TSCALX> hx(DIMSPACE,DIM_STRESS,&x(0));
    Mat<DIM,DIM_STRESS,TSCALX> tx = mip.GetJacobianInverse() * hx;
    
    y = 0;
    for (int j = 0; j < DIM; j++)   // d / dxj
      {
        IntegrationPoint ipts[4];
        
        ipts[0] = ip;
        ipts[0](j) -= eps;
        ipts[1] = ip;
        ipts[1](j) += eps;              
        ipts[2] = ip;
        ipts[2](j) -= 2*eps;
        ipts[3] = ip;
        ipts[3](j) += 2*eps;
        
        IntegrationRule ir(4, ipts);
        MappedIntegrationRule<DIM,DIMSPACE> mirl(ir, eltrans, lh);
        
        fel_u.CalcMappedShape (mirl[0], shape_ul);
        fel_u.CalcMappedShape (mirl[1], shape_ur);
        fel_u.CalcMappedShape (mirl[2], shape_ull);
        fel_u.CalcMappedShape (mirl[3], shape_urr);
        
        dshape_u_ref = (1.0/(12.0*eps)) * (8.0*shape_ur-8.0*shape_ul-shape_urr+shape_ull);
        y += dshape_u_ref * tx.Row(j);
      }
  }



  template<typename FEL, int DIMSPACE, int DIM, int DIM_STRESS>
  void CalcSIMDDShapeFE(const FEL & fel, const SIMD_MappedIntegrationRule<DIM,DIMSPACE>& mir, BareSliceMatrix<SIMD<double>> mat, double eps = 1e-4)
    {
      size_t nd_u = fel.GetNDof();

      STACK_ARRAY(SIMD<double>, mem1, 2*DIM_STRESS*nd_u);
      FlatMatrix<SIMD<double>> shape_u_tmp(nd_u*DIM_STRESS, 1, &mem1[0]);

      FlatMatrix<SIMD<double>> dshape_u_ref(nd_u*DIM_STRESS, 1, &mem1[DIM_STRESS*nd_u]);

      LocalHeapMem<10000> lh("diffopgrad-lh");

      auto & ir = mir.IR();
      for (size_t i = 0; i < mir.Size(); i++)
        {
          const SIMD<IntegrationPoint> & ip = ir[i];
          const ElementTransformation & eltrans = mir[i].GetTransformation();

          for (int j = 0; j < DIM; j++)   // d / dxj
            {
              HeapReset hr(lh);
              SIMD<IntegrationPoint> ipts[4];
              ipts[0] = ip;
              ipts[0](j) -= eps;
              ipts[1] = ip;
              ipts[1](j) += eps;              
              ipts[2] = ip;
              ipts[2](j) -= 2*eps;
              ipts[3] = ip;
              ipts[3](j) += 2*eps;

              SIMD_IntegrationRule ir(4, ipts);
              SIMD_MappedIntegrationRule<DIM,DIMSPACE> mirl(ir, eltrans, lh);

              fel.CalcMappedShape (mirl[2], shape_u_tmp);
              dshape_u_ref = 1.0/(12.0*eps) * shape_u_tmp;
              fel.CalcMappedShape (mirl[3], shape_u_tmp);
              dshape_u_ref -= 1.0/(12.0*eps) * shape_u_tmp;
              fel.CalcMappedShape (mirl[0], shape_u_tmp);
              dshape_u_ref -= 8.0/(12.0*eps) * shape_u_tmp;
              fel.CalcMappedShape (mirl[1], shape_u_tmp);
              dshape_u_ref += 8.0/(12.0*eps) * shape_u_tmp;

              // dshape_u_ref =  (8.0*shape_ur-8.0*shape_ul-shape_urr+shape_ull);
              for (size_t l = 0; l < DIM_STRESS; l++)
                for (size_t k = 0; k < nd_u; k++)
                  mat(k*DIM_STRESS*DIM+j*DIM_STRESS+l, i) = dshape_u_ref(k*DIM_STRESS+l, 0);
            }
          
          for (size_t j = 0; j < DIM_STRESS; j++)
            for (size_t k = 0; k < nd_u; k++)
              {
                Vec<DIM,SIMD<double>> dshape_u_ref, dshape_u;
                for (size_t l = 0; l < DIM; l++)
                  dshape_u_ref(l) = mat(k*DIM_STRESS*DIM+l*DIM+j, i);
                
                dshape_u = Trans(mir[i].GetJacobianInverse()) * dshape_u_ref;
                
                for (size_t l = 0; l < DIMSPACE; l++)
                  mat(k*DIM_STRESS*DIMSPACE+l*DIM_STRESS+j, i) = dshape_u(l);
              }
        }
    }
  
  template<typename FEL, int DIMSPACE, int DIM, int DIM_STRESS>
  void ApplySIMDDShapeFE (const FEL & fel_u, const SIMD_BaseMappedIntegrationRule & bmir,
                          BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y, double eps = 1e-4)
  {
    constexpr size_t BS = 64; // number of simd-points
    size_t maxnp = min2(BS, bmir.Size());
    size_t size = (maxnp+1)*SIMD<double>::Size()*500  +  5*DIM_STRESS*BS*sizeof(SIMD<double>);
    STACK_ARRAY(char, data, size);
    LocalHeap lh(data, size);
    
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
    auto & ir = mir.IR();
    const ElementTransformation & trafo = mir.GetTransformation();

    for (int k = 0; k < mir.Size(); k++)
      for (int m = 0; m < DIM_STRESS*DIMSPACE; m++)
        y(m, k) = SIMD<double> (0.0);
    
    for (size_t base = 0; base < ir.Size(); base += BS)
      {
        HeapReset hr(lh);
        size_t num = min2(BS, ir.Size()-base);
        
        FlatMatrix<SIMD<double>> hxl(DIM_STRESS, num, lh);
        FlatMatrix<SIMD<double>> hxr(DIM_STRESS, num, lh);
        FlatMatrix<SIMD<double>> hxll(DIM_STRESS, num, lh);
        FlatMatrix<SIMD<double>> hxrr(DIM_STRESS, num, lh);
        FlatMatrix<SIMD<double>> hx(DIM_STRESS, num, lh);
        
        for (int j = 0; j < DIM; j++)
          {
            // hx = (F^-1 * x).Row(j)
            {
              HeapReset hr(lh);
              SIMD_IntegrationRule irl(num*SIMD<double>::Size(), lh);
              for (int k = 0; k < irl.Size(); k++)
                {
                  irl[k] = ir[base+k];
                  irl[k](j) -= eps;
                }
              SIMD_MappedIntegrationRule<DIM,DIMSPACE> mirl(irl, trafo, lh);
              fel_u.Evaluate (mirl, x, hxl);
            }
            {
              HeapReset hr(lh);
              SIMD_IntegrationRule irr(num*SIMD<double>::Size(), lh);
              for (int k = 0; k < irr.Size(); k++)
                {
                  irr[k] = ir[base+k];              
                  irr[k](j) += eps;
                }
              SIMD_MappedIntegrationRule<DIM,DIMSPACE> mirr(irr, trafo, lh);
              fel_u.Evaluate (mirr, x, hxr);
            }
            {
              HeapReset hr(lh);
              SIMD_IntegrationRule irll(num*SIMD<double>::Size(), lh);
              for (int k = 0; k < irll.Size(); k++)
                {
                  irll[k] = ir[base+k];
                  irll[k](j) -= 2*eps;
                }
              SIMD_MappedIntegrationRule<DIM,DIMSPACE> mirll(irll, trafo, lh);
              fel_u.Evaluate (mirll, x, hxll);
            }
            {
              HeapReset hr(lh);
              SIMD_IntegrationRule irrr(num*SIMD<double>::Size(), lh);
              for (int k = 0; k < irrr.Size(); k++)
                {
                  irrr[k] = ir[base+k];              
                  irrr[k](j) += 2*eps;
                }
              SIMD_MappedIntegrationRule<DIM,DIMSPACE> mirrr(irrr, trafo, lh);
              fel_u.Evaluate (mirrr, x, hxrr);
            }
            // hx = 1.0/(2*eps) * (hxr-hxl);
            hx = 1.0/(12*eps) * (8*hxr-8*hxl-hxrr+hxll);
            for (int k = 0; k < num; k++)
              {
                auto jacinv = mir[base+k].GetJacobianInverse();
                for (int l = 0; l < DIM_STRESS; l++)
                  {
                    for (int m = 0; m < DIMSPACE; m++)
                      y(m*DIM_STRESS+l, base+k) += jacinv(j,m) * hx(l, k);
                  }
              }
          }
      }
  }




  template<typename FEL, int DIMSPACE, int DIM, int DIM_STRESS>
  void AddTransSIMDDShapeFE (const FEL & fel_u, const SIMD_BaseMappedIntegrationRule & bmir,
                             BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y, double eps = 1e-4)
  {
    constexpr size_t BS = 64; // number of simd-points
    size_t maxnp = min2(BS, bmir.Size());
    size_t size = (maxnp+1)*SIMD<double>::Size()*500;
    
    STACK_ARRAY(char, data, size);
    LocalHeap lh(data, size);
    
    auto & mir = static_cast<const SIMD_MappedIntegrationRule<DIM,DIMSPACE>&> (bmir);
    auto & ir = mir.IR();
    const ElementTransformation & trafo = mir.GetTransformation();
    
    for (size_t base = 0; base < ir.Size(); base += BS)
      {
        HeapReset hr(lh);
        size_t num = min2(BS, ir.Size()-base);
        
        FlatMatrix<SIMD<double>> hx1(DIM_STRESS, num, lh);
        FlatMatrix<SIMD<double>> hx2(DIM_STRESS, num, lh);
        
        for (size_t j = 0; j < DIM; j++)
          {
            // hx = (F^-1 * x).Row(j)
            for (size_t k = 0; k < num; k++)
              {
                auto jacinv = mir[base+k].GetJacobianInverse();
                for (int l = 0; l < DIM_STRESS; l++)
                  {
                    SIMD<double> sum = 0;
                    for (int m = 0; m < DIMSPACE; m++)
                      sum += jacinv(j,m) * x(m*DIM_STRESS+l, base+k);
                    
                    hx1(l,k) = (-(8/(12*eps)) * sum).Data();
                    hx2(l,k) = ( (1/(12*eps)) * sum).Data();
                  }
              }
            
            {
              HeapReset hr(lh);
              SIMD_IntegrationRule irl(num*SIMD<double>::Size(), lh);
              for (size_t k = 0; k < irl.Size(); k++)
                {
                  irl[k] = ir[base+k];
                  irl[k](j) -= eps;
                }
              SIMD_MappedIntegrationRule<DIM,DIMSPACE> mirl(irl, trafo, lh);
              fel_u.AddTrans (mirl, hx1, y);
              irl.NothingToDelete();
            }
            {
              HeapReset hr(lh);
              hx1 *= -1;
              SIMD_IntegrationRule irr(num*SIMD<double>::Size(), lh);
              for (int k = 0; k < irr.Size(); k++)
                {
                  irr[k] = ir[base+k];              
                  irr[k](j) += eps;
                }
              SIMD_MappedIntegrationRule<DIM,DIMSPACE> mirr(irr, trafo, lh);
              fel_u.AddTrans (mirr, hx1, y);
            }
            {
              HeapReset hr(lh);
              SIMD_IntegrationRule irl(num*SIMD<double>::Size(), lh);
              for (int k = 0; k < irl.Size(); k++)
                {
                  irl[k] = ir[base+k];
                  irl[k](j) -= 2*eps;
                }
              SIMD_MappedIntegrationRule<DIM,DIMSPACE> mirl(irl, trafo, lh);
              fel_u.AddTrans (mirl, hx2, y);
            }
            {
              HeapReset hr(lh);
              hx2 *= -1;
              SIMD_IntegrationRule irr(num*SIMD<double>::Size(), lh);
              for (int k = 0; k < irr.Size(); k++)
                {
                  irr[k] = ir[base+k];              
                  irr[k](j) += 2*eps;
                }
              SIMD_MappedIntegrationRule<DIM,DIMSPACE> mirr(irr, trafo, lh);
              fel_u.AddTrans (mirr, hx2, y);
            }
          }
      }
  }
}

#endif 
