#include <fem.hpp>
#include "l2hofetp.hpp"
#include "../fem/tscalarfe_impl.hpp"

namespace ngfem
{

  L2HighOrderFETP<ET_QUAD> :: ~L2HighOrderFETP() { ; }

  


  
  void L2HighOrderFETP<ET_QUAD> ::
  Evaluate (const SIMD_IntegrationRule & ir,
            BareSliceVector<> bcoefs,
            BareVector<SIMD<double>> values) const
  {
    static Timer t("quad evaluate");
    ThreadRegionTimer reg(t, TaskManager::GetThreadId());

    if (ir.IsTP())
      {
        static Timer tcopy("quad mult copy");

        double facx[] = { -1, 1, 1, -1 };
        double facy[] = { -1, -1, 1, 1 };
        INT<4> f = GetFaceSort (0, vnums);
        double fx = facx[f[0]];
        double fy = facy[f[0]];
        bool flip = (facx[f[0]] == facx[f[1]]);
            
        auto & irx = ir.GetIRX();
        auto & iry = ir.GetIRY();
        size_t nipx = irx.GetNIP();
        size_t nipy = iry.GetNIP();

        NgProfiler::StartThreadTimer (tcopy, TaskManager::GetThreadId());
        bool needs_copy = bcoefs.Dist() != 1;
        STACK_ARRAY(double, mem_coefs, needs_copy ? (order+1)*(order+1) : 0);
        if (needs_copy)
          {
            FlatVector<> coefs(sqr(order+1), mem_coefs);
            coefs = bcoefs;
          }
        FlatMatrix<> mat_coefs(order+1, order+1, needs_copy ? mem_coefs : &bcoefs(0));
        NgProfiler::StopThreadTimer (tcopy, TaskManager::GetThreadId());                        
          
        if (nipx == 1)
          {
            static Timer t("quad evaluate x");
            static Timer tmult("quad mult x");
            ThreadRegionTimer reg(t, TaskManager::GetThreadId());

            STACK_ARRAY(double, mem_shapex, order+1);
            FlatVector<> shapex(order+1, mem_shapex);
            LegendrePolynomial (order, fx*(2*irx[0](0)[0]-1), shapex);
              
            STACK_ARRAY(SIMD<double>, mem_shapey, (order+1)*iry.Size());
            FlatMatrix<SIMD<double>> simd_shapey(order+1, iry.Size(), mem_shapey);
            SliceMatrix<double> shapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), (double*)&mem_shapey[0]);
            for (size_t i = 0; i < iry.Size(); i++)
              LegendrePolynomial (order, fy*(2*iry[i](0)-1), simd_shapey.Col(i));
              
            FlatVector<> vec_values(nipy, (double*)&values(0));
              
            STACK_ARRAY(double, mem_tmp, order+1);
            FlatVector<> tmp(order+1, mem_tmp);

            ThreadRegionTimer regmult(tmult, TaskManager::GetThreadId());
            NgProfiler::AddThreadFlops (tmult, TaskManager::GetThreadId(), (order+1)*(order+1+nipy));
              
            if (flip)
              tmp = mat_coefs * shapex;
            else
              tmp = Trans(mat_coefs) * shapex;
            vec_values = Trans(shapey)*tmp;
          }
        else if (nipy == 1)
          {
            static Timer t("quad evaluate y");
            static Timer tmult("quad mult y");
            static Timer tleg("quad mult legendre");
            // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
            NgProfiler::StartThreadTimer (t, TaskManager::GetThreadId());
              
            NgProfiler::StartThreadTimer (tleg, TaskManager::GetThreadId());
            STACK_ARRAY(double, mem_shapey, order+1);
            FlatVector<> shapey(order+1, mem_shapey);
            LegendrePolynomial (order, fy*(2*iry[0](0)[0]-1), shapey);
            NgProfiler::StopThreadTimer (tleg, TaskManager::GetThreadId());


            FlatVector<> vec_values(nipx, (double*)&values(0));
              
            STACK_ARRAY(double, mem_tmp, order+1);
            FlatVector<> tmp(order+1, mem_tmp);

            ThreadRegionTimer regmult(tmult, TaskManager::GetThreadId());
            NgProfiler::AddThreadFlops (tmult, TaskManager::GetThreadId(), (order+1)*(order+1+nipx));
              
            if (flip)
              tmp = Trans(mat_coefs) * shapey;
            else
              tmp = mat_coefs * shapey;
            // vec_values = Trans(shapex) * tmp;
            for (size_t i = 0; i < irx.Size(); i++)
              {
                SIMD<double> sum = 0.0;
                LegendrePolynomial (order, fx*(2*irx[i](0)-1),
                                    SBLambda([tmp,&sum] (size_t nr, auto val)
                                             { sum += tmp(nr)*val; }));
                values(i) = sum;
              }
            NgProfiler::StopThreadTimer (t, TaskManager::GetThreadId());              
          }
        else
          {
            static Timer t("quad evaluate xy");
            static Timer tmult("quad mult xy");              
            ThreadRegionTimer reg(t, TaskManager::GetThreadId());

            STACK_ARRAY(SIMD<double>, mem_shapex, (order+1)*irx.Size());
            FlatMatrix<SIMD<double>> simd_shapex(order+1, irx.Size(), mem_shapex);
            SliceMatrix<double> shapex(order+1, nipx, SIMD<double>::Size()*irx.Size(), (double*)&mem_shapex[0]);
            for (size_t i = 0; i < irx.Size(); i++)
              LegendrePolynomial (order, fx*(2*irx[i](0)-1), simd_shapex.Col(i));
              
            STACK_ARRAY(SIMD<double>, mem_shapey, (order+1)*iry.Size());
            FlatMatrix<SIMD<double>> simd_shapey(order+1, iry.Size(), mem_shapey);
            SliceMatrix<double> shapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), (double*)&mem_shapey[0]);          
            for (size_t i = 0; i < iry.Size(); i++)
              LegendrePolynomial (order, fy*(2*iry[i](0)-1), simd_shapey.Col(i));
              
            FlatMatrix<> mat_values(nipx, nipy, (double*)&values(0));
              
            if (nipx <= nipy)
              {
                ThreadRegionTimer regmult(tmult, TaskManager::GetThreadId());
                NgProfiler::AddThreadFlops (tmult, TaskManager::GetThreadId(), nipx*(order+1)*((order+1)+nipy));
                  
                STACK_ARRAY(double, mem_tmp, (order+1)*nipx);
                FlatMatrix<> tmp(order+1, nipx, mem_tmp);
                  
                if (flip)
                  tmp = mat_coefs * shapex;
                else
                  tmp = Trans(mat_coefs) * shapex;
                mat_values = Trans(tmp) * shapey;
              }
            else
              {
                STACK_ARRAY(double, mem_tmp, (order+1)*nipy);
                FlatMatrix<> tmp(order+1, nipy, mem_tmp);

                ThreadRegionTimer regmult(tmult, TaskManager::GetThreadId());
                NgProfiler::AddThreadFlops (tmult, TaskManager::GetThreadId(), nipy*(order+1)*((order+1)+nipx));
                  
                if (flip)
                  tmp = Trans(mat_coefs) * shapey;
                else
                  tmp = mat_coefs * shapey;
                  
                mat_values = Trans(shapex) * tmp;
              }
          }
        return;
      }
      
    TBASE::Evaluate (ir, bcoefs, values);
  }

  void L2HighOrderFETP<ET_QUAD> ::
  AddTrans (const SIMD_IntegrationRule & ir,
            BareVector<SIMD<double>> values,
            BareSliceVector<> bcoefs) const 
  {
    static Timer t("quad AddTrans");
    ThreadRegionTimer reg(t, TaskManager::GetThreadId());

    if (ir.IsTP())
      {
        double facx[] = { -1, 1, 1, -1 };
        double facy[] = { -1, -1, 1, 1 };
        INT<4> f = GetFaceSort (0, vnums);
        double fx = facx[f[0]];
        double fy = facy[f[0]];
        bool flip = (facx[f[0]] == facx[f[1]]);
            
        auto & irx = ir.GetIRX();
        auto & iry = ir.GetIRY();
        size_t nipx = irx.GetNIP();
        size_t nipy = iry.GetNIP();

        bool needs_copy = bcoefs.Dist() != 1;

        STACK_ARRAY(double, mem_coefs, needs_copy ? (order+1)*(order+1) : 0);
        if (needs_copy)
          {
            FlatVector<> coefs(sqr(order+1), mem_coefs);
            coefs = bcoefs;
          }
        FlatMatrix<> mat_coefs(order+1, order+1, needs_copy ? mem_coefs : &bcoefs(0));

        if (nipx == 1)
          {
            static Timer t("quad add trans x");
            static Timer tmult("quad add trans mult x");
            ThreadRegionTimer reg(t, TaskManager::GetThreadId());

            STACK_ARRAY(double, mem_shapex, order+1);
            FlatVector<> shapex(order+1, mem_shapex);
            LegendrePolynomial (order, fx*(2*irx[0](0)[0]-1), shapex);
              
            STACK_ARRAY(SIMD<double>, mem_shapey, (order+1)*iry.Size());
            FlatMatrix<SIMD<double>> simd_shapey(order+1, iry.Size(), mem_shapey);
            SliceMatrix<double> shapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), (double*)&mem_shapey[0]);
            for (size_t i = 0; i < iry.Size(); i++)
              LegendrePolynomial (order, fy*(2*iry[i](0)-1), simd_shapey.Col(i));
              
            FlatVector<> vec_values(nipy, (double*)&values(0));
              
            STACK_ARRAY(double, mem_tmp, order+1);
            FlatVector<> tmp(order+1, mem_tmp);

            ThreadRegionTimer regmult(tmult, TaskManager::GetThreadId());
            NgProfiler::AddThreadFlops (tmult, TaskManager::GetThreadId(), (order+1)*(order+1+nipy));

            tmp = shapey * vec_values;
            if (!flip)
              mat_coefs += shapex * Trans(tmp);
            else
              mat_coefs += tmp * Trans(shapex);
          }
        else if (nipy == 1)
          {
            static Timer t("quad add trans y");
            static Timer tmult("quad add trans mult y");
            ThreadRegionTimer reg(t, TaskManager::GetThreadId());

            STACK_ARRAY(SIMD<double>, mem_shapex, (order+1)*irx.Size());
            FlatMatrix<SIMD<double>> simd_shapex(order+1, irx.Size(), mem_shapex);
            SliceMatrix<double> shapex(order+1, nipx, SIMD<double>::Size()*irx.Size(), (double*)&mem_shapex[0]);
            for (size_t i = 0; i < irx.Size(); i++)
              LegendrePolynomial (order, fx*(2*irx[i](0)-1), simd_shapex.Col(i));
            
            STACK_ARRAY(double, mem_shapey, order+1);
            FlatVector<> shapey(order+1, mem_shapey);
            LegendrePolynomial (order, fy*(2*iry[0](0)[0]-1), shapey);

            FlatVector<> vec_values(nipx, (double*)&values(0));
              
            STACK_ARRAY(double, mem_tmp, order+1);
            FlatVector<> tmp(order+1, mem_tmp);

            ThreadRegionTimer regmult(tmult, TaskManager::GetThreadId());
            NgProfiler::AddThreadFlops (tmult, TaskManager::GetThreadId(), (order+1)*(order+1+nipx));

            tmp = shapex * vec_values;
            if (flip)
              mat_coefs += shapey * Trans(tmp);
            else
              mat_coefs += tmp * Trans(shapey);
          }
        else
          {
            static Timer t("quad add trans xy");
            static Timer tmult("quad add trans mult xy");              
            ThreadRegionTimer reg(t, TaskManager::GetThreadId());

            STACK_ARRAY(SIMD<double>, mem_shapex, (order+1)*irx.Size());
            FlatMatrix<SIMD<double>> simd_shapex(order+1, irx.Size(), mem_shapex);
            SliceMatrix<double> shapex(order+1, nipx, SIMD<double>::Size()*irx.Size(), (double*)&mem_shapex[0]);
            for (size_t i = 0; i < irx.Size(); i++)
              LegendrePolynomial (order, fx*(2*irx[i](0)-1), simd_shapex.Col(i));
              
            STACK_ARRAY(SIMD<double>, mem_shapey, (order+1)*iry.Size());
            FlatMatrix<SIMD<double>> simd_shapey(order+1, iry.Size(), mem_shapey);
            SliceMatrix<double> shapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), (double*)&mem_shapey[0]);    
            for (size_t i = 0; i < iry.Size(); i++)
              LegendrePolynomial (order, fy*(2*iry[i](0)-1), simd_shapey.Col(i));
              
            FlatMatrix<> mat_values(nipx, nipy, (double*)&values(0));
              
            if (nipx <= nipy)
              {
                ThreadRegionTimer regmult(tmult, TaskManager::GetThreadId());
                NgProfiler::AddThreadFlops (tmult, TaskManager::GetThreadId(), nipx*(order+1)*((order+1)+nipy));
                  
                STACK_ARRAY(double, mem_tmp, (order+1)*nipx);
                FlatMatrix<> tmp(order+1, nipx, mem_tmp);

                tmp = shapey * Trans(mat_values);
                if (flip)
                  mat_coefs += tmp * Trans(shapex);
                else
                  mat_coefs += shapex * Trans(tmp);
              }
            else
              {
                STACK_ARRAY(double, mem_tmp, (order+1)*nipy);
                FlatMatrix<> tmp(order+1, nipy, mem_tmp);

                ThreadRegionTimer regmult(tmult, TaskManager::GetThreadId());
                NgProfiler::AddThreadFlops (tmult, TaskManager::GetThreadId(), nipy*(order+1)*((order+1)+nipx));

                /*
                if (flip)
                  tmp = Trans(mat_coefs) * shapey;
                else
                  tmp = mat_coefs * shapey;
                mat_values = Trans(shapex) * tmp;
                */
                tmp = shapex * mat_values;
                if (!flip)
                  mat_coefs += tmp * Trans(shapey);
                else
                  mat_coefs += shapey * Trans(tmp);
                
              }
          }

        if (needs_copy)
          {
            FlatVector<> coefs(sqr(order+1), mem_coefs);
            bcoefs.Range(0,ndof) = coefs;
          }
          
        return;
      }
    
    TBASE::AddTrans (ir, values, bcoefs);
  }


  
  
  // template class L2HighOrderFETP<ET_QUAD>;
  template class T_ScalarFiniteElement<L2HighOrderFETP<ET_QUAD>, ET_QUAD, DGFiniteElement<ET_trait<ET_QUAD>::DIM>>;


  // ********************************** HEX ****************
  
  void L2HighOrderFETP<ET_HEX> ::
  Evaluate (const SIMD_IntegrationRule & ir,
            BareSliceVector<> bcoefs,
            BareVector<SIMD<double>> values) const
  {
    static Timer t("hex evaluate");
    static Timer tmult("hex mult");
    static Timer ttrans("hex transpose");                  
    ThreadRegionTimer reg(t, TaskManager::GetThreadId());

    if (ir.IsTP())
      {
        auto & irx = ir.GetIRX();
        auto & iry = ir.GetIRY();
        auto & irz = ir.GetIRZ();
        size_t nipx = irx.GetNIP();
        size_t nipy = iry.GetNIP();
        size_t nipz = irz.GetNIP();
        size_t nip = nipx*nipy*nipz;
        size_t ndof = (order+1)*(order+1)*(order+1);
        bool needs_copy = bcoefs.Dist() != 1;
        STACK_ARRAY(double, mem_coefs, needs_copy ? ndof : 0);
        if (needs_copy)
          {
            FlatVector<> coefs(ndof, mem_coefs);
            coefs = bcoefs;
          }
        FlatMatrix<> mat_coefs(sqr(order+1), order+1, needs_copy ? mem_coefs : &bcoefs(0));

        STACK_ARRAY(SIMD<double>, mem_shapex, (order+1)*irx.Size());
        FlatMatrix<SIMD<double>> simd_shapex(order+1, irx.Size(), mem_shapex);
        SliceMatrix<double> shapex(order+1, nipx, SIMD<double>::Size()*irx.Size(), (double*)&mem_shapex[0]);
        for (size_t i = 0; i < irx.Size(); i++)
          LegendrePolynomial (order, (2*irx[i](0)-1), simd_shapex.Col(i));
              
        STACK_ARRAY(SIMD<double>, mem_shapey, (order+1)*iry.Size());
        FlatMatrix<SIMD<double>> simd_shapey(order+1, iry.Size(), mem_shapey);
        SliceMatrix<double> shapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), (double*)&mem_shapey[0]); 
        for (size_t i = 0; i < iry.Size(); i++)
          LegendrePolynomial (order, (2*iry[i](0)-1), simd_shapey.Col(i));

        STACK_ARRAY(SIMD<double>, mem_shapez, (order+1)*irz.Size());
        FlatMatrix<SIMD<double>> simd_shapez(order+1, irz.Size(), mem_shapez);
        SliceMatrix<double> shapez(order+1, nipz, SIMD<double>::Size()*irz.Size(), (double*)&mem_shapez[0]);          
        for (size_t i = 0; i < irz.Size(); i++)
          LegendrePolynomial (order, (2*irz[i](0)-1), simd_shapez.Col(i));

        NgProfiler::StartThreadTimer (ttrans, TaskManager::GetThreadId());
        STACK_ARRAY(double, memtshapez, nipz*(order+1));
        FlatMatrix<> tshapez(nipz, order+1, memtshapez);
        STACK_ARRAY(double, memtshapey, nipy*(order+1));
        FlatMatrix<> tshapey(nipy, order+1, memtshapey);
        STACK_ARRAY(double, memtshapex, nipx*(order+1));
        FlatMatrix<> tshapex(nipx, order+1, memtshapex);
        
        tshapez = Trans(shapez);
        tshapey = Trans(shapey);
        tshapex = Trans(shapex);
        NgProfiler::StopThreadTimer (ttrans, TaskManager::GetThreadId());        

        NgProfiler::AddThreadFlops (tmult, TaskManager::GetThreadId(),
                                    nipx*nipy*nipz*(order+1) + nipy*nipz*sqr(order+1) + nipz*ndof);
        NgProfiler::StartThreadTimer (tmult, TaskManager::GetThreadId());                
        
        STACK_ARRAY(double, mem1, nipz*sqr(order+1));
        FlatMatrix<> temp1(nipz, sqr(order+1), mem1);
        // temp1 = Trans(shapez)*Trans(mat_coefs);
        temp1 = tshapez*Trans(mat_coefs);

        FlatMatrix<> temp1reshape(nipz*(order+1), order+1, &temp1(0,0));
        STACK_ARRAY(double, mem2, nipy*nipz*(order+1));
        FlatMatrix<> temp2(nipy, nipz*(order+1), mem2);
        // temp2 = Trans(shapey)*Trans(temp1reshape);
        temp2 = tshapey*Trans(temp1reshape);
        
        FlatMatrix<> temp2reshape(nipz*nipy, order+1, &temp2(0,0));
        STACK_ARRAY(double, mem3, nipx*nipy*nipz);
        FlatMatrix<> temp3(nipx, nipz*nipy, mem3);
        // temp3 = Trans(shapex)*Trans(temp2reshape);
        temp3 = tshapex*Trans(temp2reshape);

        NgProfiler::StopThreadTimer (tmult, TaskManager::GetThreadId());
        
        FlatVector<> vals(nip, (double*)&values(0));
        FlatVector<> vectemp3(nip, &temp3(0,0));
        values(ir.Size()-1) = 0.0; // clear overhead
        vals = vectemp3;
        return;
      }

    TBASE::Evaluate(ir, bcoefs, values);
  }

  void L2HighOrderFETP<ET_HEX> ::
  AddTrans (const SIMD_IntegrationRule & ir,
            BareVector<SIMD<double>> values,
            BareSliceVector<> bcoefs) const 
  {
    static Timer t("hex AddTrans");
    static Timer tx("hex AddTrans x");
    static Timer ty("hex AddTrans y");
    static Timer tz("hex AddTrans z");
    static Timer txyz("hex AddTrans xyz");
    
    static Timer tmult("hex addtrans mult");
    ThreadRegionTimer reg(t, TaskManager::GetThreadId());

    if (ir.IsTP())
      {
        auto & irx = ir.GetIRX();
        auto & iry = ir.GetIRY();
        auto & irz = ir.GetIRZ();
        size_t nipx = irx.GetNIP();
        size_t nipy = iry.GetNIP();
        size_t nipz = irz.GetNIP();
        //size_t nip = nipx*nipy*nipz;
        size_t ndof = (order+1)*(order+1)*(order+1);

        /*
        bool needs_copy = bcoefs.Dist() != 1;
        STACK_ARRAY(double, mem_coefs, needs_copy ? ndof : 0);
        if (needs_copy)
          {
            FlatVector<> coefs(ndof, mem_coefs);
            coefs = bcoefs;
          }
        FlatMatrix<> mat_coefs(sqr(order+1), order+1, needs_copy ? mem_coefs : &bcoefs(0));
        */
        
        STACK_ARRAY(SIMD<double>, mem_shapex, (order+1)*irx.Size());
        FlatMatrix<SIMD<double>> simd_shapex(order+1, irx.Size(), mem_shapex);
        SliceMatrix<double> shapex(order+1, nipx, SIMD<double>::Size()*irx.Size(), (double*)&mem_shapex[0]);
        for (size_t i = 0; i < irx.Size(); i++)
          LegendrePolynomial (order, (2*irx[i](0)-1), simd_shapex.Col(i));
              
        STACK_ARRAY(SIMD<double>, mem_shapey, (order+1)*iry.Size());
        FlatMatrix<SIMD<double>> simd_shapey(order+1, iry.Size(), mem_shapey);
        SliceMatrix<double> shapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), (double*)&mem_shapey[0]);  
        for (size_t i = 0; i < iry.Size(); i++)
          LegendrePolynomial (order, (2*iry[i](0)-1), simd_shapey.Col(i));

        STACK_ARRAY(SIMD<double>, mem_shapez, (order+1)*irz.Size());
        FlatMatrix<SIMD<double>> simd_shapez(order+1, irz.Size(), mem_shapez);
        SliceMatrix<double> shapez(order+1, nipz, SIMD<double>::Size()*irz.Size(), (double*)&mem_shapez[0]);      
        for (size_t i = 0; i < irz.Size(); i++)
          LegendrePolynomial (order, (2*irz[i](0)-1), simd_shapez.Col(i));

        STACK_ARRAY(double, memtshapez, nipz*(order+1));
        FlatMatrix<> tshapez(nipz, order+1, memtshapez);
        STACK_ARRAY(double, memtshapey, nipy*(order+1));
        FlatMatrix<> tshapey(nipy, order+1, memtshapey);
        STACK_ARRAY(double, memtshapex, nipx*(order+1));
        FlatMatrix<> tshapex(nipx, order+1, memtshapex);
        
        tshapez = Trans(shapez);
        tshapey = Trans(shapey);
        tshapex = Trans(shapex);

        NgProfiler::AddThreadFlops (tmult, TaskManager::GetThreadId(),
                                    nipx*nipy*nipz*(order+1) + nipy*nipz*sqr(order+1) + nipz*ndof);
        NgProfiler::StartThreadTimer (tmult, TaskManager::GetThreadId());                

        int nr = txyz;
        if (nipx == 1) nr = tx;
        if (nipy == 1) nr = ty;
        if (nipz == 1) nr = tz;
        ThreadRegionTimer reg(nr, TaskManager::GetThreadId());        

        
        STACK_ARRAY(double, mem0, (order+1)*sqr(order+1));
        FlatMatrix<> temp0(sqr(order+1), order+1, mem0);
        STACK_ARRAY(double, mem1, nipz*sqr(order+1));
        FlatMatrix<> temp1(nipz, sqr(order+1), mem1);
        STACK_ARRAY(double, mem2, nipy*nipz*(order+1));
        FlatMatrix<> temp2(nipy, nipz*(order+1), mem2);
        // STACK_ARRAY(double, mem3, nipx*nipy*nipz);
        FlatMatrix<> temp3(nipx, nipz*nipy, (double*)&values(0));
        
        FlatMatrix<> temp2reshape(nipz*nipy, order+1, &temp2(0,0));
        temp2reshape = Trans(temp3)*tshapex;

        FlatMatrix<> temp1reshape(nipz*(order+1), order+1, &temp1(0,0));
        temp1reshape = Trans(temp2)*tshapey;

        temp0 = Trans(temp1)*tshapez;
        NgProfiler::StopThreadTimer (tmult, TaskManager::GetThreadId());

        FlatVector<> vec_coefs(sqr(order+1)*(order+1), &temp0(0));
        bcoefs.Range(0,ndof) += vec_coefs;
        
        return;
      }

    TBASE::AddTrans(ir, values, bcoefs);
  }


  void L2HighOrderFETP<ET_HEX> ::  
  AddGradTrans (const SIMD_BaseMappedIntegrationRule & mir,
                BareSliceMatrix<SIMD<double>> values,
                BareSliceVector<> bcoefs) const
  {
    static Timer t("hex AddGradTrans");
    ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    auto & ir = mir.IR();
    if (ir.IsTP())
      // if (false)
      {
        mir.TransformGradientTrans (values);
        
        auto & irx = ir.GetIRX();
        auto & iry = ir.GetIRY();
        auto & irz = ir.GetIRZ();
        size_t nipx = irx.GetNIP();
        size_t nipy = iry.GetNIP();
        size_t nipz = irz.GetNIP();
        //size_t nip = nipx*nipy*nipz;
        size_t ndof = (order+1)*(order+1)*(order+1);
        
        
        STACK_ARRAY(SIMD<double>, mem_shapex, (order+1)*irx.Size());
        FlatMatrix<SIMD<double>> simd_shapex(order+1, irx.Size(), mem_shapex);
        SliceMatrix<double> shapex(order+1, nipx, SIMD<double>::Size()*irx.Size(), (double*)&mem_shapex[0]);
        STACK_ARRAY(SIMD<double>, mem_dshapex, (order+1)*irx.Size());
        FlatMatrix<SIMD<double>> simd_dshapex(order+1, irx.Size(), mem_dshapex);
        SliceMatrix<double> dshapex(order+1, nipx, SIMD<double>::Size()*irx.Size(), (double*)&mem_dshapex[0]);
        
        for (size_t i = 0; i < irx.Size(); i++)
          {
            AutoDiff<1,SIMD<double>> adx(irx[i](0), 0);
            LegendrePolynomial (order, (2*adx-1),
                                SBLambda([&] (size_t nr, auto val)
                                         {
                                           simd_shapex(nr, i) = val.Value();
                                           simd_dshapex(nr, i) = val.DValue(0);
                                         }));
          }

        STACK_ARRAY(SIMD<double>, mem_shapey, (order+1)*iry.Size());
        FlatMatrix<SIMD<double>> simd_shapey(order+1, iry.Size(), mem_shapey);
        SliceMatrix<double> shapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), (double*)&mem_shapey[0]);
        STACK_ARRAY(SIMD<double>, mem_dshapey, (order+1)*iry.Size());
        FlatMatrix<SIMD<double>> simd_dshapey(order+1, iry.Size(), mem_dshapey);
        SliceMatrix<double> dshapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), (double*)&mem_dshapey[0]);
        
        for (size_t i = 0; i < iry.Size(); i++)
          {
            AutoDiff<1,SIMD<double>> ady(iry[i](0), 0);
            LegendrePolynomial (order, (2*ady-1),
                                SBLambda([&] (size_t nr, auto val)
                                         {
                                           simd_shapey(nr, i) = val.Value();
                                           simd_dshapey(nr, i) = val.DValue(0);
                                         }));
          }

        STACK_ARRAY(SIMD<double>, mem_shapez, (order+1)*irz.Size());
        FlatMatrix<SIMD<double>> simd_shapez(order+1, irz.Size(), mem_shapez);
        SliceMatrix<double> shapez(order+1, nipz, SIMD<double>::Size()*irz.Size(), (double*)&mem_shapez[0]);
        STACK_ARRAY(SIMD<double>, mem_dshapez, (order+1)*irz.Size());
        FlatMatrix<SIMD<double>> simd_dshapez(order+1, irz.Size(), mem_dshapez);
        SliceMatrix<double> dshapez(order+1, nipz, SIMD<double>::Size()*irz.Size(), (double*)&mem_dshapez[0]);
        
        for (size_t i = 0; i < irz.Size(); i++)
          {
            AutoDiff<1,SIMD<double>> adz(irz[i](0), 0);
            LegendrePolynomial (order, (2*adz-1),
                                SBLambda([&] (size_t nr, auto val)
                                         {
                                           simd_shapez(nr, i) = val.Value();
                                           simd_dshapez(nr, i) = val.DValue(0);
                                         }));
          }
        

        for (size_t j = 0; j < 3; j++)
          {
            STACK_ARRAY(double, memtshapez, nipz*(order+1));
            FlatMatrix<> tshapez(nipz, order+1, memtshapez);
            STACK_ARRAY(double, memtshapey, nipy*(order+1));
            FlatMatrix<> tshapey(nipy, order+1, memtshapey);
            STACK_ARRAY(double, memtshapex, nipx*(order+1));
            FlatMatrix<> tshapex(nipx, order+1, memtshapex);

            if (j == 2)
              tshapez = Trans(dshapez);
            else
              tshapez = Trans(shapez);

            if (j == 1)
              tshapey = Trans(dshapey);
            else
              tshapey = Trans(shapey);

            if (j == 0)
              tshapex = Trans(dshapex);
            else
              tshapex = Trans(shapex);
            
            
            STACK_ARRAY(double, mem0, (order+1)*sqr(order+1));
            FlatMatrix<> temp0(sqr(order+1), order+1, mem0);
            STACK_ARRAY(double, mem1, nipz*sqr(order+1));
            FlatMatrix<> temp1(nipz, sqr(order+1), mem1);
            STACK_ARRAY(double, mem2, nipy*nipz*(order+1));
            FlatMatrix<> temp2(nipy, nipz*(order+1), mem2);
            
            FlatMatrix<> temp3(nipx, nipz*nipy, (double*)&values(j,0));
            
            FlatMatrix<> temp2reshape(nipz*nipy, order+1, &temp2(0,0));
            temp2reshape = Trans(temp3)*tshapex;
            
            FlatMatrix<> temp1reshape(nipz*(order+1), order+1, &temp1(0,0));
            temp1reshape = Trans(temp2)*tshapey;
            
            temp0 = Trans(temp1)*tshapez;
            
            FlatVector<> vec_coefs(sqr(order+1)*(order+1), &temp0(0));
            bcoefs.Range(0,ndof) += vec_coefs;
          }
        
        return;
      }

    
    TBASE::AddGradTrans(mir, values, bcoefs);
  }
  

  L2HighOrderFETP<ET_HEX> :: ~L2HighOrderFETP() { ; }   
}


