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
            SliceMatrix<double> shapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), &mem_shapey[0][0]);
            for (size_t i = 0; i < iry.Size(); i++)
              LegendrePolynomial (order, fy*(2*iry[i](0)-1), simd_shapey.Col(i));
              
            FlatVector<> vec_values(nipy, &values(0)[0]);
              
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


            FlatVector<> vec_values(nipx, &values(0)[0]);
              
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
            SliceMatrix<double> shapex(order+1, nipx, SIMD<double>::Size()*irx.Size(), &mem_shapex[0][0]);
            for (size_t i = 0; i < irx.Size(); i++)
              LegendrePolynomial (order, fx*(2*irx[i](0)-1), simd_shapex.Col(i));
              
            STACK_ARRAY(SIMD<double>, mem_shapey, (order+1)*iry.Size());
            FlatMatrix<SIMD<double>> simd_shapey(order+1, iry.Size(), mem_shapey);
            SliceMatrix<double> shapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), &mem_shapey[0][0]);          
            for (size_t i = 0; i < iry.Size(); i++)
              LegendrePolynomial (order, fy*(2*iry[i](0)-1), simd_shapey.Col(i));
              
            FlatMatrix<> mat_values(nipx, nipy, &values(0)[0]);
              
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
            SliceMatrix<double> shapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), &mem_shapey[0][0]);
            for (size_t i = 0; i < iry.Size(); i++)
              LegendrePolynomial (order, fy*(2*iry[i](0)-1), simd_shapey.Col(i));
              
            FlatVector<> vec_values(nipy, &values(0)[0]);
              
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
            SliceMatrix<double> shapex(order+1, nipx, SIMD<double>::Size()*irx.Size(), &mem_shapex[0][0]);
            for (size_t i = 0; i < irx.Size(); i++)
              LegendrePolynomial (order, fx*(2*irx[i](0)-1), simd_shapex.Col(i));
            
            STACK_ARRAY(double, mem_shapey, order+1);
            FlatVector<> shapey(order+1, mem_shapey);
            LegendrePolynomial (order, fy*(2*iry[0](0)[0]-1), shapey);

            FlatVector<> vec_values(nipx, &values(0)[0]);
              
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
            SliceMatrix<double> shapex(order+1, nipx, SIMD<double>::Size()*irx.Size(), &mem_shapex[0][0]);
            for (size_t i = 0; i < irx.Size(); i++)
              LegendrePolynomial (order, fx*(2*irx[i](0)-1), simd_shapex.Col(i));
              
            STACK_ARRAY(SIMD<double>, mem_shapey, (order+1)*iry.Size());
            FlatMatrix<SIMD<double>> simd_shapey(order+1, iry.Size(), mem_shapey);
            SliceMatrix<double> shapey(order+1, nipy, SIMD<double>::Size()*iry.Size(), &mem_shapey[0][0]);          
            for (size_t i = 0; i < iry.Size(); i++)
              LegendrePolynomial (order, fy*(2*iry[i](0)-1), simd_shapey.Col(i));
              
            FlatMatrix<> mat_values(nipx, nipy, &values(0)[0]);
              
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
            bcoefs.AddSize(ndof) = coefs;
          }
          
        return;
      }
    
    TBASE::AddTrans (ir, values, bcoefs);
  }


  
  
  // template class L2HighOrderFETP<ET_QUAD>;
  template class T_ScalarFiniteElement<L2HighOrderFETP<ET_QUAD>, ET_QUAD, DGFiniteElement<ET_trait<ET_QUAD>::DIM>>;
}


