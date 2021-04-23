/*********************************************************************/
/* File:   newtonCF.cpp                                              */
/* Authors: Joachim Schoeberl, Matthias Rambausek                    */
/* Date:   Feb 2021                                                  */
/*********************************************************************/

#include <fem.hpp>

namespace ngfem
{

  // TODO: provide also a mere NewtonCF for general nonlinear systems
    
    
  class MinimizationCF : public CoefficientFunction
  {
    shared_ptr<CoefficientFunction> expression;
    shared_ptr<CoefficientFunction> startingpoint;

    std::vector<ProxyFunction *> proxies;
    Array<CoefficientFunction*> cachecf;
    
    // same parameters as scipy's newton
    // Alternatively, could one think of ParameterCFs here?
    double tol{1e-8};
    double rtol{0.0};
    int maxiter{10};
    
  public:
    MinimizationCF (shared_ptr<CoefficientFunction> aexpression,
                    shared_ptr<CoefficientFunction> astartingpoint,
                    std::optional<double> atol, 
                    std::optional<double> artol, 
                    std::optional<int> amaxiter)
      : expression(aexpression), startingpoint(astartingpoint)
    {
      expression->TraverseTree
        ( [&] (CoefficientFunction & nodecf)
          {
            auto nodeproxy = dynamic_cast<ProxyFunction*> (&nodecf);
            if (nodeproxy) 
              {
                if (!nodeproxy->IsTestFunction())
                  {
                    if (std::find(begin(proxies), end(proxies), nodeproxy) == end(proxies))
                      proxies.push_back(nodeproxy);
                  }
              }
            else if (nodecf.StoreUserData() && !cachecf.Contains(&nodecf))
              cachecf.Append (&nodecf);
          });
      if (proxies.empty()) throw Exception("MinimizedCF: don't have a proxy");
      if (proxies.size() > 1) throw Exception("MinimizedCF: only a single proxy is support ATM");
      
      // TODO: support case proxies.size() > 1
      // All proxies must originate from one FE space. However, there is no way to check this.
      // Available information for proxies:
      //  - block index: via Proxy Evaluator if this can be cast to a CompoundDifferentialOperator
      //  - block size: via Proxy Evaluator
      //  - dimension: should be the same as block size?
      //  - dims: array with dimensions for each axis
      
      // For the CF (starting point):
      //  - dim, dims (misleading for GF of CompoundFESpace)
      //  - components (CFs) for GF of CompoundFESpace
      
      // Strategy: 
      //  1. If more that one proxy has been found, sort proxies by block index; only allow each block index appear exactly once.
      //     IOW, two different proxies with for some reason same block index are forbidden.
      //  2. Determine cumulated dimensions of collected proxies.
      //  3. Compare cumulated dimension with dimension of startingpoint.
      //     '-> In case of mismatch, try whether startingpoint has components
      //     '-> In case of multiple proxies and startingpoint being a GF, the latter must have 
      //         components corresponding to the proxies
      //  4. Call SetDimensions with the appropriate information.
      
      // Note: it would be nice, if GF on product space have proper dimensions
      SetDimensions (startingpoint->Dimensions());
      
      if (atol)
        tol = *atol;
      if (artol)
        rtol = *artol;
      if (amaxiter)
        maxiter = *amaxiter;
    }
    

    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      cout << "pw eval not overloaded" << endl;
      return 0;
    }

    void Evaluate (const BaseMappedIntegrationRule & mir,
                   BareSliceMatrix<double> values) const override
    {
      // static Timer t("MinimizationCF::Eval", 2);
      // static Timer t1("MinimizationCF::Eval get Jac", 2);            
      // static Timer t2("MinimizationCF::Eval solve", 2);      
      // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
      // RegionTracer regtr(TaskManager::GetThreadId(), t);    
        
      // cout << "eval minimizatin" << endl;
      LocalHeap lh(1000000);
      
      // startingpoint -> Evaluate (mir, values);
      // cout << "starting: " << endl << values.AddSize(Dimension(), ir.Size()) << endl;

      const ElementTransformation & trafo = mir.GetTransformation();
      
      ProxyUserData ud(1, cachecf.Size(), lh);
      for (CoefficientFunction * cf : cachecf)
        ud.AssignMemory (cf, mir.Size(), cf->Dimension(), lh);
      
      const_cast<ElementTransformation&>(trafo).userdata = &ud;
      
      auto proxy = proxies[0];
      ud.AssignMemory (proxy, mir.Size(), proxy->Dimension(), lh);

      FlatMatrix xk { ud.GetMemory(proxy) };
      
      startingpoint -> Evaluate (mir, xk);
      // cout << "starting value = " << ud.GetMemory(proxy) << endl;

      
      // TODO: support multiple proxies (compound spaces)
      
      FlatMatrix<> dderiv(mir.Size(), 1,lh);
      FlatMatrix<AutoDiffDiff<1,double>> ddval(mir.Size(), 1, lh);
      
      FlatMatrix<> diags(mir.Size(), proxy->Dimension(), lh);
      FlatMatrix<> dWdB(mir.Size(), proxy->Dimension(), lh);

      FlatMatrix<> w(mir.Size(), proxy->Dimension(), lh);
      FlatMatrix<> xold(mir.Size(), proxy->Dimension(), lh);

      auto proxy1 = proxy;
      auto proxy2 = proxy;
      FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
      
      FlatVector<> rhs(proxy2->Dimension(), lh);
      FlatVector<> sol(proxy1->Dimension(), lh);
      FlatMatrix<> mat(proxy2->Dimension(), proxy1->Dimension(), lh);
      
      const auto vsemb = proxy->Evaluator()->GetVSEmbedding();
      const auto proj_dim_1 = vsemb ? vsemb->Width() : 0;
      const auto proj_dim_2 = vsemb ? vsemb->Width() : 0;
      FlatVector<> proj_rhs(proj_dim_2, lh);
      FlatVector<> proj_sol(proj_dim_1, lh);
      FlatMatrix<> proj_mat(proj_dim_2, proj_dim_1, lh);
      
      
      // TODO: use Linf norm
      const auto converged = [&](const auto& rhs_vec){
        const auto res = L2Norm(res);
        return res <= tol || (res / res_0) <= rtol;
      }
      
      // TODO: Algorithm cleanup: reorder steps
      //   1. [per quadrature point] Do not compute (off-diagonal)/invert hessian if already converged
      //   2. Unified algorithm for vsemb true/false; maybe via lambdas returning refs to the
      //      vectors and matrices that hold the respective results
      
      // TODO: exploit symmetry in hessian computation
      
      
      // How to compute norms for tolerance/convergence evaluation?
      for (int step = 0; step < maxiter; step++)
        {
          // cout << "step = " << step; // << ", starting: " << endl << xk << endl;
      
          double energy = 0;
          
          {
            // RegionTracer regtr1(TaskManager::GetThreadId(), t1); 
          for (int k = 0; k < proxy->Dimension(); k++)
            {
              ud.trialfunction = proxy;
              ud.trial_comp = k;
              ud.testfunction = proxy;
              ud.test_comp = k;
              expression -> Evaluate (mir, ddval);
              for (size_t i = 0; i < mir.Size(); i++)
                diags(i,k) = ddval(i,0).DDValue(0);
              for (size_t i = 0; i < mir.Size(); i++)
                dWdB(i,k) = ddval(i,0).DValue(0);

              if (k == 0)
                for (size_t i = 0; i < mir.Size(); i++)
                  energy += mir[i].GetWeight() * ddval(i,0).Value();
            }
          // cout << "energy old = " << energy << endl;
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                {
                  expression -> Evaluate (mir, ddval);
                  for (size_t i = 0; i < mir.Size(); i++)
                    dderiv(i,0) = ddval(i,0).DDValue(0);
                }
                proxyvalues(STAR,l,k) = dderiv.Col(0);
                
                if (proxy1 != proxy2 || k != l)  // computed mixed second derivatives
                  {
                    proxyvalues(STAR,l,k) -= diags.Col(k);
                    proxyvalues(STAR,l,k) -= diags.Col(l);
                    proxyvalues(STAR,l,k) *= 0.5;
                  }
              }
          }
          // cout << "hessions = " << proxyvalues << endl;
          // cout << "gradients = " << dWdB << endl;
          
          // newton step
          // RegionTracer regtr2(TaskManager::GetThreadId(), t2);    
          // cout << "proxy diffop = " <<  typeid(* (proxy->Evaluator())).name() << endl;
          if (vsemb)
            {
              const auto embmat = *vsemb;
              for (int i = 0; i < mir.Size(); i++)
                {
                  rhs = dWdB.Row(i);
                  if (converged(rhs))
                      continue;
                  
                  mat = proxyvalues(i, STAR,STAR);
                  
                  proj_rhs = Trans(embmat) * rhs;
                  proj_mat = Trans(embmat) * mat * embmat;
                  proj_sol(proj_rhs.Size());

                  // cout << "proj rhs = " << proj_rhs << endl;
                  // cout << "proj mat = " << proj_mat << endl;
                  
                  // cout << "mat = " << mat << ", projmat = " << proj_mat << endl;
                  CalcInverse (proj_mat);
                  proj_sol = proj_mat * proj_rhs;
                  // cout << "proj sol = " << proj_sol << endl;
                  
                  w.Row(i) = embmat * proj_sol;
                  // cout << "sol = " << w << endl;
                }
            }
          
          else
            
            for (int i = 0; i < mir.Size(); i++)
              {
                rhs = dWdB.Row(i);
                if (converged(rhs))
                    continue;
                
                mat = proxyvalues(i, STAR,STAR);
                sol(rhs.Size());
                CalcInverse (mat);
                w.Row(i) = mat * rhs;
                // cout << "sol = " << sol << endl;
              }

          xold = xk;
          double alpha = 1;
          // linesearch
          double newenergy = energy + 1;
          // cout << "w = " << endl << w << endl;
          while (newenergy > energy && alpha > 1e-10)
            {
              xk = xold - alpha * w;

              newenergy = 0;
              ud.trialfunction = proxy;
              ud.trial_comp = 0;
              ud.testfunction = proxy;
              ud.test_comp = 0;
              expression -> Evaluate (mir, ddval);
              for (size_t i = 0; i < mir.Size(); i++)
                newenergy +=  mir[i].GetWeight() * ddval(i,0).Value();

              // cout << "alpha = " << alpha << ", newen = " << newenergy << endl;
              alpha /= 2;
            }

          // cout << "energy-dec: " << newenergy-energy << endl;
        }

      // cout << "result = " << xk << endl;
      values.AddSize(mir.Size(), Dimension()) = xk;
    }
    
  };

  


  shared_ptr<CoefficientFunction>
  CreateMinimizationCF (shared_ptr<CoefficientFunction> expression,
                        shared_ptr<CoefficientFunction> startingpoint)
  {
    return make_shared<MinimizationCF> (expression, startingpoint, std::optional<double>{}, std::optional<double>{}, std::optional<int>{});
  }
}
