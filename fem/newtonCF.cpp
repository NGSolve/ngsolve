/*********************************************************************/
/* File:   newtonCF.cpp                                              */
/* Authors: Joachim Schoeberl, Matthias Rambausek                    */
/* Date:   Feb 2021                                                  */
/*********************************************************************/

#include <fem.hpp>
#include <limits>
#include <cmath>

namespace ngfem
{
    
  namespace 
  {
    template<typename T
    LInfNorm(const T& array) -> auto
    {
      std::remove_cv_t<decltype(*begin(array))> s = 0;
      for (auto item : array) 
      {
        s = std::abs(item) > s ? std::abs(item) : s;
      }
      return s;
    }
    
    
    bool has_vs_embedding(const ProxyFunction * proxy) {
      return proxy->Evaluator()->GetVSEmbedding();
    }
    
    auto get_vs_embedding(const ProxyFunction * proxy) {
      return proxy->Evaluator()->GetVSEmbedding();
    }
    
    int proxy_dof_dimension(const ProxyFunction * proxy) {
      if (const auto vsemb = get_vs_embedding(proxy))
        return vsemb->Width();
      else
        return proxy->Dimension();
    }
  }
  

  class NewtonCF : public CoefficientFunction
  {
    shared_ptr<CoefficientFunction> expression;
    shared_ptr<CoefficientFunction> startingpoint;

    FlatArray<ProxyFunction *> proxies;
    Array<CoefficientFunction*> cachecf;
    
    // The total dimension of the linear system (because of VS embeddings, this can be different from this->Dimension())
    int numeric_dim = 0;
    
    // Same parameters as scipy's newton
    // Alternatively, could one think of ParameterCFs here?
    double tol{1e-8};
    double rtol{0.0};
    int maxiter{10};
    
  public:
    NewtonCF (shared_ptr<CoefficientFunction> aexpression,
                    shared_ptr<CoefficientFunction> astartingpoint,
                    std::optional<double> atol, 
                    std::optional<double> artol, 
                    std::optional<int> amaxiter)
      : expression(aexpression), startingpoint(astartingpoint)
    {

      // NOTES: 
        
      // TODO: INTERPOLATION INTO GENERIC COMPOUND SPACES DOES NOT WORK CURRENTLY (only first component is respected)
        
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
      
      // NOTE: GFs on generic CompoundSpaces do not provide useful/usable dimension data!
        
      expression->TraverseTree
        ( [&] (CoefficientFunction & nodecf)
          {
            auto nodeproxy = dynamic_cast<ProxyFunction*> (&nodecf);
            if (nodeproxy) 
              {
                if (!nodeproxy->IsTestFunction())
                  {
                    if (std::find(proxies.begin(), proxies.end(), nodeproxy) == proxies.end())
                      proxies.Append(nodeproxy);
                  }
              }
            else if (nodecf.StoreUserData() && !cachecf.Contains(&nodecf))
              cachecf.Append (&nodecf);
          });
      if (proxies.Size() == 0) 
        throw Exception("NewtonCF: don't have a proxy");
      
      if (proxies.Size() != startingpoint.GetNComponents()) 
        throw Exception("NewtonCF: number of proxies does not match the number of components of the startingpoint");
      
      
      Array<int> dims;
      
      if (proxies.Size() == 1)
        {  
          if (proxy->Dimensions() != stpt_comp->Dimensions())
            throw Exception(std::string("NewtonCF: Dimensions of proxy and startingpoint do not agree");
          
          dims = proxy->Dimensions();
        }
      else
        { 
          // Check whether all proxies belong to a compound FE space and put them in order
          FlatArray<ProxyFunction*> sorted_proxies(proxies.Size(), nullptr);
          FlatArray<DifferentialOperator*> diffops(proxies.Size(), nullptr);
          for (const auto proxy : proxies)
            {
              const auto evaluator = dynamic_cast<const CompoundDifferentialOperator*>(proxy->Evaluator().get());
              if (!evaluator)
                {
                  throw Exception("NewtonCF: More than one proxy has been found but not all proxy evaluators are of type CompoundDifferentialOperator");
                }
              else
                {
                  if (sorted_proxies[evaluator->Component()] || std::find(diffops.begin(), diffops.end(), evaluator) != diffops.end())
                    throw Exception("NewtonCF: A proxy evaluator (component) has been detected twice");
                  diffops[evaluator->Component()] = evaluator;
                  sorted_proxies[evaluator->Component()] = proxy;
                }
            }
          // Copy over...
          std::copy(sorted_proxies.begin(), sorted_proxies.end(), proxies.begin());
          
          // Process dimensions
          for (auto comp = startingpoint.GetNComponents(); comp > 0; --comp) {
            const auto stpt_comp = startingpoint->GetComponent(comp);
            const auto proxy = proxies[comp];
            
            if (proxy->Dimensions() != stpt_comp->Dimensions())
              throw Exception(std::string("NewtonCF: Dimensions of component ") + std::to_string(comp) + " do not agree");
            
            // TODO: Does this make sense or shall we just not set dimensions in case of generic compound spaces/multiple proxies
            // Should probably be consistent with Compound GFs/CFs
            const auto pdims = proxy->Dimensions();
            for (auto dim : pdims)
               dims.Append(dim);
          }
        }
      
      SetDimensions(dims);

      for (const auto proxy : proxies)
        numeric_dim += proxy_dof_dimension(proxy);
      
      // Process options
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
      // static Timer t("NewtonCF::Eval", 2);
      // static Timer t1("NewtonCF::Eval get Jac", 2);
      // static Timer t2("NewtonCF::Eval solve", 2);
      // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
      // RegionTracer regtr(TaskManager::GetThreadId(), t);
        
      // cout << "eval minimization" << endl;
      LocalHeap lh(1000000);
      
      // startingpoint -> Evaluate (mir, values);
      // cout << "starting: " << endl << values.AddSize(Dimension(), ir.Size()) << endl;

      const ElementTransformation & trafo = mir.GetTransformation();
      
      ProxyUserData ud(1, cachecf.Size(), lh);
      for (CoefficientFunction * cf : cachecf)
        ud.AssignMemory (cf, mir.Size(), cf->Dimension(), lh);
      
      const_cast<ElementTransformation&>(trafo).userdata = &ud;
      
      for (ProxyFunction *  proxy : proxies)
        ud.AssignMemory (proxy, mir.Size(), proxy->Dimension(), lh);

      
      // Prepare data structures for blocks
      const auto nblocks = proxies.Size();
      FlatArray<FlatMatrix<double>> xk_blocks(nblocks);
      FlatArray<FlatMatrix<double>> w_blocks(nblocks);
      FlatArray<FlatMatrix<double>> xold_blocks(nblocks);
      FlatArray<FlatMatrix<double>> res_blocks(nblocks);
      FlatArray<FlatTensor<3>> lin_blocks(nblocks * nblocks);
      
      // These are only "independent" for blocks having "vsemb"; otherwise just views
      FlatArray<FlatMatrix<double>> rhs_blocks(nblocks);
      FlatArray<FlatTensor<3>> lhs_blocks(nblocks * nblocks);
      
      for (int i : Range(proxies))
        {
          const auto proxy = proxies[i];
          xk_blocks[i] = ud.GetMemory(proxy);
          xold_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
          w_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
          res_blocks[i].AssignMemory(mir.Size(), proxy->Dimension(), lh);
          
          if (has_vs_embedding(proxy))
            rhs_blocks[i].AssignMemory(mir.Size(), proxy_dof_dimension(proxy), lh);
          else
            rhs_blocks[i].AssignMemory(mir.Size(), proxy_dof_dimension(proxy), res_blocks[i].Data());
          
          for (int j = : Range(proxies))
            {
              //TODO: IMPORTANT convention: row-major/column-major?
              const auto ij = j * nblocks + i;
              lin_blocks[ij].AssignMemory(lh, mir.Size(), proxies[i]->Dimension(), proxies[j]->Dimension());
              if (has_vs_embedding(proxy))
                lhs_blocks[ij].AssignMemory(lh, mir.Size(), proxy_dof_dimension(proxies[j]), proxy_dof_dimension(proxies[i]));
              else
                lhs_blocks[ij].AssignMemory(in_blocks[ij].Data(), mir.Size(), proxy_dof_dimension(proxies[j]), proxy_dof_dimension(proxies[i]));
            }
        }
      
      // Block-agnostic data structures
      FlatMatrix<> deriv(mir.Size(), 1,lh);
      FlatMatrix<AutoDiff<1,double>> dval(mir.Size(), 1, lh);
      
      // TODO: Why? FlatMatrix has column-major layout, thus it is cheaper to take cols instead of rows to access qp data.!
      FlatMatrix<> xk(mir.Size(), Dimension(), lh);
      FlatMatrix<> xold(mir.Size(), Dimension(), lh);
      FlatMatrix<> w(mir.Size(), Dimension(), lh);
      
      FlatVector<> rhs(numeric_dim, lh);
      FlatVector<> sol(numeric_dim, lh);
      FlatMatrix<> lhs(numeric_dim, numeric_dim, lh);
      
      
      const auto converged = [&](const auto& rhs_vec) {
        const auto res = LInfNorm(res);
        return res <= tol || (res / res_0) <= rtol;
      };
      
      const auto calc_residuals = [&]() mutable -> void {
          // RegionTracer regtr1(TaskManager::GetThreadId(), t1);
          for (int block : Range(proxies))
            {
              auto proxy = proxies[block];
              auto& res = res_blocks[block];
              
              for (int k = 0; k < proxy->Dimension(); k++)
                {
                  ud.trialfunction = proxy;
                  ud.trial_comp = k;
                  ud.testfunction = proxy;
                  ud.test_comp = k;
                  expression -> Evaluate (mir, dval);
                  for (size_t qi = 0; qi < mir.Size(); qi++)
                    res(qi,k) = dval(qi,0).Value(0);
                }
              cout << "res block " << block << " = " << res << endl;
            }
      };
      
      const auto calc_linearizations = [&]() mutable -> void {
          // RegionTracer regtr2(TaskManager::GetThreadId(), t2);
          for (int block1 = 0; block1 < nblocks; ++block1)
            for (int block2 = 0; block2 < nblocks; ++block2) 
              { 
                auto proxy1 = proxies[block1];
                auto proxy2 = proxies[block2];
                auto& lin = lin_blocks[block2 * nblocks + block1];
                
                for (int k = 0; k < proxy1->Dimension(); ++k)
                  for (int l = 0; l < proxy2->Dimension(); ++l)
                    {
                      // TODO: row/col major? -> implementation is col-major, no?
                      // Is this ordering compatible with the blocks in this implementation?
                      //  -> symbolicintegrator.cpp:4656ff for block handling
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = l;
                      expression -> Evaluate (mir, dval);
                      for (size_t qi = 0; qi < mir.Size(); qi++) {
                        deriv(qi,0) = dval(qi,0).DValue(0);
                      }
                      lin(STAR,l,k) = deriv.Col(0);
                    }
                cout << "lin block (" << block2 << ", " << block1 << ") = " << lin << endl;
            }
      };
      
      const auto compute_increments = [&]() mutable -> void {
          for (size_t qi = 0; qi < mir.Size(); qi++)
          {
              
            // TODO: when to skip something because of convergence?
            //  -> the current approach assumes that evaluation of "expression" for all qpoints at once is beneficial!
            
            int offset1 = 0;
            for (int block1 : Range(proxies))
              {
                const auto proxy1 = proxies[block1];
                const auto& resb = res_blocks[block1].Row(qi);
                const auto& rhsb = rhs_blocks[block1].Row(qi);
                if (auto vsemb = get_vs_embedding(proxy1); vsemb)
                  {
                    rhsb = Trans(vsemb) * resb;
                  }
                else
                  {  
                    // nothing to do
                  }
                  
                for (int k : Range(rhsb.Size()))
                  rhs[offset1 + k] = rhsb[k];
                
                int offset2 = 0
                for (int block2 : Range(proxies))
                  {
                    const auto proxy2 = proxies[block2];
                    const auto& linb = lin_blocks[block2 * nblocks + block1](qi);
                    const auto& lhsb = lhs_blocks[block2 * nblocks + block1](qi);
                    if (auto vsemb1 = get_vs_embedding(proxy1); vsemb1)
                      if (auto vsemb2 = get_vs_embedding(proxy2); vsemb2)
                        {
                          lhsb = Trans(vsemb1) * linb * vsemb2;
                        }
                      else
                        {
                          lhsb = Trans(vsemb1) * linb;
                        }
                    else
                      if (auto vsemb2 = get_vs_embedding(proxy2); vsemb2)
                        {
                          lhsb = linb * vsemb2;
                        }
                      else
                        {
                          // nothing to do
                        }
                
                    for (int k : Range(lhsb.Width()))
                      for (int l : Range(lhsb.Height()))
                        lhs(offset2 + l, offset1 + k) = lhsb[l, k];
                        
                    offset2 += lhsb.Height();
                  }
                offset1 += rhsb.Size();
              }
            
            if(converged(rhs))
              {
                w.Row(qi) = 0;
                continue;
              }

            CalcInverse (lhs);
            sol = lhs * rhs;
            
            // handle vs embedding
            // TODO: Better go for a layout that uses Col instead of Row?
            auto& wi = w.Row(qi);
            int offset_w = 0;
            int offset_sol = 0;
            for (int block : Range(proxies))
              {
                const auto proxy = proxies[block];
                if (const auto vsemb = get_vs_embedding(proxy); vsemb)
                  wi.Range(offset_w + proxy->Dimension()) = vsemb * sol.Range(offset_sol);
                offset_w += proxy->Dimension();
                offset_sol += proxy_dof_dimension(proxy);
              }
            }

          }
      };
      
      const auto distribute_xk_to_blocks = [&]() mutable -> void {
        // TODO: Better got for a layout that uses Col instead of Row?
        for (size_t qi = 0; qi < mir.Size(); ++qi) 
          {
            auto xk_qi = xk.Row(qi);
            int offset = 0;
            for (int block : Range(proxies))
              {
                auto xkb_qi = xk_blocks[block].Row(qi);
                // Q: Does this assign values?
                xkb_qi = xk_qi.Range(offset, offset + xkb_qi.Size());
                offset += xkb_qi.Size();
              }
          }
      };
    
      
      // Evaluate starting point
      startingpoint.GetComponent(i) -> Evaluate (mir, xk);
      xold = xk;
      
      // Distribute starting point to proxies
      distribute_xk_to_blocks();
        
      cout << "starting value = " << xk << endl;
      cout << "blocks:" << endl;
      for (int i : Range(proxies))
          cout << i << " -> " << ud.GetMemory(proxies[i]) << endl;
      cout << endl;
      
      // The actual algorithm
      for (int step = 0; step < maxiter; step++)
        {            
          calc_residuals();
          
          if(converged(res))
              break;
          
          calc_linearizations();
          compute_increments();

          xk = xold - w;
          distribute_xk_to_blocks();
        }
        
      if (!converged(res))
          xk = std::numeric_limits<double>::quiet_NaN();

      // cout << "result = " << xk << endl;
      values.AddSize(mir.Size(), Dimension()) = xk;
    }
    
  };
    
    
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
      

      double energy = 0;
      
      const auto converged = [&](const auto& rhs_vec){
        const auto res = L2Norm(res);
        return res <= tol || (res / res_0) <= rtol;
      };
      
      
      const auto calc_energy_rhs_and_diags = [&]() mutable -> void {
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
      };
      
      const auto calc_off_diagonals = [&]() mutable -> void {
          // TODO: exploit symmetry
          // RegionTracer regtr2(TaskManager::GetThreadId(), t2);
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                {
                  expression -> Evaluate (mir, ddval);
                  for (size_t i = 0; i < mir.Size(); i++) {
                    dderiv(i,0) = ddval(i,0).DDValue(0);
                  }
                }
                proxyvalues(STAR,l,k) = dderiv.Col(0);

                if (proxy1 != proxy2 || k != l)  // computed mixed second derivatives
                  {
                    proxyvalues(STAR,l,k) -= diags.Col(k);
                    proxyvalues(STAR,l,k) -= diags.Col(l);
                    proxyvalues(STAR,l,k) *= 0.5;
                  }
              }
      };
      
      const auto compute_increments = [&]() mutable -> void {
        for (size_t i = 0; i < mir.Size(); i++) 
          {
            
            
            if(converged(dWdB.Row(i)))
              {
                w.Row(i) = 0;
                continue;
              }
              
            rhs = dWdB.Row(i);
            mat = proxyvalues(i, STAR,STAR);
            if (vsemb) {
              proj_rhs = Trans(vsemb) * rhs;
              proj_mat = Trans(vsemb) * mat * vsemb;
              CalcInverse (proj_mat);
              proj_sol = proj_mat * proj_rhs;
              w.Row(i) = vsemb * proj_sol;
            }
            else {
              CalcInverse (mat);
              w.Row(i) = mat * rhs;
            }
          }
      };
      
      
      const auto linesearch = [&]() mutable -> void {// linesearch
          xold = xk;
          double alpha = 1;
          double newenergy = energy + 1;
          
          ud.trialfunction = proxy;
          ud.trial_comp = 0;
          ud.testfunction = proxy;
          ud.test_comp = 0;
          
          // cout << "w = " << endl << w << endl;
          while (newenergy > energy && alpha > 1e-10)
            {
              xk = xold - alpha * w;

              newenergy = 0;
              expression -> Evaluate (mir, ddval);
              for (size_t i = 0; i < mir.Size(); i++) {
                  newenergy +=  mir[i].GetWeight() * ddval(i,0).Value();
              }

              // cout << "alpha = " << alpha << ", newen = " << newenergy << endl;
              alpha /= 2;
            }
      };
    
      // The actual algorithm
      for (int step = 0; step < maxiter; step++)
        {            
          calc_energy_rhs_and_diags();
          if(converged(dWdB))
              break;
          
          calc_off_diagonals();
          compute_increments();
          linesearch();
        }
        
      if (!converged(dWdB))
          xk = std::numeric_limits<double>::quiet_NaN();

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
  
  shared_ptr<CoefficientFunction>
  CreateNewtonCF (shared_ptr<CoefficientFunction> expression,
                  shared_ptr<CoefficientFunction> startingpoint,
                  std::optional<double> atol,
                  std::optional<double> rtol,
                  std::optional<int> maxiter,
                 )
  {
    return make_shared<NewtonCF> (expression, startingpoint, atol, rtol, maxiter);
  }
}
