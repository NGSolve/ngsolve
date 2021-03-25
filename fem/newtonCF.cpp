/*********************************************************************/
/* File:   newtonCF.cpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   Feb 2021                                                  */
/*********************************************************************/

#include <fem.hpp>

namespace ngfem
{

  class MinimizationCF : public CoefficientFunction
  {
    shared_ptr<CoefficientFunction> expression;
    shared_ptr<CoefficientFunction> startingpoint;

    ProxyFunction * proxy = nullptr;
    Array<CoefficientFunction*> cachecf;
  public:
    MinimizationCF (shared_ptr<CoefficientFunction> aexpression,
                    shared_ptr<CoefficientFunction> astartingpoint)
      : expression(aexpression), startingpoint(astartingpoint)
    {
      SetDimensions (startingpoint->Dimensions());
      
      expression->TraverseTree
        ( [&] (CoefficientFunction & nodecf)
          {
            auto nodeproxy = dynamic_cast<ProxyFunction*> (&nodecf);
            if (nodeproxy) 
              {
                if (!nodeproxy->IsTestFunction())
                  {
                    if (proxy && proxy != nodeproxy) throw Exception("MinimizeCF: only one proxy allowed");
                    proxy = nodeproxy;
                  }
              }
            else if (nodecf.StoreUserData() && !cachecf.Contains(&nodecf))
              cachecf.Append (&nodecf);
          });
      if (!proxy) throw Exception("MinimizedCF: don't have a proxy");
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
      // ud.fel = &fel;   // how do I get that ? looks like we don't need it anymore ???
      ud.AssignMemory (proxy, mir.Size(), proxy->Dimension(), lh);

      FlatMatrix xk { ud.GetMemory(proxy) };
      
      startingpoint -> Evaluate (mir, xk);
      // cout << "starting value = " << ud.GetMemory(proxy) << endl;

      
      FlatMatrix<> dderiv(mir.Size(), 1,lh);
      FlatMatrix<AutoDiffDiff<1,double>> ddval(mir.Size(), 1, lh);
      
      FlatMatrix<> diags(mir.Size(), proxy->Dimension(), lh);
      FlatMatrix<> dWdB(mir.Size(), proxy->Dimension(), lh);

      FlatMatrix<> w(mir.Size(), proxy->Dimension(), lh);
      FlatMatrix<> xold(mir.Size(), proxy->Dimension(), lh);

      for (int step = 0; step < 5; step++)
        {
          // cout << "step = " << step; // << ", starting: " << endl << xk << endl;
      
          double energy = 0;

          auto proxy1 = proxy;
          auto proxy2 = proxy;
          FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());

          
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
          if (auto vsemb = proxy->Evaluator()->GetVSEmbedding(); vsemb)
            {
              for (int i = 0; i < mir.Size(); i++)
                {
                  Matrix mat = proxyvalues(i, STAR,STAR);
                  Vector rhs = dWdB.Row(i);

                  Vector proj_rhs = Trans(*vsemb) * rhs;
                  Matrix proj_mat = Trans(*vsemb) * mat * *vsemb;
                  Vector proj_sol(proj_rhs.Size());

                  // cout << "proj rhs = " << proj_rhs << endl;
                  // cout << "proj mat = " << proj_mat << endl;
                  
                  // cout << "mat = " << mat << ", projmat = " << proj_mat << endl;
                  CalcInverse (proj_mat);
                  proj_sol = proj_mat * proj_rhs;
                  // cout << "proj sol = " << proj_sol << endl;
                  
                  w.Row(i) = *vsemb * proj_sol;
                  // cout << "sol = " << w << endl;
                }
            }
          
          else
            
            for (int i = 0; i < mir.Size(); i++)
              {
                Matrix mat = proxyvalues(i, STAR,STAR);
                Vector rhs = dWdB.Row(i);
                Vector sol(rhs.Size());
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
    return make_shared<MinimizationCF> (expression, startingpoint);
  }
}
