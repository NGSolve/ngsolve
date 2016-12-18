/*********************************************************************/
/* File:   symbolicintegrator.cpp                                    */
/* Author: Joachim Schoeberl                                         */
/* Date:   August 2015                                               */
/*********************************************************************/
/* 
   Symbolic integrators
*/

#include <fem.hpp>

namespace ngfem
{

  /*
  Array<int> ProxyFunction ::
  Dimensions() const
  {
    int dim = evaluator->Dim();
    int blockdim = evaluator->BlockDim();
    if (blockdim == 1)
      return Array<int> ({dim});
    else
      return Array<int> ({dim/blockdim, blockdim});
  }
  */

  void ProxyFunction ::
  GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    auto dims = Dimensions();

    string header = "\n\
    {flatmatrix} {values};\n\
    ProxyUserData * {ud} = (ProxyUserData*)mir.GetTransformation().userdata;\n\
    {\n\
      if (!{ud})\n\
        throw Exception (\"cannot evaluate ProxyFunction without userdata\");\n\
          ";

    if(!testfunction) {
      header+=
"      if ({ud}->fel) {\n\
          if ({ud}->HasMemory ({this})) {\n";
      if(code.is_simd) {
        header += "auto x = {ud}->GetAMemory ({this});\n";
        header += "{values}.AssignMemory(x.Height(), x.Width(), &x(0,0));\n";
      } else {
        header += "auto x = {ud}->GetMemory ({this});\n";
        header += "{values}.AssignMemory(x.Height(), x.Width(), &x(0,0));\n";
      }
      header+=
"          }\n\
          else\n\
            throw Exception(\"userdata has no memory!\");\n\
      }\n";
    }
    header += "}\n";

    TraverseDimensions( dims, [&](int ind, int i, int j) {
        header += Var("comp", index,i,j).Declare("{scal_type}", 0.0);
        if(!testfunction && code.deriv==2)
        {
          header += "if(( ({ud}->trialfunction == {this}) && ({ud}->trial_comp=="+ToString(ind)+"))\n"+
          " ||  (({ud}->testfunction == {this}) && ({ud}->test_comp=="+ToString(ind)+")))\n";
        }
        else
          header += "if({ud}->{comp_string}=="+ToString(ind)+")\n";
        header += Var("comp", index,i,j).S() + string("{get_component}") + " = 1.0;\n";
    });
    string body = "";

    TraverseDimensions( dims, [&](int ind, int i, int j) {
        body += Var(index, i,j).Declare("{scal_type}", 0.0);
        body += Var(index, i,j).Assign(CodeExpr("0.0"), false);
    });

    if(!testfunction) {
      body += "if ({ud}->fel) {\n";
      TraverseDimensions( dims, [&](int ind, int i, int j) {
          string var = Var(index, i,j);
          if(code.deriv)
            var += ".Value()";
          string values = "{values}";
          if(code.is_simd)
            values += "(" + ToString(ind) + ",i)";
          else
            values += "(i," + ToString(ind) + ")";

          body += var + " = " + values + ";\n";
      });
      if(code.deriv)
        body += "}\n";
      else
        body += "} else ";
    }
    if(!testfunction && code.deriv==2)
      body += "{\n";
    else
      body += "if({ud}->{func_string} == {this}) {\n";
    if(testfunction)
      TraverseDimensions( dims, [&](int ind, int i, int j) {
        if(code.deriv==0) body += Var(index,i,j).Assign( Var("comp", index,i,j), false );
        if(code.deriv==1) body += Var(index,i,j).Assign( Var("comp", index,i,j), false );
        if(code.deriv==2) body += Var(index,i,j).Call("DValue","0").Assign( Var("comp", index,i,j).Call("DValue","0"), false );
      });
    else
      TraverseDimensions( dims, [&](int ind, int i, int j) {
        if(code.deriv==0) body += Var(index,i,j).Assign( Var("comp", index,i,j), false );
        if(code.deriv>=1) body += Var(index,i,j).Call("DValue","0").Assign( Var("comp", index,i,j).Call("DValue","0"), false );
      });
    body += "}\n";

    string func_string = testfunction ? "testfunction" : "trialfunction";
    string comp_string = testfunction ? "test_comp" : "trial_comp";
    std::map<string,string> variables;
    variables["ud"] = "tmp_"+ToString(index)+"_0";
    variables["this"] = "reinterpret_cast<ProxyFunction*>("+ToString(this)+")";
    variables["func_string"] = testfunction ? "testfunction" : "trialfunction";
    variables["comp_string"] = testfunction ? "test_comp" : "trial_comp";
    variables["testfunction"] = ToString(testfunction);

    variables["flatmatrix"] = code.is_simd ? "FlatMatrix<SIMD<double>>" : "FlatMatrix<double>";

    variables["values"] = Var("values", index).S();

    variables["get_component"] = "";
    if(code.deriv==1)
      variables["get_component"] = testfunction ? ".Value()" : ".DValue(0)";
    if(code.deriv==2)
      variables["get_component"] = ".DValue(0)";

    string scal_type = "double";
    if(code.is_simd) scal_type = "SIMD<" + scal_type + ">";
    if(code.deriv==1) scal_type = "AutoDiff<1,"+scal_type+">";
    if(code.deriv==2) scal_type = "AutoDiffDiff<1,"+scal_type+">";
    variables["scal_type"] = scal_type;
    code.header += Code::Map(header, variables);
    code.body += Code::Map(body, variables);
  }
  
  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationPoint & mip,
            FlatVector<> result) const
  {
    ProxyUserData * ud = (ProxyUserData*)mip.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");

    if (!testfunction && ud->fel)
      {
        static bool first = true;
        if (first) cerr << "ProxyFunction::Evaluate (mip) ... should not be here" << endl;
        first = false;
        // if (ud->HasMemory (this))
        // result = ud->GetMemory (this);
        // else
        evaluator->Apply (*ud->fel, mip, *ud->elx, result, *ud->lh);
        return;
      }

    result = 0;
    if (ud->testfunction == this)
      result (ud->test_comp) = 1;
    if (ud->trialfunction == this)
      result (ud->trial_comp) = 1;
  }

  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationRule & mir,
            FlatMatrix<> result) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");
    
    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          result = ud->GetMemory (this);
        else
          evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
        return;
      }

    result = 0;
    if (ud->testfunction == this)
      result.Col(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      result.Col(ud->trial_comp) = 1;
  }

  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationRule & ir,
            FlatMatrix<Complex> result) const
  {
    size_t dim = Dimension();
    STACK_ARRAY(double, hmem, ir.Size()*dim);
    FlatMatrix<> temp(ir.Size(), dim, &hmem[0]);
    Evaluate (ir, temp);
    result = temp;
  }
  

  void ProxyFunction ::
  Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<SIMD<double>> result) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");
    
    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          result.AddSize(Dimension(), mir.Size()) = BareSliceMatrix<SIMD<double>> (ud->GetAMemory (this));
        else
          {
            static bool first = true;
            if (first) cerr << "ProxyFunction::Evaluate(simd_mir) ... should not be here" << endl;
            first = false;
            evaluator->Apply (*ud->fel, mir, *ud->elx, result); // , *ud->lh);
          }
        return;
      }

    result.AddSize(Dimension(), mir.Size()) = 0;
    if (ud->testfunction == this)
      result.Row(ud->test_comp).AddSize(mir.Size()) = 1;
    if (ud->trialfunction == this)
      result.Row(ud->trial_comp).AddSize(mir.Size()) = 1;
  }

  void ProxyFunction ::
  Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<SIMD<Complex>> result) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");
    
    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          result.AddSize(Dimension(), mir.Size()) = BareSliceMatrix<SIMD<double>> (ud->GetAMemory (this));
        else
          {
            throw ExceptionNOSIMD ("ProxyFunction :: Evaluate(SIMD<Complex>) without precomputed trial values");
          }
        return;
      }

    result.AddSize(Dimension(), mir.Size()) = 0;
    if (ud->testfunction == this)
      result.Row(ud->test_comp).AddSize(mir.Size()) = 1;
    if (ud->trialfunction == this)
      result.Row(ud->trial_comp).AddSize(mir.Size()) = 1;
  }

  
  void ProxyFunction ::
  Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
            FlatArray<AFlatMatrix<double>*> input,
            AFlatMatrix<double> result) const
  {
    ProxyFunction::Evaluate (mir, result);
  }

  
  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationPoint & ip,
            FlatVector<Complex> result) const
  {
    STACK_ARRAY(double, hmem, result.Size());
    FlatVector<> result_double(result.Size(), &hmem[0]);
    Evaluate (ip, result_double);
    result = result_double;
  }

  void ProxyFunction ::
  EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                 FlatMatrix<> result,
                 FlatMatrix<> deriv) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    deriv = 0;
    result = 0;

    // static Timer t("ProxyFunction EvaluateDeriv");
    // t.Start();
    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory(this))
          result = ud->GetMemory (this);          
        else
          {
            static bool first = true;
            if (first) cerr << "ProxyFunction::EvaluateDeriv ... should not be here" << endl;
            first = false;
            evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
          }
      }
    // t.Stop();
    
    if (ud->testfunction == this)
      result.Col(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv.Col(ud->trial_comp) = 1;
  }


  void ProxyFunction ::
  EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                  FlatMatrix<> result,
                  FlatMatrix<> deriv,
                  FlatMatrix<> dderiv) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    result = 0;
    deriv = 0;
    dderiv = 0;

    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          result = ud->GetMemory (this);
        else
          {
            static bool first = true;
            if (first) cerr << "ProxyFunction::EvaluateDDeriv ... should not be here" << endl;
            first = false;
            evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
          }
      }
    if (ud->testfunction == this)
      deriv.Col(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv.Col(ud->trial_comp) = 1;
  }
  
  void ProxyFunction ::
  EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & mir, 
                 AFlatMatrix<double> result, AFlatMatrix<double> deriv) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    result = 0;
    deriv = 0;

    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          result = ABareSliceMatrix<> (ud->GetAMemory (this));
        else
          {
            static bool first = true;
            if (first) cerr << "ProxyFunction::EvaluateDDeriv ... should not be here" << endl;
            first = false;
            evaluator->Apply (*ud->fel, mir, *ud->elx, result);
          }
      }
    if (ud->testfunction == this)
      result.Row(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv.Row(ud->trial_comp) = 1;
  }
  
  void ProxyFunction ::
  EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & mir, 
                  AFlatMatrix<double> result, AFlatMatrix<double> deriv,
                  AFlatMatrix<double> dderiv) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction");

    result = 0;
    deriv = 0;
    dderiv = 0;

    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          result = ABareSliceMatrix<> (ud->GetAMemory (this));
        else
          {
            static bool first = true;
            if (first) cerr << "ProxyFunction::EvaluateDDeriv ... should not be here" << endl;
            first = false;
            evaluator->Apply (*ud->fel, mir, *ud->elx, result);
          }
      }
    if (ud->testfunction == this)
      deriv.Row(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv.Row(ud->trial_comp) = 1;
  }
  
  void ProxyFunction ::  
  NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const
  {
    if (!testfunction && ud.fel)
      {
        nonzero = true;
        return;
      }

    nonzero = false;
    if (ud.testfunction == this)
      nonzero(ud.test_comp) = true;
    if (ud.trialfunction == this)
      nonzero(ud.trial_comp) = true;
  }


  SymbolicLinearFormIntegrator ::
  SymbolicLinearFormIntegrator(shared_ptr<CoefficientFunction> acf, VorB avb,
                               bool aelement_boundary)
    : cf(acf), vb(avb), element_boundary(aelement_boundary)
  {
    simd_evaluate = true;
    
    if (cf->Dimension() != 1)
      throw Exception ("SymbolicLFI needs scalar-valued CoefficientFunction");
    cf->TraverseTree
      ([&] (CoefficientFunction & nodecf)
       {
         auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
         if (proxy && !proxies.Contains(proxy))
           proxies.Append (proxy);
       });
  }

  /*
  template <typename SCAL> 
  void SymbolicLinearFormIntegrator ::
  T_CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatVector<SCAL> elvec,
                       LocalHeap & lh) const
  {
    // static Timer t("symbolicLFI - CalcElementVector", 2); RegionTimer reg(t);
    
    HeapReset hr(lh);
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    
    FlatVector<SCAL> elvec1(elvec.Size(), lh);
    
    FlatMatrix<SCAL> values(ir.Size(), 1, lh);
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    
    elvec = 0;
    for (auto proxy : proxies)
      {
        // td.Start();
        FlatMatrix<SCAL> proxyvalues(ir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            
            cf -> Evaluate (mir, values);
            for (int i = 0; i < mir.Size(); i++)
              proxyvalues(i,k) = mir[i].GetWeight() * values(i,0);
          }
        // td.Stop();
        // tb.Start();
        proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
        // tb.Stop();
        elvec += elvec1;
      }
  }
  */

  template <typename SCAL>   
  void SymbolicLinearFormIntegrator ::
  T_CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatVector<SCAL> elvec,
                       LocalHeap & lh) const
  {
    if (element_boundary)
      { // not yet simded
        elvec = 0;
    
        auto eltype = trafo.GetElementType();
        int nfacet = ElementTopology::GetNFacets(eltype);
        
        Facet2ElementTrafo transform(eltype); 
        
        for (int k = 0; k < nfacet; k++)
          {
            HeapReset hr(lh);
            ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
            
            IntegrationRule ir_facet(etfacet, 2*fel.Order());
            IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
            BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
            
            ProxyUserData ud;
            const_cast<ElementTransformation&>(trafo).userdata = &ud;

            mir.ComputeNormalsAndMeasure (eltype, k);
            
            FlatVector<SCAL> elvec1(elvec.Size(), lh);
            FlatMatrix<SCAL> val(mir.Size(), 1,lh);
            for (auto proxy : proxies)
              {
                HeapReset hr(lh);
                FlatMatrix<SCAL> proxyvalues(mir.Size(), proxy->Dimension(), lh);
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    ud.testfunction = proxy;
                    ud.test_comp = k;
                    cf -> Evaluate (mir, val);
                    proxyvalues.Col(k) = val.Col(0);
                  }
                
                for (size_t i = 0; i < mir.Size(); i++)
                  // proxyvalues.Row(i) *= ir_facet[i].Weight() * measure(i);
                  proxyvalues.Row(i) *= ir_facet[i].Weight() * mir[i].GetMeasure();
                
                proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
                elvec += elvec1;
              }
          }
        return;
      }
    if (simd_evaluate)
      {
        try
          {
            // static Timer t("symbolicLFI - CalcElementVector (SIMD)", 2); RegionTimer reg(t);    
            HeapReset hr(lh);
            SIMD_IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
            auto & mir = trafo(ir, lh);
            
            ProxyUserData ud;
            const_cast<ElementTransformation&>(trafo).userdata = &ud;
            
            elvec = 0;
            for (auto proxy : proxies)
              {
                FlatMatrix<SIMD<SCAL>> proxyvalues(proxy->Dimension(), ir.Size(), lh);
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    ud.testfunction = proxy;
                    ud.test_comp = k;
                    
                    cf -> Evaluate (mir, proxyvalues.Rows(k,k+1));
                    for (size_t i = 0; i < mir.Size(); i++)
                      proxyvalues(k,i) *= mir[i].GetWeight();
                  }
                
                proxy->Evaluator()->AddTrans(fel, mir, proxyvalues, elvec);
              }
          }
        catch (ExceptionNOSIMD e)
          {
            cout << IM(4) << e.What() << endl
                 << "switching back to standard evaluation" << endl;
            simd_evaluate = false;
            T_CalcElementVector (fel, trafo, elvec, lh);
          }
      }
    else
      {
        // static Timer t("symbolicLFI - CalcElementVector", 2); RegionTimer reg(t);
        HeapReset hr(lh);
        // IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
        IntegrationRule ir = fel.GetIR(2*fel.Order());
        BaseMappedIntegrationRule & mir = trafo(ir, lh);
        
        FlatVector<SCAL> elvec1(elvec.Size(), lh);
        
        FlatMatrix<SCAL> values(ir.Size(), 1, lh);
        ProxyUserData ud;
        const_cast<ElementTransformation&>(trafo).userdata = &ud;
        
        elvec = 0;
        for (auto proxy : proxies)
          {
            // td.Start();
            FlatMatrix<SCAL> proxyvalues(ir.Size(), proxy->Dimension(), lh);
            for (int k = 0; k < proxy->Dimension(); k++)
              {
                ud.testfunction = proxy;
                ud.test_comp = k;
                
                cf -> Evaluate (mir, values);
                for (int i = 0; i < mir.Size(); i++)
                  proxyvalues(i,k) = mir[i].GetWeight() * values(i,0);
              }
            // td.Stop();
            // tb.Start();
            proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
            // tb.Stop();
            elvec += elvec1;
          }
      }
  }



  void 
  SymbolicLinearFormIntegrator ::
  CalcElementVector (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatVector<double> elvec,
                     LocalHeap & lh) const
  {
    T_CalcElementVector (fel, trafo, elvec, lh);
  }
  
  void 
  SymbolicLinearFormIntegrator ::
  CalcElementVector (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatVector<Complex> elvec,
                     LocalHeap & lh) const
  {
    T_CalcElementVector (fel, trafo, elvec, lh);
  }
  

  

  SymbolicBilinearFormIntegrator ::
  SymbolicBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb,
                                  bool aelement_boundary)
    : cf(acf), vb(avb), element_boundary(aelement_boundary)
  {
    simd_evaluate = true;
    
    if (cf->Dimension() != 1)
        throw Exception ("SymblicBFI needs scalar-valued CoefficientFunction");
    trial_cum.Append(0);
    test_cum.Append(0);    
    cf->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
          if (proxy) 
            {
              if (proxy->IsTestFunction())
                {
                  if (!test_proxies.Contains(proxy))
                    {
                      test_proxies.Append (proxy);
                      test_cum.Append(test_cum.Last()+proxy->Dimension());
                    }
                }
              else
                {                                         
                  if (!trial_proxies.Contains(proxy))
                    {
                      trial_proxies.Append (proxy);
                      trial_cum.Append(trial_cum.Last()+proxy->Dimension());
                    }
                }
            }
        });
    cout << IM(3) << "num test_proxies " << test_proxies.Size() << endl;
    cout << IM(3) << "num trial_proxies " << trial_proxies.Size() << endl;
    cout << IM(3) << "cumulated test_proxy dims  " << test_cum << endl;
    cout << IM(3) << "cumulated trial_proxy dims " << trial_cum << endl;

    elementwise_constant = cf -> ElementwiseConstant();
    cout << IM(3) << "element-wise constant = " << elementwise_constant << endl;

    // find non-zeros
    int cnttest = 0, cnttrial = 0;
    for (auto proxy : trial_proxies)
      cnttrial += proxy->Dimension();
    for (auto proxy : test_proxies)
      cnttest += proxy->Dimension();
    nonzeros = Matrix<bool>(cnttest, cnttrial);

    ProxyUserData ud;
    Vector<bool> nzvec(1);
    int k = 0;
    for (int k1 : test_proxies.Range())
      for (int k2 : Range(0,test_proxies[k1]->Dimension()))
        {
          int l = 0;
          for (int l1 : trial_proxies.Range())
            for (int l2 : Range(0,trial_proxies[l1]->Dimension()))
              {
                ud.trialfunction = trial_proxies[l1];
                ud.trial_comp = l2;
                ud.testfunction = test_proxies[k1];
                ud.test_comp = k2;
                cf -> NonZeroPattern (ud, nzvec);
                nonzeros(k,l) = nzvec(0);
                l++;
              }
          k++;
        }
    
    nonzeros_proxies = Matrix<bool>(test_proxies.Size(), trial_proxies.Size());
    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          
          bool is_nonzero = false;
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              if (nonzeros(test_cum[l1]+l, trial_cum[k1]+k))
                is_nonzero = true;
          nonzeros_proxies(l1, k1) = is_nonzero;
        }
    
    cout << IM(3) << "nonzeros: " << endl << nonzeros << endl;
    cout << IM(3) << "nonzeros_proxies: " << endl << nonzeros_proxies << endl;
  }

  /*
  template <typename T> 
  struct HMatType {
    typedef FlatMatrix<T> tmat;
    typedef FlatVector<T> tvec;
  };
  template <>
  struct HMatType<double> {
    typedef AFlatMatrix<double> tmat;
    typedef AFlatVector<double> tvec;
  };
  */

  IntegrationRule SymbolicBilinearFormIntegrator :: GetIntegrationRule (const FiniteElement & fel) const
  {
    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;

    int trial_difforder = 99, test_difforder = 99;
    for (auto proxy : trial_proxies)
      trial_difforder = min2(trial_difforder, proxy->Evaluator()->DiffOrder());
    for (auto proxy : test_proxies)
      test_difforder = min2(test_difforder, proxy->Evaluator()->DiffOrder());

    int intorder = fel_trial.Order()+fel_test.Order();
    auto et = fel.ElementType();
    if (et == ET_TRIG || et == ET_TET)
      intorder -= test_difforder+trial_difforder;
    return IntegrationRule (et, intorder);
  }
    
  SIMD_IntegrationRule SymbolicBilinearFormIntegrator :: Get_SIMD_IntegrationRule (const FiniteElement & fel) const
  {
    /*
    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;
    */
    bool is_mixed = typeid(fel) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe = static_cast<const MixedFiniteElement*> (&fel);    
    const FiniteElement & fel_trial = is_mixed ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = is_mixed ? mixedfe->FETest() : fel;

    
    int trial_difforder = 99, test_difforder = 99;
    for (auto proxy : trial_proxies)
      trial_difforder = min2(trial_difforder, proxy->Evaluator()->DiffOrder());
    for (auto proxy : test_proxies)
      test_difforder = min2(test_difforder, proxy->Evaluator()->DiffOrder());

    int intorder = fel_trial.Order()+fel_test.Order();
    auto et = fel.ElementType();
    if (et == ET_TRIG || et == ET_TET)
      intorder -= test_difforder+trial_difforder;
    return SIMD_IntegrationRule (et, intorder);
  }



  
  template <typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatMatrix<SCAL> elmat,
                       LocalHeap & lh) const
    
  { 
    static Timer t("symbolicBFI - CalcElementMatrix", 2);
    // static Timer tstart("symboliBFI - CalcElementMatrix startup", 2);
    // static Timer tstart1("symboliBFI - CalcElementMatrix startup 1", 2);
    // static Timer tmain("symboliBFI - CalcElementMatrix main", 2);

    /*
    static Timer td("symboliBFI - CalcElementMatrix dmats", 2);
    static Timer tb("symboliBFI - CalcElementMatrix diffops", 2);
    static Timer tdb("symboliBFI - CalcElementMatrix D * B", 2);
    static Timer tlapack("symboliBFI - CalcElementMatrix lapack", 2);
    */
    
    // tstart.Start();
    
    if (element_boundary)
      {
        switch (trafo.SpaceDim())
          {
          case 1:
            T_CalcElementMatrixEB<1,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          case 2:
            T_CalcElementMatrixEB<2,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          case 3:
            T_CalcElementMatrixEB<3,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          default:
            throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
          }
      }
    
    RegionTimer reg(t);

    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;

    elmat = 0;

    IntegrationRule ir = GetIntegrationRule (fel);
    // IntegrationRule ir = fel_trial.GetIR(intorder);
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    
    
    // tstart.Stop();
    bool symmetric_so_far = true;
    int k1 = 0;
    for (auto proxy1 : trial_proxies)
      {
        int l1 = 0;
        for (auto proxy2 : test_proxies)
          {
            bool is_diagonal = proxy1->Dimension() == proxy2->Dimension();
            bool is_nonzero = false;

            for (int k = 0; k < proxy1->Dimension(); k++)
              for (int l = 0; l < proxy2->Dimension(); l++)
                if (nonzeros(l1+l, k1+k))
                  {
                    if (k != l) is_diagonal = false;
                    is_nonzero = true;
                  }


            if (is_nonzero)
              {
                HeapReset hr(lh);
                bool samediffop = (*(proxy1->Evaluator()) == *(proxy2->Evaluator())) && !mixedfe;
                // td.Start();
                FlatTensor<3,SCAL> proxyvalues(lh, mir.Size(), proxy1->Dimension(), proxy2->Dimension());
                FlatVector<SCAL> diagproxyvalues(mir.Size()*proxy1->Dimension(), lh);
                FlatMatrix<SCAL> val(mir.Size(), 1, lh);
                
                
                if (!is_diagonal)
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    for (int l = 0; l < proxy2->Dimension(); l++)
                      {
                        if (nonzeros(l1+l, k1+k))
                          {
                            if (k != l) is_diagonal = false;
                            is_nonzero = true;
                            ud.trialfunction = proxy1;
                            ud.trial_comp = k;
                            ud.testfunction = proxy2;
                            ud.test_comp = l;
                            
                            cf -> Evaluate (mir, val);
                            proxyvalues(STAR,k,l) = val.Col(0);
                          }
                        else
                          proxyvalues(STAR,k,l) = 0.0;
                      }
                else
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    {
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = k;

                      if (!elementwise_constant)
                        {
                          cf -> Evaluate (mir, val);
                          diagproxyvalues.Slice(k, proxy1->Dimension()) = val.Col(0);
                        }
                      else
                        {
                          cf -> Evaluate (mir[0], val.Row(0));
                          diagproxyvalues.Slice(k, proxy1->Dimension()) = val(0,0);
                        }
                    }
            
                // td.Stop();

                if (!mir.IsComplex())
                  {
                    if (!is_diagonal)
                      for (int i = 0; i < mir.Size(); i++)
                        proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
                    else
                      for (int i = 0; i < mir.Size(); i++)
                        diagproxyvalues.Range(proxy1->Dimension()*IntRange(i,i+1)) *= mir[i].GetWeight();
                  }
                else
                  { // pml
                    if (!is_diagonal)
                      for (int i = 0; i < mir.Size(); i++)
                        proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
                    else
                      for (int i = 0; i < mir.Size(); i++)
                        diagproxyvalues.Range(proxy1->Dimension()*IntRange(i,i+1)) *=
                          static_cast<const ScalMappedIntegrationPoint<SCAL>&> (mir[i]).GetJacobiDet()*ir[i].Weight();
                  }
                IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                SliceMatrix<SCAL> part_elmat = elmat.Rows(r2).Cols(r1);
                // FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                // FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

                
                // enum { BS = 16 };
                constexpr size_t BS = 16;
                for (size_t i = 0; i < mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    size_t bs = min2(size_t(BS), mir.Size()-i);
                    
                    AFlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);
                    AFlatMatrix<SCAL> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                    AFlatMatrix<SCAL_SHAPES> bbmat2 = samediffop ?
                      bbmat1 : AFlatMatrix<SCAL_SHAPES>(elmat.Height(), bs*proxy2->Dimension(), lh);

                    // tb.Start();
                    BaseMappedIntegrationRule & bmir = mir.Range(i, i+bs, lh);
                    proxy1->Evaluator()->CalcMatrix(fel_trial, bmir, Trans(bbmat1), lh);
                    
                    if (!samediffop)
                      proxy2->Evaluator()->CalcMatrix(fel_test, bmir, Trans(bbmat2), lh);
                    // tb.Stop();

                    // tdb.Start();
                    if (is_diagonal)
                      {
                        AFlatVector<SCAL> diagd(bs*proxy1->Dimension(), lh);
                        diagd = diagproxyvalues.Range(i*proxy1->Dimension(),
                                                      (i+bs)*proxy1->Dimension());
                        
                        /*
                        for (int i = 0; i < diagd.Size(); i++)
                          bdbmat1.Col(i) = diagd(i) * bbmat1.Col(i);
                        */
                        MultMatDiagMat(bbmat1, diagd, bdbmat1);
                        // tdb.AddFlops (bbmat1.Height()*bbmat1.Width());
                      }
                    else
                      {
                        for (int j = 0; j < bs; j++)
                          {
                            int ii = i+j;
                            IntRange r1 = proxy1->Dimension() * IntRange(j,j+1);
                            IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                            // bdbmat1.Cols(r2) = bbmat1.Cols(r1) * proxyvalues(ii,STAR,STAR);
                            MultMatMat (bbmat1.Cols(r1), proxyvalues(ii,STAR,STAR), bdbmat1.Cols(r2));
                          }
                        // tdb.AddFlops (proxy1->Dimension()*proxy2->Dimension()*bs*bbmat1.Height());
                      }
                    //  tdb.Stop();
                    
                    // tlapack.Start();
                    // elmat.Rows(r2).Cols(r1) += bbmat2.Rows(r2) * Trans(bdbmat1.Rows(r1));
                    // AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), elmat.Rows(r2).Cols(r1));

                    symmetric_so_far &= samediffop && is_diagonal;

                    if (symmetric_so_far)
                      AddABtSym (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);
                    else
                      AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);

                    // tlapack.Stop();
                    // tlapack.AddFlops (r2.Size()*r1.Size()*bdbmat1.Width());
                  }

                if (symmetric_so_far)
                  for (int i = 0; i < part_elmat.Height(); i++)
                    for (int j = i+1; j < part_elmat.Width(); j++)
                      part_elmat(i,j) = part_elmat(j,i);
              }
            
            l1 += proxy2->Dimension();  
          }
        k1 += proxy1->Dimension();
      }
  }

#define SIMD_CALCMATRIX
#ifdef SIMD_CALCMATRIX
  template <>
  void SymbolicBilinearFormIntegrator ::
  T_CalcElementMatrix<double,double> (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatMatrix<> elmat,
                       LocalHeap & lh) const
    
  {
    typedef double SCAL;
    typedef double SCAL_SHAPES;
    // static Timer t("symbolicBFI - CalcElementMatrix", 2);
    // static Timer tsimd("symbolicBFI - CalcElementMatrix - simd", 2);
    // static Timer tstart("symboliBFI - CalcElementMatrix startup", 2);
    // static Timer tstart1("symboliBFI - CalcElementMatrix startup 1", 2);
    // static Timer tmain("symboliBFI - CalcElementMatrix main", 2);

    // static Timer td("symboliBFI - CalcElementMatrix dmats", 2);
    // static Timer tb("symboliBFI - CalcElementMatrix diffops", 2);
    // static Timer tdb("symboliBFI - CalcElementMatrix D * B", 2);
    // static Timer tlapack("symboliBFI - CalcElementMatrix lapack", 2);
    
    // tstart.Start();
    
    if (element_boundary)
      {
        switch (trafo.SpaceDim())
          {
          case 1:
            T_CalcElementMatrixEB<1,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          case 2:
            T_CalcElementMatrixEB<2,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          case 3:
            T_CalcElementMatrixEB<3,SCAL, SCAL_SHAPES> (fel, trafo, elmat, lh);
            return;
          default:
            throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
          }
      }
    
    // RegionTimer reg(t);

    // const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    bool is_mixedfe = typeid(fel) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe = static_cast<const MixedFiniteElement*> (&fel);
    // if (is_mixedfe != ( mixedfe != nullptr) ) cout << "different" << endl;
    const FiniteElement & fel_trial = is_mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = is_mixedfe ? mixedfe->FETest() : fel;

    
    elmat = 0;
    if (simd_evaluate)
      try
        {
          // RegionTimer reg(tsimd);          
          SIMD_IntegrationRule ir = Get_SIMD_IntegrationRule (fel);
          SIMD_BaseMappedIntegrationRule & mir = trafo(ir, lh);

          ProxyUserData ud;
          const_cast<ElementTransformation&>(trafo).userdata = &ud;

          bool symmetric_so_far = true;
          int k1 = 0;
          for (auto proxy1 : trial_proxies)
            {
              int l1 = 0;
              for (auto proxy2 : test_proxies)
                {
                  bool is_diagonal = proxy1->Dimension() == proxy2->Dimension();
                  bool is_nonzero = false;
                  
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    for (int l = 0; l < proxy2->Dimension(); l++)
                      if (nonzeros(l1+l, k1+k))
                        {
                          if (k != l) is_diagonal = false;
                          is_nonzero = true;
                        }

                  if (is_nonzero)
                    {
                      HeapReset hr(lh);
                      bool samediffop = (*(proxy1->Evaluator()) == *(proxy2->Evaluator())) && !is_mixedfe;
                      // td.Start();
                      AFlatMatrix<SCAL> proxyvalues(proxy1->Dimension()*proxy2->Dimension(), ir.GetNIP(), lh);
                      AFlatMatrix<SCAL> diagproxyvalues(proxy1->Dimension(), ir.GetNIP(), lh);
                      AFlatMatrix<SCAL> val(1, ir.GetNIP(), lh);

                      
                      if (!is_diagonal)
                        for (size_t k = 0, kk = 0; k < proxy1->Dimension(); k++)
                          for (size_t l = 0; l < proxy2->Dimension(); l++, kk++)
                            {
                              if (nonzeros(l1+l, k1+k))
                                {
                                  if (k != l) is_diagonal = false;
                                  is_nonzero = true;
                                  ud.trialfunction = proxy1;
                                  ud.trial_comp = k;
                                  ud.testfunction = proxy2;
                                  ud.test_comp = l;
                                  
                                  cf -> Evaluate (mir, proxyvalues.Rows(kk,kk+1));
                                }
                              else
                                proxyvalues.Row(kk) = 0.0;
                            }
                      else
                        for (int k = 0; k < proxy1->Dimension(); k++)
                          {
                            ud.trialfunction = proxy1;
                            ud.trial_comp = k;
                            ud.testfunction = proxy2;
                            ud.test_comp = k;
                            
                            // if (!elementwise_constant)
                            cf -> Evaluate (mir, diagproxyvalues.Rows(k,k+1));
                            /*
                              else
                              {
                                cf -> Evaluate (mir[0], val);
                                diagproxyvalues.Row(k) = val(0,0);
                                }
                              */
                          }
                      // td.Stop();
                      
                      if (!is_diagonal)
                        for (size_t i = 0; i < mir.Size(); i++)
                          {
                            auto fac = mir[i].GetWeight();
                            for (size_t j = 0; j < proxyvalues.Height(); j++)
                              proxyvalues.Get(j,i) *= fac;
                          }
                      else
                        for (size_t i = 0; i < mir.Size(); i++)
                          {
                            auto fac = mir[i].GetWeight();
                            for (size_t j = 0; j < proxy1->Dimension(); j++)
                              diagproxyvalues.Get(j,i) *= fac;
                          }

                      IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                      IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                      SliceMatrix<SCAL> part_elmat = elmat.Rows(r2).Cols(r1);
                      FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                      FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

                      // enum { BS = 16 };
                      // for (int i = 0; i < mir.Size(); i+=BS)
                        {
                          HeapReset hr(lh);
                          // int bs = min2(int(BS), mir.Size()-i);
                          AFlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width()*proxy1->Dimension(), ir.GetNIP(), lh);
                          AFlatMatrix<SCAL> bdbmat1(elmat.Width()*proxy2->Dimension(), ir.GetNIP(), lh);
                          AFlatMatrix<SCAL_SHAPES> bbmat2 = samediffop ?
                            bbmat1 : AFlatMatrix<SCAL_SHAPES>(elmat.Height()*proxy2->Dimension(), ir.GetNIP(), lh);
                          
                          AFlatMatrix<SCAL_SHAPES> hbdbmat1(elmat.Width(), proxy2->Dimension()*SIMD<double>::Size()*ir.Size(),
                                                            &bdbmat1.Get(0,0));
                          AFlatMatrix<SCAL_SHAPES> hbbmat2(elmat.Height(), proxy2->Dimension()*SIMD<double>::Size()*ir.Size(),
                                                           &bbmat2.Get(0,0));

                          // tb.Start();
                          // BaseMappedIntegrationRule & bmir = mir.Range(i, i+bs, lh);
                          bbmat1 = 0.0;
                          bbmat2 = 0.0;
                          proxy1->Evaluator()->CalcMatrix(fel_trial, mir, bbmat1);
                          
                          if (!samediffop)
                            proxy2->Evaluator()->CalcMatrix(fel_test, mir, bbmat2);

                          // tb.Stop();
                          // tdb.Start();
                          if (is_diagonal)
                            { // too much work, use r1 ... 
                              for (size_t i = 0, ii = 0; i < elmat.Width(); i++)
                                for (size_t j = 0; j < proxy1->Dimension(); j++, ii++)
                                  for (size_t k = 0; k < mir.Size(); k++)
                                    bdbmat1.Get(ii,k) = bbmat1.Get(ii,k) * diagproxyvalues.Get(j, k);

                              // AFlatVector<SCAL> diagd(bs*proxy1->Dimension(), lh);
                              // diagd = diagproxyvalues.Range(i*proxy1->Dimension(),
                              // (i+bs)*proxy1->Dimension());
                              // MultMatDiagMat(bbmat1, diagd, bdbmat1);
                            }
                          else
                            {
                              bdbmat1 = 0.0;
                              // for (size_t i = 0; i < elmat.Width(); i++)
                              for (auto i : r1)
                                for (size_t j = 0; j < proxy2->Dimension(); j++)
                                  for (size_t k = 0; k < proxy1->Dimension(); k++)
                                    {
                                      auto res = bdbmat1.Row(i*proxy2->Dimension()+j);
                                      auto a = bbmat1.Row(i*proxy1->Dimension()+k);
                                      auto b = proxyvalues.Row(k*proxy2->Dimension()+j);
                                      for (size_t l = 0; l < mir.Size(); l++)
                                        res.Get(l) += a.Get(l) * b.Get(l);
                                    }
                            }

                          // tdb.Stop();
                          
                          // tlapack.Start();
                          // elmat.Rows(r2).Cols(r1) += bbmat2.Rows(r2) * Trans(bdbmat1.Rows(r1));
                          // AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), elmat.Rows(r2).Cols(r1));

                          symmetric_so_far &= samediffop && is_diagonal;
                          if (symmetric_so_far)
                            AddABtSym (hbbmat2.Rows(r2), hbdbmat1.Rows(r1), part_elmat);
                          else
                            AddABt (hbbmat2.Rows(r2), hbdbmat1.Rows(r1), part_elmat);
                          
                          // tlapack.Stop();
                          // tlapack.AddFlops (r2.Size()*r1.Size()*bdbmat1.Width());
                        }
                      
                      if (symmetric_so_far)
                        for (size_t i = 0; i < part_elmat.Height(); i++)
                          for (size_t j = i+1; j < part_elmat.Width(); j++)
                            part_elmat(i,j) = part_elmat(j,i);
                    }
                  
                  l1 += proxy2->Dimension();  
                }
              k1 += proxy1->Dimension();
            }
          // throw ExceptionNOSIMD("test1");
          // *testout << "elmat = " << endl << elmat << endl;
          return;
        }
      catch (ExceptionNOSIMD e)
        {
          cout << IM(4) << e.What() << endl
               << "switching to scalar evaluation" << endl;
          simd_evaluate = false;
          T_CalcElementMatrix (fel, trafo, elmat, lh);
          return;
        }
    

    // IntegrationRule ir(trafo.GetElementType(), intorder);
    IntegrationRule ir = GetIntegrationRule (fel);    
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    
    
    // tstart.Stop();
    bool symmetric_so_far = true;
    int k1 = 0;
    for (auto proxy1 : trial_proxies)
      {
        int l1 = 0;
        for (auto proxy2 : test_proxies)
          {
            bool is_diagonal = proxy1->Dimension() == proxy2->Dimension();
            bool is_nonzero = false;

            for (int k = 0; k < proxy1->Dimension(); k++)
              for (int l = 0; l < proxy2->Dimension(); l++)
                if (nonzeros(l1+l, k1+k))
                  {
                    if (k != l) is_diagonal = false;
                    is_nonzero = true;
                  }


            if (is_nonzero)
              {
                HeapReset hr(lh);
                bool samediffop = *(proxy1->Evaluator()) == *(proxy2->Evaluator());
                // td.Start();
                FlatTensor<3,SCAL> proxyvalues(lh, mir.Size(), proxy1->Dimension(), proxy2->Dimension());
                FlatVector<SCAL> diagproxyvalues(mir.Size()*proxy1->Dimension(), lh);
                FlatMatrix<SCAL> val(mir.Size(), 1, lh);
                
                
                if (!is_diagonal)
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    for (int l = 0; l < proxy2->Dimension(); l++)
                      {
                        if (nonzeros(l1+l, k1+k))
                          {
                            if (k != l) is_diagonal = false;
                            is_nonzero = true;
                            ud.trialfunction = proxy1;
                            ud.trial_comp = k;
                            ud.testfunction = proxy2;
                            ud.test_comp = l;
                            
                            cf -> Evaluate (mir, val);
                            proxyvalues(STAR,k,l) = val.Col(0);
                          }
                        else
                          proxyvalues(STAR,k,l) = 0.0;
                      }
                else
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    {
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = k;

                      if (!elementwise_constant)
                        {
                          cf -> Evaluate (mir, val);
                          diagproxyvalues.Slice(k, proxy1->Dimension()) = val.Col(0);
                        }
                      else
                        {
                          cf -> Evaluate (mir[0], val.Row(0));
                          diagproxyvalues.Slice(k, proxy1->Dimension()) = val(0,0);
                        }
                    }
            
                // td.Stop();

                if (!mir.IsComplex())
                  {
                    if (!is_diagonal)
                      for (int i = 0; i < mir.Size(); i++)
                        proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
                    else
                      for (int i = 0; i < mir.Size(); i++)
                        diagproxyvalues.Range(proxy1->Dimension()*IntRange(i,i+1)) *= mir[i].GetWeight();
                  }
                else
                  { // pml
                    if (!is_diagonal)
                      for (int i = 0; i < mir.Size(); i++)
                        proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
                    else
                      for (int i = 0; i < mir.Size(); i++)
                        diagproxyvalues.Range(proxy1->Dimension()*IntRange(i,i+1)) *=
                          static_cast<const ScalMappedIntegrationPoint<SCAL>&> (mir[i]).GetJacobiDet()*ir[i].Weight();
                  }
                IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                SliceMatrix<SCAL> part_elmat = elmat.Rows(r2).Cols(r1);
                FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

                
                // enum { BS = 16 };
                constexpr size_t BS = 16;
                for (int i = 0; i < mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    int bs = min2(size_t(BS), mir.Size()-i);
                    
                    AFlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);
                    AFlatMatrix<SCAL> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                    AFlatMatrix<SCAL_SHAPES> bbmat2 = samediffop ?
                      bbmat1 : AFlatMatrix<SCAL_SHAPES>(elmat.Height(), bs*proxy2->Dimension(), lh);

                    // tb.Start();
                    BaseMappedIntegrationRule & bmir = mir.Range(i, i+bs, lh);
                    proxy1->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat1), lh);
                    
                    if (!samediffop)
                      proxy2->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat2), lh);
                    // tb.Stop();

                    // tdb.Start();
                    if (is_diagonal)
                      {
                        AFlatVector<SCAL> diagd(bs*proxy1->Dimension(), lh);
                        diagd = diagproxyvalues.Range(i*proxy1->Dimension(),
                                                      (i+bs)*proxy1->Dimension());
                        /*
                        for (int i = 0; i < diagd.Size(); i++)
                          bdbmat1.Col(i) = diagd(i) * bbmat1.Col(i);
                        */
                        MultMatDiagMat(bbmat1, diagd, bdbmat1);
                        // tdb.AddFlops (bbmat1.Height()*bbmat1.Width());
                      }
                    else
                      {
                        for (int j = 0; j < bs; j++)
                          {
                            int ii = i+j;
                            IntRange r1 = proxy1->Dimension() * IntRange(j,j+1);
                            IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                            // bdbmat1.Cols(r2) = bbmat1.Cols(r1) * proxyvalues(ii,STAR,STAR);
                            MultMatMat (bbmat1.Cols(r1), proxyvalues(ii,STAR,STAR), bdbmat1.Cols(r2));
                          }
                        // tdb.AddFlops (proxy1->Dimension()*proxy2->Dimension()*bs*bbmat1.Height());
                      }
                    // tdb.Stop();
                    // tlapack.Start();
                    // elmat.Rows(r2).Cols(r1) += bbmat2.Rows(r2) * Trans(bdbmat1.Rows(r1));
                    // AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), elmat.Rows(r2).Cols(r1));

                    symmetric_so_far &= samediffop && is_diagonal;
                    if (symmetric_so_far)
                      AddABtSym (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);
                    else
                      AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);
                    // tlapack.Stop();
                    // tlapack.AddFlops (r2.Size()*r1.Size()*bdbmat1.Width());
                  }

                if (symmetric_so_far)
                  for (int i = 0; i < part_elmat.Height(); i++)
                    for (int j = i+1; j < part_elmat.Width(); j++)
                      part_elmat(i,j) = part_elmat(j,i);
              }
            
            l1 += proxy2->Dimension();  
          }
        k1 += proxy1->Dimension();
      }
  }
#endif  

  void 
  SymbolicBilinearFormIntegrator ::
  CalcElementMatrix (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    T_CalcElementMatrix<double> (fel, trafo, elmat, lh);
  }
  
  void 
  SymbolicBilinearFormIntegrator ::
  CalcElementMatrix (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatMatrix<Complex> elmat,
                     LocalHeap & lh) const
  {
    if (fel.ComplexShapes() || trafo.IsComplex())
      T_CalcElementMatrix<Complex,Complex> (fel, trafo, elmat, lh);
    else
      T_CalcElementMatrix<Complex,double> (fel, trafo, elmat, lh);
  }


  

  template <int D, typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcElementMatrixEB (const FiniteElement & fel,
                           const ElementTransformation & trafo, 
                           FlatMatrix<SCAL> elmat,
                           LocalHeap & lh) const
      
    {
      static Timer t("symbolicBFI - CalcElementMatrix EB", 2);
      /*
      static Timer tir("symbolicBFI - CalcElementMatrix EB - intrules", 2);
      static Timer td("symbolicBFI - CalcElementMatrix EB - dmats", 2);
      static Timer tdb("symbolicBFI - CalcElementMatrix EB - b*d", 2);
      static Timer tb("symbolicBFI - CalcElementMatrix EB - bmats", 2);
      static Timer tmult("symbolicBFI - CalcElementMatrix EB - mult", 2);
      */
      // RegionTimer reg(t);

      elmat = 0;

      const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
      const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
      const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;
      
      auto eltype = trafo.GetElementType();
      int nfacet = ElementTopology::GetNFacets(eltype);

      Facet2ElementTrafo transform(eltype); 

      for (int k = 0; k < nfacet; k++)
        {
          // tir.Start();
          HeapReset hr(lh);
          ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
        
          IntegrationRule ir_facet(etfacet, fel_trial.Order()+fel_test.Order());
          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
          BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
          
          ProxyUserData ud;
          const_cast<ElementTransformation&>(trafo).userdata = &ud;

          mir.ComputeNormalsAndMeasure(eltype, k);
          
          for (int k1 : Range(trial_proxies))
            for (int l1 : Range(test_proxies))
              {
                if (!nonzeros_proxies(l1, k1)) continue;
                
                auto proxy1 = trial_proxies[k1];
                auto proxy2 = test_proxies[l1];
                
                HeapReset hr(lh);
                FlatTensor<3,SCAL> proxyvalues(lh, mir.Size(), proxy1->Dimension(), proxy2->Dimension());
                FlatMatrix<SCAL> val(mir.Size(), 1, lh);
                
                // td.Start();
                for (int k = 0; k < proxy1->Dimension(); k++)
                  for (int l = 0; l < proxy2->Dimension(); l++)
                    {
                      ud.trialfunction = proxy1;
                      ud.trial_comp = k;
                      ud.testfunction = proxy2;
                      ud.test_comp = l;
                      
                      cf->Evaluate (mir, val);
                      for (int i = 0; i < mir.Size(); i++)
                        val(i) *= ir_facet[i].Weight() * mir[i].GetMeasure(); 

                      proxyvalues(STAR,k,l) = val.Col(0);
                    }
                // td.Stop();
                /*
                for (int i = 0; i < mir.Size(); i++)
                  {
                    tb.Start();
                    FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                    FlatMatrix<SCAL,ColMajor> dbmat1(proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
                    
                    proxy1->Evaluator()->CalcMatrix(fel, mir[i], bmat1, lh);
                    proxy2->Evaluator()->CalcMatrix(fel, mir[i], bmat2, lh);
                    tb.Stop();
                    tmult.Start();
                    IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                    IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                    
                    dbmat1 = proxyvalues(i,STAR,STAR) * bmat1;                    
                    elmat.Rows(r2).Cols(r1) += Trans (bmat2.Cols(r2)) * dbmat1.Cols(r1);
                    tmult.Stop();
                  }
                */
                
                IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                SliceMatrix<SCAL> part_elmat = elmat.Rows(r2).Cols(r1);
                // enum { BS = 16 };
                constexpr size_t BS = 16;
                for (size_t i = 0; i < mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    int bs = min2(BS, mir.Size()-i);
                    
                    AFlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);
                    AFlatMatrix<SCAL> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                    AFlatMatrix<SCAL_SHAPES> bbmat2(elmat.Height(), bs*proxy2->Dimension(), lh);

                    // tb.Start();
                    BaseMappedIntegrationRule & bmir = mir.Range(i, i+bs, lh);
                    proxy1->Evaluator()->CalcMatrix(fel_trial, bmir, Trans(bbmat1), lh);
                    proxy2->Evaluator()->CalcMatrix(fel_test, bmir, Trans(bbmat2), lh);
                    // tb.Stop();
                    bdbmat1 = 0.0;
                    // tdb.Start();

                    auto part_bbmat1 = bbmat1.Rows(r1);
                    auto part_bdbmat1 = bdbmat1.Rows(r1);
                    auto part_bbmat2 = bbmat2.Rows(r2);
                    
                    for (int j = 0; j < bs; j++)
                      {
                        IntRange rj1 = proxy1->Dimension() * IntRange(j,j+1);
                        IntRange rj2 = proxy2->Dimension() * IntRange(j,j+1);
                        MultMatMat (part_bbmat1.Cols(rj1), proxyvalues(i+j,STAR,STAR), part_bdbmat1.Cols(rj2));
                      }

                    // tdb.Stop();
                    
                    // tmult.Start();                    
                    AddABt (part_bbmat2, part_bdbmat1, part_elmat);
                    // part_elmat += part_bbmat2 * Trans(part_bdbmat1);
                    // tmult.Stop();
                    // tmult.AddFlops (r2.Size() * r1.Size() * bbmat2.Width());
                  }                
              }
        }
    }

  
  void 
  SymbolicBilinearFormIntegrator ::
  CalcLinearizedElementMatrix (const FiniteElement & fel,
                               const ElementTransformation & trafo, 
                               FlatVector<double> elveclin,
                               FlatMatrix<double> elmat,
                               LocalHeap & lh) const
  {
    /*
      CalcElementMatrix(fel, trafo, elmat, lh);
      return;
    */

      
    if (element_boundary)
      {
        switch (trafo.SpaceDim())
          {
          case 1:
            T_CalcLinearizedElementMatrixEB<1,double,double> (fel, trafo, elveclin, elmat, lh);
            return;
          case 2:
            T_CalcLinearizedElementMatrixEB<2,double,double> (fel, trafo, elveclin, elmat, lh);            
            return;
          case 3:
            T_CalcLinearizedElementMatrixEB<3,double,double> (fel, trafo, elveclin, elmat, lh);
            return;
          default:
            throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
          }
      }

    
    // static Timer t("symbolicbfi - calclinearized", 2);
    // static Timer td("symbolicbfi - calclinearized dmats", 2);
    // RegionTimer reg(t);



    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;

    /*
    int trial_difforder = 99, test_difforder = 99;
    for (auto proxy : trial_proxies)
      trial_difforder = min2(trial_difforder, proxy->Evaluator()->DiffOrder());
    for (auto proxy : test_proxies)
      test_difforder = min2(test_difforder, proxy->Evaluator()->DiffOrder());

    int intorder = fel_trial.Order()+fel_test.Order();
    auto et = trafo.GetElementType();
    if (et == ET_TRIG || et == ET_TET)
      intorder -= test_difforder+trial_difforder;
    */
    
    // IntegrationRule ir(trafo.GetElementType(), intorder);
    IntegrationRule ir = GetIntegrationRule (fel);
    // IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elveclin;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      {
        ud.AssignMemory (proxy, ir.Size(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel_trial, mir, elveclin, ud.GetMemory(proxy), lh);
      }
    
    FlatMatrix<> val(mir.Size(), 1, lh), deriv(mir.Size(), 1, lh);
    elmat = 0;
    
    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          // td.Start(); 
          FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              // if (nonzeros(test_cum[l1]+l, trial_cum[k1]+k))
              if (true)
                {
                  ud.trialfunction = proxy1;
                  ud.trial_comp = k;
                  ud.testfunction = proxy2;
                  ud.test_comp = l;
                  
                  cf -> EvaluateDeriv (mir, val, deriv);
                  proxyvalues(STAR,l,k) = deriv.Col(0);
                }
              else
                proxyvalues(STAR,l,k) = 0;
          // td.Stop();

          for (int i = 0; i < mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();

          // t.AddFlops (double (mir.Size()) * proxy1->Dimension()*elmat.Width()*elmat.Height());

          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

          // enum { BS = 16 };
          constexpr size_t BS = 16;
          for (int i = 0; i < mir.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  proxy1->Evaluator()->CalcMatrix(fel_trial, mir[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel_test, mir[ii], bmat2, lh);
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
              elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }
  }



  template <int D, typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcLinearizedElementMatrixEB (const FiniteElement & fel,
                                   const ElementTransformation & trafo, 
                                   FlatVector<double> elveclin,
                                   FlatMatrix<double> elmat,
                                   LocalHeap & lh) const
  {
    static Timer t("symbolicbfi - calclinearized EB", 2);
    static Timer td("symbolicbfi - calclinearized EB dmats", 2);
    RegionTimer reg(t);
    
    // IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    // BaseMappedIntegrationRule & mir = trafo(ir, lh);
    /*
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elveclin;
    ud.lh = &lh;
    */
    elmat = 0;

      auto eltype = trafo.GetElementType();
      int nfacet = ElementTopology::GetNFacets(eltype);

      Facet2ElementTrafo transform(eltype); 

      for (int k = 0; k < nfacet; k++)
        {
          HeapReset hr(lh);
          ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
        
          IntegrationRule ir_facet(etfacet, 2*fel.Order());
          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
          BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
          mir.ComputeNormalsAndMeasure(eltype, k);
          
          ProxyUserData ud(trial_proxies.Size(), lh);          
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          ud.fel = &fel;
          ud.elx = &elveclin;
          ud.lh = &lh;

          for (ProxyFunction * proxy : trial_proxies)
            {
              ud.AssignMemory (proxy, mir.Size(), proxy->Dimension(), lh);
              proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetMemory(proxy), lh);
            }
    
    
          FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh);
          for (int k1 : Range(trial_proxies))
            for (int l1 : Range(test_proxies))
              {
                HeapReset hr(lh);
                auto proxy1 = trial_proxies[k1];
                auto proxy2 = test_proxies[l1];
                td.Start();
                FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
                
                for (int k = 0; k < proxy1->Dimension(); k++)
                  for (int l = 0; l < proxy2->Dimension(); l++)
                    // if (nonzeros(test_cum[l1]+l, trial_cum[k1]+k))
                    if (true)
                      {
                        ud.trialfunction = proxy1;
                        ud.trial_comp = k;
                        ud.testfunction = proxy2;
                        ud.test_comp = l;
                        
                        cf -> EvaluateDeriv (mir, val, deriv);
                        proxyvalues(STAR,l,k) = deriv.Col(0);
                      }
                    else
                      proxyvalues(STAR,l,k) = 0.0;
                        
                td.Stop();

                for (int i = 0; i < mir.Size(); i++)
                  proxyvalues(i,STAR,STAR) *= ir_facet[i].Weight() * mir[i].GetMeasure(); 
                
                t.AddFlops (double (mir.Size()) * proxy1->Dimension()*elmat.Width()*elmat.Height());
                
                FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
                
                // enum { BS = 16 };
                constexpr size_t BS = 16;
                for (int i = 0; i < mir.Size(); i+=BS)
                  {
                    int rest = min2(BS, mir.Size()-i);
                    HeapReset hr(lh);
                    FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);
                    
                    for (int j = 0; j < rest; j++)
                      {
                        int ii = i+j;
                        IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                        proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                        proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                        bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                        bbmat2.Rows(r2) = bmat2;
                      }
                    
                    IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                    IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                    elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
                  }


                /*
                int i = 0;
                
                enum { BS = 16 };
                for ( ; i+BS <= mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    FlatMatrix<double,ColMajor> bdbmat1(BS*proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<double,ColMajor> bbmat2(BS*proxy2->Dimension(), elmat.Height(), lh);
                    
                    for (int j = 0; j < BS; j++)
                      {
                        int ii = i+j;
                        IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                        proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                        proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                        bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                        bbmat2.Rows(r2) = bmat2;
                      }
                    elmat += Trans (bbmat2) * bdbmat1 | Lapack;
                  }
                
                
                if (i < mir.Size())
                  {
                    HeapReset hr(lh);
                    int rest = mir.Size()-i;
                    FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
                    FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);
                    
                    for (int j = 0; j < rest; j++)
                      {
                        int ii = i+j;
                        IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                        proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                        proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                        bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                        bbmat2.Rows(r2) = bmat2;
                      }
                    elmat += Trans (bbmat2) * bdbmat1 | Lapack;
                  }
                */
              }
        }
  }









  
  
  void
  SymbolicBilinearFormIntegrator :: ApplyElementMatrix (const FiniteElement & fel, 
                                                        const ElementTransformation & trafo, 
                                                        const FlatVector<double> elx, 
                                                        FlatVector<double> ely,
                                                        void * precomputed,
                                                        LocalHeap & lh) const
  {
    /*
    if (element_boundary)
      {
        FlatMatrix<> elmat(elx.Size(), lh);
        CalcElementMatrix(fel, trafo, elmat, lh);
        ely = elmat * elx;
        return;
      }
    */
    // const TPHighOrderFE * tpfel = dynamic_cast<const TPHighOrderFE *>(&fel);
    // if(tpfel)
    // {
      // ApplyElementMatrixTP(fel,trafo,elx,ely,precomputed,lh);
      // return;
    // }
    if (element_boundary)
      {
        switch (trafo.SpaceDim())
          {
          case 1:
            T_ApplyElementMatrixEB<1,double,double> (fel, trafo, elx, ely, precomputed, lh);
            return;
          case 2:
            T_ApplyElementMatrixEB<2,double,double> (fel, trafo, elx, ely, precomputed, lh);            
            return;
          case 3:
            T_ApplyElementMatrixEB<3,double,double> (fel, trafo, elx, ely, precomputed, lh);            
            return;
          default:
            throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
          }
      }


    if (simd_evaluate)
      try
        {
          static Timer tall("SymbolicBFI::Apply - all", 4); RegionTimer rall(tall);

          const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
          const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
          const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;

          HeapReset hr(lh);

          SIMD_IntegrationRule simd_ir = Get_SIMD_IntegrationRule (fel);
          auto & simd_mir = trafo(simd_ir, lh);
          
          ProxyUserData ud(trial_proxies.Size(), lh);
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          ud.fel = &fel;
          ud.elx = &elx;
          ud.lh = &lh;
          for (ProxyFunction * proxy : trial_proxies)
            ud.AssignMemory (proxy, simd_ir.GetNIP(), proxy->Dimension(), lh);
          
          for (ProxyFunction * proxy : trial_proxies)
            proxy->Evaluator()->Apply(fel_trial, simd_mir, elx, ud.GetAMemory(proxy)); 
          
          ely = 0;
          for (auto proxy : test_proxies)
            {
              HeapReset hr(lh);

              AFlatMatrix<double> simd_proxyvalues(proxy->Dimension(), simd_ir.GetNIP(), lh);
              for (int k = 0; k < proxy->Dimension(); k++)
                {
                  ud.testfunction = proxy;
                  ud.test_comp = k;
                  cf -> Evaluate (simd_mir, simd_proxyvalues.Rows(k,k+1));
                }
              
              for (size_t i = 0; i < simd_proxyvalues.Height(); i++)
                {
                  auto row = simd_proxyvalues.Row(i);
                  for (size_t j = 0; j < row.VSize(); j++)
                    row.Get(j) *= simd_mir[j].GetMeasure().Data() * simd_ir[j].Weight().Data();
                }

              proxy->Evaluator()->AddTrans(fel_test, simd_mir, simd_proxyvalues, ely); 
            }
          return;
        }
      catch (ExceptionNOSIMD e)
        {
          cout << IM(4) << e.What() << endl
               << "switching to scalar evaluation" << endl;
          simd_evaluate = false;
          ApplyElementMatrix (fel, trafo, elx, ely, precomputed, lh);
          return;
        }

    
    // no simd
    HeapReset hr(lh);
    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;
 
    ProxyUserData ud(trial_proxies.Size(), lh);    
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;

    // IntegrationRule ir(trafo.GetElementType(), fel_trial.Order()+fel_test.Order());
    IntegrationRule ir = GetIntegrationRule (fel);

    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);

    for (ProxyFunction * proxy : trial_proxies)
      proxy->Evaluator()->Apply(fel_trial, mir, elx, ud.GetMemory(proxy), lh);
    
    ely = 0;
    FlatVector<> ely1(ely.Size(), lh);

    FlatMatrix<> val(mir.Size(), 1,lh);
    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir, val);
            proxyvalues.Col(k) = val.Col(0);
          }
        
        for (int i = 0; i < mir.Size(); i++)
          proxyvalues.Row(i) *= mir[i].GetWeight();

        proxy->Evaluator()->ApplyTrans(fel_test, mir, proxyvalues, ely1, lh);
        ely += ely1;
      }
  }

 void
   TensorProductBilinearFormIntegrator :: ApplyElementMatrix (const FiniteElement & fel, 
                                                         const ElementTransformation & trafo, 
                                                         const FlatVector<double> elx, 
                                                         FlatVector<double> ely,
                                                         void * precomputed,
                                                         LocalHeap & lh) const
   {
     static Timer tallnosimd("SymbolicBFI::Apply"); RegionTimer rallnosimd(tallnosimd);
     const TPHighOrderFE & tpfel = dynamic_cast<const TPHighOrderFE &> (fel);
     // if (element_boundary)
       // {
         // switch (trafo.SpaceDim())
           // {
           // case 1:
             // T_ApplyElementMatrixEB<1,double,double> (fel, trafo, elx, ely, precomputed, lh);
             // return;
           // case 2:
             // T_ApplyElementMatrixEB<2,double,double> (fel, trafo, elx, ely, precomputed, lh);            
             // return;
           // case 3:
             // T_ApplyElementMatrixEB<3,double,double> (fel, trafo, elx, ely, precomputed, lh);            
             // return;
           // default:
             // throw Exception ("Illegal space dimension" + ToString(trafo.SpaceDim()));
           // }
       // }
     HeapReset hr(lh);
     ely = 0;
     ArrayMem<const IntegrationRule *,2> irs(2);
     const IntegrationRule ir0 = SelectIntegrationRule(tpfel.elements[0]->ElementType(),2*tpfel.elements[0]->Order());
     const IntegrationRule ir1 = SelectIntegrationRule(tpfel.elements[1]->ElementType(),2*tpfel.elements[1]->Order());
     irs[0] = &ir0;
     irs[1] = &ir1;
     TPIntegrationRule ir(irs);
     BaseMappedIntegrationRule & mir = trafo(ir, lh);
     int niptp = ir0.Size()*ir1.Size();
     ProxyUserData ud(trial_proxies.Size(), lh);    
     const_cast<ElementTransformation&>(trafo).userdata = &ud;
     ud.fel = &fel;
     ud.elx = &elx;
     ud.lh = &lh;
     for (ProxyFunction * proxy : trial_proxies)
       ud.AssignMemory (proxy, niptp, proxy->Dimension(), lh);
     for (ProxyFunction * proxy : trial_proxies)
         proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
     FlatVector<> ely1(ely.Size(), lh);
     FlatMatrix<> val(niptp, 1,lh);
     for (auto proxy : test_proxies)
       {
         HeapReset hr(lh);
         FlatMatrix<> proxyvalues(niptp, proxy->Dimension(), lh);
         for (int k = 0; k < proxy->Dimension(); k++)
           {
             ud.testfunction = proxy;
             ud.test_comp = k;
             cf -> Evaluate (mir, val);
             proxyvalues.Col(k) = val.Col(0);
           }
         TPMappedIntegrationRule & tpmir = dynamic_cast<TPMappedIntegrationRule& >(mir);
         for(int i=0,ii=0;i<tpmir.GetIRs()[0]->Size();i++)
           for(int j=0;j<tpmir.GetIRs()[1]->Size();j++,ii++)
             proxyvalues.Row(ii) *= tpmir.GetIRs()[0]->operator[](i).GetWeight()*tpmir.GetIRs()[1]->operator[](j).GetWeight();
         
         proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, ely1, lh);
         ely += ely1;
       }
   }  
   
  

  template <int D, typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_ApplyElementMatrixEB (const FiniteElement & fel, 
                          const ElementTransformation & trafo, 
                          const FlatVector<double> elx, 
                          FlatVector<double> ely,
                          void * precomputed,
                          LocalHeap & lh) const
  {
    /*
    ProxyUserData ud(trial_proxies.Size(), lh);    
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;

    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    cout << "ir.GetNIP() = " << ir.GetNIP() << endl;
    */
    
    ely = 0;
    
    auto eltype = trafo.GetElementType();
    int nfacet = ElementTopology::GetNFacets(eltype);

    Facet2ElementTrafo transform(eltype); 

    for (int k = 0; k < nfacet; k++)
      {
        HeapReset hr(lh);
        ngfem::ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype, k);
          
        IntegrationRule ir_facet(etfacet, 2*fel.Order());
        IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
        BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);

        ProxyUserData ud(trial_proxies.Size(), lh);    
        const_cast<ElementTransformation&>(trafo).userdata = &ud;
        ud.fel = &fel;
        ud.elx = &elx;
        ud.lh = &lh;
    
        for (ProxyFunction * proxy : trial_proxies)
          ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
        
        for (ProxyFunction * proxy : trial_proxies)
          proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
        
        mir.ComputeNormalsAndMeasure (eltype, k);
        
        FlatVector<> ely1(ely.Size(), lh);
        FlatMatrix<> val(mir.Size(), 1,lh);
        for (auto proxy : test_proxies)
          {
            HeapReset hr(lh);
            FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
            for (int k = 0; k < proxy->Dimension(); k++)
              {
                ud.testfunction = proxy;
                ud.test_comp = k;
                cf -> Evaluate (mir, val);
                proxyvalues.Col(k) = val.Col(0);
              }
            
            for (size_t i = 0; i < mir.Size(); i++)
              proxyvalues.Row(i) *= ir_facet[i].Weight() * mir[i].GetMeasure();
            
            proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, ely1, lh);
            ely += ely1;
          }
      }
  }


  SymbolicFacetLinearFormIntegrator ::
  SymbolicFacetLinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb /* , bool eb */)
    : cf(acf), vb(avb) // , element_boundary(eb)
  {
    if (cf->Dimension() != 1)
      throw Exception ("SymblicLFI needs scalar-valued CoefficientFunction");
    test_cum.Append(0);    
    cf->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
          if (proxy) 
            if (proxy->IsTestFunction())
              if (!proxies.Contains(proxy))
                {
                  proxies.Append (proxy);
                  test_cum.Append(test_cum.Last()+proxy->Dimension());
                }
        });
  }


  void SymbolicFacetLinearFormIntegrator ::
  CalcFacetVector (const FiniteElement & fel1, int LocalFacetNr,
                   const ElementTransformation & trafo1, FlatArray<int> & ElVertices,
                   const ElementTransformation & strafo,  
                   FlatVector<double> elvec,
                   LocalHeap & lh) const
  {
    static Timer t("SymbolicFacetLFI::CalcFacetVector - boundary", 2);
    HeapReset hr(lh);


    elvec = 0;
    
    FlatVector<> elvec1(elvec.Size(), lh);

    int maxorder = fel1.Order();

    auto eltype1 = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr, ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    // auto & smir = strafo(ir_facet, lh);
    
    // evaluate proxy-values
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    const_cast<ElementTransformation&>(strafo).userdata = &ud;

    RegionTimer reg(t);
    
    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr);
    
    FlatMatrix<> val(ir_facet.Size(), 1,lh);
    for (auto proxy : proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(ir_facet.Size(), proxy->Dimension(), lh);
        
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir1, val);  // needed for grad(u), mesh_size, but: index = material index
            // cf -> Evaluate (smir, val);
            proxyvalues.Col(k) = val.Col(0);
          }

        for (int i = 0; i < mir1.Size(); i++)
          proxyvalues.Row(i) *= mir1[i].GetMeasure() * ir_facet[i].Weight(); // use also smir here ? 

        elvec1 = 0.0;
        proxy->Evaluator()->ApplyTrans(fel1, mir1, proxyvalues, elvec1, lh);
        elvec += elvec1;
      }
  }



  
  
  SymbolicFacetBilinearFormIntegrator ::
  SymbolicFacetBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb, bool eb)
    : cf(acf), vb(avb), element_boundary(eb)
  {
    simd_evaluate = true;
    
    if (cf->Dimension() != 1)
        throw Exception ("SymblicBFI needs scalar-valued CoefficientFunction");
    trial_cum.Append(0);
    test_cum.Append(0);    
    cf->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
          if (proxy) 
            {
              if (proxy->IsTestFunction())
                {
                  if (!test_proxies.Contains(proxy))
                    {
                      test_proxies.Append (proxy);
                      test_cum.Append(test_cum.Last()+proxy->Dimension());
                    }
                }
              else
                {                                         
                  if (!trial_proxies.Contains(proxy))
                    {
                      trial_proxies.Append (proxy);
                      trial_cum.Append(trial_cum.Last()+proxy->Dimension());
                    }
                }
            }
        });

    neighbor_testfunction = false;
    for (auto proxy : test_proxies)
      if (proxy->IsOther())
        neighbor_testfunction = true;
    
    cout << IM(3) << "num test_proxies " << test_proxies.Size() << endl;
    cout << IM(3) << "num trial_proxies " << trial_proxies.Size() << endl;
    cout << IM(3) << "cumulated test_proxy dims  " << test_cum << endl;
    cout << IM(3) << "cumulated trial_proxy dims " << trial_cum << endl;
  }


  void SymbolicFacetBilinearFormIntegrator ::
  CalcFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                   const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                   const FiniteElement & fel2, int LocalFacetNr2,
                   const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                   FlatMatrix<double> elmat,
                   LocalHeap & lh) const
  {
    elmat = 0.0;

    if (LocalFacetNr2==-1) throw Exception ("SymbolicFacetBFI: LocalFacetNr2==-1");

    int maxorder = max2 (fel1.Order(), fel2.Order());

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    
    Facet2ElementTrafo transform2(eltype2, ElVertices2); 
    IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, ir_facet, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(ir_facet_vol2, lh);

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());

          mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> Evaluate (mir1, val);
                proxyvalues(STAR,l,k) = val.Col(0);
              }

          for (int i = 0; i < mir1.Size(); i++)
            // proxyvalues(i,STAR,STAR) *= measure(i) * ir_facet[i].Weight();
            proxyvalues(i,STAR,STAR) *= mir1[i].GetMeasure() * ir_facet[i].Weight();

          IntRange trial_range = proxy1->IsOther() ? IntRange(fel1.GetNDof(), elmat.Width()) : IntRange(0, fel1.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(fel1.GetNDof(), elmat.Height()) : IntRange(0, fel1.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);

          // enum { BS = 16 };
          constexpr size_t BS = 16;
          for (int i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), loc_elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), loc_elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  if (proxy1->IsOther())
                    proxy1->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat1, lh);
                  else
                    proxy1->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat1, lh);
                  
                  if (proxy2->IsOther())
                    proxy2->Evaluator()->CalcMatrix(fel2, mir2[ii], bmat2, lh);
                  else
                    proxy2->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat2, lh);
                  
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(proxy1->IsOther() ? fel2 : fel1);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(proxy2->IsOther() ? fel2 : fel1);
              loc_elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }
  }



  void SymbolicFacetBilinearFormIntegrator ::
  CalcFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                   const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                   const ElementTransformation & strafo, FlatArray<int> & SElVertices1,
                   FlatMatrix<double> elmat,
                   LocalHeap & lh) const
  {
    // cout << "calc boundary facet matrix (DG)" << endl;
    elmat = 0.0;

    int maxorder = fel1.Order();

    auto eltype1 = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices1);
    Facet2SurfaceElementTrafo stransform(strafo.GetElementType(), SElVertices1); 
    
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    IntegrationRule & ir_facet_surf = stransform(ir_facet, lh);  // not yet used ???
    
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);          
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          if (proxy1->IsOther() || proxy2->IsOther()) continue;

          FlatTensor<3> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> Evaluate (mir1, val);
                proxyvalues(STAR,l,k) = val.Col(0);
              }

          for (int i = 0; i < mir1.Size(); i++)
            proxyvalues(i,STAR,STAR) *=  mir1[i].GetMeasure() * ir_facet[i].Weight();

          // IntRange trial_range = proxy1->IsOther() ? IntRange(fel1.GetNDof(), elmat.Width()) : IntRange(0, fel1.GetNDof());
          // IntRange test_range  = proxy2->IsOther() ? IntRange(fel1.GetNDof(), elmat.Height()) : IntRange(0, fel1.GetNDof());

          // auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

          // enum { BS = 16 };
          constexpr size_t BS = 16;          
          for (int i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  proxy1->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat2, lh);
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(fel1);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(fel1);
              elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }

  }


  void SymbolicFacetBilinearFormIntegrator ::
  CalcLinearizedFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                             const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                             const ElementTransformation & strafo, FlatArray<int> & SElVertices1,  
                             FlatVector<double> elveclin, FlatMatrix<double> elmat,
                             LocalHeap & lh) const
  {
    elmat = 0.0;

    int maxorder = fel1.Order();

    auto eltype1 = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices1);
    Facet2SurfaceElementTrafo stransform(strafo.GetElementType(), SElVertices1); 
    
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    IntegrationRule & ir_facet_surf = stransform(ir_facet, lh);  // not yet used ???
    
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);          
    
    // new 
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    ud.fel = &fel1;
    ud.elx = &elveclin;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      {
        ud.AssignMemory (proxy, ir_facet_vol1.Size(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel1, mir1, elveclin, ud.GetMemory(proxy), lh);
      }
    
    FlatMatrix<> val(mir1.Size(), 1, lh), deriv(mir1.Size(), 1, lh);
    elmat = 0;
    // endnew


    
    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          // FlatMatrix<> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          if (proxy1->IsOther() || proxy2->IsOther()) continue;

          FlatTensor<3> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> EvaluateDeriv (mir1, val, deriv);
                proxyvalues(STAR,l,k) = deriv.Col(0);
              }

          for (int i = 0; i < mir1.Size(); i++)
            proxyvalues(i,STAR,STAR) *=  mir1[i].GetMeasure() * ir_facet[i].Weight();

          // auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

          // enum { BS = 16 };
          constexpr size_t BS = 16;          
          for (int i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  proxy1->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel1, mir1[ii], bmat2, lh);
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(fel1);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(fel1);
              elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }
  }
  


  void SymbolicFacetBilinearFormIntegrator ::
  ApplyFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                    const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                    const FiniteElement & fel2, int LocalFacetNr2,
                    const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                    FlatVector<double> elx, FlatVector<double> ely,
                    LocalHeap & lh) const
  {
    // const TPHighOrderFE * tpfel = dynamic_cast<const TPHighOrderFE *>(&fel1);
    // if(tpfel)
    // {
      // int facet_x_y = 0;
      // if(LocalFacetNr1>=10)
      // {
        // facet_x_y = 1;
        // LocalFacetNr1-=10;
      // }
      // ApplyFacetMatrixTP(fel1,LocalFacetNr1,trafo1,ElVertices1,fel2,LocalFacetNr2,trafo2,ElVertices2,elx,ely,facet_x_y,lh);
      // return;
    // }
    if (simd_evaluate)
      {
        try
          {
            static Timer tall("SymbolicFacetBFI::Apply - all", 2); // RegionTimer rall(tall);
            static Timer tstart("SymbolicFacetBFI::Apply - startup", 2);
            static Timer tapply("SymbolicFacetBFI::Apply - apply", 2);
            static Timer tcoef("SymbolicFacetBFI::Apply - coef", 2);
            static Timer tapplyt("SymbolicFacetBFI::Apply - apply-trans", 2); 
            
            HeapReset hr(lh);
            // tall.Start();
            
            // tstart.Start();
            /*
              Matrix<> elmat(elx.Size());
              CalcFacetMatrix(fel1, LocalFacetNr1, trafo1, ElVertices1,
              fel2, LocalFacetNr2, trafo2, ElVertices2, elmat, lh);
              ely = elmat * elx;
              return;
            */
            
            ely = 0;
            
            int maxorder = max2 (fel1.Order(), fel2.Order());
            
            auto eltype1 = trafo1.GetElementType();
            auto eltype2 = trafo2.GetElementType();
            auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);
            
            Facet2ElementTrafo transform1(eltype1, ElVertices1); 
            Facet2ElementTrafo transform2(eltype2, ElVertices2); 
            
            SIMD_IntegrationRule simd_ir_facet(etfacet, 2*maxorder);
            
            auto & simd_ir_facet_vol1 = transform1(LocalFacetNr1, simd_ir_facet, lh);
            auto & simd_mir1 = trafo1(simd_ir_facet_vol1, lh);
            
            auto & simd_ir_facet_vol2 = transform2(LocalFacetNr2, simd_ir_facet, lh);
            auto & simd_mir2 = trafo2(simd_ir_facet_vol2, lh);
            
            simd_mir1.ComputeNormalsAndMeasure(eltype1, LocalFacetNr1);
            
            // evaluate proxy-values
            ProxyUserData ud(trial_proxies.Size(), lh);
            const_cast<ElementTransformation&>(trafo1).userdata = &ud;
            ud.fel = &fel1;   // necessary to check remember-map
            // ud.elx = &elx;
            ud.lh = &lh;
            for (ProxyFunction * proxy : trial_proxies)
              ud.AssignMemory (proxy, simd_ir_facet.GetNIP(), proxy->Dimension(), lh);
            // tstart.Stop();
            // tapply.Start();
            
            for (ProxyFunction * proxy : trial_proxies)
              {
                IntRange trial_range  = proxy->IsOther() ? IntRange(fel1.GetNDof(), elx.Size()) : IntRange(0, fel1.GetNDof());
                trial_range = proxy->Evaluator()->BlockDim() * trial_range;
                
                if (proxy->IsOther())
                  proxy->Evaluator()->Apply(fel2, simd_mir2, elx.Range(trial_range), ud.GetAMemory(proxy)); // , lh);
                else
                  proxy->Evaluator()->Apply(fel1, simd_mir1, elx.Range(trial_range), ud.GetAMemory(proxy)); // , lh);
                // tapply.AddFlops (trial_range.Size() * simd_ir_facet_vol1.GetNIP());
              }
            // tapply.Stop();
            
            for (auto proxy : test_proxies)
              {
                HeapReset hr(lh);
                // tcoef.Start();
                AFlatMatrix<double> simd_proxyvalues(proxy->Dimension(), simd_ir_facet.GetNIP(), lh);        
                
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    ud.testfunction = proxy;
                    ud.test_comp = k;
                    cf -> Evaluate (simd_mir1, simd_proxyvalues.Rows(k,k+1));
                  }
                
                for (int i = 0; i < simd_proxyvalues.Height(); i++)
                  {
                    auto row = simd_proxyvalues.Row(i);
                    for (int j = 0; j < row.VSize(); j++)
                      row.Get(j) *= simd_mir1[j].GetMeasure().Data() * simd_ir_facet[j].Weight().Data();
                  }
                // tcoef.Stop();
                // tapplyt.Start();
                IntRange test_range  = proxy->IsOther() ? IntRange(fel1.GetNDof(), elx.Size()) : IntRange(0, fel1.GetNDof());
                int blockdim = proxy->Evaluator()->BlockDim();
                test_range = blockdim * test_range;
                
                if (proxy->IsOther())
                  proxy->Evaluator()->AddTrans(fel2, simd_mir2, simd_proxyvalues, ely.Range(test_range));
                else
                  proxy->Evaluator()->AddTrans(fel1, simd_mir1, simd_proxyvalues, ely.Range(test_range));
                // tapplyt.AddFlops (test_range.Size() * simd_ir_facet_vol1.GetNIP());                
                // tapplyt.Stop();
              }
            // tall.Stop();
          }
        catch (ExceptionNOSIMD e)
          {
            cout << IM(4) << "caught in SymbolicFacetInegtrator::Apply: " << endl
                 << e.What() << endl;
            simd_evaluate = false;
            ApplyFacetMatrix (fel1, LocalFacetNr1, trafo1, ElVertices1,
                              fel2, LocalFacetNr2, trafo2, ElVertices2,
                              elx, ely, lh);
          }
        return;
      }
        
    
    static Timer tall("SymbolicFacetBFI::Apply - all", 2); RegionTimer rall(tall);
    /*
    static Timer t("SymbolicFacetBFI::Apply", 2);
    static Timer ts1("SymbolicFacetBFI::Apply start 1", 2);
    static Timer ts2("SymbolicFacetBFI::Apply start 2", 2);
    static Timer t1("SymbolicFacetBFI::Apply 1", 2);
    static Timer t2("SymbolicFacetBFI::Apply 2", 2);
    static Timer t3("SymbolicFacetBFI::Apply 3", 2);
    */
    
    HeapReset hr(lh);
    // ts1.Start();
    /*
    Matrix<> elmat(elx.Size());
    CalcFacetMatrix(fel1, LocalFacetNr1, trafo1, ElVertices1,
                    fel2, LocalFacetNr2, trafo2, ElVertices2, elmat, lh);
    ely = elmat * elx;
    return;
    */
    
    ely = 0;
    
    FlatVector<> ely1(ely.Size(), lh);

    int maxorder = max2 (fel1.Order(), fel2.Order());

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    
    Facet2ElementTrafo transform2(eltype2, ElVertices2); 
    IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, ir_facet, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(ir_facet_vol2, lh);


    // ts1.Stop();
    // ts2.Start();

    
    // evaluate proxy-values
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    ud.fel = &fel1;   // necessary to check remember-map
    // ud.elx = &elx;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir_facet.Size(), proxy->Dimension(), lh);

    for (ProxyFunction * proxy : trial_proxies)
      {
	IntRange trial_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*fel1.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim()*fel1.GetNDof());
	if (proxy->IsOther()) 
          proxy->Evaluator()->Apply(fel2, mir2, elx.Range(trial_range), ud.GetMemory(proxy), lh);
        else
          proxy->Evaluator()->Apply(fel1, mir1, elx.Range(trial_range), ud.GetMemory(proxy), lh);
      }

    // ts2.Stop();
    // RegionTimer reg(t);
    // t.Start();


    /*
        FlatVector<> measure(mir1.Size(), lh);
        switch (trafo1.SpaceDim())
          {
	  case 1:
            {
              Vec<1> normal_ref = ElementTopology::GetNormals<1>(eltype1)[LocalFacetNr1];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<1,1>&> (mir1[i]);
                  Mat<1> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<1> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<1,1>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
	      break;
            }
          case 2:
            {
              Vec<2> normal_ref = ElementTopology::GetNormals<2>(eltype1)[LocalFacetNr1];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<2,2>&> (mir1[i]);
                  Mat<2> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<2> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<2,2>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
              break;
            }
          default:
            cout << "Symbolic DG in " << trafo1.SpaceDim() << " not available" << endl;
          }
    */

    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);


    
    FlatMatrix<> val(ir_facet.Size(), 1,lh);

    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        // t1.Start();
        FlatMatrix<> proxyvalues(ir_facet.Size(), proxy->Dimension(), lh);
        
        // t1.Stop();
        // t2.Start();
        
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir1, val);
            proxyvalues.Col(k) = val.Col(0);
          }

        // t2.Stop();
        // t3.Start();

        for (int i = 0; i < mir1.Size(); i++)
          // proxyvalues.Row(i) *= measure(i) * ir_facet[i].Weight();
          proxyvalues.Row(i) *= mir1[i].GetMeasure() * ir_facet[i].Weight();
        
        ely1 = 0.0;
	IntRange test_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*fel1.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim()*fel1.GetNDof());
	if (proxy->IsOther())
	  proxy->Evaluator()->ApplyTrans(fel2, mir2, proxyvalues, ely1.Range(test_range), lh);
        else
          proxy->Evaluator()->ApplyTrans(fel1, mir1, proxyvalues, ely1.Range(test_range), lh);
        
        ely += ely1;
        // t3.Stop();
      }
    // t.Stop();
  }




  void SymbolicFacetBilinearFormIntegrator ::
  ApplyFacetMatrix (const FiniteElement & fel1, int LocalFacetNr,
                    const ElementTransformation & trafo1, FlatArray<int> & ElVertices,
                    const ElementTransformation & strafo, FlatArray<int> & SElVertices,
                    FlatVector<double> elx, FlatVector<double> ely,
                    LocalHeap & lh) const
  {
    // const TPHighOrderFE * tpfel = dynamic_cast<const TPHighOrderFE *>(&fel1);
    // if(tpfel)
    // {
      // int facet_x_y = 0;
      // if(LocalFacetNr>=10)
      // {
        // facet_x_y = 1;
        // LocalFacetNr-=10;
      // }
      // ApplyFacetMatrixTP(fel1,LocalFacetNr,trafo1,ElVertices,strafo,elx,ely,facet_x_y,lh);
      // return;
    // }    
    if (simd_evaluate)
      {
        try
          {
            static Timer t("SymbolicFacetBFI::ApplyFacetMatrix - boundary", 2);
            
            HeapReset hr(lh);
            
            ely = 0;
            
            int maxorder = fel1.Order();
            
            auto eltype1 = trafo1.GetElementType();
            auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);
            
            SIMD_IntegrationRule ir_facet(etfacet, 2*maxorder);
            Facet2ElementTrafo transform1(eltype1, ElVertices);
            Facet2SurfaceElementTrafo stransform(strafo.GetElementType(), SElVertices); 
            auto & ir_facet_vol1 = transform1(LocalFacetNr, ir_facet, lh);
            auto & ir_facet_surf = stransform(ir_facet, lh);
            auto & mir1 = trafo1(ir_facet_vol1, lh);
            auto & smir = strafo(ir_facet_surf, lh);

            mir1.ComputeNormalsAndMeasure(eltype1, LocalFacetNr);
            
            ProxyUserData ud(trial_proxies.Size(), lh);
            const_cast<ElementTransformation&>(trafo1).userdata = &ud;
            const_cast<ElementTransformation&>(strafo).userdata = &ud;
            ud.fel = &fel1;   // necessary to check remember-map
            // ud.elx = &elx;
            ud.lh = &lh;
            
            
            for (ProxyFunction * proxy : trial_proxies)
              ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
            
            for (ProxyFunction * proxy : trial_proxies)
              if (! (proxy->IsOther() && proxy->BoundaryValues()) )
                proxy->Evaluator()->Apply(fel1, mir1, elx, ud.GetAMemory(proxy));
            
            for (ProxyFunction * proxy : trial_proxies)
              if (proxy->IsOther() && proxy->BoundaryValues())
                proxy->BoundaryValues()->Evaluate (smir, ud.GetAMemory(proxy));
            
            RegionTimer reg(t);
            
            
            for (auto proxy : test_proxies)
              {
                HeapReset hr(lh);
                AFlatMatrix<double> proxyvalues(proxy->Dimension(), ir_facet.GetNIP(), lh);
                
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    ud.testfunction = proxy;
                    ud.test_comp = k;
                    cf -> Evaluate (mir1, proxyvalues.Rows(k,k+1));            
                  }
                
                for (int i = 0; i < proxyvalues.Height(); i++)
                  {
                    auto row = proxyvalues.Row(i);
                    for (int j = 0; j < row.VSize(); j++)
                      row.Get(j) *= mir1[j].GetMeasure().Data() * ir_facet[j].Weight().Data();
                  }
                
                if (proxy->IsOther() && proxy->BoundaryValues())
                  ; // nothing to do 
                else
                  proxy->Evaluator()->AddTrans(fel1, mir1, proxyvalues, ely);
              }
            
          }
        catch (ExceptionNOSIMD e)
          {
            cout << IM(4) << "caught in SymbolicFacetInegtrator::ApplyBnd: " << endl
                 << e.What() << endl;
            simd_evaluate = false;
            ApplyFacetMatrix (fel1, LocalFacetNr, trafo1, ElVertices,
                              strafo, SElVertices, 
                              elx, ely, lh);
          }
        return;
        
      }

    
    static Timer t("SymbolicFacetBFI::ApplyFacetMatrix - boundary", 2);
    
    HeapReset hr(lh);

    /*
    Matrix<> elmat(elx.Size());
    CalcFacetMatrix(fel1, LocalFacetNr, trafo1, ElVertices, strafo, elmat, lh);
    ely = elmat * elx;
    return;
    */

    ely = 0;
    
    FlatVector<> ely1(ely.Size(), lh);

    int maxorder = fel1.Order();

    auto eltype1 = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

    IntegrationRule ir_facet(etfacet, 2*maxorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices);
    Facet2SurfaceElementTrafo stransform(strafo.GetElementType(), SElVertices); 
    
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr, ir_facet, lh);
    IntegrationRule & ir_facet_surf = stransform(ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    auto & smir = strafo(ir_facet_surf, lh);
    
    // evaluate proxy-values
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    const_cast<ElementTransformation&>(strafo).userdata = &ud;
    ud.fel = &fel1;   // necessary to check remember-map
    // ud.elx = &elx;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir_facet.Size(), proxy->Dimension(), lh);
    
    for (ProxyFunction * proxy : trial_proxies)
      if (! (proxy->IsOther() && proxy->BoundaryValues()))
        proxy->Evaluator()->Apply(fel1, mir1, elx, ud.GetMemory(proxy), lh);

    for (ProxyFunction * proxy : trial_proxies)    
      if (proxy->IsOther() && proxy->BoundaryValues())
        proxy->BoundaryValues()->Evaluate (smir, ud.GetMemory(proxy));

    RegionTimer reg(t);
    
    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr);
    
    FlatMatrix<> val(ir_facet.Size(), 1,lh);
    for (auto proxy : test_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(ir_facet.Size(), proxy->Dimension(), lh);
        /*
        FlatVector<> measure(mir1.Size(), lh);
        switch (trafo1.SpaceDim())
          {
          case 1:
            {
              Vec<1> normal_ref = ElementTopology::GetNormals<1>(eltype1)[LocalFacetNr];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<1,1>&> (mir1[i]);
                  Mat<1> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<1> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<1,1>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
                break;
            }
          case 2:
            {
              Vec<2> normal_ref = ElementTopology::GetNormals<2>(eltype1)[LocalFacetNr];
              for (int i = 0; i < mir1.Size(); i++)
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<2,2>&> (mir1[i]);
                  Mat<2> inv_jac = mip.GetJacobianInverse();
                  double det = mip.GetMeasure();
                  Vec<2> normal = det * Trans (inv_jac) * normal_ref;       
                  double len = L2Norm (normal);    // that's the surface measure 
                  normal /= len;                   // normal vector on physical element
                  const_cast<MappedIntegrationPoint<2,2>&> (mip).SetNV(normal);
                  measure(i) = len;
                }
                break;
            }
          default:
            cout << "Symbolic DG in " << trafo1.SpaceDim() << " not available" << endl;
          }
        */
        
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir1, val);
            proxyvalues.Col(k) = val.Col(0);
          }

        for (int i = 0; i < mir1.Size(); i++)
          // proxyvalues.Row(i) *= measure(i) * ir_facet[i].Weight();
          proxyvalues.Row(i) *= mir1[i].GetMeasure() * ir_facet[i].Weight();

        ely1 = 0.0;
        if (proxy->IsOther() && proxy->BoundaryValues())
          ;  // nothing to do 
        else
          proxy->Evaluator()->ApplyTrans(fel1, mir1, proxyvalues, ely1, lh);
        ely += ely1;
      }
  }
  

  
void TensorProductFacetBilinearFormIntegrator :: 
   ApplyFacetMatrix (const FiniteElement & avolumefel1, int LocalFacetNr1,
                       const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                       const FiniteElement & avolumefel2, int LocalFacetNr2,
                       const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                       FlatVector<double> elx, FlatVector<double> ely, 
                       LocalHeap & lh) const
   {
     int facet_x_y = 0;
     if(LocalFacetNr1>=10)
     {
       facet_x_y = 1;
       LocalFacetNr1-=10;
     }
     static Timer tall("SymbolicFacetBFI::ApplyFacetMatrixTP, inner facet total"); RegionTimer rall(tall);
     static Timer ttrialproxy("SymbolicFacetBFI::ApplyFacetMatrixTP, inner facet eval trial proxies");
     static Timer ttestproxy("SymbolicFacetBFI::ApplyFacetMatrixTP, inner facet eval test proxies");
     static Timer tcoeff("SymbolicFacetBFI::ApplyFacetMatrixTP, inner facet eval coeff");
     static Timer tstartup("SymbolicFacetBFI::ApplyFacetMatrixTP, inner facet startup");
     //RegionTracer(TaskManager::GetThreadId(), tall);
     tstartup.Start();
     HeapReset hr(lh);
     const TPElementTransformation & eltrans1 = dynamic_cast<const TPElementTransformation &>(trafo1);
     const TPElementTransformation & eltrans2 = dynamic_cast<const TPElementTransformation &>(trafo2);
     const TPHighOrderFE & volumefel1 = dynamic_cast<const TPHighOrderFE &>(avolumefel1);
     const TPHighOrderFE & volumefel2 = dynamic_cast<const TPHighOrderFE &>(avolumefel2);
     ely = 0;
     ArrayMem<const IntegrationRule *,2> irs(2);
     FlatVector<> ely1(ely.Size(), lh);
     int maxorder1 = max2 (volumefel1.elements[facet_x_y]->Order(), volumefel2.elements[facet_x_y]->Order());
     int maxorder2 = max2 (volumefel1.elements[1-facet_x_y]->Order(), volumefel2.elements[1-facet_x_y]->Order());
     auto eltype1 = eltrans1.GetTrafo(facet_x_y).GetElementType();
     auto eltype2 = eltrans2.GetTrafo(facet_x_y).GetElementType();

     auto eltypevol = eltrans1.GetTrafo(1-facet_x_y).GetElementType();
     auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

     const IntegrationRule & ir_facet = SelectIntegrationRule(etfacet, 2*maxorder1);
     const IntegrationRule & ir_vol1 = SelectIntegrationRule(eltypevol, 2*maxorder2);
     Facet2ElementTrafo transform1(eltype1, ElVertices1);
     const IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
     Facet2ElementTrafo transform2(eltype2, ElVertices2);
     const IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, ir_facet, lh);
     irs[facet_x_y] = &ir_facet_vol1;
     irs[1-facet_x_y] = &ir_vol1;    
     TPIntegrationRule tpir1(irs);
     BaseMappedIntegrationRule & tpmir1 = eltrans1(tpir1,lh);
     irs[facet_x_y] = &ir_facet_vol2;
     TPIntegrationRule tpir2(irs);
     BaseMappedIntegrationRule & tpmir2 = eltrans2(tpir2,lh);
     // int tpnip = ir_facet_vol1.Size()*ir_vol1.Size();
     dynamic_cast<TPMappedIntegrationRule &>(tpmir1).SetFacet(facet_x_y);
     dynamic_cast<TPMappedIntegrationRule &>(tpmir2).SetFacet(facet_x_y);
     // evaluate proxy-values
     ProxyUserData ud(trial_proxies.Size(), lh);
     const_cast<ElementTransformation&>(trafo1).userdata = &ud;
     ud.fel = &volumefel1;   // necessary to check remember-map
     // ud.elx = &elx;
     ud.lh = &lh;
     FlatVector<> coef_temp(volumefel2.GetNDof(),lh),elx1(volumefel2.GetNDof(),lh);
     for (ProxyFunction * proxy : trial_proxies)
       ud.AssignMemory (proxy, tpmir1.Size(), proxy->Dimension(), lh);
     tstartup.Stop();
     ttrialproxy.Start();
     for (ProxyFunction * proxy : trial_proxies)
     {
       IntRange trial_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*volumefel1.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim()*volumefel1.GetNDof());       
       if (proxy->IsOther()) 
         proxy->Evaluator()->Apply(volumefel2, tpmir2, elx.Range(trial_range), ud.GetMemory(proxy), lh);
       else
         proxy->Evaluator()->Apply(volumefel1, tpmir1, elx.Range(trial_range), ud.GetMemory(proxy), lh);       
     }
     ttrialproxy.Stop();
     dynamic_cast<TPMappedIntegrationRule &>(tpmir1).GetIRs()[facet_x_y]->ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);
     FlatMatrixFixWidth<1> val(tpmir1.Size(),lh);
     for (auto proxy : test_proxies)
     {
       HeapReset hr(lh);
       FlatMatrix<> proxyvalues(tpmir1.Size(), proxy->Dimension(), lh);
       for (int k = 0; k < proxy->Dimension(); k++)
       {
         tcoeff.Start();
         ud.testfunction = proxy;
         ud.test_comp = k;
         cf -> Evaluate (tpmir1, val);
         proxyvalues.Col(k) = val.Col(0);
       }
       TPMappedIntegrationRule & ttpmir = dynamic_cast<TPMappedIntegrationRule &>(tpmir1);
       int ii=0;
       for(int i=0;i<ttpmir.GetIRs()[0]->Size();i++)
         for(int j=0;j<ttpmir.GetIRs()[1]->Size();j++)
           proxyvalues.Row(ii++) *= (*ttpmir.GetIRs()[0])[i].GetWeight()*(*ttpmir.GetIRs()[1])[j].GetWeight();
       tcoeff.Stop();
       ttestproxy.Start();
       ely1 = 0.0;
       IntRange test_range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*volumefel1.GetNDof(), elx.Size()) : IntRange(0, proxy->Evaluator()->BlockDim()*volumefel1.GetNDof());
       if (proxy->IsOther())
         proxy->Evaluator()->ApplyTrans(volumefel2, tpmir2, proxyvalues, ely1.Range(test_range), lh);     
       else
         proxy->Evaluator()->ApplyTrans(volumefel1, tpmir1, proxyvalues, ely1.Range(test_range), lh);
       ely += ely1;
       ttestproxy.Stop();
     }
   }
 
   void TensorProductFacetBilinearFormIntegrator :: 
   ApplyFacetMatrix (const FiniteElement & fel, int LocalFacetNr,
                       const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                       const ElementTransformation & seltrans, FlatArray<int> & SElVertices, 
                       FlatVector<double> elx, FlatVector<double> ely, 
                       LocalHeap & lh) const
   {
     int xfacet = 0;
     if(LocalFacetNr>=10)
     {
       xfacet = 1;
       LocalFacetNr-=10;
     }   
     static Timer tall("SymbolicFacetBFI::ApplyFacetMatrixTP, boundary facet"); RegionTimer rall(tall);  
     HeapReset hr(lh);
     ely = 0;
     FlatVector<> ely1(ely.Size(), lh);
     const TPHighOrderFE & volumefel = dynamic_cast<const TPHighOrderFE &>(fel);
     int order_fac = volumefel.elements[xfacet]->Order();
     int order_vol = volumefel.elements[1-xfacet]->Order();
     auto eltype1 = dynamic_cast<const TPElementTransformation &>(eltrans).GetTrafo(xfacet).GetElementType();
     
     auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);
     auto etvol = dynamic_cast<const TPElementTransformation &>(eltrans).GetTrafo(1-xfacet).GetElementType();
     ArrayMem<const IntegrationRule *,2> irs(2);
     const IntegrationRule & ir_facet = SelectIntegrationRule(etfacet, 2*order_fac);
     Facet2ElementTrafo transform(eltype1, ElVertices); 
     const IntegrationRule & ir_facet_vol = transform(LocalFacetNr, ir_facet, lh);
     const IntegrationRule &ir_vol = SelectIntegrationRule(etvol,2*order_vol);
     irs[xfacet] = &ir_facet_vol;
     irs[1-xfacet] = &ir_vol;
     TPIntegrationRule tpir(irs);
     BaseMappedIntegrationRule & tpmir = eltrans(tpir, lh);
     auto & smir = seltrans(ir_facet, lh);
     TPMappedIntegrationRule tpsmir(tpir,eltrans);
     //tpsmir.GetIRs().SetSize(2);
     tpsmir.GetIRs()[xfacet] = &smir;
     tpsmir.GetIRs()[1-xfacet] = dynamic_cast<TPMappedIntegrationRule &>(tpmir).GetIRs()[1-xfacet];
     // evaluate proxy-values
     ProxyUserData ud(trial_proxies.Size(), lh);
     // int niptp=irs[0]->Size()*irs[1]->Size();
     dynamic_cast<TPMappedIntegrationRule &>(tpmir).SetFacet(xfacet);
     dynamic_cast<TPMappedIntegrationRule &>(tpsmir).SetFacet(xfacet);
     const_cast<ElementTransformation&>(eltrans).userdata = &ud;
     const_cast<ElementTransformation&>(seltrans).userdata = &ud;
     ud.fel = &volumefel;   // necessary to check remember-map
     // ud.elx = &elx;
     ud.lh = &lh;
     for (ProxyFunction * proxy : trial_proxies)
       ud.AssignMemory (proxy, tpmir.Size(), proxy->Dimension(), lh);
     for (ProxyFunction * proxy : trial_proxies)
       if (! (proxy->IsOther() && proxy->BoundaryValues()))
         proxy->Evaluator()->Apply(volumefel, tpmir, elx, ud.GetMemory(proxy), lh);
     for (ProxyFunction * proxy : trial_proxies)    
       if (proxy->IsOther() && proxy->BoundaryValues())
         proxy->BoundaryValues()->Evaluate (tpsmir, ud.GetMemory(proxy));
 
     // RegionTimer reg(t);
     dynamic_cast<TPMappedIntegrationRule &>(tpmir).GetIRs()[xfacet]->ComputeNormalsAndMeasure (eltype1, LocalFacetNr);
     FlatMatrixFixWidth<1> val(tpmir.Size(),lh);
     for (auto proxy : test_proxies)
       {
         HeapReset hr(lh);
         FlatMatrix<> proxyvalues(tpmir.Size(), proxy->Dimension(), lh);
         for (int k = 0; k < proxy->Dimension(); k++)
           {
             ud.testfunction = proxy;
             ud.test_comp = k;
             cf -> Evaluate (tpmir, val);
             proxyvalues.Col(k) = val.Col(0);
           }
         TPMappedIntegrationRule & ttpmir = dynamic_cast<TPMappedIntegrationRule &>(tpmir);
         int ii=0;
         for(int i=0;i<ttpmir.GetIRs()[0]->Size();i++)
           for(int j=0;j<ttpmir.GetIRs()[1]->Size();j++)
             proxyvalues.Row(ii++) *= ttpmir.GetIRs()[0]->operator[](i).GetWeight()*ttpmir.GetIRs()[1]->operator[](j).GetWeight();
 
         ely1 = 0.0;
         if (proxy->IsOther() && proxy->BoundaryValues())
           ;  // nothing to do 
         else
           proxy->Evaluator()->ApplyTrans(volumefel, tpmir, proxyvalues, ely1, lh);
         ely += ely1;
       }
  }  
  
  SymbolicEnergy :: SymbolicEnergy (shared_ptr<CoefficientFunction> acf,
                                    VorB avb)
    : cf(acf), vb(avb)
  {
    simd_evaluate = true;
    
    if (cf->Dimension() != 1)
      throw Exception ("SymblicEnergy needs scalar-valued CoefficientFunction");
    
    
    cf->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
          if (proxy) 
            {
              if (!proxy->IsTestFunction())
                {                                         
                  if (!trial_proxies.Contains(proxy))
                    trial_proxies.Append (proxy);
                }
            }
        });
  }
  

  void 
  SymbolicEnergy :: CalcLinearizedElementMatrix (const FiniteElement & fel,
                                                 const ElementTransformation & trafo, 
                                                 FlatVector<double> elveclin,
                                                 FlatMatrix<double> elmat,
                                                 LocalHeap & lh) const
  {
    // static Timer t("symbolicenergy - calclinearized", 2);
    // static Timer tint("symbolicenergy - calclinearized intrules", 2);
    // static Timer tapply("symbolicenergy - calclinearized apply", 2);
    // static Timer td("symbolicenergy - calclinearized dmats", 2);
    // static Timer tb("symbolicenergy - calclinearized bmats", 2);
    // static Timer tbd("symbolicenergy - calclinearized bd", 2);
    // static Timer tmult("symbolicenergy - calclinearized mult", 2);
    // RegionTimer reg(t);

    if (simd_evaluate)
      //if (false)
      {
        try
          {
            // tint.Start();
            IntegrationRule std_ir(trafo.GetElementType(), 2*fel.Order());
            auto & std_mir = trafo(std_ir, lh);

            SIMD_IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
            auto & mir = trafo(ir, lh);
            // tint.Stop();
            // tapply.Start();
            ProxyUserData ud(trial_proxies.Size(), lh);
            const_cast<ElementTransformation&>(trafo).userdata = &ud;
            ud.fel = &fel;
            ud.elx = &elveclin;
            ud.lh = &lh;
            for (ProxyFunction * proxy : trial_proxies)
              {
                ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);
                proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetAMemory(proxy));
              }
            // tapply.Stop();
    
            AFlatMatrix<double> val(1, ir.GetNIP(),lh), deriv(1, ir.GetNIP(),lh), dderiv(1, ir.GetNIP(), lh);
            elmat = 0;
            // td.Start();            
            FlatArray<AFlatMatrix<double>> diags(trial_proxies.Size(), lh);
            for (int k1 : Range(trial_proxies))
              {
                auto proxy = trial_proxies[k1];
                new(&diags[k1]) AFlatMatrix<double>(proxy->Dimension(), ir.GetNIP(), lh);
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    ud.trialfunction = proxy;
                    ud.trial_comp = k;
                    ud.testfunction = proxy;
                    ud.test_comp = k;
                    cf -> EvaluateDDeriv (mir, val, deriv, diags[k1].Rows(k,k+1));
                  }
              }
            // td.Stop();
            for (int k1 : Range(trial_proxies))
              for (int l1 : Range(trial_proxies))
                {
                  HeapReset hr(lh);
                  auto proxy1 = trial_proxies[k1];
                  auto proxy2 = trial_proxies[l1];
                  // td.Start();

                  FlatTensor<3> proxyvalues(lh, ir.GetNIP(), proxy1->Dimension(), proxy2->Dimension());
                  
                  for (int k = 0; k < proxy1->Dimension(); k++)
                    for (int l = 0; l < proxy2->Dimension(); l++)
                      {
                        ud.trialfunction = proxy1;
                        ud.trial_comp = k;
                        ud.testfunction = proxy2;
                        ud.test_comp = l;
                        
                        cf -> EvaluateDDeriv (mir, val, deriv, dderiv);
                        proxyvalues(STAR,k,l) = dderiv.Row(0);
                        
                        if (proxy1 != proxy2 || k != l)  // computed mixed second derivatives
                          {
                            proxyvalues(STAR,k,l) -= diags[k1].Row(k);
                            proxyvalues(STAR,k,l) -= diags[l1].Row(l);
                            proxyvalues(STAR,k,l) *= 0.5;
                          }
                      }
                  // td.Stop();

                  for (int i = 0; i < ir.GetNIP(); i++)
                    proxyvalues(i,STAR,STAR) *= std_mir[i].GetWeight();
                  
                  // t.AddFlops (double (ir.GetNIP()) * proxy1->Dimension()*elmat.Width()*elmat.Height());
                  
                  FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                  FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

                  int i = 0;
                  IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                  IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
                  SliceMatrix<> part_elmat = elmat.Rows(r2).Cols(r1);
                  
                  // enum { BS = 16 };
                  constexpr size_t BS = 16;
                  for ( ; i < ir.GetNIP(); i+=BS)
                    {
                      HeapReset hr(lh);
                      int bs = min2(size_t(BS), ir.GetNIP()-i);

                      AFlatMatrix<double> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);                      
                      AFlatMatrix<double> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                      AFlatMatrix<double> bbmat2(elmat.Height(), bs*proxy2->Dimension(), lh);

                      // tb.Start();
                      BaseMappedIntegrationRule & bmir = std_mir.Range(i, i+bs, lh);
                      proxy1->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat1), lh);
                      proxy2->Evaluator()->CalcMatrix(fel, bmir, Trans(bbmat2), lh);
                      // tb.Stop();
                      
                      // tbd.Start();
                      for (int j = 0; j < bs; j++)
                        {
                          IntRange rj2 = proxy2->Dimension() * IntRange(j,j+1);
                          IntRange rj1 = proxy1->Dimension() * IntRange(j,j+1);
                          // bdbmat1.Rows(r1).Cols(rj2) = bbmat1.Rows(r1).Cols(rj1) * Trans (proxyvalues(i+j,STAR,STAR));
                          MultMatMat (bbmat1.Rows(r1).Cols(rj1), proxyvalues(i+j,STAR,STAR), bdbmat1.Rows(r1).Cols(rj2));
                        }
                      // tbd.Stop();
                      
                      // tmult.Start();
                      AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), part_elmat);
                      // tmult.Stop();
                      // tmult.AddFlops (double(r1.Size())*r2.Size()*bbmat2.Width());
                    }
                }
          }
        
        catch (ExceptionNOSIMD e)
          {
            cout << IM(4) << e.What() << endl
                 << "switching back to standard evaluation (in SymbolicEnergy::CalcLinearized)" << endl;
            simd_evaluate = false;
            CalcLinearizedElementMatrix (fel, trafo, elveclin, elmat, lh);
          }
        return;
      }
    
    
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elveclin;
    ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      {
        ud.AssignMemory (proxy, ir.Size(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetMemory(proxy), lh);        
      }
    
    FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh), dderiv(mir.Size(), 1,lh);
    
    elmat = 0;
    


    FlatArray<FlatMatrix<>> diags(trial_proxies.Size(), lh);
    for (int k1 : Range(trial_proxies))
      {
        auto proxy = trial_proxies[k1];
        diags[k1].AssignMemory(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> EvaluateDDeriv (mir, val, deriv, dderiv);

            diags[k1].Col(k) = dderiv.Col(0);
          }
      }
           
    
    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(trial_proxies))
        {
          HeapReset hr(lh);
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = trial_proxies[l1];
          // td.Start();
          // Tensor<3> proxyvalues(mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> EvaluateDDeriv (mir, val, deriv, dderiv);
                proxyvalues(STAR,l,k) = dderiv.Col(0);
                
                if (proxy1 != proxy2 || k != l)  // computed mixed second derivatives
                  {
                    proxyvalues(STAR,l,k) -= diags[k1].Col(k);
                    proxyvalues(STAR,l,k) -= diags[l1].Col(l);
                    proxyvalues(STAR,l,k) *= 0.5;
                  }
              }
          // td.Stop();

          /*
          for (int i = 0; i < mir.Size(); i++)
            {
              HeapReset hr(lh);
              proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
              
              FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> dbmat1(proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
              
              proxy1->Evaluator()->CalcMatrix(fel, mir[i], bmat1, lh);
              proxy2->Evaluator()->CalcMatrix(fel, mir[i], bmat2, lh);
              dbmat1 = proxyvalues(i,STAR,STAR) * bmat1;
              elmat += Trans (bmat2) * dbmat1 | Lapack;
            }
          */
          
          for (int i = 0; i < mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();

          // t.AddFlops (double (mir.Size()) * proxy1->Dimension()*elmat.Width()*elmat.Height());

          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
          int i = 0;

          // enum { BS = 16 };
          constexpr size_t BS = 16;
          for ( ; i+BS <= mir.Size(); i+=BS)
            {
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bdbmat1(BS*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(BS*proxy2->Dimension(), elmat.Height(), lh);

              for (int j = 0; j < BS; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  // tb.Start();
                  proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                  // tb.Stop();
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }
              // tmult.Start();
              elmat += Trans (bbmat2) * bdbmat1 | Lapack;
              // tmult.Stop();
            }


          if (i < mir.Size())
            {
              HeapReset hr(lh);
              int rest = mir.Size()-i;
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);
              
              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  // tb.Start();
                  proxy1->Evaluator()->CalcMatrix(fel, mir[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel, mir[ii], bmat2, lh);
                  // tb.Stop();
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }
              // tmult.Start();
              elmat += Trans (bbmat2) * bdbmat1 | Lapack;
              // tmult.Stop();
            }
          /*
          for ( ; i < mir.Size(); i++)
            {
              HeapReset hr(lh);
              // proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();
              
              FlatMatrix<double,ColMajor> dbmat1(proxy2->Dimension(), elmat.Width(), lh);
              
              proxy1->Evaluator()->CalcMatrix(fel, mir[i], bmat1, lh);
              proxy2->Evaluator()->CalcMatrix(fel, mir[i], bmat2, lh);
              dbmat1 = proxyvalues(i,STAR,STAR) * bmat1;
              elmat += Trans (bmat2) * dbmat1 | Lapack;
            }
          */
        }
  }
  
  
  
  double SymbolicEnergy :: Energy (const FiniteElement & fel, 
                                   const ElementTransformation & trafo, 
                                   FlatVector<double> elx, 
                                   LocalHeap & lh) const
  {
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;

    for (ProxyFunction * proxy : trial_proxies)
      {
        // ud.remember[proxy] = Matrix<> (ir.Size(), proxy->Dimension());
        // proxy->Evaluator()->Apply(fel, mir, elveclin, ud.remember[proxy], lh);
        ud.AssignMemory (proxy, ir.Size(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
      }


    FlatMatrix<> values(mir.Size(), 1, lh);
    cf -> Evaluate(mir, values);

    double sum = 0;
    for (int i = 0; i < mir.Size(); i++)
      sum += mir[i].GetWeight() * values(i,0);
    return sum;
  }

  void
  SymbolicEnergy :: ApplyElementMatrix (const FiniteElement & fel, 
                                        const ElementTransformation & trafo, 
                                        const FlatVector<double> elx, 
                                        FlatVector<double> ely,
                                        void * precomputed,
                                        LocalHeap & lh) const
  {
    static Timer t("SymbolicEnergy::ApplyElementMatrix", 2); 
    // static Timer ts("SymbolicEnergy::ApplyElementMatrix start", 2);
    static Timer ta("SymbolicEnergy::ApplyElementMatrix apply", 2);
    static Timer tc("SymbolicEnergy::ApplyElementMatrix coef", 2);
    static Timer tt("SymbolicEnergy::ApplyElementMatrix applyT", 2); 
    
    ProxyUserData ud(trial_proxies.Size(), lh);        
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    ud.elx = &elx;
    ud.lh = &lh;


    if (simd_evaluate)
      {
        try
          {
            RegionTimer reg(t);            
            // ts.Start();
            HeapReset hr(lh);
            SIMD_IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
            auto & mir = trafo(ir, lh);
            
            for (ProxyFunction * proxy : trial_proxies)
              ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);
            // ts.Stop();
            // ta.Start();
            for (ProxyFunction * proxy : trial_proxies)
              proxy->Evaluator()->Apply(fel, mir, elx, ud.GetAMemory(proxy));
            // ta.Stop();
            
            ely = 0;
            AFlatMatrix<double> val(1, ir.GetNIP(), lh);
            for (auto proxy : trial_proxies)
              {
                HeapReset hr(lh);
                AFlatMatrix<double> proxyvalues(proxy->Dimension(), ir.GetNIP(), lh);
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    // RegionTimer reg(tc);
                    ud.trialfunction = proxy;
                    ud.trial_comp = k;
                    cf -> EvaluateDeriv (mir, val, proxyvalues.Rows(k,k+1));
                  }

                for (int i = 0; i < proxyvalues.Height(); i++)
                  {
                    auto row = proxyvalues.Row(i);
                    for (int j = 0; j < row.VSize(); j++)
                      row.Get(j) *= mir[j].GetMeasure().Data() * ir[j].Weight().Data();
                  }

                // tt.Start();
                
                proxy->Evaluator()->AddTrans(fel, mir, proxyvalues, ely);
                // tt.Stop();
              }
          }
        catch (ExceptionNOSIMD e)
          {
            cout << IM(4) << e.What() << endl
                 << "switching back to standard evaluation (in SymbolicEnergy::CalcLinearized)" << endl;              
            simd_evaluate = false;
            ApplyElementMatrix (fel, trafo, elx, ely, precomputed, lh);
          }
        return;
      }
    
    HeapReset hr(lh);
    IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);

    for (ProxyFunction * proxy : trial_proxies)
      proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
    
    ely = 0;
    FlatVector<> ely1(ely.Size(), lh);

    FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh);
      
    for (auto proxy : trial_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            cf -> EvaluateDeriv (mir, val, deriv);
            proxyvalues.Col(k) = deriv.Col(0);
          }
        
        for (int i = 0; i < mir.Size(); i++)
          proxyvalues.Row(i) *= mir[i].GetWeight();
        
        proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, ely1, lh);
        ely += ely1;
      }
  }
  
  
}

