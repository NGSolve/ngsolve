/*********************************************************************/
/* File:   symbolicintegrator.cpp                                    */
/* Author: Joachim Schoeberl                                         */
/* Date:   August 2015                                               */
/*********************************************************************/
/* 
   Symbolic integrators
*/

#include <cassert>
#include <variant>
// #include <fem.hpp>
#include "integratorcf.hpp"
#include "symbolicintegrator.hpp"
#include "diffop.hpp"

namespace ngfem
{
  bool symbolic_integrator_uses_diff = false;
  
  ProxyFunction ::
  ProxyFunction (shared_ptr<ngcomp::FESpace> afes,
                 bool atestfunction, bool ais_complex,
                 shared_ptr<DifferentialOperator> aevaluator, 
                 shared_ptr<DifferentialOperator> aderiv_evaluator,
                 shared_ptr<DifferentialOperator> atrace_evaluator,
                 shared_ptr<DifferentialOperator> atrace_deriv_evaluator,
		 shared_ptr<DifferentialOperator> attrace_evaluator,
		 shared_ptr<DifferentialOperator> attrace_deriv_evaluator)
    
    : CoefficientFunction(aevaluator ? aevaluator->Dim() : 1,  /* ais_complex */ false),
      fes(afes),
      testfunction(atestfunction), is_other(false),
      evaluator(aevaluator), 
      deriv_evaluator(aderiv_evaluator),
      trace_evaluator(atrace_evaluator), 
      trace_deriv_evaluator(atrace_deriv_evaluator),
      ttrace_evaluator(attrace_evaluator),
      ttrace_deriv_evaluator(attrace_deriv_evaluator)
  {
    /*
    if (!aevaluator)
      {
        cerr << "don't have a primal evaluator" << endl;
        if (trace_evaluator) cout << "trace-eval = " << typeid(*trace_evaluator).name() << endl;
        if (deriv_evaluator) cout << "deriv-eval = " << typeid(*deriv_evaluator).name() << endl;
        throw Exception("don't have a primal evaluator");
      }
    */
    /*
    if (deriv_evaluator || trace_deriv_evaluator)
      deriv_proxy = make_shared<ProxyFunction> (fes, testfunction, is_complex,
                                                deriv_evaluator, nullptr,
                                                trace_deriv_evaluator, nullptr,
						ttrace_deriv_evaluator, nullptr);
    */
    if (evaluator)
      SetDimensions (evaluator->Dimensions());
    else if (trace_evaluator)
      SetDimensions (trace_evaluator->Dimensions());
    else if (ttrace_evaluator)
      SetDimensions (ttrace_evaluator->Dimensions());
    else
      throw Exception("a proxy needs at least one evaluator");
    elementwise_constant = true;

    SetVariable(true);
  }

  string ProxyFunction :: GetDescription () const
  {
    return string(testfunction ? "test-function" : "trial-function")
      + string(" diffop = ")
      + (evaluator ? evaluator->Name() : (trace_evaluator ? trace_evaluator->Name() : string("???")))
      + ((OrderDt() > 0) ? " order-dt="+ToString(OrderDt()) : string());
  }

  shared_ptr<ProxyFunction> ProxyFunction :: Deriv() const
  {
    if (auto sp = deriv_proxy.lock())
      return sp;
    
    if (deriv_evaluator || trace_deriv_evaluator)
      {
        auto sp = make_shared<ProxyFunction> (fes, testfunction, is_complex,
                                              deriv_evaluator, nullptr,
                                              trace_deriv_evaluator, nullptr,
                                              ttrace_deriv_evaluator, nullptr);
        sp -> primaryproxy = dynamic_pointer_cast<ProxyFunction>(const_cast<ProxyFunction*>(this)->shared_from_this());
        sp->is_other = is_other;
        const_cast<ProxyFunction*>(this)->deriv_proxy = sp;
        return sp;
      }
        
    return nullptr;
  }
  
  shared_ptr<ProxyFunction> ProxyFunction :: Trace() const
  {
    shared_ptr<ProxyFunction> trace_proxy = nullptr;
    
    if (trace_evaluator)
      trace_proxy = make_shared<ProxyFunction> (fes, testfunction, is_complex,
                                                trace_evaluator, trace_deriv_evaluator,
                                                ttrace_evaluator, ttrace_deriv_evaluator,
                                                nullptr, nullptr);
    else if (auto trace_from_diffop = evaluator->GetTrace())
      trace_proxy = make_shared<ProxyFunction>(fes, testfunction, is_complex,
                                               trace_from_diffop, nullptr,
                                               nullptr, nullptr,
                                               nullptr, nullptr);
    else
      throw Exception (string("don't have a trace, primal evaluator = ")+
                       evaluator->Name());
    
    for (int i = 0; i < additional_diffops.Size(); i++)
      if (auto add_diffop_trace = additional_diffops[i]->GetTrace())
        trace_proxy->SetAdditionalEvaluator (additional_diffops.GetName(i), add_diffop_trace);
    trace_proxy->is_other = is_other;
    return trace_proxy;
  }

  shared_ptr<ProxyFunction> ProxyFunction :: Dt() const
  {
    if (dt.lock())
      return dt.lock();

    auto newdt = make_shared<ProxyFunction> (*this);
    newdt -> anti_dt = const_pointer_cast<ProxyFunction>(dynamic_pointer_cast<const ProxyFunction>(this->shared_from_this()));
    
    dt = newdt;
    return newdt;
  }

  shared_ptr<ProxyFunction> ProxyFunction :: AntiDt() const
  {
    return anti_dt;
  }

  int ProxyFunction :: OrderDt() const
  {
    int cnt = 0;
    const ProxyFunction * ptr = this;
    while ( (ptr = ptr->AntiDt().get()) ) cnt++;
    return cnt;
  }

  
  shared_ptr<ProxyFunction> ProxyFunction :: Other(shared_ptr<CoefficientFunction> _boundary_values) const
  {
    auto other = make_shared<ProxyFunction> (fes, testfunction, is_complex, evaluator, deriv_evaluator, trace_evaluator, trace_deriv_evaluator,ttrace_evaluator, ttrace_deriv_evaluator);
    other->is_other = true;
    // if (other->deriv_proxy)
    // other->deriv_proxy->is_other = true;
    other->boundary_values = _boundary_values;

    for (int i = 0; i < additional_diffops.Size(); i++)
      other->SetAdditionalEvaluator (additional_diffops.GetName(i), additional_diffops[i]);
    
    return other;
  }


  shared_ptr<ProxyFunction> ProxyFunction :: GetAdditionalProxy (string name) const
  {
    if (additional_proxies.Used(name))
      if (auto sp = additional_proxies[name].lock())
        return sp;
    
    if (additional_diffops.Used(name))
      {
        auto adddiffop = make_shared<ProxyFunction> (fes, testfunction, is_complex, additional_diffops[name], nullptr, nullptr, nullptr, nullptr, nullptr);
        if (is_other)
          adddiffop->is_other = true;
        adddiffop -> primaryproxy = dynamic_pointer_cast<ProxyFunction>(const_cast<ProxyFunction*>(this)->shared_from_this());
        additional_proxies.Set(name, adddiffop);
        return adddiffop;
      }
    return shared_ptr<ProxyFunction>();
  }

  void ProxyFunction ::
  GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    // auto dims = Dimensions();
    
    string header = "\n\
    {flatmatrix} {values};\n\
    ProxyUserData * {ud} = (ProxyUserData*)mir.GetTransformation().userdata;\n\
    {\n\
      // if (!{ud})\n                                                      \
      // throw Exception (\"cannot evaluate ProxyFunction without userdata\");\n \
          ";

    if(!testfunction) {
      header+=
"      if ({ud}->fel) {\n\
     //     if ({ud}->HasMemory ({this})) {\n";
      if(code.is_simd) {
        header += "auto x = {ud}->GetAMemory ({this});\n";
        header += "{values}.AssignMemory(x.Height(), x.Width(), &x(0,0));\n";
      } else {
        header += "auto x = {ud}->GetMemory ({this});\n";
        header += "{values}.AssignMemory(x.Height(), x.Width(), &x(0,0));\n";
      }
      header+=
"      //    }\n\
       //   else\n\
       //     throw Exception(\"userdata has no memory!\");\n\
      }\n";
    }
    header += "}\n";
    if(code.deriv || !testfunction)
      header += "const bool {has_values} = {ud}->HasMemory({this});\n";
    else
      header += "constexpr bool {has_values} = true;\n"; // always need values in case code.deriv==0 (i.e. evaluation of proxy)

    for (int i = 0; i < this->Dimension(); i++) {
      header += Var("comp", index,i,this->Dimensions()).Declare("{scal_type}", 0.0);
        if(!testfunction && code.deriv==2)
          {
          header += "if(( ({ud}->trialfunction == {this}) && ({ud}->trial_comp=="+ToLiteral(i)+"))\n"+
          " ||  (({ud}->testfunction == {this}) && ({ud}->test_comp=="+ToLiteral(i)+")))\n";
        }
        else
          header += "if({ud}->{comp_string}=="+ToLiteral(i)+" && {ud}->{func_string} == {this})\n";
        header += Var("comp", index,i,this->Dimensions()).S() + string("{get_component}") + " = 1.0;\n";
    }
    
    string body = "";

    if (code_uses_tensors)
      {
        body += "Tens<" + code.res_type;
        for (auto d : this->Dimensions())
          body += ',' + ToLiteral(d);
        body += "> var_" + ToLiteral(index) + ";\n";
      }
    
    for (int i = 0; i < this->Dimension(); i++) {
      if (!code_uses_tensors)
        body += Var(index, i, this->Dimensions()).Declare("{scal_type}", 0.0);   // why do we initizlize ? 
      body += Var(index, i, this->Dimensions()).Assign(CodeExpr("0.0"), false);
    }

  
    if(!testfunction) {
      body += "// if ({ud}->fel) {\n";
      body += "if ({has_values}) {\n";
      for (int i = 0; i < Dimension(); i++) {
        string var = Var(index, i, this->Dimensions());
        if(code.deriv)
          var += ".Value()";
        string values = "{values}";
        if(code.is_simd)
          values += "(" + ToLiteral(i) + ",i)";
        else
          values += "(i," + ToLiteral(i) + ")";
        
        body += var + " = " + values + ";\n";
      }
      
      if(code.deriv)
        body += "}\n";
      else
        body += "} else ";
    }
    body += "{\n";
    if(testfunction)
      for (int i = 0; i < Dimension(); i++) {
        if(code.deriv==0) body += Var(index,i,Dimensions()).Assign( Var("comp", index,i,Dimensions()), false );
        if(code.deriv==1) body += Var(index,i,Dimensions()).Assign( Var("comp", index,i,Dimensions()), false );
        if(code.deriv==2) body += Var(index,i,Dimensions()).Call("DValue","0")
                            .Assign( Var("comp", index,i, Dimensions()).Call("DValue","0"), false );
      }
    else
      for (int i = 0; i < this->Dimension(); i++) {
        if(code.deriv==0) body += Var(index,i,Dimensions()).Assign( Var("comp", index,i,Dimensions()), false );
        if(code.deriv>=1) body += Var(index,i,Dimensions()).Call("DValue","0").Assign( Var("comp", index,i,Dimensions()).Call("DValue","0"), false );
      };
    body += "}\n";
    if (!testfunction && code.deriv == 0)
      body += "\n";
  

    string func_string = testfunction ? "testfunction" : "trialfunction";
    string comp_string = testfunction ? "test_comp" : "trial_comp";
    std::map<string,string> variables;
    variables["ud"] = "tmp_"+ToLiteral(index)+"_0";
    variables["this"] = "reinterpret_cast<ProxyFunction*>("+code.AddPointer(this)+")";
    variables["func_string"] = testfunction ? "testfunction" : "trialfunction";
    variables["comp_string"] = testfunction ? "test_comp" : "trial_comp";
    variables["testfunction"] = ToLiteral(testfunction);

    variables["flatmatrix"] = code.is_simd ? "FlatMatrix<SIMD<double>>" : "FlatMatrix<double>";

    variables["has_values"] = Var("has_values", index).S();
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
        // evaluator->Apply (*ud->fel, mip, *ud->elx, result, *ud->lh);
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
            BareSliceMatrix<> hresult) const
  {
    auto result = hresult.AddSize(mir.Size(), Dimension());
    
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");

    // for shape-wise evaluate
    if (ud->trial_elvec && !testfunction)
      {
        const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (ud->fel);
        const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : *ud->fel;
        evaluator->Apply (fel_trial, mir, *ud->trial_elvec, result, *ud->lh);
        return;
      }
    if (ud->test_elvec && testfunction)
      {
        const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (ud->fel);
        const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : *ud->fel;        
        evaluator->Apply (fel_test, mir, *ud->test_elvec, result, *ud->lh);
        return;
      }


    /*
    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory (this))
          result = ud->GetMemory (this);
        else
          cerr << "Proxy evaluate, but no userdata" << endl;
          // evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
        return;
      }
    */
    if (ud->HasMemory (this))
      {
        result = ud->GetMemory (this);
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
            BareSliceMatrix<Complex> result) const
  {
    size_t dim = Dimension();
    STACK_ARRAY(double, hmem, ir.Size()*dim);
    FlatMatrix<> temp(ir.Size(), dim, &hmem[0]);
    Evaluate (ir, temp);
    result.AddSize(ir.Size(), dim) = temp;
  }
  

  void ProxyFunction ::
  Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<SIMD<double>> result) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    if (!ud) 
      throw Exception ("cannot evaluate ProxyFunction without userdata");

    /*
    if (!testfunction && ud->fel)
      {
        result.AddSize(Dimension(), mir.Size()) = BareSliceMatrix<SIMD<double>> (ud->GetAMemory (this));
        return;
      }
    */
    
    if (ud->HasMemory(this))
      {
        result.AddSize(Dimension(), mir.Size()) = BareSliceMatrix<SIMD<double>> (ud->GetAMemory (this));
        return;
      }
    
    result.AddSize(Dimension(), mir.Size()) = 0;
    if (ud->testfunction == this)
      result.Row(ud->test_comp).Range(0,mir.Size()) = 1;
    if (ud->trialfunction == this)
      result.Row(ud->trial_comp).Range(0,mir.Size()) = 1;
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
      result.Row(ud->test_comp).Range(0,mir.Size()) = 1;
    if (ud->trialfunction == this)
      result.Row(ud->trial_comp).Range(0,mir.Size()) = 1;
  }

  /*
  void ProxyFunction ::
  Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
            FlatArray<AFlatMatrix<double>*> input,
            AFlatMatrix<double> result) const
  {
    ProxyFunction::Evaluate (mir, result);
  }
  */
  
  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationPoint & ip,
            FlatVector<Complex> result) const
  {
    STACK_ARRAY(double, hmem, result.Size());
    FlatVector<> result_double(result.Size(), &hmem[0]);
    Evaluate (ip, result_double);
    result = result_double;
  }

  /*
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

    if (!testfunction && ud->fel)
      {
        if (ud->HasMemory(this))
          result = ud->GetMemory (this);          
        else
          {
            static bool first = true;
            if (first) cerr << "ProxyFunction::EvaluateDeriv ... should not be here" << endl;
            first = false;
            // evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
          }
      }
    
    if (ud->testfunction == this)
      result.Col(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv.Col(ud->trial_comp) = 1;
  }
  */

  /*
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
            // evaluator->Apply (*ud->fel, mir, *ud->elx, result, *ud->lh);
          }
      }
    if (ud->testfunction == this)
      deriv.Col(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv.Col(ud->trial_comp) = 1;
  }
  */

  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationRule & mir,
            BareSliceMatrix<AutoDiff<1,double>> values) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    assert (ud);
    // assert (ud->fel);

    size_t np = mir.Size();
    size_t dim = Dimension();
    
    values.AddSize(np, dim) = AutoDiff<1,double> (0.0);

    if (IsTrialFunction())
      {
        auto val = ud->GetMemory(this);
        for (size_t i = 0; i < dim; i++)
          for (size_t j = 0; j < np; j++)
            values(j,i).Value() = val(j,i);
      }

    if (ud->testfunction == this)
      {
        auto col = values.Col(ud->test_comp);        
        for (size_t i = 0; i < np; i++)
          col(i).Value() = 1;
      }
    if (ud->trialfunction == this)
      {
        auto col = values.Col(ud->trial_comp);
        for (size_t i = 0; i < np; i++)
          col(i).DValue(0) = 1;
      }
  }
  

  void ProxyFunction ::
  Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;
    assert (ud);
    // assert (ud->fel);

    size_t np = mir.Size();
    size_t dim = Dimension();
    
    values.AddSize(dim, np) = AutoDiff<1,SIMD<double>> (0.0);

    if (IsTrialFunction())
      {
        auto val = ud->GetAMemory(this);
        for (size_t i = 0; i < dim; i++)
          for (size_t j = 0; j < np; j++)
            values(i,j).Value() = val(i,j);
      }

    if (ud->testfunction == this)
      {
        auto row = values.Row(ud->test_comp);        
        for (size_t i = 0; i < np; i++)
          row(i).Value() = 1;
      }
    if (ud->trialfunction == this)
      {
        auto row = values.Row(ud->trial_comp);
        for (size_t i = 0; i < np; i++)
          row(i).DValue(0) = 1;
      }
  }

  void ProxyFunction ::
  Evaluate (const BaseMappedIntegrationRule & mir, 
            BareSliceMatrix<AutoDiffDiff<1,double>> values) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;

    assert (ud);
    // assert (ud->fel);

    size_t np = mir.Size();
    size_t dim = Dimension();

    values.AddSize(np, dim) = AutoDiffDiff<1,double> (0.0);

    if (!testfunction) 
      {
        auto val = ud->GetMemory(this);
        for (size_t i = 0; i < dim; i++)
          for (size_t j = 0; j < np; j++)
            values(j,i).Value() = val(j,i);
      }

    if (ud->testfunction == this)
      for (size_t i = 0; i < np; i++)
        values(i, ud->test_comp).DValue(0) = 1;
    if (ud->trialfunction == this)
      for (size_t i = 0; i < np; i++)
        values(i, ud->trial_comp).DValue(0) = 1;
  }


  void ProxyFunction ::
  Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
            BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const
  {
    ProxyUserData * ud = (ProxyUserData*)mir.GetTransformation().userdata;

    assert (ud);
    // assert (ud->fel);

    size_t np = mir.Size();
    size_t dim = Dimension();

    values.AddSize(dim, np) = AutoDiffDiff<1,SIMD<double>> (0.0);

    if (!testfunction) //  && ud->fel)
      {
        auto val = ud->GetAMemory(this);
        for (size_t i = 0; i < dim; i++)
          for (size_t j = 0; j < np; j++)
            values(i,j).Value() = val(i,j);
      }

    if (ud->testfunction == this)
      for (size_t i = 0; i < np; i++)
        values(ud->test_comp, i).DValue(0) = 1;
    if (ud->trialfunction == this)
      for (size_t i = 0; i < np; i++)
        values(ud->trial_comp, i).DValue(0) = 1;
  }

  
  /*
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
            // evaluator->Apply (*ud->fel, mir, *ud->elx, result);
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
            // evaluator->Apply (*ud->fel, mir, *ud->elx, result);
          }
      }
    if (ud->testfunction == this)
      deriv.Row(ud->test_comp) = 1;
    if (ud->trialfunction == this)
      deriv.Row(ud->trial_comp) = 1;
  }
  */
  
  void ProxyFunction ::  
  NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                  FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
  {
    nonzero = false;
    nonzero_deriv = false;
    nonzero_dderiv = false;

    if (ud.eval_deriv == 1)
      {
        if (!testfunction)
          nonzero = true;
        if (ud.testfunction == this)
          nonzero(ud.test_comp) = true;
        if (ud.trialfunction == this)
          nonzero_deriv(ud.trial_comp) = true;
        return;
      }
    
    if (ud.fel)
      {
        if (!testfunction)
          nonzero = true;
        if (ud.testfunction == this)
          nonzero_deriv(ud.test_comp) = true;
        if (ud.trialfunction == this)
          nonzero_deriv(ud.trial_comp) = true;
        return;
      }

    if (ud.testfunction == this)
      nonzero(ud.test_comp) = true;
    if (ud.trialfunction == this)
      nonzero(ud.trial_comp) = true;
  }

  void ProxyFunction ::
  NonZeroPattern (const class ProxyUserData & ud,
                  FlatVector<AutoDiffDiff<1,NonZero>> values) const
  {
    Vector<bool> nz(values.Size()), nzd(values.Size()), nzdd(values.Size());
    NonZeroPattern (ud, nz, nzd, nzdd);
    for (size_t i = 0; i < values.Size(); i++)
      {
        values(i).Value() = nz(i);
        values(i).DValue(0) = nzd(i);
        values(i).DDValue(0) = nzdd(i);
      }
  }



  
  shared_ptr<CoefficientFunction> ProxyFunction :: 
  Operator (const string & name) const
  {
    if (deriv_evaluator && deriv_evaluator->Name() == name)
      return Deriv();
    return GetAdditionalProxy (name); 
  }

  
  
  shared_ptr<CoefficientFunction> ProxyFunction :: 
  Operator (shared_ptr<DifferentialOperator> diffop) const
  {
    // whatever component is the diffop, take the component we are
    
    // first strip components
    while (auto compdiffop = dynamic_pointer_cast<CompoundDifferentialOperator> (diffop))
      diffop = compdiffop->BaseDiffOp();
    
    // figure out the compound-component from our base diffop
    Array<int> comps;
    auto hp = evaluator;
    while (auto comphp = dynamic_pointer_cast<CompoundDifferentialOperator> (hp))
      {
        comps.Append (comphp->Component());
        hp = comphp->BaseDiffOp();
      }
    // and now wrap ...
    for (int i = comps.Size()-1; i >= 0; i--)
      diffop = make_shared<CompoundDifferentialOperator> (diffop, comps[i]);
    
    auto proxy = make_shared<ProxyFunction> (fes, testfunction, is_complex, diffop, nullptr, nullptr, nullptr, nullptr, nullptr);
    proxy -> primaryproxy = dynamic_pointer_cast<ProxyFunction>(const_cast<ProxyFunction*>(this)->shared_from_this());    
    if (is_other)
      proxy->is_other = true;
    return proxy;
  }

  shared_ptr<CoefficientFunction>
  ProxyFunction :: Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const
  {
    //if (var == shape.get())
    if (dynamic_cast<const DiffShapeCF*>(var))                
      return evaluator->DiffShape (const_cast<ProxyFunction*>(this)->shared_from_this(), dir);
    else if (var == this)
      return dir;

    if (primaryproxy.get() == var)
      return dir -> Operator(evaluator);
    /*
    if (auto proxyvar = dynamic_cast<const ProxyFunction*> (var))
      {
        // cout << "I am a proxy" << endl;
        if (this == proxyvar->Deriv().get())
          return dynamic_pointer_cast<ProxyFunction> (dir) -> Deriv();
        if (this == proxyvar->Trace().get())
          {
            // cout << "I am the trace" << endl;
            return dynamic_pointer_cast<ProxyFunction> (dir) -> Trace();
          }
      }
    */

    /*if (Dimension() == 1)
      return make_shared<ConstantCoefficientFunction>(0);
    else
      {
        auto zero1 = make_shared<ConstantCoefficientFunction>(0);
        Array<shared_ptr<CoefficientFunction>> zero_array(Dimension());
        for (auto & z : zero_array)
          z = zero1;
        auto zerovec = MakeVectorialCoefficientFunction(move(zero_array));
        zerovec->SetDimensions(Dimensions());
        return zerovec;
        }*/
    return ZeroCF(Dimensions());

  }

  
  SymbolicLinearFormIntegrator ::
  SymbolicLinearFormIntegrator(shared_ptr<CoefficientFunction> acf, VorB avb,
                               VorB aelement_vb)
    : cf(acf), vb(avb), element_vb(aelement_vb)
  {
    simd_evaluate = true;
    
    if (cf->Dimension() != 1)
      throw Exception ("SymbolicLFI needs scalar-valued CoefficientFunction");
    cf->TraverseTree
      ([&] (CoefficientFunction & nodecf)
       {
         if (auto proxy = dynamic_cast<ProxyFunction*> (&nodecf))
           {
             if (!proxies.Contains(proxy))
               proxies.Append (proxy);
           }
         else
           if (nodecf.StoreUserData() && !gridfunction_cfs.Contains(&nodecf))
             gridfunction_cfs.Append (&nodecf);
       });

    for (auto proxy : proxies)
      if (!proxy->Evaluator()->SupportsVB(vb))
        throw Exception ("Testfunction does not support "+ToString(vb)+"-forms, maybe a Trace() operator is missing");

    cache_cfs = FindCacheCF(*cf);
    
    dcf_dtest.SetSize(proxies.Size());

    if (symbolic_integrator_uses_diff)
      for (auto i : Range(proxies))
      {
        try
          {
            CoefficientFunction::T_DJC cache;
            dcf_dtest[i] = cf->DiffJacobi(proxies[i], cache);
            // cout << "dcf_dtest = " << *dcf_dtest[i] << endl;
          }
        catch (const Exception& e)
          {
              cout << IM(5) << "dcf_dtest has thrown exception " << e.What() << endl;
          }
      }
  }

  /*
  template <typename SCAL> 
  void SymbolicLinearFormIntegrator ::
  T_CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatVector<SCAL> elvec,
                       LocalHeap & lh) const
  {
    // static Timer t("symbolicLFI - CalcElementVector", NoTracing); RegionTimer reg(t);
    
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

  /*
  static Timer telvec("SymbolicLFI::CalcElVec");
  static Timer telvec_mapping("SymbolicLFI::mapping");
  static Timer telvec_zero("SymbolicLFI::zero");
  static Timer telvec_applytrans("SymbolicLFI::applytrans");
  static Timer telvec_dvec("SymbolicLFI::dvec");
  */
  
  template <typename SCAL>   
  void SymbolicLinearFormIntegrator ::
  T_CalcElementVector (const FiniteElement & fel,
                       const ElementTransformation & trafo, 
                       FlatVector<SCAL> elvec,
                       LocalHeap & lh) const
  {
    if (element_vb != VOL)
      { // not yet simded
        elvec = 0;
    
        auto eltype = trafo.GetElementType();
        Facet2ElementTrafo transform(eltype, element_vb); 
        int nfacet = transform.GetNFacets();
        
        for (int k = 0; k < nfacet; k++)
          {
            HeapReset hr(lh);
            ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
            
            const IntegrationRule& ir_facet = GetIntegrationRule(etfacet, 2*fel.Order()+bonus_intorder);
            IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
            BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
            
            ProxyUserData ud(0, gridfunction_cfs.Size(), lh);
            for (CoefficientFunction * cf : gridfunction_cfs)
              ud.AssignMemory (cf, ir_facet.GetNIP(), cf->Dimension(), lh);
            
            const_cast<ElementTransformation&>(trafo).userdata = &ud;
            PrecomputeCacheCF(cache_cfs, mir, lh);

            // mir.ComputeNormalsAndMeasure (eltype, k);
            
            FlatVector<SCAL> elvec1(elvec.Size(), lh);
            FlatMatrix<SCAL> val(mir.Size(), 1,lh);
            for (auto j : Range(proxies))
              {
                HeapReset hr(lh);
                auto proxy = proxies[j];
                FlatMatrix<SCAL> proxyvalues(mir.Size(), proxy->Dimension(), lh);

                if (dcf_dtest[j])
                    dcf_dtest[j]->Evaluate (mir, proxyvalues);
                else
                  for (int k = 0; k < proxy->Dimension(); k++)
                    {
                      ud.testfunction = proxy;
                      ud.test_comp = k;
                      cf -> Evaluate (mir, val);
                      proxyvalues.Col(k) = val.Col(0);
                    }
                
                for (auto i : Range(mir))
                    proxyvalues.Row(i) *= ir_facet[i].Weight() * mir[i].GetMeasure();
                  // proxyvalues.Row(i) *= ir_facet[i].Weight() * measure(i);
                
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
            // static Timer t("symbolicLFI - CalcElementVector (SIMD)", NoTracing); RegionTimer reg(t);
            // size_t tid = TaskManager::GetThreadId();
            // NgProfiler::StartThreadTimer(telvec, tid);
            
            HeapReset hr(lh);
            // NgProfiler::StartThreadTimer(telvec_mapping, tid);
            const SIMD_IntegrationRule& ir = GetSIMDIntegrationRule(trafo.GetElementType(), 2*fel.Order()+bonus_intorder);
            auto & mir = trafo(ir, lh);
            // NgProfiler::StopThreadTimer(telvec_mapping, tid);
            
            // NgProfiler::StartThreadTimer(telvec_zero, tid);            
            // ProxyUserData ud;
            ProxyUserData ud(0, gridfunction_cfs.Size(), lh);
            const_cast<ElementTransformation&>(trafo).userdata = &ud;
            for (CoefficientFunction * cf : gridfunction_cfs)
              ud.AssignMemory (cf, ir.GetNIP(), cf->Dimension(), lh);
            
            PrecomputeCacheCF(cache_cfs, mir, lh);
            
            elvec = 0;
            // NgProfiler::StopThreadTimer(telvec_zero, tid);                        
            for (auto j : Range(proxies))
              {
                // NgProfiler::StartThreadTimer(telvec_dvec, tid);
                auto proxy = proxies[j];
                FlatMatrix<SIMD<SCAL>> proxyvalues(proxy->Dimension(), ir.Size(), lh);
                if (dcf_dtest[j])
                  dcf_dtest[j]->Evaluate (mir, proxyvalues);
                else
                  for (size_t k = 0; k < proxy->Dimension(); k++)
                    {
                      ud.testfunction = proxy;
                      ud.test_comp = k;

                      cf -> Evaluate (mir, proxyvalues.Rows(k,k+1));
                    }
                    
                for (auto i : Range(proxyvalues.Height()))
                  {
                    auto row = proxyvalues.Row(i);
                    for (auto j : Range(row.Size()))
                      row(j) *= mir[j].GetWeight();
                  }
                    
                // NgProfiler::StopThreadTimer(telvec_dvec, tid);
                
                // NgProfiler::StartThreadTimer(telvec_applytrans, tid);                                
                proxy->Evaluator()->AddTrans(fel, mir, proxyvalues, elvec);
                // NgProfiler::StopThreadTimer(telvec_applytrans, tid);
              }
            // NgProfiler::StopThreadTimer(telvec, tid);
          }
        catch (const ExceptionNOSIMD& e)
          {
            cout << IM(6) << e.What() << endl
                 << "switching back to standard evaluation" << endl;
            simd_evaluate = false;
            T_CalcElementVector (fel, trafo, elvec, lh);
          }
      }
    else
      {
        // static Timer t("symbolicLFI - CalcElementVector", NoTracing); RegionTimer reg(t);
        HeapReset hr(lh);
        // IntegrationRule ir(trafo.GetElementType(), 2*fel.Order());
        const IntegrationRule& ir = GetIntegrationRule(trafo.GetElementType(),2*fel.Order()+bonus_intorder);
        BaseMappedIntegrationRule & mir = trafo(ir, lh);
        
        FlatVector<SCAL> elvec1(elvec.Size(), lh);
        
        FlatMatrix<SCAL> values(ir.Size(), 1, lh);
        ProxyUserData ud(0, gridfunction_cfs.Size(), lh);
        const_cast<ElementTransformation&>(trafo).userdata = &ud;
        PrecomputeCacheCF(cache_cfs, mir, lh);

        for (CoefficientFunction * cf : gridfunction_cfs)
          ud.AssignMemory (cf, ir.GetNIP(), cf->Dimension(), lh);
        
        elvec = 0;
        for (auto j : Range(proxies))
          {
            auto proxy = proxies[j];
            FlatMatrix<SCAL> proxyvalues(ir.Size(), proxy->Dimension(), lh);
            if (dcf_dtest[j])
              dcf_dtest[j]->Evaluate (mir, proxyvalues);
            else
              for (int k = 0; k < proxy->Dimension(); k++)
                {
                  ud.testfunction = proxy;
                  ud.test_comp = k;
                  cf -> Evaluate (mir, values);
                  proxyvalues.Col(k) = values.Col(0);
                }
                
            for (int i = 0; i < mir.Size(); i++)
              proxyvalues.Row(i) *= mir[i].GetWeight();
                
            proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, elvec1, lh);
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
                                  VorB aelement_vb)
    : cf(acf), vb(avb), element_vb(aelement_vb)
  {
    simd_evaluate = true;
    
    if (cf->Dimension() != 1)
        throw Exception ("SymbolicBFI needs scalar-valued CoefficientFunction");
    trial_cum.Append(0);
    test_cum.Append(0);
    has_interpolate = false;
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
          else
            if (nodecf.StoreUserData() && !gridfunction_cfs.Contains(&nodecf))
              gridfunction_cfs.Append (&nodecf);

          if (nodecf.GetDescription() == "InterpolationCF")
            {
              has_interpolate = true;
              cout << IM(3) << "integrand has an Interpolation Operator" << endl;
            }
        });

    cache_cfs = FindCacheCF(*cf);

    for (auto proxy : trial_proxies)
      if (proxy->OrderDt() != 0)
        throw Exception("time-derivatives not allowed in spatial integrators");
    for (auto proxy : trial_proxies)
      if (!proxy->Evaluator()->SupportsVB(vb))
        throw Exception ("Trialfunction does not support "+ToString(vb)+"-forms, maybe a Trace() operator is missing, type = "+proxy->Evaluator()->Name());
    for (auto proxy : test_proxies)
      if (!proxy->Evaluator()->SupportsVB(vb))
        throw Exception ("Testfunction does not support "+ToString(vb)+"-forms, maybe a Trace() operator is missing");
    
    cout << IM(6) << "num test_proxies " << test_proxies.Size() << endl;
    cout << IM(6) << "num trial_proxies " << trial_proxies.Size() << endl;
    cout << IM(6) << "cumulated test_proxy dims  " << test_cum << endl;
    cout << IM(6) << "cumulated trial_proxy dims " << trial_cum << endl;

    elementwise_constant = cf -> ElementwiseConstant();
    cout << IM(6) << "element-wise constant = " << elementwise_constant << endl;

    // find non-zeros
    int cnttest = 0, cnttrial = 0;
    for (auto proxy : trial_proxies)
      cnttrial += proxy->Dimension();
    for (auto proxy : test_proxies)
      cnttest += proxy->Dimension();
    nonzeros = Matrix<bool>(cnttest, cnttrial);
    nonzeros_deriv = Matrix<bool>(cnttest, cnttrial);

    ProxyUserData ud;
    // Vector<bool> nzvec(1), nzdvec(1), nzddvec(1);
    Vector<AutoDiffDiff<1,NonZero>> nzvec(1);
    int k = 0;
    for (int k1 : test_proxies.Range())
      for (int k2 : Range(test_proxies[k1]->Dimension()))
        {
          int l = 0;
          for (int l1 : trial_proxies.Range())
            for (int l2 : Range(trial_proxies[l1]->Dimension()))
              {
                ud.trialfunction = trial_proxies[l1];
                ud.trial_comp = l2;
                ud.testfunction = test_proxies[k1];
                ud.test_comp = k2;
                cf -> NonZeroPattern (ud, nzvec);
                nonzeros(k,l) = nzvec(0).Value();
                l++;
              }
          k++;
        }

    // derivative ...
    DummyFE<ET_TRIG> dummyfe;
    ud.fel = &dummyfe;
    ud.eval_deriv = 1;
    k = 0;
    for (int k1 : test_proxies.Range())
      for (int k2 : Range(test_proxies[k1]->Dimension()))
        {
          int l = 0;
          for (int l1 : trial_proxies.Range())
            for (int l2 : Range(trial_proxies[l1]->Dimension()))
              {
                ud.trialfunction = trial_proxies[l1];
                ud.trial_comp = l2;
                ud.testfunction = test_proxies[k1];
                ud.test_comp = k2;
                cf -> NonZeroPattern (ud, nzvec); 
                // nonzeros_deriv(k,l) = nzdvec(0);
                nonzeros_deriv(k,l) = nzvec(0).DValue(0);
                l++;
              }
          k++;
        }


    
    nonzeros_proxies = Matrix<bool>(test_proxies.Size(), trial_proxies.Size());
    diagonal_proxies = Matrix<bool>(test_proxies.Size(), trial_proxies.Size());
    same_diffops = Matrix<bool>(test_proxies.Size(), trial_proxies.Size());
    is_symmetric = true;
    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          
          bool is_nonzero = false;
          bool is_diagonal = proxy1->Dimension() == proxy2->Dimension();
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              if (nonzeros(test_cum[l1]+l, trial_cum[k1]+k))
                {
                  is_nonzero = true;
                  if (k != l) is_diagonal = false;
                }
          nonzeros_proxies(l1, k1) = is_nonzero;
          diagonal_proxies(l1, k1) = is_diagonal;
          same_diffops(l1,k1) = *(proxy2->Evaluator()) == *(proxy1->Evaluator());
          // same_diffops(l1,k1) = proxy2->Evaluator() == proxy1->Evaluator();
          // same objects, not only equal objects, what implies also the same space
          // should be ok now, since we check also for MixedFE
          // important to recognize equatl CompoundDiffOps (which are never the same)
          if (nonzeros_proxies(l1,k1) && (!diagonal_proxies(l1,k1) || !same_diffops(l1,k1))) is_symmetric = false;
        }
    
    cout << IM(6) << "nonzeros: " << endl << nonzeros << endl;
    cout << IM(6) << "nonzeros_deriv: " << endl << nonzeros_deriv << endl;
    cout << IM(6) << "nonzeros_proxies: " << endl << nonzeros_proxies << endl;
    cout << IM(6) << "symmetric: " << endl << is_symmetric << endl;
    cout << IM(6) << "same diffops: " << endl << same_diffops << endl;
    trial_difforder = 99, test_difforder = 99;
    for (auto proxy : trial_proxies)
      trial_difforder = min2(trial_difforder, proxy->Evaluator()->DiffOrder());
    for (auto proxy : test_proxies)
      test_difforder = min2(test_difforder, proxy->Evaluator()->DiffOrder());
    if (trial_proxies.Size() == 0) trial_difforder = 0;

    dcf_dtest.SetSize(test_proxies.Size());
    ddcf_dtest_dtrial.SetSize(test_proxies.Size(), trial_proxies.Size());

    if (symbolic_integrator_uses_diff)
      {
        for (auto i : Range(test_proxies))
          {
            try
              {
                CoefficientFunction::T_DJC cache;
                dcf_dtest[i] = cf->DiffJacobi(test_proxies[i], cache);
              }
            catch (const Exception& e)
              {
                cout << IM(5) << "dcf_dtest has thrown exception " << e.What() << endl;
              }

            for (auto j : Range(trial_proxies))
              {
                if (dcf_dtest[i])
                  try
                    {
                      CoefficientFunction::T_DJC cache2;
                        ddcf_dtest_dtrial(i, j) = dcf_dtest[i]->DiffJacobi(trial_proxies[j], cache2);
                    }
                  catch (const Exception& e)
                    {
                      cout << IM(5) << "ddcf_dtest_dtrial has thrown exception " << e.What() << endl;
                    }
              }
          }
      }
  }


  const IntegrationRule& SymbolicBilinearFormIntegrator ::
  GetIntegrationRule (const FiniteElement & fel, LocalHeap & /* lh */) const
  {
    if (userdefined_intrules[fel.ElementType()]) return *userdefined_intrules[fel.ElementType()];
    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;

    /*
    int trial_difforder = 99, test_difforder = 99;
    for (auto proxy : trial_proxies)
      trial_difforder = min2(trial_difforder, proxy->Evaluator()->DiffOrder());
    for (auto proxy : test_proxies)
      test_difforder = min2(test_difforder, proxy->Evaluator()->DiffOrder());
    if (trial_proxies.Size() == 0) trial_difforder = 0;
    */
    
    int intorder = fel_trial.Order()+fel_test.Order()+bonus_intorder;
    auto et = fel.ElementType();
    if (et == ET_TRIG || et == ET_TET)
      intorder -= test_difforder+trial_difforder;
    return SelectIntegrationRule (et, intorder);
  }

  const SIMD_IntegrationRule& SymbolicBilinearFormIntegrator ::
  Get_SIMD_IntegrationRule (const FiniteElement & fel, LocalHeap & lh) const
  {
    if (userdefined_simd_intrules[fel.ElementType()] ) return *userdefined_simd_intrules[fel.ElementType()];

    bool is_mixed = typeid(fel) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe = static_cast<const MixedFiniteElement*> (&fel);    
    const FiniteElement & fel_trial = is_mixed ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = is_mixed ? mixedfe->FETest() : fel;

    /*
    int trial_difforder = 99, test_difforder = 99;
    for (auto proxy : trial_proxies)
      trial_difforder = min2(trial_difforder, proxy->Evaluator()->DiffOrder());
    for (auto proxy : test_proxies)
      test_difforder = min2(test_difforder, proxy->Evaluator()->DiffOrder());
    if (trial_proxies.Size() == 0) trial_difforder = 0;
    */
    
    int intorder = fel_trial.Order()+fel_test.Order()+bonus_intorder;
    auto et = fel.ElementType();
    if (et == ET_TRIG || et == ET_TET)
      intorder -= test_difforder+trial_difforder;
    return SIMD_SelectIntegrationRule (et, intorder);
  }



  template <typename SCAL>
  void ExtendSymmetric (SliceMatrix<SCAL> elmat)
  {
    /*
    size_t h = elmat.Height();
    for (size_t i = 0; i+1 < h; i++)
      for (size_t j = i+1; j < h; j++)
        elmat(i,j) = elmat(j,i);
    */

    size_t h = elmat.Height();
    size_t d = elmat.Dist();
    size_t i = 0, di = 0;
    for ( ; i+2 < h; i+=2, di+=2*d)
      {
        elmat(di+i+1) = elmat(di+d+i);
        size_t j = i+2, dj = di+2*d;
        for ( ; j+1 < h; j+=2, dj+=2*d)
          {
            SCAL tmp00 = elmat(dj+i);
            SCAL tmp01 = elmat(dj+i+1);
            SCAL tmp10 = elmat(dj+d+i);
            SCAL tmp11 = elmat(dj+d+i+1);
            elmat(di+j) = tmp00;
            elmat(di+d+j) = tmp01;
            elmat(di+j+1) = tmp10;
            elmat(di+d+j+1) = tmp11;
          }
        if (j < h)
          {
            SCAL tmp0 = elmat(dj+i);
            SCAL tmp1 = elmat(dj+i+1);
            elmat(di+j) = tmp0;
            elmat(di+d+j) = tmp1;
          }
      }
    /*
    for ( ; i+1 < h; i++)
      for (size_t j = i+1; j < h; j++)
        elmat(i,j) = elmat(j,i);
    */
    if (i+1 < h)
      elmat(di+i+1) = elmat(di+d+i);
  }

  void ExtendSymmetric1 (SliceMatrix<double> elmat)
  {
    ExtendSymmetric (elmat);
  }
  

  /*
  Timer timer_SymbBFI("SymbolicBFI");
  Timer timer_SymbBFIstart("SymbolicBFI start");
  Timer timer_SymbBFIscale("SymbolicBFI scale");
  Timer timer_SymbBFIbd("SymbolicBFI bd");
  Timer timer_SymbBFIbmat("SymbolicBFI bmat");
  Timer timer_SymbBFIdmat("SymbolicBFI dmat");
  Timer timer_SymbBFImult("SymbolicBFI mult");
  Timer timer_SymbBFImultsym("SymbolicBFI multsym");
  */

  template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcElementMatrixAdd (const FiniteElement & fel,
                          const ElementTransformation & trafo, 
                          FlatMatrix<SCAL_RES> elmat,
                          bool & symmetric_so_far,
                          LocalHeap & lh) const
    
  {
    static Timer t(string("SymbolicBFI::CalcElementMatrixAdd")+typeid(SCAL).name()+typeid(SCAL_SHAPES).name()+typeid(SCAL_RES).name(), NoTracing);
//    static Timer tdmat("SymbolicBFI::CalcDMat - simd", NoTracing);
//    static Timer tmult("SymbolicBFI::mult - simd", NoTracing);
    RegionTimer reg(t);
    // RegionTracer regtr(TaskManager::GetThreadId(), t);    

    auto save_userdata = trafo.PushUserData(); 
    
    if (element_vb != VOL)
      {
        // static Timer t(string("SymbolicBFI::EB ")+typeid(SCAL).name()+typeid(SCAL_SHAPES).name()+typeid(SCAL_RES).name(), NoTracing);
        // RegionTimer reg(t);
        // RegionTracer regtr(TaskManager::GetThreadId(), t);    
        
        T_CalcElementMatrixEBAdd<SCAL, SCAL_SHAPES, SCAL_RES> (fel, trafo, elmat, lh);
        if (!IsSymmetric().IsTrue()) symmetric_so_far = false;
        return;
      }

    if (has_interpolate)
      {
        T_CalcElementMatrixAddShapeWise<SCAL, SCAL_SHAPES, SCAL_RES> (fel, trafo, elmat, lh);
        if (!IsSymmetric().IsTrue()) symmetric_so_far = false;
        return;
      }
    

    bool is_mixedfe = typeid(fel) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe = static_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = is_mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = is_mixedfe ? mixedfe->FETest() : fel;
    // size_t first_std_eval = 0;
    if (simd_evaluate)
      try
        {
          // static Timer tsimd(string("SymbolicBFI::CalcElementMatrixAddSIMD")+typeid(SCAL).name()+typeid(SCAL_SHAPES).name()+typeid(SCAL_RES).name(), NoTracing);          
          // RegionTracer regsimd(TaskManager::GetThreadId(), tsimd);
 
          const SIMD_IntegrationRule& bigir = Get_SIMD_IntegrationRule (fel, lh);

          constexpr size_t BS = 64;
          for (size_t ii = 0; ii < bigir.Size(); ii+=BS)
            {
              HeapReset hr(lh);
              int bs = min2(BS, bigir.Size()-ii);
              auto ir = bigir.Range(ii, ii+bs);
          
          SIMD_BaseMappedIntegrationRule & mir = trafo(ir, lh);
          // NgProfiler::StopThreadTimer (timer_SymbBFIstart, TaskManager::GetThreadId());

          ProxyUserData ud;
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          PrecomputeCacheCF(cache_cfs, mir, lh);

          // bool symmetric_so_far = true;
          int k1 = 0;
          int k1nr = 0;
          
          for (auto proxy1 : trial_proxies)
            {
              int l1 = 0;
              int l1nr = 0;
              for (auto proxy2 : test_proxies)
                {
                  size_t dim_proxy1 = proxy1->Dimension();
                  size_t dim_proxy2 = proxy2->Dimension();
                  
                  size_t tt_pair = l1nr*trial_proxies.Size()+k1nr;
                  // first_std_eval = k1nr*test_proxies.Size()+l1nr;  // in case of SIMDException
                  bool is_nonzero = nonzeros_proxies(tt_pair);
                  bool is_diagonal = diagonal_proxies(tt_pair);

                  if (is_nonzero)
                    {
                      HeapReset hr(lh);
                      bool samediffop = same_diffops(tt_pair) && !is_mixedfe;
                      // td.Start();

                      FlatMatrix<SIMD<SCAL>> proxyvalues(dim_proxy1*dim_proxy2, ir.Size(), lh);
                      FlatMatrix<SIMD<SCAL>> diagproxyvalues(dim_proxy1, ir.Size(), lh);
                      FlatMatrix<SIMD<SCAL>> val(1, ir.Size(), lh);

                      IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                      IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                      SliceMatrix<SCAL_RES> part_elmat = elmat.Rows(r2).Cols(r1);

                      FlatMatrix<SIMD<SCAL_SHAPES>> bbmat1(elmat.Width() * dim_proxy1, ir.Size(), lh);
                      FlatMatrix<SIMD<SCAL>> bdbmat1(elmat.Width() * dim_proxy2, ir.Size(), lh);
                      FlatMatrix<SIMD<SCAL_SHAPES>> bbmat2 =
                              samediffop ?
                                bbmat1
                                :
                                FlatMatrix<SIMD<SCAL_SHAPES>>(elmat.Height() * dim_proxy2, ir.Size(), lh);

                      FlatMatrix<SIMD<SCAL>> hbdbmat1(elmat.Width(), dim_proxy2 * ir.Size(), bdbmat1.Data());
                      FlatMatrix<SIMD<SCAL_SHAPES>> hbbmat2(elmat.Height(), dim_proxy2 * ir.Size(), bbmat2.Data());

                      if (ddcf_dtest_dtrial(l1nr, k1nr))
                        {
//                          RegionTimer regdmat(tdmat);
//                          cout << "use ddcf_dtest_dtrial" << endl;
                          ddcf_dtest_dtrial(l1nr, k1nr)->Evaluate(mir, proxyvalues);

                          if (is_diagonal)
                            for (auto k : Range(dim_proxy1))
                              diagproxyvalues.Row(k) = proxyvalues.Row(k*(dim_proxy1 + 1));
                        }
                      else
                        {
//                          RegionTimer regdmat(tdmat);
                          if (!is_diagonal)
                            {
                              for (size_t k = 0, kk = 0; k < dim_proxy1; k++)
                                for (size_t l = 0; l < dim_proxy2; l++, kk++)
                                  {
                                    if (nonzeros(l1+l, k1+k))
                                      {
                                        ud.trialfunction = proxy1;
                                        ud.trial_comp = k;
                                        ud.testfunction = proxy2;
                                        ud.test_comp = l;

                                        cf -> Evaluate(mir, proxyvalues.Rows(kk,kk+1));
                                      }
                                    else;
                                      // proxyvalues.Row(kk) = 0.0;
                                  }
                            }
                          else
                            {
                              for (size_t k = 0; k < dim_proxy1; k++)
                                {
                                  ud.trialfunction = proxy1;
                                  ud.trial_comp = k;
                                  ud.testfunction = proxy2;
                                  ud.test_comp = k;

                                  cf -> Evaluate (mir, diagproxyvalues.Rows(k,k+1));
                                }
                            }
                          // td.Stop();
                        }



                      // NgProfiler::StartThreadTimer (timer_SymbBFIscale, TaskManager::GetThreadId());
                      FlatVector<SIMD<double>> weights(ir.Size(), lh);
                      if (!is_diagonal)
                        for (size_t i = 0; i < ir.Size(); i++)
                          // proxyvalues.Col(i) *= mir[i].GetWeight();
                          weights(i) = mir[i].GetWeight();
                      else
                        for (size_t i = 0; i < ir.Size(); i++)
                          diagproxyvalues.Col(i) *= mir[i].GetWeight();

                      // NgProfiler::StopThreadTimer (timer_SymbBFIscale, TaskManager::GetThreadId());
                      // bbmat1 = 0.0;
                      // bbmat2 = 0.0;
                      {
                        // RegionTimer regbmat(timer_SymbBFIbmat);
                        proxy1->Evaluator()->CalcMatrix(fel_trial, mir, bbmat1);
                        if (!samediffop)
                          proxy2->Evaluator()->CalcMatrix(fel_test, mir, bbmat2);
                      }

                      if (is_diagonal)
                        {
                          // static Timer t("diag DB", NoTracing);
                          // RegionTracer reg(TaskManager::GetThreadId(), t);
                          // NgProfiler::StartThreadTimer (timer_SymbBFIbd, TaskManager::GetThreadId());                      
                          
                          /*
                          size_t ii = r1.First()*dim_proxy1;
                          for (size_t i : r1)
                            for (size_t j = 0; j < dim_proxy1; j++, ii++)
                              bdbmat1.Row(ii) = pw_mult(bbmat1.Row(ii), diagproxyvalues.Row(j));
                          */
                          
                          // size_t sr1 = r1.Size();
                          for (size_t j = 0; j < dim_proxy1; j++)
                            {
                              auto hbbmat1 = bbmat1.RowSlice(j,dim_proxy1).Rows(r1);
                              auto hbdbmat1 = bdbmat1.RowSlice(j,dim_proxy1).Rows(r1);
                              
                              for (size_t k = 0; k < bdbmat1.Width(); k++)
                                hbdbmat1.Col(k).Range(0,r1.Size()) = diagproxyvalues(j,k) * hbbmat1.Col(k);
                            }
                        }
                      else
                        {
                          // static Timer t("DB", NoTracing);
                          // RegionTracer reg(TaskManager::GetThreadId(), t);
                          
                          // bdbmat1 = 0.0;
                          hbdbmat1.Rows(r1) = 0.0; 
                          /*
                          for (auto i : r1)
                            for (size_t j = 0; j < dim_proxy2; j++)
                              for (size_t k = 0; k < dim_proxy1; k++)
                                {
                                  auto res = bdbmat1.Row(i*dim_proxy2+j);
                                  auto a = bbmat1.Row(i*dim_proxy1+k);
                                  auto b = proxyvalues.Row(k*dim_proxy2+j);
                                  res += pw_mult(a,b);
                                }
                          */
                          /*
                          for (size_t j = 0; j < dim_proxy2; j++)
                            for (size_t k = 0; k < dim_proxy1; k++)
                              if (nonzeros(l1+j, k1+k))
                                {
                                  auto hproxyvalues = proxyvalues.Row(k*dim_proxy2+j);
                                  auto hbbmat1 = bbmat1.RowSlice(k, dim_proxy1);
                                  auto hbdbmat1 = bdbmat1.RowSlice(j, dim_proxy2);
                                    // for (auto i : r1)
                                    // hbdbmat1.Row(i).AddSize(ir.Size()) += pw_mult(hbbmat1.Row(i), hproxyvalues);
                                  for (size_t i = 0; i < ir.Size(); i++)
                                    hbdbmat1.Col(i).Range(r1) += hproxyvalues(i) * hbbmat1.Col(i).Range(r1);
                                }
                          */
//                          RegionTimer regmult(tmult);
                          for (size_t j = 0; j < dim_proxy2; j++)
                            for (size_t k = 0; k < dim_proxy1; k++)
                              if (nonzeros(l1+j, k1+k))
                                {
                                  auto proxyvalues_jk = symbolic_integrator_uses_diff ?
                                          proxyvalues.Row(j*dim_proxy1+k) : proxyvalues.Row(k*dim_proxy2+j);
                                  auto bbmat1_k = bbmat1.RowSlice(k, dim_proxy1).Rows(r1);
                                  auto bdbmat1_j = bdbmat1.RowSlice(j, dim_proxy2).Rows(r1);

                                  for (size_t i = 0; i < ir.Size(); i++)
                                    bdbmat1_j.Col(i).Range(0,r1.Size()) += proxyvalues_jk(i)*weights(i) * bbmat1_k.Col(i);
                                }
                        }
                      
                      // elmat.Rows(r2).Cols(r1) += bbmat2.Rows(r2) * Trans(bdbmat1.Rows(r1));
                      // AddABt (bbmat2.Rows(r2), bdbmat1.Rows(r1), elmat.Rows(r2).Cols(r1));

                      symmetric_so_far &= samediffop && is_diagonal;
                      /*
                      if (symmetric_so_far)
                        AddABtSym (AFlatMatrix<double>(hbbmat2.Rows(r2)),
                                   AFlatMatrix<double> (hbdbmat1.Rows(r1)), part_elmat);
                      else
                        AddABt (AFlatMatrix<double> (hbbmat2.Rows(r2)),
                                AFlatMatrix<double> (hbdbmat1.Rows(r1)), part_elmat);
                      */

                      {
                        // static Timer t("AddABt", NoTracing);
                        // RegionTracer reg(TaskManager::GetThreadId(), t);
                        
                        if (symmetric_so_far)
                        {
                          /*
                            RegionTimer regdmult(timer_SymbBFImultsym);
                            NgProfiler::AddThreadFlops(timer_SymbBFImultsym, TaskManager::GetThreadId(),
                            SIMD<double>::Size()*2*r2.Size()*(r1.Size()+1)*hbbmat2.Width() / 2);
                          */
                          AddABtSym (hbbmat2.Rows(r2), hbdbmat1.Rows(r1), part_elmat);
                        }
                      else
                        {
                          /*
                          RegionTimer regdmult(timer_SymbBFImult);
                          NgProfiler::AddThreadFlops(timer_SymbBFImult, TaskManager::GetThreadId(),
                                                     SIMD<double>::Size()*2*r2.Size()*r1.Size()*hbbmat2.Width());
                          */
                          AddABt (hbbmat2.Rows(r2), hbdbmat1.Rows(r1), part_elmat);
                        }
                      }
                      if (symmetric_so_far)
                        {
                          ExtendSymmetric (part_elmat);
                          /*
                          size_t h = part_elmat.Height();
                          for (size_t i = 0; i+1 < h; i++)
                            for (size_t j = i+1; j < h; j++)
                              part_elmat(i,j) = part_elmat(j,i);
                          */
                        }
                    }
              
                  l1 += proxy2->Dimension();
                  l1nr++;
                }
              k1 += proxy1->Dimension();
              k1nr++;
            }
            }
          // ir.NothingToDelete();
          return;
        }
      catch (const ExceptionNOSIMD& e)
        {
          cout << IM(6) << e.What() << endl
               << "switching to scalar evaluation" << endl;
          simd_evaluate = false;
          throw ExceptionNOSIMD("in TCalcElementMatrixAdd");
          // T_CalcElementMatrixAdd<SCAL, SCAL_SHAPES, SCAL_RES> (fel, trafo, elmat, lh);
          // return;
        }
    

    // IntegrationRule ir(trafo.GetElementType(), intorder);
    const IntegrationRule& ir = GetIntegrationRule (fel, lh);
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    PrecomputeCacheCF(cache_cfs, mir, lh);
    
    // tstart.Stop();
    // bool symmetric_so_far = true;
    int k1 = 0;
    int k1nr = 0;
    for (auto proxy1 : trial_proxies)
      {
        int l1 = 0;
        int l1nr = 0;
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

            if (is_nonzero) //   && k1nr*test_proxies.Size()+l1nr >= first_std_eval)
              {
                HeapReset hr(lh);
                bool samediffop = (*(proxy1->Evaluator()) == *(proxy2->Evaluator())) && !is_mixedfe;
                // td.Start();
                FlatTensor<3,SCAL> proxyvalues(lh, mir.Size(), proxy1->Dimension(), proxy2->Dimension());
                FlatVector<SCAL> diagproxyvalues(mir.Size()*proxy1->Dimension(), lh);
                FlatMatrix<SCAL> val(mir.Size(), 1, lh);

                IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                SliceMatrix<SCAL_RES> part_elmat = elmat.Rows(r2).Cols(r1);
                FlatMatrix<SCAL_SHAPES, ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                FlatMatrix<SCAL_SHAPES, ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

                if (ddcf_dtest_dtrial(l1nr, k1nr))
                  {
//                    cout << "use ddcf_dtest_dtrial (NO SIMD)" << endl;
                    // TODO: optimize for element-wise constant case?
                    FlatMatrix<SCAL> mproxyvalues(mir.Size(), proxy1->Dimension() * proxy2->Dimension(),
                                                  proxyvalues.Data());
                    ddcf_dtest_dtrial(l1nr, k1nr)->Evaluate(mir, mproxyvalues);
                    if (is_diagonal)
                        for (auto k: Range(proxy1->Dimension()))
                            diagproxyvalues.Slice(k, proxy1->Dimension()) = proxyvalues(STAR, k, k);
                  }
                else
                  {
                    if (!is_diagonal)
                      for (int k = 0; k < proxy1->Dimension(); k++)
                        for (int l = 0; l < proxy2->Dimension(); l++)
                          {
                            if (nonzeros(l1+l, k1+k))
                              {
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

                constexpr size_t BS = 16;
                for (size_t i = 0; i < mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    int bs = min2(size_t(BS), mir.Size()-i);
                    
                    FlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);
                    FlatMatrix<SCAL> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                    FlatMatrix<SCAL_SHAPES> bbmat2 = samediffop ?
                      bbmat1 : FlatMatrix<SCAL_SHAPES>(elmat.Height(), bs*proxy2->Dimension(), lh);

                    // tb.Start();
                    BaseMappedIntegrationRule & bmir = mir.Range(i, i+bs, lh);

                    proxy1->Evaluator()->CalcMatrix(fel_trial, bmir, Trans(bbmat1), lh);

                    if (!samediffop)
                      proxy2->Evaluator()->CalcMatrix(fel_test, bmir, Trans(bbmat2), lh);
                    // tb.Stop();

                    // tdb.Start();
                    if (is_diagonal)
                      {
                        FlatVector<SCAL> diagd(bs*proxy1->Dimension(), lh);
                        diagd = diagproxyvalues.Range(i*proxy1->Dimension(),
                                                      (i+bs)*proxy1->Dimension());
                        for (size_t i = 0; i < diagd.Size(); i++)
                          bdbmat1.Col(i) = diagd(i) * bbmat1.Col(i);
                        // MultMatDiagMat(bbmat1, diagd, bdbmat1);
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
            l1nr++;
          }
        k1 += proxy1->Dimension();
        k1nr++;
      }
  }


  void 
  SymbolicBilinearFormIntegrator ::
  CalcElementMatrix (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    elmat = 0.0;
    bool symmetric_so_far = true;
    try
      {
        T_CalcElementMatrixAdd<double,double,double> (fel, trafo, elmat, symmetric_so_far, lh);
      }
    catch (ExceptionNOSIMD & e)
      {
        elmat = 0.0;        
        T_CalcElementMatrixAdd<double,double,double> (fel, trafo, elmat, symmetric_so_far, lh);        
      }
  }
  
  void 
  SymbolicBilinearFormIntegrator ::
  CalcElementMatrix (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatMatrix<Complex> elmat,
                     LocalHeap & lh) const
  {
    try
      {
        elmat = 0.0;
        bool symmetric_so_far = true;
        if (fel.ComplexShapes() || trafo.IsComplex())
          T_CalcElementMatrixAdd<Complex,Complex,Complex> (fel, trafo, elmat, symmetric_so_far, lh);
        else
          {
            if (cf->IsComplex())
              T_CalcElementMatrixAdd<Complex,double,Complex> (fel, trafo, elmat, symmetric_so_far, lh);
            else
              T_CalcElementMatrixAdd<double,double,Complex> (fel, trafo, elmat, symmetric_so_far, lh);
          }
      }
    catch (const ExceptionNOSIMD& e)  // retry with simd_evaluate is off
      {
        elmat = 0.0;
        bool symmetric_so_far = true;        
        if (fel.ComplexShapes() || trafo.IsComplex())
          T_CalcElementMatrixAdd<Complex,Complex,Complex> (fel, trafo, elmat, symmetric_so_far, lh);
        else
          {
            if (cf->IsComplex())
              T_CalcElementMatrixAdd<Complex,double,Complex> (fel, trafo, elmat, symmetric_so_far, lh);
            else
              T_CalcElementMatrixAdd<double,double,Complex> (fel, trafo, elmat, symmetric_so_far, lh);
          }
      }    
  }

  void 
  SymbolicBilinearFormIntegrator ::
  CalcElementMatrixAdd (const FiniteElement & fel,
                        const ElementTransformation & trafo, 
                        FlatMatrix<double> elmat,
                        bool & symmetric_so_far,
                        LocalHeap & lh) const
  {
    T_CalcElementMatrixAdd<double,double,double> (fel, trafo, elmat, symmetric_so_far, lh);
  }
  
  void 
  SymbolicBilinearFormIntegrator ::
  CalcElementMatrixAdd (const FiniteElement & fel,
                        const ElementTransformation & trafo, 
                        FlatMatrix<Complex> elmat,
                        bool & symmetric_so_far,
                        LocalHeap & lh) const
  {
    if (fel.ComplexShapes() || trafo.IsComplex())
      T_CalcElementMatrixAdd<Complex,Complex,Complex> (fel, trafo, elmat, symmetric_so_far, lh);
    else
      if (cf->IsComplex())
        T_CalcElementMatrixAdd<Complex,double,Complex> (fel, trafo, elmat, symmetric_so_far, lh);
      else
        T_CalcElementMatrixAdd<double,double,Complex> (fel, trafo, elmat, symmetric_so_far, lh);
  }


  

  template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcElementMatrixEBAdd (const FiniteElement & fel,
                            const ElementTransformation & trafo, 
                            FlatMatrix<SCAL_RES> elmat,
                            LocalHeap & lh) const
      
    {
      static Timer t("symbolicBFI - CalcElementMatrix EB", NoTracing);
      /*
      static Timer tir("symbolicBFI - CalcElementMatrix EB - intrules", NoTracing);
      static Timer td("symbolicBFI - CalcElementMatrix EB - dmats", NoTracing);
      static Timer tdb("symbolicBFI - CalcElementMatrix EB - b*d", NoTracing);
      static Timer tb("symbolicBFI - CalcElementMatrix EB - bmats", NoTracing);
      static Timer tmult("symbolicBFI - CalcElementMatrix EB - mult", NoTracing);
      */
      // RegionTimer reg(t);

      // elmat = 0;

      const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
      const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
      const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;
      
      auto eltype = trafo.GetElementType();

      Facet2ElementTrafo transform(eltype, element_vb); 
      int nfacet = transform.GetNFacets();

      if (simd_evaluate)
        // if (false)  // throwing the no-simd exception after some terms already added is still a problem 
        {
          try
            {
              for (int k = 0; k < nfacet; k++)
                {
                  HeapReset hr(lh);
                  ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
                  SIMD_IntegrationRule ir_facet1(etfacet, fel_trial.Order()+fel_test.Order()+bonus_intorder);
                  SIMD_IntegrationRule & ir_facet(userdefined_simd_intrules[etfacet] ?
                                                  *userdefined_simd_intrules[etfacet]
                                                  : ir_facet1);
                  auto & ir_facet_vol = transform(k, ir_facet, lh);
                  
                  auto & mir = trafo(ir_facet_vol, lh);
          
                  ProxyUserData ud;
                  const_cast<ElementTransformation&>(trafo).userdata = &ud;

                  PrecomputeCacheCF(cache_cfs, mir, lh);

                  // mir.ComputeNormalsAndMeasure(eltype, k);
                  
                  for (int k1 : Range(trial_proxies))
                    for (int l1 : Range(test_proxies))
                      {
                        if (!nonzeros_proxies(l1, k1)) continue;
                        
                        auto proxy1 = trial_proxies[k1];
                        auto proxy2 = test_proxies[l1];
                        size_t dim_proxy1 = proxy1->Dimension();
                        size_t dim_proxy2 = proxy2->Dimension();
                        HeapReset hr(lh);
                        FlatMatrix<SIMD<SCAL>> proxyvalues(dim_proxy1*dim_proxy2, ir_facet.Size(), lh);
                
                        // td.Start();
                        for (int k = 0; k < dim_proxy1; k++)
                          for (int l = 0; l < dim_proxy2; l++)
                            {
                              ud.trialfunction = proxy1;
                              ud.trial_comp = k;
                              ud.testfunction = proxy2;
                              ud.test_comp = l;

                              auto kk = l + k*dim_proxy2;
                              cf->Evaluate (mir, proxyvalues.Rows(kk, kk+1));
                              for (size_t i = 0; i < mir.Size(); i++)
                                proxyvalues(kk, i) *= mir[i].GetWeight(); 
                            }

                        
                        IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                        IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                        SliceMatrix<SCAL_RES> part_elmat = elmat.Rows(r2).Cols(r1);
                        
                        FlatMatrix<SIMD<SCAL_SHAPES>> bbmat1(elmat.Width()*dim_proxy1, mir.Size(), lh);
                        bool samediffop = false; // not yet available
                        FlatMatrix<SIMD<SCAL_SHAPES>> bbmat2 = samediffop ?
                          bbmat1 : FlatMatrix<SIMD<SCAL_SHAPES>>(elmat.Height()*dim_proxy2, mir.Size(), lh);

                        {
                          // RegionTimer regbmat(timer_SymbBFIbmat);
                          proxy1->Evaluator()->CalcMatrix(fel_trial, mir, bbmat1);
                          if (!samediffop)
                            proxy2->Evaluator()->CalcMatrix(fel_test, mir, bbmat2);
                        }

                        


                        if (dim_proxy2 < dim_proxy1)
                          {
                            FlatMatrix<SIMD<SCAL>> bdbmat1(elmat.Width()*dim_proxy2, mir.Size(), lh);
                            FlatMatrix<SIMD<SCAL>> hbdbmat1(elmat.Width(), dim_proxy2*mir.Size(),
                                                            &bdbmat1(0,0));
                            
                            {
                              // static Timer t("SymbolicBFI::EB - DB  V1", NoTracing);
                              // RegionTracer regtr(TaskManager::GetThreadId(), t);    
                              hbdbmat1.Rows(r1) = 0.0;
                              for (size_t j = 0; j < dim_proxy2; j++)
                                for (size_t k = 0; k < dim_proxy1; k++)
                                  // if (nonzeros(l1+j, k1+k))
                                  {
                                    auto proxyvalues_jk = proxyvalues.Row(k*dim_proxy2+j);
                                    auto bbmat1_k = bbmat1.RowSlice(k, dim_proxy1).Rows(r1);
                                    auto bdbmat1_j = bdbmat1.RowSlice(j, dim_proxy2).Rows(r1);
                                    
                                    for (size_t i = 0; i < mir.Size(); i++)
                                      bdbmat1_j.Col(i).Range(0,r1.Size()) += proxyvalues_jk(i)*bbmat1_k.Col(i);
                                  }
                            }
                            
                            // static Timer t("SymbolicBFI::EB - AddABt V1", NoTracing);
                            // RegionTracer regtr(TaskManager::GetThreadId(), t);    
                            
                            FlatMatrix<SIMD<SCAL_SHAPES>> hbbmat2(elmat.Height(), dim_proxy2*mir.Size(),
                                                                  &bbmat2(0,0));
                            AddABt (hbbmat2.Rows(r2), hbdbmat1.Rows(r1), part_elmat);
                          }
                        else
                          {
                            FlatMatrix<SIMD<SCAL>> bdbmat2(elmat.Height()*dim_proxy1, mir.Size(), lh);
                            FlatMatrix<SIMD<SCAL>> hbdbmat2(elmat.Height(), dim_proxy1*mir.Size(),
                                                            &bdbmat2(0,0));
                            
                            {
                              // static Timer t("SymbolicBFI::EB - DB V2", NoTracing);
                              // RegionTracer regtr(TaskManager::GetThreadId(), t);    
                              hbdbmat2.Rows(r2) = 0.0;
                              for (size_t j = 0; j < dim_proxy2; j++)
                                for (size_t k = 0; k < dim_proxy1; k++)
                                  // if (nonzeros(l1+j, k1+k))
                                  {
                                    auto proxyvalues_jk = proxyvalues.Row(k*dim_proxy2+j);
                                    auto bbmat2_j = bbmat2.RowSlice(j, dim_proxy2).Rows(r2);
                                    auto bdbmat2_k = bdbmat2.RowSlice(k, dim_proxy1).Rows(r2);
                                    
                                    for (size_t i = 0; i < mir.Size(); i++)
                                      bdbmat2_k.Col(i).Range(0,r2.Size()) += proxyvalues_jk(i)*bbmat2_j.Col(i);
                                  }
                            }
                            
                            // static Timer t("SymbolicBFI::EB - AddABt V2", NoTracing);
                            // RegionTracer regtr(TaskManager::GetThreadId(), t);    
                            
                            FlatMatrix<SIMD<SCAL_SHAPES>> hbbmat1(elmat.Width(), dim_proxy1*mir.Size(),
                                                                  &bbmat1(0,0));
                            AddABt (hbdbmat2.Rows(r2), hbbmat1.Rows(r1), part_elmat);
                          }
                      }
                }
              return;
            }
          
          catch (const ExceptionNOSIMD& e)
            {
              cout << IM(6) << e.What() << endl
                   << "switching to scalar evaluation, may be a problem with Add" << endl;
              simd_evaluate = false;
              throw ExceptionNOSIMD("disabled simd-evaluate in AddElementMatrixEB");
            }
        }
      
      for (int k = 0; k < nfacet; k++)
        {
          // tir.Start();
          HeapReset hr(lh);
          ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
          const IntegrationRule& ir_facet = GetIntegrationRule(etfacet, fel_trial.Order()+fel_test.Order()+bonus_intorder);
          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
          
          BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
          
          ProxyUserData ud;
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          // mir.ComputeNormalsAndMeasure(eltype, k);

          PrecomputeCacheCF(cache_cfs, mir, lh);

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
                SliceMatrix<SCAL_RES> part_elmat = elmat.Rows(r2).Cols(r1);

                constexpr size_t BS = 16;
                for (size_t i = 0; i < mir.Size(); i+=BS)
                  {
                    HeapReset hr(lh);
                    int bs = min2(BS, mir.Size()-i);
                    
                    FlatMatrix<SCAL_SHAPES> bbmat1(elmat.Width(), bs*proxy1->Dimension(), lh);
                    FlatMatrix<SCAL> bdbmat1(elmat.Width(), bs*proxy2->Dimension(), lh);
                    FlatMatrix<SCAL_SHAPES> bbmat2(elmat.Height(), bs*proxy2->Dimension(), lh);

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


  
  template <typename SCAL, typename SCAL_SHAPES, typename SCAL_RES>
  void SymbolicBilinearFormIntegrator ::  
  T_CalcElementMatrixAddShapeWise (const FiniteElement & fel,
                                   const ElementTransformation & trafo, 
                                   FlatMatrix<SCAL_RES> elmat,
                                   LocalHeap & lh) const
  {
    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;

    HeapReset hr(lh);
    
    const IntegrationRule& ir = GetIntegrationRule (fel, lh);
    BaseMappedIntegrationRule & mir = trafo(ir, lh);
    
    ProxyUserData ud(trial_proxies.Size()+test_proxies.Size(), 0, lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;

    PrecomputeCacheCF(cache_cfs, mir, lh);

    ud.fel = &fel;

    FlatVector<> vtrial(fel_trial.GetNDof(), lh);
    FlatVector<> vtest(fel_test.GetNDof(), lh);
    ud.trial_elvec = &vtrial;
    ud.test_elvec = &vtest;
    ud.lh = &lh;
    
    FlatMatrix<SCAL> val(mir.Size(), 1, lh);
    
    for (int k = 0; k < fel_trial.GetNDof(); k++)
      for (int l = 0; l < fel_test.GetNDof(); l++)
        {
          vtrial = 0.0;
          vtest = 0.0;
          vtrial(k) = 1;
          vtest(l) = 1;
          
          cf -> Evaluate (mir, val);
          
          SCAL sum = 0.0;
          for (int i = 0; i < mir.Size(); i++)
            sum += val(i) *mir[i].GetWeight();
          elmat(l,k) += sum;                
        }
  }

  void
  SymbolicBilinearFormIntegrator ::
  CalcLinearizedElementMatrix (const FiniteElement & fel,
                               const ElementTransformation & trafo,
                               FlatVector<Complex> elveclin,
                               FlatMatrix<Complex> elmat,
                               LocalHeap & lh) const
  {
    if(linearization)
      {
        T_CalcLinearizedElementMatrixFrozen(fel, trafo, elveclin, elmat, lh);
        return;
      }
    CalcElementMatrix(fel, trafo, elmat, lh);
  }
  
  void 
  SymbolicBilinearFormIntegrator ::
  CalcLinearizedElementMatrix (const FiniteElement & fel,
                               const ElementTransformation & trafo, 
                               FlatVector<double> elveclin,
                               FlatMatrix<double> elmat,
                               LocalHeap & lh) const
  {
    if(linearization)
      {
        T_CalcLinearizedElementMatrixFrozen(fel, trafo, elveclin, elmat, lh);
        return;
      }
    auto save_userdata = trafo.PushUserData();
    
    
    /*
      CalcElementMatrix(fel, trafo, elmat, lh);
      return;
    */

      
    if (element_vb != VOL)
      {
        T_CalcLinearizedElementMatrixEB<double,double> (fel, trafo, elveclin, elmat, lh);
        return;
      }

    
    static Timer t("symbolicbfi - calclinearized", NoTracing);
    RegionTimer reg(t);
    
    // static Timer td("symbolicbfi - calclinearized dmats", NoTracing);
    // RegionTimer reg(t);

    if (simd_evaluate)
      // if (false)
      try
        {
          const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
          const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
          const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;

          const SIMD_IntegrationRule& ir = Get_SIMD_IntegrationRule (fel, lh);
          SIMD_BaseMappedIntegrationRule & mir = trafo(ir, lh);

          ProxyUserData ud(trial_proxies.Size(), gridfunction_cfs.Size(), lh);
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          ud.fel = &fel;
          // ud.elx = &elveclin;
          // ud.lh = &lh;

          for (ProxyFunction * proxy : trial_proxies)
            {
              ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);
              proxy->Evaluator()->Apply(fel_trial, mir, elveclin, ud.GetAMemory(proxy));
            }
          for (CoefficientFunction * cf : gridfunction_cfs)
            ud.AssignMemory (cf, ir.GetNIP(), cf->Dimension(), lh);
    
          // AFlatMatrix<> val(1, mir.IR().GetNIP(), lh);
          FlatMatrix<AutoDiff<1,SIMD<double>>> val(1, mir.Size(), lh);
          elmat = 0;

          IntRange unified_r1(0, 0);
          for (int l1 : Range(test_proxies))
            {
              HeapReset hr(lh);              
              auto proxy2 = test_proxies[l1];
              FlatMatrix<SIMD<double>> bdbmat1(elmat.Width()*proxy2->Dimension(), ir.Size(), lh);
              FlatMatrix<SIMD<double>> hbdbmat1(elmat.Width(), proxy2->Dimension()*ir.Size(),
                                                &bdbmat1(0,0));
              bdbmat1 = 0.0;
              
              for (int k1 : Range(trial_proxies))
                {
                  HeapReset hr(lh);
                  auto proxy1 = trial_proxies[k1];
                  
                  FlatMatrix<SIMD<double>> proxyvalues(proxy1->Dimension()*proxy2->Dimension(), ir.Size(), lh);
                  
                  for (size_t k = 0, kk = 0; k < proxy1->Dimension(); k++)
                    for (size_t l = 0; l < proxy2->Dimension(); l++, kk++)
                      {
                        ud.trialfunction = proxy1;
                        ud.trial_comp = k;
                        ud.testfunction = proxy2;
                        ud.test_comp = l;
                        // cf -> EvaluateDeriv (mir, val, proxyvalues.Rows(kk,kk+1));
                        cf -> Evaluate (mir, val);
                        auto row = proxyvalues.Row(kk);
                        for (auto j : Range(mir.Size()))
                          row(j) = val(j).DValue(0);
                        /*
                        if (HSum(L2Norm2(row)) > 0 && !nonzeros_deriv(test_cum[l1]+l, trial_cum[k1]+k))
                          {
                            cout << "wrong nonzero deriv" << endl;
                            cout << "nonzeros_deriv = " << endl << nonzeros_deriv << endl;
                          }
                        */
                      }
                  
                  for (size_t i = 0; i < mir.Size(); i++)
                    proxyvalues.Col(i) *= mir[i].GetWeight();

                  IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                  
                  FlatMatrix<SIMD<double>> bbmat1(elmat.Width()*proxy1->Dimension(), ir.Size(), lh);
                  
                  // bbmat1 = 0.0;
                  proxy1->Evaluator()->CalcMatrix(fel_trial, mir, bbmat1);
                  for (auto i : r1)
                    for (size_t j = 0; j < proxy2->Dimension(); j++)
                      for (size_t k = 0; k < proxy1->Dimension(); k++)
                        {
                          bdbmat1.Row(i*proxy2->Dimension()+j) +=
                            pw_mult (bbmat1.Row(i*proxy1->Dimension()+k),
                                     proxyvalues.Row(k*proxy2->Dimension()+j));
                        }

                  if (unified_r1.Next() == 0)
                    unified_r1 = r1;
                  else
                    unified_r1 = IntRange (min2(r1.First(), unified_r1.First()),
                                           max2(r1.Next(), unified_r1.Next()));
                }

              IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
              
              FlatMatrix<SIMD<double>> bbmat2(elmat.Height()*proxy2->Dimension(), ir.Size(), lh);
              FlatMatrix<SIMD<double>> hbbmat2(elmat.Height(), proxy2->Dimension()*ir.Size(),
                                              &bbmat2(0,0));
              // bbmat2 = 0.0;
              proxy2->Evaluator()->CalcMatrix(fel_test, mir, bbmat2);
              
              AddABt (hbbmat2.Rows(r2), hbdbmat1.Rows(unified_r1), elmat.Rows(r2).Cols(unified_r1));
              // AddABt (hbbmat2.Rows(r2), hbdbmat1, elmat.Rows(r2));
            }
          /*
          Matrix<> helmat = elmat;
          simd_evaluate = false;
          CalcLinearizedElementMatrix (fel, trafo, elveclin, elmat, lh);
          simd_evaluate = true;          
          double err = L2Norm(helmat-elmat);
          if (err > 1e-5) cout << "err = " << err << endl;
          */
          return;
        }
      catch (const ExceptionNOSIMD& e)
        {
          cout << IM(6) << e.What() << endl
               << "switching to scalar evaluation in CalcLinearized" << endl;
          simd_evaluate = false;
          CalcLinearizedElementMatrix (fel, trafo, elveclin, elmat, lh);
          return;
        }
    
    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;

    const IntegrationRule& ir = GetIntegrationRule (fel, lh);
    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    // ud.elx = &elveclin;
    // ud.lh = &lh;

    for (ProxyFunction * proxy : trial_proxies)
      {
        ud.AssignMemory (proxy, ir.Size(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel_trial, mir, elveclin, ud.GetMemory(proxy), lh);
      }
    
    FlatMatrix<> val(mir.Size(), 1, lh), deriv(mir.Size(), 1, lh);
    FlatMatrix<AutoDiff<1>> dval(mir.Size(), 1, lh);
    elmat = 0;
    
    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              // if (nonzeros(test_cum[l1]+l, trial_cum[k1]+k)) // does no work for non-linear 
              if (true)
                {
                  ud.trialfunction = proxy1;
                  ud.trial_comp = k;
                  ud.testfunction = proxy2;
                  ud.test_comp = l;
                  
                  // cf -> EvaluateDeriv (mir, val, deriv);
                  // proxyvalues(STAR,l,k) = deriv.Col(0);
                  cf -> Evaluate (mir, dval);
                  for (size_t i = 0; i < mir.Size(); i++)
                    proxyvalues(i,l,k) = dval(i,0).DValue(0);
                }
              else
                proxyvalues(STAR,l,k) = 0;

          for (int i = 0; i < mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();

          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

          constexpr size_t BS = 16;
          for (size_t i = 0; i < mir.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<double,ColMajor> bbmat1(rest*proxy1->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);

              proxy1->Evaluator()->CalcLinearizedMatrix(fel_trial, mir.Range(i,i+rest,lh), elveclin, bbmat1, lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r1 = proxy1->Dimension() * IntRange(j,j+1);                  
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  // proxy1->Evaluator()->CalcMatrix(fel_trial, mir[ii], bmat1, lh);
                  bmat1 = bbmat1.Rows(r1);
                  proxy2->Evaluator()->CalcMatrix(fel_test, mir[ii], bmat2, lh);
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
              elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
            }
        }
  }



  template <typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_CalcLinearizedElementMatrixEB (const FiniteElement & fel1,
                                   const ElementTransformation & trafo, 
                                   FlatVector<double> elveclin,
                                   FlatMatrix<double> elmat,
                                   LocalHeap & lh) const
  {
    // size_t tid = TaskManager::GetThreadId();    
    static Timer t("symbolicbfi - calclinearized EB", NoTracing);
    static Timer tnosimd("symbolicbfi - calclinearized EB nosimd", NoTracing);
    static Timer td("symbolicbfi - calclinearized EB dmats", NoTracing);
    static Timer tir("symbolicbfi - calclinearized EB intrule", NoTracing);
    RegionTimer reg(t);
    
    elmat = 0;
    
    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel1);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel1;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel1;

    auto eltype = trafo.GetElementType();
    
    Facet2ElementTrafo transform(eltype, element_vb); 
    int nfacet = transform.GetNFacets();
    
    if (simd_evaluate)
      try
        {
          for (int k = 0; k < nfacet; k++)
            {
              HeapReset hr(lh);
              ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
              // NgProfiler::StartThreadTimer(tir, tid);
              const SIMD_IntegrationRule& ir_facet = GetSIMDIntegrationRule(etfacet, fel_trial.Order()+fel_test.Order()+bonus_intorder);
              const SIMD_IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
              SIMD_BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
              // mir.ComputeNormalsAndMeasure(eltype, k);
              // NgProfiler::StopThreadTimer(tir, tid);
              ProxyUserData ud(trial_proxies.Size(), gridfunction_cfs.Size(), lh);
              const_cast<ElementTransformation&>(trafo).userdata = &ud;
              ud.fel = &fel1;

              for (ProxyFunction * proxy : trial_proxies)
                {
                  ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
                  proxy->Evaluator()->Apply(fel_trial, mir, elveclin, ud.GetAMemory(proxy));
                }
              for (CoefficientFunction * cf : gridfunction_cfs)
                ud.AssignMemory (cf, ir_facet.GetNIP(), cf->Dimension(), lh);

              FlatMatrix<AutoDiff<1,SIMD<double>>> val(1, mir.Size(), lh);
              
              for (int l1 : Range(test_proxies))
                {
                  HeapReset hr(lh);              
                  auto proxy2 = test_proxies[l1];
                  FlatMatrix<SIMD<double>> bdbmat1(elmat.Width()*proxy2->Dimension(), ir_facet.Size(), lh);
                  FlatMatrix<SIMD<double>> hbdbmat1(elmat.Width(), proxy2->Dimension()*ir_facet.Size(),
                                                    &bdbmat1(0,0));
                  bdbmat1 = 0.0;
                  
                  for (int k1 : Range(trial_proxies))
                    {
                      HeapReset hr(lh);
                      auto proxy1 = trial_proxies[k1];
                      
                      FlatMatrix<SIMD<double>> proxyvalues(proxy1->Dimension()*proxy2->Dimension(), ir_facet.Size(), lh);
                      
                      for (size_t k = 0, kk = 0; k < proxy1->Dimension(); k++)
                        for (size_t l = 0; l < proxy2->Dimension(); l++, kk++)
                          {
                            // RegionTimer reg(td);                            
                            ud.trialfunction = proxy1;
                            ud.trial_comp = k;
                            ud.testfunction = proxy2;
                            ud.test_comp = l;
                            
                            cf -> Evaluate (mir, val);
                            auto row = proxyvalues.Row(kk);
                            for (auto j : Range(mir.Size()))
                              row(j) = val(j).DValue(0);
                          }

                      for (size_t i = 0; i < mir.Size(); i++)
                        proxyvalues.Col(i) *= mir[i].GetWeight();                        
                      
                      IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);

                      FlatMatrix<SIMD<double>> bbmat1(elmat.Width()*proxy1->Dimension(), ir_facet.Size(), lh);
                      proxy1->Evaluator()->CalcMatrix(fel_trial, mir, bbmat1);

                      for (auto i : r1)
                        for (size_t j = 0; j < proxy2->Dimension(); j++)
                          for (size_t k = 0; k < proxy1->Dimension(); k++)
                            {
                              bdbmat1.Row(i*proxy2->Dimension()+j) +=
                                pw_mult (bbmat1.Row(i*proxy1->Dimension()+k),
                                         proxyvalues.Row(k*proxy2->Dimension()+j));
                            }
                    }
                  
                  IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                  FlatMatrix<SIMD<double>> bbmat2(elmat.Height()*proxy2->Dimension(), ir_facet.Size(), lh);
                  FlatMatrix<SIMD<double>> hbbmat2(elmat.Height(), proxy2->Dimension()*ir_facet.Size(),
                                                   &bbmat2(0,0));

                  proxy2->Evaluator()->CalcMatrix(fel_test, mir, bbmat2);
                  AddABt (hbbmat2.Rows(r2), hbdbmat1, elmat.Rows(r2));
                }
            }

          return;
        }
      catch (const ExceptionNOSIMD& e)
        {
          cout << IM(6) << e.What() << endl
               << "switching to scalar evaluation in CalcLinearizedEB" << endl;
          simd_evaluate = false;
          T_CalcLinearizedElementMatrixEB<SCAL,SCAL_SHAPES> (fel1, trafo, elveclin, elmat, lh);
          return;
        }


    RegionTimer regnosimd(tnosimd);

      
    for (int k = 0; k < nfacet; k++)
      {
          HeapReset hr(lh);
          ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
        
          const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, fel_trial.Order()+fel_test.Order()+bonus_intorder);
          IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
          BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
          //  mir.ComputeNormalsAndMeasure(eltype, k);
          
          ProxyUserData ud(trial_proxies.Size(), lh);          
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          ud.fel = &fel1;

          for (ProxyFunction * proxy : trial_proxies)
            {
              ud.AssignMemory (proxy, mir.Size(), proxy->Dimension(), lh);
              proxy->Evaluator()->Apply(fel_trial, mir, elveclin, ud.GetMemory(proxy), lh);
            }
    
          FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh);
          FlatMatrix<AutoDiff<1>> dval(mir.Size(), 1,lh);
          
          for (int k1 : Range(trial_proxies))
            for (int l1 : Range(test_proxies))
              {
                HeapReset hr(lh);
                auto proxy1 = trial_proxies[k1];
                auto proxy2 = test_proxies[l1];
                // td.Start(tid);
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
                        
                        cf -> Evaluate (mir, dval);
                        for (size_t i = 0; i < mir.Size(); i++)
                          proxyvalues(i,l,k) = dval(i,0).DValue(0);
                      }
                    else
                      proxyvalues(STAR,l,k) = 0.0;
                        
                // td.Stop(tid);

                for (int i = 0; i < mir.Size(); i++)
                  proxyvalues(i,STAR,STAR) *= ir_facet[i].Weight() * mir[i].GetMeasure(); 
                
                t.AddFlops (double (mir.Size()) * proxy1->Dimension()*elmat.Width()*elmat.Height());
                
                FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
                FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);
                
                constexpr size_t BS = 16;
                for (size_t i = 0; i < mir.Size(); i+=BS)
                  {
                    int rest = min2(BS, mir.Size()-i);
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
                    
                    IntRange r1 = proxy1->Evaluator()->UsedDofs(fel_trial);
                    IntRange r2 = proxy2->Evaluator()->UsedDofs(fel_test);
                    // elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1) | Lapack;
                    AddABt (Trans(bbmat2).Rows(r2), Trans(bdbmat1).Rows(r1),
                            SliceMatrix<> (elmat.Rows(r2).Cols(r1)));
                  }
              }
        }
  }

  template<typename SCAL>
  void SymbolicBilinearFormIntegrator ::
  T_CalcLinearizedElementMatrixFrozen (const FiniteElement & fel,
                                       const ElementTransformation & trafo, 
                                       FlatVector<SCAL> elveclin,
                                       FlatMatrix<SCAL> elmat,
                                       LocalHeap & lh) const
    {
      linearization->CalcElementMatrix(fel, trafo, elmat, lh);
    }
  
  void
  SymbolicBilinearFormIntegrator :: ApplyElementMatrix (const FiniteElement & fel, 
                                                        const ElementTransformation & trafo, 
                                                        const FlatVector<double> elx, 
                                                        FlatVector<double> ely,
                                                        void * precomputed,
                                                        LocalHeap & lh) const
  {
    auto save_userdata = trafo.PushUserData();
    
    if (element_vb != VOL)
      {
        T_ApplyElementMatrixEB<double,double> (fel, trafo, elx, ely, precomputed, lh);
        return;
      }

    if (simd_evaluate)
      try
        {
          static Timer tall("SymbolicBFI::Apply - all SIMD");
          static Timer tpre("SymbolicBFI::Apply - precomput SIMD");
          static Timer teval("SymbolicBFI::Apply - evaluate SIMD");
          static Timer taddBt("SymbolicBFI::Apply - addBt SIMD");          
          // RegionTimer rall(tall);
 
          bool is_mixed = typeid(fel) == typeid(const MixedFiniteElement&);
          const MixedFiniteElement * mixedfe = static_cast<const MixedFiniteElement*> (&fel);    
          const FiniteElement & fel_trial = is_mixed ? mixedfe->FETrial() : fel;
          const FiniteElement & fel_test = is_mixed ? mixedfe->FETest() : fel;

          HeapReset hr(lh);

          const SIMD_IntegrationRule& simd_ir = Get_SIMD_IntegrationRule (fel, lh);
          auto & simd_mir = trafo(simd_ir, lh);
          
          ProxyUserData ud(trial_proxies.Size(), gridfunction_cfs.Size(), lh);
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          ud.fel = &fel;

          {
            // RegionTimer rpre(tpre);          
          PrecomputeCacheCF(cache_cfs, simd_mir, lh);

          for (ProxyFunction * proxy : trial_proxies)
            ud.AssignMemory (proxy, simd_ir.GetNIP(), proxy->Dimension(), lh);
          for (CoefficientFunction * cf : gridfunction_cfs)
            ud.AssignMemory (cf, simd_ir.GetNIP(), cf->Dimension(), lh);
          
          for (ProxyFunction * proxy : trial_proxies)
            proxy->Evaluator()->Apply(fel_trial, simd_mir, elx, ud.GetAMemory(proxy)); 
          }
          
          ely = 0;
          // for (auto proxy : test_proxies)
          for (auto i : Range(test_proxies))
            {
              auto proxy = test_proxies[i];
              
              HeapReset hr(lh);

              FlatMatrix<SIMD<double>> simd_proxyvalues(proxy->Dimension(), simd_ir.Size(), lh);
              {
                // RegionTimer reval(teval);

              if (dcf_dtest[i])
                dcf_dtest[i]->Evaluate (simd_mir, simd_proxyvalues);
              else
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    ud.testfunction = proxy;
                    ud.test_comp = k;
                    cf -> Evaluate (simd_mir, simd_proxyvalues.Rows(k,k+1));
                  }
              }
              
              for (size_t i = 0; i < simd_proxyvalues.Height(); i++)
                {
                  auto row = simd_proxyvalues.Row(i);
                  for (size_t j = 0; j < row.Size(); j++)
                    row(j) *= simd_mir[j].GetWeight(); //  * simd_ir[j].Weight();
                }

              // RegionTimer rBt(taddBt);          
              proxy->Evaluator()->AddTrans(fel_test, simd_mir, simd_proxyvalues, ely); 
            }
          return;
        }
    
      catch (const ExceptionNOSIMD& e)
        {
          cout << IM(6) << e.What() << endl
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

    IntegrationRule ir = GetIntegrationRule (fel, lh);

    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    PrecomputeCacheCF(cache_cfs, mir, lh);

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


 
  
  
  template <typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_ApplyElementMatrixEB (const FiniteElement & fel, 
                          const ElementTransformation & trafo, 
                          const FlatVector<double> elx, 
                          FlatVector<double> ely,
                          void * precomputed,
                          LocalHeap & lh) const
  {
    // static Timer t("symbolicbfi - Apply EB", NoTracing);
    // static Timer tir("symbolicbfi - Apply EB, intrule", NoTracing);
    // static Timer teval("symbolicbfi - Apply EB, evaluate", NoTracing);
    // static Timer td("symbolicbfi - Apply EB, evaluate D", NoTracing);
    // static Timer ttrans("symbolicbfi - Apply EB, trans", NoTracing);
    
    // RegionTimer reg(t);
    
    ely = 0;

    auto eltype = trafo.GetElementType();

    Facet2ElementTrafo transform(eltype, element_vb); 
    int nfacet = transform.GetNFacets();

    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;
    
    if (simd_evaluate)
      try
        {
          for (int k = 0; k < nfacet; k++)
            {
              HeapReset hr(lh);
              // NgProfiler::StartThreadTimer(tir, tid);                              
              ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
              const SIMD_IntegrationRule & ir_facet =
                GetSIMDIntegrationRule(etfacet,fel_trial.Order()+fel_test.Order()+bonus_intorder);
              auto & ir_facet_vol = transform(k, ir_facet, lh);
              auto & mir = trafo(ir_facet_vol, lh);
              // mir.ComputeNormalsAndMeasure (eltype, k);
              // NgProfiler::StopThreadTimer(tir, tid);

              // NgProfiler::StartThreadTimer(teval, tid);                                            
              ProxyUserData ud(trial_proxies.Size(), gridfunction_cfs.Size(), lh);              
              const_cast<ElementTransformation&>(trafo).userdata = &ud;
              ud.fel = &fel_trial;

              PrecomputeCacheCF(cache_cfs, mir, lh);
          
              for (ProxyFunction * proxy : trial_proxies)
                ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
              for (CoefficientFunction * cf : gridfunction_cfs)
                ud.AssignMemory (cf, ir_facet.GetNIP(), cf->Dimension(), lh);
          
              for (ProxyFunction * proxy : trial_proxies)
                proxy->Evaluator()->Apply(fel_trial, mir, elx, ud.GetAMemory(proxy)); 
          
              // NgProfiler::StopThreadTimer(teval, tid);
              // for (auto proxy : test_proxies)
              for (auto i : Range(test_proxies))
                {
                  auto proxy = test_proxies[i];
                  HeapReset hr(lh);
                  // NgProfiler::StartThreadTimer(td, tid);
                  FlatMatrix<SIMD<double>> simd_proxyvalues(proxy->Dimension(), ir_facet.Size(), lh);
                  for (int k = 0; k < proxy->Dimension(); k++)
                    {
                      ud.testfunction = proxy;
                      ud.test_comp = k;
                      cf -> Evaluate (mir, simd_proxyvalues.Rows(k,k+1));
                    }
                  
                  for (size_t i = 0; i < simd_proxyvalues.Height(); i++)
                    {
                      auto row = simd_proxyvalues.Row(i);
                      for (size_t j = 0; j < row.Size(); j++)
                        row(j) *= mir[j].GetWeight(); //  * simd_ir[j].Weight();
                    }
                  // NgProfiler::StopThreadTimer(td, tid);
                  
                  // NgProfiler::StartThreadTimer(ttrans, tid);
                  proxy->Evaluator()->AddTrans(fel_test, mir, simd_proxyvalues, ely);
                  // NgProfiler::StopThreadTimer(ttrans, tid);                                                  
                }
            }
          return;
        }
      catch (const ExceptionNOSIMD& e)
        {
          cout << IM(6) << e.What() << endl
               << "switching to scalar evaluation" << endl;
          simd_evaluate = false;
          T_ApplyElementMatrixEB<SCAL,SCAL_SHAPES> (fel, trafo, elx, ely, precomputed, lh);
          return;
        }
    
    
    
    for (int k = 0; k < nfacet; k++)
      {
        HeapReset hr(lh);
        ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
        const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, fel_trial.Order()+fel_test.Order()+bonus_intorder);
        IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
        BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
        // mir.ComputeNormalsAndMeasure (eltype, k);
        

        ProxyUserData ud(trial_proxies.Size(), lh);    
        const_cast<ElementTransformation&>(trafo).userdata = &ud;
        ud.fel = &fel_trial;

        PrecomputeCacheCF(cache_cfs, mir, lh);

        for (ProxyFunction * proxy : trial_proxies)
          ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
        
        for (ProxyFunction * proxy : trial_proxies)
          proxy->Evaluator()->Apply(fel_trial, mir, elx, ud.GetMemory(proxy), lh);
        
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
            
            proxy->Evaluator()->ApplyTrans(fel_test, mir, proxyvalues, ely1, lh);
            ely += ely1;
          }
      }
  }







  
  void
  SymbolicBilinearFormIntegrator :: ApplyElementMatrixTrans (const FiniteElement & fel, 
                                                             const ElementTransformation & trafo, 
                                                             const FlatVector<double> elx, 
                                                             FlatVector<double> ely,
                                                             void * precomputed,
                                                             LocalHeap & lh) const
  {
    auto save_userdata = trafo.PushUserData();
    
    if (element_vb != VOL)
      {
        T_ApplyElementMatrixTransEB<double,double> (fel, trafo, elx, ely, precomputed, lh);
        return;
      }

    if (simd_evaluate)
      try
        {
          // static Timer tall("SymbolicBFI::Apply - all", NoTracing, NoTiming); RegionTimer rall(tall);

          bool is_mixed = typeid(fel) == typeid(const MixedFiniteElement&);
          const MixedFiniteElement * mixedfe = static_cast<const MixedFiniteElement*> (&fel);    
          const FiniteElement & fel_trial = is_mixed ? mixedfe->FETrial() : fel;
          const FiniteElement & fel_test = is_mixed ? mixedfe->FETest() : fel;

          HeapReset hr(lh);

          const SIMD_IntegrationRule& simd_ir = Get_SIMD_IntegrationRule (fel, lh);
          auto & simd_mir = trafo(simd_ir, lh);
          
          ProxyUserData ud(test_proxies.Size(), gridfunction_cfs.Size(), lh);
          const_cast<ElementTransformation&>(trafo).userdata = &ud;
          ud.fel = &fel;

          PrecomputeCacheCF(cache_cfs, simd_mir, lh);

          for (ProxyFunction * proxy : test_proxies)
            ud.AssignMemory (proxy, simd_ir.GetNIP(), proxy->Dimension(), lh);
          for (CoefficientFunction * cf : gridfunction_cfs)
            ud.AssignMemory (cf, simd_ir.GetNIP(), cf->Dimension(), lh);
          
          for (ProxyFunction * proxy : test_proxies)
            proxy->Evaluator()->Apply(fel_test, simd_mir, elx, ud.GetAMemory(proxy)); 
          
          ely = 0;
          for (auto proxy : trial_proxies)
            {
              HeapReset hr(lh);

              FlatMatrix<SIMD<double>> simd_proxyvalues(proxy->Dimension(), simd_ir.Size(), lh);
              for (int k = 0; k < proxy->Dimension(); k++)
                {
                  ud.trialfunction = proxy;
                  ud.trial_comp = k;
                  cf -> Evaluate (simd_mir, simd_proxyvalues.Rows(k,k+1));
                }

              for (size_t i = 0; i < simd_proxyvalues.Height(); i++)
                {
                  auto row = simd_proxyvalues.Row(i);
                  for (size_t j = 0; j < row.Size(); j++)
                    row(j) *= simd_mir[j].GetWeight(); //  * simd_ir[j].Weight();
                }
              
              proxy->Evaluator()->AddTrans(fel_trial, simd_mir, simd_proxyvalues, ely); 
            }
          return;
        }
    
      catch (const ExceptionNOSIMD& e)
        {
          cout << IM(6) << e.What() << endl
               << "switching to scalar evaluation" << endl;
          simd_evaluate = false;
          ApplyElementMatrix (fel, trafo, elx, ely, precomputed, lh);
          return;
        }

    HeapReset hr(lh);
    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;
 
    ProxyUserData ud(test_proxies.Size(), lh);    
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;

    IntegrationRule ir = GetIntegrationRule (fel, lh);

    BaseMappedIntegrationRule & mir = trafo(ir, lh);

    PrecomputeCacheCF(cache_cfs, mir, lh);

    for (ProxyFunction * proxy : test_proxies)
      ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);

    for (ProxyFunction * proxy : test_proxies)
      proxy->Evaluator()->Apply(fel_test, mir, elx, ud.GetMemory(proxy), lh);
    
    ely = 0;
    FlatVector<> ely1(ely.Size(), lh);

    FlatMatrix<> val(mir.Size(), 1,lh);
    for (auto proxy : trial_proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            cf -> Evaluate (mir, val);
            proxyvalues.Col(k) = val.Col(0);
          }
        
        for (int i = 0; i < mir.Size(); i++)
          proxyvalues.Row(i) *= mir[i].GetWeight();

        proxy->Evaluator()->ApplyTrans(fel_trial, mir, proxyvalues, ely1, lh);
        ely += ely1;
      }
  }


 
  
  
  template <typename SCAL, typename SCAL_SHAPES>
  void SymbolicBilinearFormIntegrator ::
  T_ApplyElementMatrixTransEB (const FiniteElement & fel, 
                               const ElementTransformation & trafo, 
                               const FlatVector<double> elx, 
                               FlatVector<double> ely,
                               void * precomputed,
                               LocalHeap & lh) const
  {
    // static Timer t("symbolicbfi - Apply EB", NoTracing);
    // static Timer tir("symbolicbfi - Apply EB, intrule", NoTracing);
    // static Timer teval("symbolicbfi - Apply EB, evaluate", NoTracing);
    // static Timer td("symbolicbfi - Apply EB, evaluate D", NoTracing);
    // static Timer ttrans("symbolicbfi - Apply EB, trans", NoTracing);
    
    // RegionTimer reg(t);
    
    ely = 0;

    auto eltype = trafo.GetElementType();

    Facet2ElementTrafo transform(eltype, element_vb); 
    int nfacet = transform.GetNFacets();

    const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fel);
    const FiniteElement & fel_trial = mixedfe ? mixedfe->FETrial() : fel;
    const FiniteElement & fel_test = mixedfe ? mixedfe->FETest() : fel;
    
    if (simd_evaluate)
      try
        {
          for (int k = 0; k < nfacet; k++)
            {
              HeapReset hr(lh);
              // NgProfiler::StartThreadTimer(tir, tid);                              
              ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
              const SIMD_IntegrationRule & ir_facet =
                GetSIMDIntegrationRule(etfacet,fel_trial.Order()+fel_test.Order()+bonus_intorder);
              auto & ir_facet_vol = transform(k, ir_facet, lh);
              auto & mir = trafo(ir_facet_vol, lh);
              // mir.ComputeNormalsAndMeasure (eltype, k);
              // NgProfiler::StopThreadTimer(tir, tid);

              // NgProfiler::StartThreadTimer(teval, tid);                                            
              ProxyUserData ud(test_proxies.Size(), gridfunction_cfs.Size(), lh);              
              const_cast<ElementTransformation&>(trafo).userdata = &ud;
              ud.fel = &fel_trial;

              PrecomputeCacheCF(cache_cfs, mir, lh);
          
              for (ProxyFunction * proxy : test_proxies)
                ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
              for (CoefficientFunction * cf : gridfunction_cfs)
                ud.AssignMemory (cf, ir_facet.GetNIP(), cf->Dimension(), lh);
          
              for (ProxyFunction * proxy : test_proxies)
                proxy->Evaluator()->Apply(fel_test, mir, elx, ud.GetAMemory(proxy)); 
          
              // NgProfiler::StopThreadTimer(teval, tid);
              for (auto proxy : trial_proxies)
                {
                  HeapReset hr(lh);
                  // NgProfiler::StartThreadTimer(td, tid);
                  FlatMatrix<SIMD<double>> simd_proxyvalues(proxy->Dimension(), ir_facet.Size(), lh);
                  for (int k = 0; k < proxy->Dimension(); k++)
                    {
                      ud.trialfunction = proxy;
                      ud.trial_comp = k;
                      cf -> Evaluate (mir, simd_proxyvalues.Rows(k,k+1));
                    }
                  
                  for (size_t i = 0; i < simd_proxyvalues.Height(); i++)
                    {
                      auto row = simd_proxyvalues.Row(i);
                      for (size_t j = 0; j < row.Size(); j++)
                        row(j) *= mir[j].GetWeight(); //  * simd_ir[j].Weight();
                    }
                  // NgProfiler::StopThreadTimer(td, tid);
                  
                  // NgProfiler::StartThreadTimer(ttrans, tid);
                  proxy->Evaluator()->AddTrans(fel_trial, mir, simd_proxyvalues, ely);
                  // NgProfiler::StopThreadTimer(ttrans, tid);                                                  
                }
            }
          return;
        }
      catch (const ExceptionNOSIMD& e)
        {
          cout << IM(6) << e.What() << endl
               << "switching to scalar evaluation" << endl;
          simd_evaluate = false;
          T_ApplyElementMatrixEB<SCAL,SCAL_SHAPES> (fel, trafo, elx, ely, precomputed, lh);
          return;
        }
    
    
    for (int k = 0; k < nfacet; k++)
      {
        HeapReset hr(lh);
        ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
        const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, fel_trial.Order()+fel_test.Order()+bonus_intorder);
        IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
        BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
        // mir.ComputeNormalsAndMeasure (eltype, k);
        

        ProxyUserData ud(test_proxies.Size(), lh);    
        const_cast<ElementTransformation&>(trafo).userdata = &ud;
        ud.fel = &fel_test;

        PrecomputeCacheCF(cache_cfs, mir, lh);
    
        for (ProxyFunction * proxy : test_proxies)
          ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
        
        for (ProxyFunction * proxy : test_proxies)
          proxy->Evaluator()->Apply(fel_test, mir, elx, ud.GetMemory(proxy), lh);
        
        FlatVector<> ely1(ely.Size(), lh);
        FlatMatrix<> val(mir.Size(), 1,lh);
        for (auto proxy : trial_proxies)
          {
            HeapReset hr(lh);
            FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
            for (int k = 0; k < proxy->Dimension(); k++)
              {
                ud.trialfunction = proxy;
                ud.trial_comp = k;
                cf -> Evaluate (mir, val);
                proxyvalues.Col(k) = val.Col(0);
              }
            
            for (size_t i = 0; i < mir.Size(); i++)
              proxyvalues.Row(i) *= ir_facet[i].Weight() * mir[i].GetMeasure();
            
            proxy->Evaluator()->ApplyTrans(fel_trial, mir, proxyvalues, ely1, lh);
            ely += ely1;
          }
      }
  }







  
  SymbolicFacetLinearFormIntegrator ::
  SymbolicFacetLinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb /* , bool eb */)
    : cf(acf), vb(avb) // , element_boundary(eb)
  {
    if (cf->Dimension() != 1)
      throw Exception ("SymbolicLFI needs scalar-valued CoefficientFunction");
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
    cache_cfs = FindCacheCF(*cf);
  }


  template<typename TSCAL>
  void SymbolicFacetLinearFormIntegrator ::
  T_CalcFacetVector (const FiniteElement & fel1, int LocalFacetNr1,
                     const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                     const FiniteElement & fel2, int LocalFacetNr2,
                     const ElementTransformation & trafo2, FlatArray<int> & ElVertices2, 
                     FlatVector<TSCAL> elvec,
                     LocalHeap &lh) const 
  {
    static Timer t("SymbolicFacetLFI::CalcFacetVector - inner", NoTracing);
    HeapReset hr(lh);

    elvec = 0;
    
    int maxorder = max2 (fel1.Order(), fel2.Order());
    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, 2*maxorder+bonus_intorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    Facet2ElementTrafo transform2(eltype2, ElVertices2); 
    IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, ir_facet, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(ir_facet_vol2, lh);

    mir1.SetOtherMIR (&mir2);
    mir2.SetOtherMIR (&mir1);

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    PrecomputeCacheCF(cache_cfs, mir1, lh);

    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);
    mir2.ComputeNormalsAndMeasure (eltype2, LocalFacetNr2);
    
    FlatMatrix<TSCAL> val(ir_facet.Size(), 1,lh);
    for (auto proxy : proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<TSCAL> proxyvalues(ir_facet.Size(), proxy->Dimension(), lh);

        
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

        IntRange range  = proxy->IsOther() ? IntRange(proxy->Evaluator()->BlockDim()*fel1.GetNDof(), elvec.Size()) : IntRange(0, proxy->Evaluator()->BlockDim()*fel1.GetNDof());

        FlatVector<TSCAL> localvec(range.Size(), lh);
        localvec = 0.0;

        if (proxy->IsOther())
          proxy->Evaluator()->ApplyTrans(fel2, mir2, proxyvalues, localvec, lh);
        else
          proxy->Evaluator()->ApplyTrans(fel1, mir1, proxyvalues, localvec, lh);
        elvec.Range(range) += localvec;
      }
  }



  template<typename TSCAL>
  void SymbolicFacetLinearFormIntegrator ::
  T_CalcFacetVector (const FiniteElement & fel1, int LocalFacetNr,
                     const ElementTransformation & trafo1, FlatArray<int> & ElVertices,
                     const ElementTransformation & strafo,
                     FlatVector<TSCAL> elvec,
                     LocalHeap & lh) const
  {
    static Timer t("SymbolicFacetLFI::CalcFacetVector - boundary", NoTracing);
    HeapReset hr(lh);

    elvec = 0;
    
    FlatVector<TSCAL> elvec1(elvec.Size(), lh);

    int maxorder = fel1.Order();

    auto eltype1 = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

    const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, 2*maxorder+bonus_intorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr, ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    auto & smir = strafo(ir_facet, lh);

    mir1.SetOtherMIR (&smir);
    smir.SetOtherMIR (&mir1);

    // evaluate proxy-values
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    const_cast<ElementTransformation&>(strafo).userdata = &ud;

    PrecomputeCacheCF(cache_cfs, mir1, lh);

    RegionTimer reg(t);
    
    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr);
    
    FlatMatrix<TSCAL> val(ir_facet.Size(), 1,lh);
    for (auto proxy : proxies)
      {
        HeapReset hr(lh);
        FlatMatrix<TSCAL> proxyvalues(ir_facet.Size(), proxy->Dimension(), lh);
        
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

  void SymbolicFacetLinearFormIntegrator::
  CalcFacetVector (const FiniteElement & volumefel1, int LocalFacetNr1,
          		   const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
    		       const FiniteElement & volumefel2, int LocalFacetNr2,
    		       const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
    		       FlatVector<double> elvec,
    		       LocalHeap & lh) const
  {
    T_CalcFacetVector(volumefel1, LocalFacetNr1, eltrans1, ElVertices1, 
                      volumefel2, LocalFacetNr2, eltrans2, ElVertices2, 
                      elvec, lh);
  }
  
  void SymbolicFacetLinearFormIntegrator::
  CalcFacetVector (const FiniteElement & volumefel1, int LocalFacetNr1,
    		       const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
    		       const FiniteElement & volumefel2, int LocalFacetNr2,
    		       const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,	 
    		       FlatVector<Complex> elvec,
    		       LocalHeap & lh) const
  {
    T_CalcFacetVector(volumefel1, LocalFacetNr1, eltrans1, ElVertices1, 
                      volumefel2, LocalFacetNr2, eltrans2, ElVertices2, 
                      elvec, lh);
  }

  void SymbolicFacetLinearFormIntegrator::
  CalcFacetVector(const FiniteElement & volumefel, int LocalFacetNr,
                  const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                  const ElementTransformation & seltrans,
                  FlatVector<double> elvec,
                  LocalHeap & lh) const
  {
    T_CalcFacetVector(volumefel, LocalFacetNr, eltrans, ElVertices, seltrans,
                      elvec, lh);
  }

  void SymbolicFacetLinearFormIntegrator::
  CalcFacetVector(const FiniteElement & volumefel, int LocalFacetNr,
                  const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                  const ElementTransformation & seltrans,
                  FlatVector<Complex> elvec,
                  LocalHeap & lh) const
  {
    T_CalcFacetVector(volumefel, LocalFacetNr, eltrans, ElVertices, seltrans,
                      elvec, lh);
  }
  
  SymbolicFacetBilinearFormIntegrator ::
  SymbolicFacetBilinearFormIntegrator (shared_ptr<CoefficientFunction> acf, VorB avb, bool eb)
    : cf(acf), vb(avb), element_boundary(eb)
  {
    simd_evaluate = true;
    
    if (cf->Dimension() != 1)
        throw Exception ("SymbolicBFI needs scalar-valued CoefficientFunction");
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
          else
            if (nodecf.StoreUserData() && !gridfunction_cfs.Contains(&nodecf))
              gridfunction_cfs.Append (&nodecf);
        });

    neighbor_testfunction = false;
    for (auto proxy : test_proxies)
      if (proxy->IsOther())
        neighbor_testfunction = true;

    cache_cfs = FindCacheCF(*cf);

    // comment in for experimental new Apply
    dcf_dtest.SetSize(test_proxies.Size());
    if (symbolic_integrator_uses_diff)    
      for (int i = 0; i < test_proxies.Size(); i++)
        {
          try
            {
              CoefficientFunction::T_DJC cache;              
              dcf_dtest[i] = cf->DiffJacobi(test_proxies[i], cache);
              // cout << "dcf_dtest = " << *dcf_dtest[i] << endl;
            }
          catch (const Exception& e)
            {
            cout << IM(5) << "dcf_dtest has thrown exception " << e.What() << endl;
            }
        }

    
    cout << IM(6) << "num test_proxies " << test_proxies.Size() << endl;
    cout << IM(6) << "num trial_proxies " << trial_proxies.Size() << endl;
    cout << IM(6) << "cumulated test_proxy dims  " << test_cum << endl;
    cout << IM(6) << "cumulated trial_proxy dims " << trial_cum << endl;
  }

  void SymbolicFacetBilinearFormIntegrator::
  CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                   const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                   const FiniteElement & volumefel2, int LocalFacetNr2,
                   const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                   FlatMatrix<double> elmat,
                   LocalHeap & lh) const
  {
    T_CalcFacetMatrix(volumefel1, LocalFacetNr1, eltrans1, ElVertices1,
                      volumefel2, LocalFacetNr2, eltrans2, ElVertices2,
                      elmat, lh);
  }

  void SymbolicFacetBilinearFormIntegrator::
  CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                   const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                   const ElementTransformation & seltrans, FlatArray<int> & SElVertices,  
                   FlatMatrix<double> elmat,
                   LocalHeap & lh) const
  {
    T_CalcFacetMatrix(volumefel, LocalFacetNr, eltrans, ElVertices,
                      seltrans, SElVertices, elmat, lh);
  }

  void SymbolicFacetBilinearFormIntegrator::
  CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                   const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                   const FiniteElement & volumefel2, int LocalFacetNr2,
                   const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                   FlatMatrix<Complex> elmat,
                   LocalHeap & lh) const
  {
      if (volumefel1.ComplexShapes() && volumefel2.ComplexShapes())
    T_CalcFacetMatrix<Complex,Complex>(volumefel1, LocalFacetNr1, eltrans1, ElVertices1,
                      volumefel2, LocalFacetNr2, eltrans2, ElVertices2,
                      elmat, lh);
      else
    T_CalcFacetMatrix<Complex,double>(volumefel1, LocalFacetNr1, eltrans1, ElVertices1,
                      volumefel2, LocalFacetNr2, eltrans2, ElVertices2,
                      elmat, lh);
  }

  void SymbolicFacetBilinearFormIntegrator::
  CalcFacetMatrix (const FiniteElement & volumefel, int LocalFacetNr,
                   const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
                   const ElementTransformation & seltrans, FlatArray<int> & SElVertices,  
                   FlatMatrix<Complex> elmat,
                   LocalHeap & lh) const
  {
      if (volumefel.ComplexShapes())
    T_CalcFacetMatrix<Complex,Complex>(volumefel, LocalFacetNr, eltrans, ElVertices,
                      seltrans, SElVertices, elmat, lh);
      else
    T_CalcFacetMatrix<Complex,double>(volumefel, LocalFacetNr, eltrans, ElVertices,
                      seltrans, SElVertices, elmat, lh);
  }

  template<typename TSCAL, typename SCAL_SHAPES>
  void SymbolicFacetBilinearFormIntegrator ::
  T_CalcFacetMatrix (const FiniteElement & fel1, int LocalFacetNr1,
                     const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                     const FiniteElement & fel2, int LocalFacetNr2,
                     const ElementTransformation & trafo2, FlatArray<int> & ElVertices2,
                     FlatMatrix<TSCAL> elmat,
                     LocalHeap & lh) const
  {
    elmat = 0.0;

    if (LocalFacetNr2==-1) throw Exception ("SymbolicFacetBFI: LocalFacetNr2==-1");

    bool is_mixedfe1 = typeid(fel1) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe1 = static_cast<const MixedFiniteElement*> (&fel1);
    const FiniteElement & fel1_trial = is_mixedfe1 ? mixedfe1->FETrial() : fel1;
    const FiniteElement & fel1_test = is_mixedfe1 ? mixedfe1->FETest() : fel1;
    bool is_mixedfe2 = typeid(fel2) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe2 = static_cast<const MixedFiniteElement*> (&fel2);
    const FiniteElement & fel2_trial = is_mixedfe2 ? mixedfe2->FETrial() : fel2;
    const FiniteElement & fel2_test = is_mixedfe2 ? mixedfe2->FETest() : fel2;

    if (is_mixedfe1 != is_mixedfe2) throw Exception("both sides should have the same type of mixity");
    
    int maxorder = max2(max2 (fel1_trial.Order(), fel1_test.Order()), max2 (fel2_trial.Order(), fel2_test.Order()));

    auto eltype1 = trafo1.GetElementType();
    auto eltype2 = trafo2.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

    const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, 2*maxorder+bonus_intorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    
    Facet2ElementTrafo transform2(eltype2, ElVertices2); 
    IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, ir_facet, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(ir_facet_vol2, lh);

    mir1.SetOtherMIR (&mir2);
    mir2.SetOtherMIR (&mir1);

    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    PrecomputeCacheCF(cache_cfs, mir1, lh);

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<TSCAL> val(mir1.Size(), 1,lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];

          FlatTensor<3,TSCAL> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());

          mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);
          mir2.ComputeNormalsAndMeasure (eltype2, LocalFacetNr2);
          
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

          IntRange trial_range  = proxy1->IsOther() ? IntRange(proxy1->Evaluator()->BlockDim()*fel1_trial.GetNDof(), elmat.Width()) : IntRange(0, proxy1->Evaluator()->BlockDim()*fel1_trial.GetNDof());
          IntRange test_range  = proxy2->IsOther() ? IntRange(proxy2->Evaluator()->BlockDim()*fel1_test.GetNDof(), elmat.Height()) : IntRange(0, proxy2->Evaluator()->BlockDim()*fel1_test.GetNDof());

          auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), loc_elmat.Width(), lh);
          FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), loc_elmat.Height(), lh);

          constexpr size_t BS = 16;
          for (size_t i = 0; i < mir1.Size(); i+=BS)
            {
              int rest = min2(size_t(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<TSCAL,ColMajor> bdbmat1(rest*proxy2->Dimension(), loc_elmat.Width(), lh);
              FlatMatrix<TSCAL,ColMajor> bbmat2(rest*proxy2->Dimension(), loc_elmat.Height(), lh);

              for (int j = 0; j < rest; j++)
                {
                  int ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  if (proxy1->IsOther())
                    proxy1->Evaluator()->CalcMatrix(fel2_trial, mir2[ii], bmat1, lh);
                  else
                    proxy1->Evaluator()->CalcMatrix(fel1_trial, mir1[ii], bmat1, lh);
                  
                  if (proxy2->IsOther())
                    proxy2->Evaluator()->CalcMatrix(fel2_test, mir2[ii], bmat2, lh);
                  else
                    proxy2->Evaluator()->CalcMatrix(fel1_test, mir1[ii], bmat2, lh);
                  
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(proxy1->IsOther() ? fel2_trial : fel1_trial);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(proxy2->IsOther() ? fel2_test : fel1_test);
              loc_elmat.Rows(r2).Cols(r1) += Trans (bbmat2.Cols(r2)) * bdbmat1.Cols(r1);
            }
        }
  }

  template<typename TSCAL, typename SCAL_SHAPES>
  void SymbolicFacetBilinearFormIntegrator ::
  T_CalcFacetMatrix(const FiniteElement & fel1, int LocalFacetNr1,
                    const ElementTransformation & trafo1, FlatArray<int> & ElVertices1,
                    const ElementTransformation & strafo, FlatArray<int> & SElVertices1,
                    FlatMatrix<TSCAL> elmat,
                    LocalHeap & lh) const
  {
    bool is_mixedfe1 = typeid(fel1) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe1 = static_cast<const MixedFiniteElement*> (&fel1);
    const FiniteElement & fel1_trial = is_mixedfe1 ? mixedfe1->FETrial() : fel1;
    const FiniteElement & fel1_test = is_mixedfe1 ? mixedfe1->FETest() : fel1;

    elmat = 0.0;

    //int maxorder = fel1.Order();
    int maxorder = max2 (fel1_trial.Order(), fel1_test.Order());

    auto etvol = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (etvol, LocalFacetNr1);

    const IntegrationRule& ir_facet = GetIntegrationRule(etfacet, 2*maxorder+bonus_intorder);
    Facet2ElementTrafo transform1(etvol, ElVertices1);
    Facet2SurfaceElementTrafo stransform(strafo.GetElementType(), SElVertices1); 
    
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    // IntegrationRule & ir_facet_surf = stransform(ir_facet, lh);  // not yet used ???
    
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    mir1.ComputeNormalsAndMeasure (etvol, LocalFacetNr1);          
    
    ProxyUserData ud;
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;

    PrecomputeCacheCF(cache_cfs, mir1, lh);

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(test_proxies))
        {
          HeapReset hr(lh);
          FlatMatrix<TSCAL> val(mir1.Size(), 1, lh);
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = test_proxies[l1];
          if (proxy1->IsOther() || proxy2->IsOther()) continue;

          FlatTensor<3, TSCAL> proxyvalues(lh, mir1.Size(), proxy2->Dimension(), proxy1->Dimension());
          
          for (size_t k = 0; k < proxy1->Dimension(); k++)
            for (size_t l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;
                
                cf -> Evaluate (mir1, val);
                proxyvalues(STAR,l,k) = val.Col(0);
              }

          for (size_t i = 0; i < mir1.Size(); i++)
            proxyvalues(i,STAR,STAR) *=  mir1[i].GetMeasure() * ir_facet[i].Weight();

          FlatMatrix<SCAL_SHAPES,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<SCAL_SHAPES,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

          constexpr size_t BS = 16;          
          for (size_t i = 0; i < mir1.Size(); i+=BS)
            {
              size_t rest = min2(size_t(BS), mir1.Size()-i);
              HeapReset hr(lh);
              FlatMatrix<TSCAL,ColMajor> bdbmat1(rest*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<TSCAL,ColMajor> bbmat2(rest*proxy2->Dimension(), elmat.Height(), lh);

              for (size_t j = 0; j < rest; j++)
                {
                  size_t ii = i+j;
                  IntRange r2 = proxy2->Dimension() * IntRange(j,j+1);
                  proxy1->Evaluator()->CalcMatrix(fel1_trial, mir1[ii], bmat1, lh);
                  proxy2->Evaluator()->CalcMatrix(fel1_test, mir1[ii], bmat2, lh);
                  bdbmat1.Rows(r2) = proxyvalues(ii,STAR,STAR) * bmat1;
                  bbmat2.Rows(r2) = bmat2;
                }

              IntRange r1 = proxy1->Evaluator()->UsedDofs(fel1_trial);
              IntRange r2 = proxy2->Evaluator()->UsedDofs(fel1_test);
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

    const IntegrationRule& ir_facet = GetIntegrationRule(etfacet, 2*maxorder+bonus_intorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices1);
    Facet2SurfaceElementTrafo stransform(strafo.GetElementType(), SElVertices1); 
    
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    // IntegrationRule & ir_facet_surf = stransform(ir_facet, lh);  // not yet used ???
    
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);          
    
    // new 
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    ud.fel = &fel1;

    // ud.elx = &elveclin;
    // ud.lh = &lh;
    for (ProxyFunction * proxy : trial_proxies)
      {
        ud.AssignMemory (proxy, ir_facet_vol1.Size(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel1, mir1, elveclin, ud.GetMemory(proxy), lh);
      }
    
    // FlatMatrix<> val(mir1.Size(), 1, lh), deriv(mir1.Size(), 1, lh);
    FlatMatrix<AutoDiff<1>> dval(mir1.Size(), 1, lh);
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
                
                // cf -> EvaluateDeriv (mir1, val, deriv);
                // proxyvalues(STAR,l,k) = deriv.Col(0);

                cf -> Evaluate (mir1, dval);
                for (size_t i = 0; i < mir1.Size(); i++)
                  proxyvalues(i,l,k) = dval(i,0).DValue(0);
              }

          for (int i = 0; i < mir1.Size(); i++)
            proxyvalues(i,STAR,STAR) *=  mir1[i].GetMeasure() * ir_facet[i].Weight();

          // auto loc_elmat = elmat.Rows(test_range).Cols(trial_range);
          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

          constexpr size_t BS = 16;          
          for (size_t i = 0; i < mir1.Size(); i+=BS)
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
    bool is_mixedfe1 = typeid(fel1) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe1 = static_cast<const MixedFiniteElement*> (&fel1);
    const FiniteElement & fel1_trial = is_mixedfe1 ? mixedfe1->FETrial() : fel1;
    const FiniteElement & fel1_test = is_mixedfe1 ? mixedfe1->FETest() : fel1;
    bool is_mixedfe2 = typeid(fel2) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe2 = static_cast<const MixedFiniteElement*> (&fel2);
    const FiniteElement & fel2_trial = is_mixedfe2 ? mixedfe2->FETrial() : fel2;
    const FiniteElement & fel2_test = is_mixedfe2 ? mixedfe2->FETest() : fel2;

    if (simd_evaluate)
      {
        try
          {
            static Timer tall("SymbolicFacetBFI::Apply - all", NoTracing); // RegionTimer rall(tall);
            static Timer tstart("SymbolicFacetBFI::Apply - startup", NoTracing);
            static Timer tapply("SymbolicFacetBFI::Apply - apply", NoTracing);
            static Timer tcoef("SymbolicFacetBFI::Apply - coef", NoTracing);
            static Timer tapplyt("SymbolicFacetBFI::Apply - apply-trans", NoTracing); 

            RegionTimer reg(tall);
            // RegionTracer rt(TaskManager::GetThreadId(), tall);
            HeapReset hr(lh);
            // tall.Start();
            
            // tstart.Start();
            
            ely = 0;
            
            int maxorder = max2(max2 (fel1_trial.Order(), fel1_test.Order()), max2 (fel2_trial.Order(), fel2_test.Order()));
            
            auto eltype1 = trafo1.GetElementType();
            auto eltype2 = trafo2.GetElementType();
            auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);
            
            Facet2ElementTrafo transform1(eltype1, ElVertices1); 
            Facet2ElementTrafo transform2(eltype2, ElVertices2); 
            
            const SIMD_IntegrationRule& simd_ir_facet = GetSIMDIntegrationRule(etfacet, 2*maxorder+bonus_intorder);
            
            auto & simd_ir_facet_vol1 = transform1(LocalFacetNr1, simd_ir_facet, lh);
            auto & simd_mir1 = trafo1(simd_ir_facet_vol1, lh);
            
            auto & simd_ir_facet_vol2 = transform2(LocalFacetNr2, simd_ir_facet, lh);
            auto & simd_mir2 = trafo2(simd_ir_facet_vol2, lh);

            simd_mir1.SetOtherMIR(&simd_mir2);
            simd_mir2.SetOtherMIR(&simd_mir1);

            
            simd_mir1.ComputeNormalsAndMeasure(eltype1, LocalFacetNr1);
            simd_mir2.ComputeNormalsAndMeasure(eltype2, LocalFacetNr2);
            
            // evaluate proxy-values
            ProxyUserData ud(trial_proxies.Size(), gridfunction_cfs.Size(), lh);
            const_cast<ElementTransformation&>(trafo1).userdata = &ud;
            ud.fel = &fel1_trial;   // necessary to check remember-map

            PrecomputeCacheCF(cache_cfs, simd_mir1, lh);

            for (ProxyFunction * proxy : trial_proxies)
              ud.AssignMemory (proxy, simd_ir_facet.GetNIP(), proxy->Dimension(), lh);
            for (CoefficientFunction * cf : gridfunction_cfs)
              ud.AssignMemory (cf, simd_ir_facet.GetNIP(), cf->Dimension(), lh);
            // tstart.Stop();
            // tapply.Start();

            // RegionTracer rtapply(TaskManager::GetThreadId(), tapply);                        
            for (ProxyFunction * proxy : trial_proxies)
              {
                IntRange trial_range  = proxy->IsOther() ?
                  IntRange(proxy->Evaluator()->BlockDim()*fel1_trial.GetNDof(), elx.Size()) :
                  IntRange(0, proxy->Evaluator()->BlockDim()*fel1_trial.GetNDof());
                
                if (proxy->IsOther())
                  proxy->Evaluator()->Apply(fel2_trial, simd_mir2, elx.Range(trial_range), ud.GetAMemory(proxy));
                else
                  proxy->Evaluator()->Apply(fel1_trial, simd_mir1, elx.Range(trial_range), ud.GetAMemory(proxy));
                // tapply.AddFlops (trial_range.Size() * simd_ir_facet_vol1.GetNIP());
              }
            // tapply.Stop();

            // RegionTracer rtapplyt(TaskManager::GetThreadId(), tapplyt);            
            // for (auto proxy : test_proxies)
            for (auto i : Range(test_proxies))
              {
                auto proxy = test_proxies[i];
                HeapReset hr(lh);
                // tcoef.Start();
                FlatMatrix<SIMD<double>> simd_proxyvalues(proxy->Dimension(), simd_ir_facet.Size(), lh);        

                if (dcf_dtest[i])
                  dcf_dtest[i]->Evaluate (simd_mir1, simd_proxyvalues);
                else
                  for (int k = 0; k < proxy->Dimension(); k++)
                    {
                      ud.testfunction = proxy;
                      ud.test_comp = k;
                      cf -> Evaluate (simd_mir1, simd_proxyvalues.Rows(k,k+1));
                    }
                
                for (int i = 0; i < simd_proxyvalues.Height(); i++)
                  {
                    auto row = simd_proxyvalues.Row(i);
                    for (int j = 0; j < row.Size(); j++)
                      row(j) *= simd_mir1[j].GetMeasure() * simd_ir_facet[j].Weight();
                  }
                // tcoef.Stop();
                // tapplyt.Start();
                IntRange test_range  = proxy->IsOther() ? IntRange(fel1_test.GetNDof(), elx.Size()) : IntRange(0, fel1_test.GetNDof());
                int blockdim = proxy->Evaluator()->BlockDim();
                test_range = blockdim * test_range;
                
                if (proxy->IsOther())
                  proxy->Evaluator()->AddTrans(fel2_test, simd_mir2, simd_proxyvalues, ely.Range(test_range));
                else
                  proxy->Evaluator()->AddTrans(fel1_test, simd_mir1, simd_proxyvalues, ely.Range(test_range));
                // tapplyt.AddFlops (test_range.Size() * simd_ir_facet_vol1.GetNIP());                
                // tapplyt.Stop();
              }
            // tall.Stop();
          }
        catch (const ExceptionNOSIMD& e)
          {
            cout << IM(6) << "caught in SymbolicFacetInegtrator::Apply: " << endl
                 << e.What() << endl;
            simd_evaluate = false;
            ApplyFacetMatrix (fel1, LocalFacetNr1, trafo1, ElVertices1,
                              fel2, LocalFacetNr2, trafo2, ElVertices2,
                              elx, ely, lh);
          }
        return;
      }
        
    
    static Timer tall("SymbolicFacetBFI::Apply - all", NoTracing); RegionTimer rall(tall);
    /*
    static Timer t("SymbolicFacetBFI::Apply", NoTracing);
    static Timer ts1("SymbolicFacetBFI::Apply start 1", NoTracing);
    static Timer ts2("SymbolicFacetBFI::Apply start 2", NoTracing);
    static Timer t1("SymbolicFacetBFI::Apply 1", NoTracing);
    static Timer t2("SymbolicFacetBFI::Apply 2", NoTracing);
    static Timer t3("SymbolicFacetBFI::Apply 3", NoTracing);
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

    const IntegrationRule& ir_facet = GetIntegrationRule(etfacet, 2*maxorder+bonus_intorder);
    
    Facet2ElementTrafo transform1(eltype1, ElVertices1); 
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr1, ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    
    Facet2ElementTrafo transform2(eltype2, ElVertices2); 
    IntegrationRule & ir_facet_vol2 = transform2(LocalFacetNr2, ir_facet, lh);
    BaseMappedIntegrationRule & mir2 = trafo2(ir_facet_vol2, lh);

    mir1.SetOtherMIR (&mir2);
    mir2.SetOtherMIR (&mir1);

    // ts1.Stop();
    // ts2.Start();

    
    // evaluate proxy-values
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    ud.fel = &fel1;   // necessary to check remember-map

    PrecomputeCacheCF(cache_cfs, mir1, lh);

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

    mir1.ComputeNormalsAndMeasure (eltype1, LocalFacetNr1);
    mir2.ComputeNormalsAndMeasure (eltype2, LocalFacetNr2);

    
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
  CalcTraceValues (const FiniteElement & volumefel, int LocalFacetNr,
		   const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
		   FlatVector<double> & trace, FlatVector<double> elx, LocalHeap & lh) const
  {

    if (simd_evaluate)
      {
        try
          {
	    int maxorder = volumefel.Order();
            
            auto eltype = eltrans.GetElementType();
            auto etfacet = ElementTopology::GetFacetType (eltype, LocalFacetNr);
            
            Facet2ElementTrafo transform(eltype, ElVertices); 
            
            const SIMD_IntegrationRule & simd_ir_facet = GetSIMDIntegrationRule(etfacet, 2*maxorder+bonus_intorder);
	    
            auto & simd_ir_facet_vol = transform(LocalFacetNr, simd_ir_facet, lh);
            auto & simd_mir = eltrans(simd_ir_facet_vol, lh);
            simd_mir.ComputeNormalsAndMeasure(eltype, LocalFacetNr);

            ProxyUserData ud(trial_proxies.Size(), lh);
            const_cast<ElementTransformation&>(eltrans).userdata = &ud;
            ud.fel = &volumefel;   // necessary to check remember-map
            // ud.lh = &lh;
            for (ProxyFunction * proxy : trial_proxies)
	      if(proxy->IsOther())
		{
		  ud.AssignMemory (proxy, simd_ir_facet.GetNIP(), proxy->Dimension(), lh);
		  IntRange trial_range  = IntRange(0, volumefel.GetNDof());
		  trial_range = proxy->Evaluator()->BlockDim() * trial_range;
		  proxy->Evaluator()->Apply(volumefel, simd_mir, elx.Range(trial_range), ud.GetAMemory(proxy)); // , lh);
		}
	    
	    size_t st = 0;
	    for (ProxyFunction * proxy : trial_proxies)
              if(proxy->IsOther())
		st += ud.GetAMemory(proxy).AsVector().Size()*SIMD<double>::Size();
		
	    
	    trace.AssignMemory(st, lh);
	    st = 0;
	    
	    for (ProxyFunction * proxy : trial_proxies)
              if(proxy->IsOther())
		for(auto k:Range(ud.GetAMemory(proxy).AsVector().Size()))
		  for(auto j:Range(SIMD<double>::Size()))
		    trace [st++] = ud.GetAMemory(proxy)(k)[j];

	    return;
	  }
        catch (const ExceptionNOSIMD& e)
          {
            cout << IM(6) << "caught in SymbolicFacetInegtrator::CalcTraceValues: " << endl
                 << e.What() << endl;
            simd_evaluate = false;
	    CalcTraceValues (volumefel, LocalFacetNr, eltrans, ElVertices, trace, elx, lh);
          }
	return;
      }  

    int maxorder = volumefel.Order();

    auto eltype = eltrans.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype, LocalFacetNr);

    const IntegrationRule& ir_facet = GetIntegrationRule(etfacet, 2*maxorder + bonus_intorder);
    
    Facet2ElementTrafo transform(eltype, ElVertices); 
    IntegrationRule & ir_facet_vol = transform(LocalFacetNr, ir_facet, lh);
    BaseMappedIntegrationRule & mir = eltrans(ir_facet_vol, lh);
    mir.ComputeNormalsAndMeasure(eltype, LocalFacetNr);
    
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(eltrans).userdata = &ud;
    ud.fel = &volumefel;   // necessary to check remember-map
    // ud.lh = &lh;

    for (ProxyFunction * proxy : trial_proxies)
      if(proxy->IsOther())
      {
	ud.AssignMemory (proxy, ir_facet.Size(), proxy->Dimension(), lh);
	IntRange trial_range  = IntRange(0, proxy->Evaluator()->BlockDim()*volumefel.GetNDof());
	proxy->Evaluator()->Apply(volumefel, mir, elx.Range(trial_range), ud.GetMemory(proxy), lh);
      }
    
    size_t st = 0;
    for (ProxyFunction * proxy : trial_proxies)
      if(proxy->IsOther())
	{
	  st += ud.GetMemory(proxy).AsVector().Size();
	}

    trace.AssignMemory(st, lh);
    st = 0;
    
    for (ProxyFunction * proxy : trial_proxies)
      if(proxy->IsOther())
	for(auto k:Range(ud.GetMemory(proxy).AsVector().Size()))
	  trace[st++] = ud.GetMemory(proxy).AsVector()[k];
    
  }//end CalcTraceValues (double)

  void SymbolicFacetBilinearFormIntegrator ::
  ApplyFromTraceValues (const FiniteElement & volumefel, int LocalFacetNr,
			const ElementTransformation & eltrans, FlatArray<int> & ElVertices,
			FlatVector<double> trace,
			FlatVector<double> elx, FlatVector<double> ely, 
			LocalHeap & lh) const
  {
    if (simd_evaluate)
      {
        try
          {
            ely = 0.0;

	    int maxorder = volumefel.Order();
	    auto eltype = eltrans.GetElementType();
            auto etfacet = ElementTopology::GetFacetType (eltype, LocalFacetNr);

            Facet2ElementTrafo transform(eltype, ElVertices); 
            const SIMD_IntegrationRule & simd_ir_facet = GetSIMDIntegrationRule(etfacet, 2*maxorder+bonus_intorder);

            auto & simd_ir_facet_vol = transform(LocalFacetNr, simd_ir_facet, lh);
            auto & simd_mir = eltrans(simd_ir_facet_vol, lh);
            // simd_mir.ComputeNormalsAndMeasure(eltype, LocalFacetNr);

            ProxyUserData ud(trial_proxies.Size(), lh);
            const_cast<ElementTransformation&>(eltrans).userdata = &ud;
            ud.fel = &volumefel;   // necessary to check remember-map
            // ud.elx = &elx;
            // ud.lh = &lh;

            PrecomputeCacheCF(cache_cfs, simd_mir, lh);
	    
	    size_t ctrace = 0;	    
	    for (ProxyFunction * proxy : trial_proxies)
	      {
		ud.AssignMemory (proxy, simd_ir_facet.GetNIP(), proxy->Dimension(), lh);
		
		if(proxy->IsOther())
		  {
		    auto mem = ud.GetAMemory(proxy);
 		    for(auto k:Range(mem.AsVector().Size()))
		      {
			mem(k) = SIMD<double>(&trace[ctrace]);
			ctrace+=SIMD<double>::Size();
		      }
		  }
		else
		  {
		    IntRange trial_range  = IntRange(0, volumefel.GetNDof());
		    trial_range = proxy->Evaluator()->BlockDim() * trial_range;
		    proxy->Evaluator()->Apply(volumefel, simd_mir, elx.Range(trial_range), ud.GetAMemory(proxy));
		  }
		
	      }

	    for (auto proxy : test_proxies)
              if(!proxy->IsOther())
		{
		  HeapReset hr(lh);
		  FlatMatrix<SIMD<double>> simd_proxyvalues(proxy->Dimension(), simd_ir_facet.Size(), lh); 
                
		  for (int k = 0; k < proxy->Dimension(); k++)
		    {
		      ud.testfunction = proxy;
		      ud.test_comp = k;
		      cf -> Evaluate (simd_mir, simd_proxyvalues.Rows(k,k+1));
		    }
                
		  for (int i = 0; i < simd_proxyvalues.Height(); i++)
		    {
		      auto row = simd_proxyvalues.Row(i);
		      for (int j = 0; j < row.Size(); j++)
			row(j) *= simd_mir[j].GetMeasure() * simd_ir_facet[j].Weight();
		    }
		  IntRange test_range  = IntRange(0, volumefel.GetNDof());
		  int blockdim = proxy->Evaluator()->BlockDim();
		  test_range = blockdim * test_range;
                
		  proxy->Evaluator()->AddTrans(volumefel, simd_mir, simd_proxyvalues, ely.Range(test_range));
		}
	  }
        catch (const ExceptionNOSIMD& e)
          {
            cout << IM(6) << "caught in SymbolicFacetInegtrator::CalcTraceValues: " << endl
                 << e.What() << endl;
            simd_evaluate = false;
	    ApplyFromTraceValues(volumefel, LocalFacetNr, eltrans, ElVertices,
				 trace, elx, ely, lh);
          }
	return;
      }

    ely = 0.0;
    FlatVector<> ely1(ely.Size(), lh);

    int maxorder = volumefel.Order();
    auto eltype = eltrans.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype, LocalFacetNr);
    const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, 2*maxorder+bonus_intorder);

    Facet2ElementTrafo transform(eltype, ElVertices); 
    IntegrationRule & ir_facet_vol = transform(LocalFacetNr, ir_facet, lh);
    BaseMappedIntegrationRule & mir = eltrans(ir_facet_vol, lh);
    mir.ComputeNormalsAndMeasure(eltype, LocalFacetNr);
        
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(eltrans).userdata = &ud;
    ud.fel = &volumefel;   // necessary to check remember-map
    // ud.lh = &lh;

    PrecomputeCacheCF(cache_cfs, mir, lh);

    size_t ctrace = 0;
    //copy traces to user memory
    for (ProxyFunction * proxy : trial_proxies)
      {
	ud.AssignMemory (proxy, ir_facet.Size(), proxy->Dimension(), lh);
	if(proxy->IsOther())
	  {
	    auto mem = ud.GetMemory(proxy).AsVector();
	    for(auto k:Range(mem.Size()))
	      mem[k] = trace[ctrace++];
	  }
	else
	  {
	    IntRange trial_range  = IntRange(0, proxy->Evaluator()->BlockDim()*volumefel.GetNDof());
	    proxy->Evaluator()->Apply(volumefel, mir, elx.Range(trial_range), ud.GetMemory(proxy), lh);
	  }
      }

    FlatMatrix<> val(ir_facet.Size(), 1,lh);
    for (auto proxy : test_proxies)
      if(!proxy->IsOther())
	{
	  HeapReset hr(lh);

	  FlatMatrix<> proxyvalues(ir_facet.Size(), proxy->Dimension(), lh);

	  for (int k = 0; k < proxy->Dimension(); k++)
	    {
	      ud.testfunction = proxy;
	      ud.test_comp = k;
	      cf -> Evaluate (mir, val);
	      proxyvalues.Col(k) = val.Col(0);
	    }


	  for (int i = 0; i < mir.Size(); i++)
	    proxyvalues.Row(i) *= mir[i].GetMeasure() * ir_facet[i].Weight();

	  ely1 = 0.0;
	  IntRange test_range  = IntRange(0, proxy->Evaluator()->BlockDim()*volumefel.GetNDof());
	  
	  proxy->Evaluator()->ApplyTrans(volumefel, mir, proxyvalues, ely1.Range(test_range), lh);

	  ely += ely1;
	}
  }//end ApplyFromTraceValues

  void SymbolicFacetBilinearFormIntegrator ::
  ApplyFacetMatrix (const FiniteElement & fel1, int LocalFacetNr,
                    const ElementTransformation & trafo1, FlatArray<int> & ElVertices,
                    const ElementTransformation & strafo, FlatArray<int> & SElVertices,
                    FlatVector<double> elx, FlatVector<double> ely,
                    LocalHeap & lh) const
  {
    bool is_mixedfe1 = typeid(fel1) == typeid(const MixedFiniteElement&);
    const MixedFiniteElement * mixedfe1 = static_cast<const MixedFiniteElement*> (&fel1);
    const FiniteElement & fel1_trial = is_mixedfe1 ? mixedfe1->FETrial() : fel1;
    const FiniteElement & fel1_test = is_mixedfe1 ? mixedfe1->FETest() : fel1;
    if (simd_evaluate)
      {
        try
          {
            static Timer t("SymbolicFacetBFI::ApplyFacetMatrix - boundary", NoTracing);
            
            HeapReset hr(lh);
            
            ely = 0;
            
            int maxorder = max2 (fel1_trial.Order(), fel1_test.Order());
            
            auto eltype1 = trafo1.GetElementType();
            auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);
            
            const SIMD_IntegrationRule & ir_facet = GetSIMDIntegrationRule(etfacet, 2*maxorder + bonus_intorder);
            Facet2ElementTrafo transform1(eltype1, ElVertices);
            Facet2SurfaceElementTrafo stransform(strafo.GetElementType(), SElVertices); 
            auto & ir_facet_vol1 = transform1(LocalFacetNr, ir_facet, lh);
            auto & ir_facet_surf = stransform(ir_facet, lh);
            auto & mir1 = trafo1(ir_facet_vol1, lh);
            auto & smir = strafo(ir_facet_surf, lh);

            mir1.SetOtherMIR(&smir);
            smir.SetOtherMIR(&mir1);

            mir1.ComputeNormalsAndMeasure(eltype1, LocalFacetNr);
            
            ProxyUserData ud(trial_proxies.Size(), lh);
            const_cast<ElementTransformation&>(trafo1).userdata = &ud;
            const_cast<ElementTransformation&>(strafo).userdata = &ud;
            ud.fel = &fel1_trial;   // necessary to check remember-map
            // ud.elx = &elx;
            // ud.lh = &lh;

            PrecomputeCacheCF(cache_cfs, mir1, lh);

            for (ProxyFunction * proxy : trial_proxies)
              ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
            
            for (ProxyFunction * proxy : trial_proxies)
              if (! (proxy->IsOther() && proxy->BoundaryValues()) )
                proxy->Evaluator()->Apply(fel1_trial, mir1, elx, ud.GetAMemory(proxy));
            
            for (ProxyFunction * proxy : trial_proxies)
              if (proxy->IsOther() && proxy->BoundaryValues())
                proxy->BoundaryValues()->Evaluate (smir, ud.GetAMemory(proxy));
            
            // RegionTimer reg(t);
            
            
            for (auto proxy : test_proxies)
              {
                HeapReset hr(lh);
                FlatMatrix<SIMD<double>> proxyvalues(proxy->Dimension(), ir_facet.Size(), lh);
                
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    ud.testfunction = proxy;
                    ud.test_comp = k;
                    cf -> Evaluate (mir1, proxyvalues.Rows(k,k+1));            
                  }
                
                for (int i = 0; i < proxyvalues.Height(); i++)
                  {
                    auto row = proxyvalues.Row(i);
                    for (int j = 0; j < row.Size(); j++)
                      row(j) *= mir1[j].GetMeasure() * ir_facet[j].Weight();
                  }
                
                if (proxy->IsOther() && proxy->BoundaryValues())
                  ; // nothing to do 
                else
                  proxy->Evaluator()->AddTrans(fel1_test, mir1, proxyvalues, ely);
              }
            
          }
        catch (const ExceptionNOSIMD& e)
          {
            cout << IM(6) << "caught in SymbolicFacetInegtrator::ApplyBnd: " << endl
                 << e.What() << endl;
            simd_evaluate = false;
            ApplyFacetMatrix (fel1, LocalFacetNr, trafo1, ElVertices,
                              strafo, SElVertices, 
                              elx, ely, lh);
          }
        return;
        
      }

    
    static Timer t("SymbolicFacetBFI::ApplyFacetMatrix - boundary", NoTracing);
    
    HeapReset hr(lh);

    /*
    Matrix<> elmat(elx.Size());
    CalcFacetMatrix(fel1, LocalFacetNr, trafo1, ElVertices, strafo, elmat, lh);
    ely = elmat * elx;
    return;
    */

    ely = 0;
    
    FlatVector<> ely1(ely.Size(), lh);

    int maxorder = max2 (fel1_trial.Order(), fel1_test.Order());

    auto eltype1 = trafo1.GetElementType();
    auto etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr);

    const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, 2*maxorder + bonus_intorder);
    Facet2ElementTrafo transform1(eltype1, ElVertices);
    Facet2SurfaceElementTrafo stransform(strafo.GetElementType(), SElVertices); 
    
    IntegrationRule & ir_facet_vol1 = transform1(LocalFacetNr, ir_facet, lh);
    IntegrationRule & ir_facet_surf = stransform(ir_facet, lh);
    BaseMappedIntegrationRule & mir1 = trafo1(ir_facet_vol1, lh);
    auto & smir = strafo(ir_facet_surf, lh);

    mir1.SetOtherMIR (&smir);
    smir.SetOtherMIR (&mir1);
    
    // evaluate proxy-values
    ProxyUserData ud(trial_proxies.Size(), lh);
    const_cast<ElementTransformation&>(trafo1).userdata = &ud;
    const_cast<ElementTransformation&>(strafo).userdata = &ud;
    ud.fel = &fel1_trial;   // necessary to check remember-map
    // ud.elx = &elx;
    // ud.lh = &lh;

    PrecomputeCacheCF(cache_cfs, mir1, lh);

    for (ProxyFunction * proxy : trial_proxies)
      ud.AssignMemory (proxy, ir_facet.Size(), proxy->Dimension(), lh);
    
    for (ProxyFunction * proxy : trial_proxies)
      if (! (proxy->IsOther() && proxy->BoundaryValues()))
        proxy->Evaluator()->Apply(fel1_trial, mir1, elx, ud.GetMemory(proxy), lh);

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
        
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir1, val);
            proxyvalues.Col(k) = val.Col(0);
          }

        for (size_t i = 0; i < mir1.Size(); i++)
          proxyvalues.Row(i) *= mir1[i].GetMeasure() * ir_facet[i].Weight();

        ely1 = 0.0;
        if (proxy->IsOther() && proxy->BoundaryValues())
          ;  // nothing to do 
        else
          proxy->Evaluator()->ApplyTrans(fel1_trial, mir1, proxyvalues, ely1, lh);
        ely += ely1;
      }
  }
  
  
  SymbolicEnergy :: SymbolicEnergy (shared_ptr<CoefficientFunction> acf,
                                    VorB avb, VorB aelement_vb)
    : cf(acf), vb(avb), element_vb(aelement_vb)
  {
    simd_evaluate = true;
    // if (element_boundary) simd_evaluate = false;
    
    if (cf->Dimension() != 1)
      throw Exception ("SymbolicEnergy needs scalar-valued CoefficientFunction");
    
    trial_cum.Append(0);
    cf->TraverseTree
      ( [&] (CoefficientFunction & nodecf)
        {
          auto proxy = dynamic_cast<ProxyFunction*> (&nodecf);
          if (proxy) 
            {
              if (!proxy->IsTestFunction())
                {                                         
                  if (!trial_proxies.Contains(proxy))
                    {
                      trial_proxies.Append (proxy);
                      trial_cum.Append(trial_cum.Last()+proxy->Dimension());
                    }             
                }
            }
        });

    for (auto proxy : trial_proxies)
      if (!proxy->Evaluator()->SupportsVB(vb))
        throw Exception ("Trialfunction does not support "+ToString(vb)+"-forms, maybe a Trace() operator is missing, type = "+proxy->Evaluator()->Name());

    nonzeros = Matrix<bool>(trial_cum.Last(), trial_cum.Last());
    nonzeros_proxies = Matrix<bool>(trial_proxies.Size(), trial_proxies.Size());
    nonzeros_proxies = false;

    ProxyUserData ud;
    DummyFE<ET_TRIG> dummyfe;
    ud.fel = &dummyfe;
    Vector<AutoDiffDiff<1,NonZero>> nzvec(1);
    int k = 0;
    for (int k1 : trial_proxies.Range())
      for (int k2 : Range(trial_proxies[k1]->Dimension()))
        {
          int l = 0;
          for (int l1 : trial_proxies.Range())
            for (int l2 : Range(trial_proxies[l1]->Dimension()))
              {
                ud.trialfunction = trial_proxies[l1];
                ud.trial_comp = l2;
                ud.testfunction = trial_proxies[k1];
                ud.test_comp = k2;
                cf -> NonZeroPattern (ud, nzvec);
                // nzddvec(0) = true;
                // nonzeros(k,l) = nzddvec(0);
                nonzeros(k,l) = nzvec(0).DDValue(0);
                if (nzvec(0).DDValue(0))
                  nonzeros_proxies(k1,l1) = true;
                l++;
              }
          k++;
        }
    int cnt = 0;
    for (auto i : Range(nonzeros.Height()))
      {
        for (auto j : Range(nonzeros.Width()))
          {
            if (nonzeros(i,j)) cnt++;
            cout << IM(6) << (nonzeros(i,j) ? "1" : "0");
          }
        cout << IM(6) << endl;
      }



    
    dcf.SetSize(trial_proxies.Size());
    ddcf.SetSize(sqr(trial_proxies.Size()));    
    if (symbolic_integrator_uses_diff)
      {
        try
          {
            cout << IM(5) << "cf = " << *cf << endl;
            for (int i = 0; i < trial_proxies.Size(); i++)
              {
                CoefficientFunction::T_DJC cache;                
                auto diffi = cf->DiffJacobi(trial_proxies[i], cache);
                dcf[i] = diffi;
                cout << IM(5) << "diffi = " << *diffi << endl;
                /*
                  // compile time too long, at the moment
                for (int j = 0; j < trial_proxies.Size(); j++)
                  {
                    // cout << "diff_" << i << "," << j << " = " << endl;
                    CoefficientFunction::T_DJC cache;                                    
                    ddcf[i*trial_proxies.Size()+j] = diffi->DiffJacobi(trial_proxies[j], cache);
                    cout << IM(5) <<  "ddcf = " << *ddcf[i*trial_proxies.Size()+j] << endl;
                  }
                */
              }
          }
        catch (const Exception& e)
          {
            cout << IM(5) << "dcf_dtest has thrown exception " << e.What() << endl;
          }
      }
    
    cout << IM(6) << "nonzero: " << cnt << "/" << sqr(nonzeros.Height()) << endl;
    cout << IM(6) << "nonzero-proxies: " << endl << nonzeros_proxies << endl;
  }


  const IntegrationRule& SymbolicEnergy ::
  GetIntegrationRule (const FiniteElement & fel, LocalHeap & /* lh */) const
  {
    if (userdefined_intrules[fel.ElementType()]) return *userdefined_intrules[fel.ElementType()];

    int trial_difforder = 99; 
    for (auto proxy : trial_proxies)
      trial_difforder = min2(trial_difforder, proxy->Evaluator()->DiffOrder());
    if (trial_proxies.Size() == 0) trial_difforder = 0;
    
    int intorder = 2*fel.Order()+bonus_intorder;
    auto et = fel.ElementType();
    if (et == ET_TRIG || et == ET_TET)
      intorder -= 2*trial_difforder;
    return SelectIntegrationRule (et, intorder);
  }

  const SIMD_IntegrationRule& SymbolicEnergy ::
  Get_SIMD_IntegrationRule (const FiniteElement & fel, LocalHeap & lh) const
  {
    if (userdefined_simd_intrules[fel.ElementType()] ) return *userdefined_simd_intrules[fel.ElementType()];

    int trial_difforder = 99; 
    for (auto proxy : trial_proxies)
      trial_difforder = min2(trial_difforder, proxy->Evaluator()->DiffOrder());
    if (trial_proxies.Size() == 0) trial_difforder = 0;
    
    int intorder = 2*fel.Order()+bonus_intorder;
    auto et = fel.ElementType();
    if (et == ET_TRIG || et == ET_TET)
      intorder -= 2*trial_difforder;
    return SIMD_SelectIntegrationRule (et, intorder);
  }


  void SymbolicEnergy :: 
  CalcElementMatrix (const FiniteElement & fel,
                     const ElementTransformation & trafo, 
                     FlatMatrix<double> elmat,
                     LocalHeap & lh) const
  {
    FlatVector<double> elveclin(elmat.Height(), lh);
    elveclin = 0.0;
    CalcLinearizedElementMatrix (fel, trafo, elveclin, elmat, lh);
  }
  

  void 
  SymbolicEnergy :: CalcLinearizedElementMatrix (const FiniteElement & fel,
                                                 const ElementTransformation & trafo, 
                                                 FlatVector<double> elveclin,
                                                 FlatMatrix<double> elmat,
                                                 LocalHeap & lh) const
  {
    auto save_userdata = trafo.PushUserData();
    
    RegionTimer reg(timer);

    if (simd_evaluate)
      {
        try
          {
            elmat = 0;
            if (element_vb == VOL)
              {
                const SIMD_IntegrationRule& ir = Get_SIMD_IntegrationRule(fel, lh);
                auto & mir = trafo(ir, lh);
                AddLinearizedElementMatrix (fel, trafo, mir, elveclin, elmat, lh);
              }
            else
              {
                auto eltype = trafo.GetElementType();
        
                Facet2ElementTrafo transform(eltype, element_vb); 
                int nfacet = transform.GetNFacets();
                
                for (int k = 0; k < nfacet; k++)
                  {
                    HeapReset hr(lh);
                    ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
                    
                    auto & ir_facet = GetSIMDIntegrationRule(etfacet, 2*fel.Order()+bonus_intorder);
                    auto & ir_facet_vol = transform(k, ir_facet, lh);
                    auto & mir = trafo(ir_facet_vol, lh);
                    // mir.ComputeNormalsAndMeasure (eltype, k);
                    AddLinearizedElementMatrix (fel, trafo, mir, elveclin, elmat, lh);            
                  }
              }
          }

        catch (const ExceptionNOSIMD& e)
          {
            cout << IM(6) << e.What() << endl
                 << "switching back to standard evaluation (in SymbolicEnergy::CalcLinearized)" << endl;
            simd_evaluate = false;
            CalcLinearizedElementMatrix (fel, trafo, elveclin, elmat, lh);
          }
        return;
      }



    if (element_vb == VOL)
      {
        // const IntegrationRule& ir = GetIntegrationRule(trafo.GetElementType(), 2*fel.Order());
        const IntegrationRule& ir = GetIntegrationRule(fel, lh);
        BaseMappedIntegrationRule & mir = trafo(ir, lh);
        
        ProxyUserData ud(trial_proxies.Size(), lh);
        const_cast<ElementTransformation&>(trafo).userdata = &ud;
        ud.fel = &fel;
        // ud.elx = &elveclin;
        // ud.lh = &lh;
        for (ProxyFunction * proxy : trial_proxies)
          {
            ud.AssignMemory (proxy, ir.Size(), proxy->Dimension(), lh);
            proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetMemory(proxy), lh);        
          }
        
        elmat = 0;
        
        AddLinearizedElementMatrix (fel, ud, mir, elveclin, elmat, lh);
      }
    else // element_boundary
      {
        elmat = 0;
        auto eltype = trafo.GetElementType();

        Facet2ElementTrafo transform(eltype, element_vb); 
        int nfacet = transform.GetNFacets();
        
        for (int k = 0; k < nfacet; k++)
          {
            HeapReset hr(lh);
            ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
            
            const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, 2*fel.Order()+bonus_intorder);
            IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
            BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
            // mir.ComputeNormalsAndMeasure (eltype, k);
            
            ProxyUserData ud(trial_proxies.Size(), lh);    
            const_cast<ElementTransformation&>(trafo).userdata = &ud;
            ud.fel = &fel;
            // ud.elx = &elveclin;
            // ud.lh = &lh;
            
            for (ProxyFunction * proxy : trial_proxies)
              {
                ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
                proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetMemory(proxy), lh);
              }

            AddLinearizedElementMatrix (fel, ud, mir, elveclin, elmat, lh);            
          }
      }
  }

  void SymbolicEnergy ::
  AddLinearizedElementMatrix (const FiniteElement & fel,
                              ProxyUserData & ud,
                              const BaseMappedIntegrationRule & mir, 
                              FlatVector<double> elveclin,
                              FlatMatrix<double> elmat,
                              LocalHeap & lh) const
  {
    static Timer t("SymbolicEnergy::AddLinearizedElementMatrix (nosimd)", NoTracing);
    static Timer tdmat("SymbolicEnergy::CalcDMat (nosimd)", NoTracing);
    // Timer & tdmat = const_cast<Timer&> (timer);
    static Timer tbmat("SymbolicEnergy::CalcBMat (nosimd)", NoTracing);
    static Timer tmult("SymbolicEnergy::mult", NoTracing);
    RegionTimer reg(t);

    HeapReset hr(lh);
    
    FlatMatrix<> dderiv(mir.Size(), 1,lh);
    FlatMatrix<AutoDiffDiff<1,double>> ddval(mir.Size(), 1, lh);
    
    FlatArray<FlatMatrix<>> diags(trial_proxies.Size(), lh);
    FlatArray<FlatMatrix<>> dWdB(trial_proxies.Size(), lh);
    for (int k1 : Range(trial_proxies))
      {
        auto proxy = trial_proxies[k1];
        diags[k1].AssignMemory(mir.Size(), proxy->Dimension(), lh);
        dWdB[k1].AssignMemory(mir.Size(), proxy->Dimension(), lh);
        for (int k = 0; k < proxy->Dimension(); k++)
          {
            ud.trialfunction = proxy;
            ud.trial_comp = k;
            ud.testfunction = proxy;
            ud.test_comp = k;
            cf -> Evaluate (mir, ddval);
            for (size_t i = 0; i < mir.Size(); i++)
              diags[k1](i,k) = ddval(i,0).DDValue(0);
            for (size_t i = 0; i < mir.Size(); i++)
              dWdB[k1](i,k) = ddval(i,0).DValue(0);
          }
      }

    for (int k1 : Range(trial_proxies))
      for (int l1 : Range(trial_proxies))
        {
          HeapReset hr(lh);
          if (!nonzeros_proxies(k1,l1)) continue;
          
          auto proxy1 = trial_proxies[k1];
          auto proxy2 = trial_proxies[l1];

          FlatTensor<3> proxyvalues(lh, mir.Size(), proxy2->Dimension(), proxy1->Dimension());
          for (int k = 0; k < proxy1->Dimension(); k++)
            for (int l = 0; l < proxy2->Dimension(); l++)
              {
                ud.trialfunction = proxy1;
                ud.trial_comp = k;
                ud.testfunction = proxy2;
                ud.test_comp = l;

                {
                  RegionTimer reg(tdmat);
                  NgProfiler::AddThreadFlops (tdmat, TaskManager::GetThreadId(), 1);
                  if (nonzeros(trial_cum[k1]+k, trial_cum[l1]+l))
                    {
                      cf -> Evaluate (mir, ddval);
                      for (size_t i = 0; i < mir.Size(); i++)
                        dderiv(i,0) = ddval(i,0).DDValue(0);
                    }
                  else
                    dderiv = 0.0;
                }
                proxyvalues(STAR,l,k) = dderiv.Col(0);

                if (proxy1 != proxy2 || k != l)  // computed mixed second derivatives
                  {
                    proxyvalues(STAR,l,k) -= diags[k1].Col(k);
                    proxyvalues(STAR,l,k) -= diags[l1].Col(l);
                    proxyvalues(STAR,l,k) *= 0.5;
                  }
              }

          for (int i = 0; i < mir.Size(); i++)
            proxyvalues(i,STAR,STAR) *= mir[i].GetWeight();

          FlatMatrix<double,ColMajor> bmat1(proxy1->Dimension(), elmat.Width(), lh);
          FlatMatrix<double,ColMajor> bmat2(proxy2->Dimension(), elmat.Height(), lh);

          IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
          IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);
          SliceMatrix<> part_elmat = elmat.Rows(r2).Cols(r1);

          constexpr size_t BS = 16;
          size_t i = 0;
          for ( ; i < mir.Size(); i+=BS)
            {
              HeapReset hr(lh);
              int bs = min2(size_t(BS), mir.IR().GetNIP()-i);
              FlatMatrix<double,ColMajor> bbmat1(bs*proxy1->Dimension(), elmat.Width(), lh);              
              FlatMatrix<double,ColMajor> bdbmat1(bs*proxy2->Dimension(), elmat.Width(), lh);
              FlatMatrix<double,ColMajor> bbmat2(bs*proxy2->Dimension(), elmat.Height(), lh);

              proxy1->Evaluator()->CalcLinearizedMatrix(fel, mir.Range(i,i+bs,lh), elveclin, bbmat1, lh);
              
              for (size_t j = 0; j < bs; j++)
                {
                  IntRange r3 = proxy2->Dimension() * IntRange(j,j+1);
                  IntRange rb1 = proxy1->Dimension() * IntRange(j,j+1);                  
                  bdbmat1.Rows(r3).Cols(r1) = proxyvalues(i+j,STAR,STAR) * bbmat1.Rows(rb1).Cols(r1);                  
                }
              
              // elmat += Trans (bbmat2) * bdbmat1 | Lapack;
              // AddABt (Trans(bbmat2), Trans(bdbmat1), elmat);
              RegionTimer reg(tmult);                                
              // AddABt (Trans(bbmat2).Rows(r2), Trans(bdbmat1).Rows(r1), part_elmat);
              // part_elmat += Trans(bbmat2).Rows(r2) * bdbmat1.Cols(r1);

              if (k1 == l1)
                part_elmat += Trans(bbmat1).Rows(r1) * bdbmat1.Cols(r1);
              else
                {
                  proxy2->Evaluator()->CalcLinearizedMatrix(fel, mir.Range(i,i+bs,lh), elveclin, bbmat2, lh);
                  part_elmat += Trans(bbmat2).Rows(r2) * bdbmat1.Cols(r1);                  
                }

            }
        }



    // contributions due to non-linear diffops:
    // for (auto & proxy : trial_proxies)
    for (int k : Range(trial_proxies))
      {
        auto & proxy = trial_proxies[k];
        if (proxy->Evaluator()->IsNonlinear())
          {
            //cout << "have a nonlienar proxy" << endl;
            //cout << "dWdB = " << dWdB[k] << endl;
            //cout << "elmat h/w =" << elmat.Height() << " | " << elmat.Width() << endl;
            proxy->Evaluator()->CalcHessianAdd (fel, mir, dWdB[k], elveclin, elmat, lh);
            //cout << "elmat = " << endl << elmat << endl;
          }
      }
  }
  
  void SymbolicEnergy :: AddLinearizedElementMatrix (const FiniteElement & fel,
                                                     const ElementTransformation & trafo, 
                                                     // ProxyUserData & ud, 
                                                     const SIMD_BaseMappedIntegrationRule & mir, 
                                                     FlatVector<double> elveclin,
                                                     FlatMatrix<double> elmat,
                                                     LocalHeap & lh) const
  {
    static Timer t("SymbolicEnergy::AddLinearizedElementMatrix - simd", NoTracing);
    static Timer tdmat("SymbolicEnergy::CalcDMat - simd", NoTracing);
    static Timer tdmat2("SymbolicEnergy::CalcDMat2 - simd", NoTracing);
    static Timer tbmat("SymbolicEnergy::CalcBMat - simd", NoTracing);
    static Timer tmult("SymbolicEnergy::mult - simd", NoTracing);
    
    RegionTimer reg(t);

    ProxyUserData ud(trial_proxies.Size(), lh);    
    const_cast<ElementTransformation&>(trafo).userdata = &ud;
    ud.fel = &fel;
    for (ProxyFunction * proxy : trial_proxies)
      {
        ud.AssignMemory (proxy, mir.IR().GetNIP(), proxy->Dimension(), lh);
        proxy->Evaluator()->Apply(fel, mir, elveclin, ud.GetAMemory(proxy));
      }



    
            FlatMatrix<AutoDiffDiff<1,SIMD<double>>> ddval(1, mir.Size(), lh);
            FlatArray<FlatMatrix<SIMD<double>>> diags(trial_proxies.Size(), lh);
            for (int k1 : Range(trial_proxies))
              {
                auto proxy = trial_proxies[k1];
                new(&diags[k1]) FlatMatrix<SIMD<double>>(proxy->Dimension(), mir.Size(), lh);
                if (nonzeros_proxies(k1,k1))
                  for (int k = 0; k < proxy->Dimension(); k++)
                    {
                      ud.trialfunction = proxy;
                      ud.trial_comp = k;
                      ud.testfunction = proxy;
                      ud.test_comp = k;
                      // cf -> EvaluateDDeriv (mir, AFlatMatrix<>(val), AFlatMatrix<>(deriv), AFlatMatrix<> (diags[k1].Rows(k,k+1)));
                      auto diagrow = diags[k1].Row(k);
                      cf -> Evaluate (mir, ddval);
                      for (auto i : Range(diagrow))
                        diagrow(i) = ddval(i).DDValue(0);
                    }
                else
                  diags[k1] = 0.0;
              }

            for (int k1 : Range(trial_proxies))
              for (int l1 : Range(trial_proxies))
                {
                  HeapReset hr(lh);
                  if (!nonzeros_proxies(k1,l1)) continue;
                  if (k1 < l1) continue;
                  auto proxy1 = trial_proxies[k1];
                  auto proxy2 = trial_proxies[l1];

                  size_t dim_proxy1 = proxy1->Dimension();
                  size_t dim_proxy2 = proxy2->Dimension();

                  FlatMatrix<SIMD<double>> proxyvalues2(dim_proxy1*dim_proxy2, mir.Size(), lh);

                  {
                  RegionTimer reg(tdmat);

                  /*
                  if (ddcf[k1*trial_proxies.Size()+l1])
                    {
                      ddcf[k1*trial_proxies.Size()+l1]->Evaluate(mir, proxyvalues2);
                    }
                  */
                  if (dcf[k1])
                    {
                      HeapReset hr(lh);
                      FlatMatrix<AutoDiff<1,SIMD<double>>> dval(dim_proxy1, mir.Size(), lh);
                      for (int l = 0; l < dim_proxy2; l++)
                        {
                          ud.trialfunction = proxy2;
                          ud.trial_comp = l;
                        
                          dcf[k1] -> Evaluate (mir, dval);
                          
                          for (int k = 0; k < dim_proxy1; k++)
                            for (auto i : Range(mir.Size()))
                              proxyvalues2(k*dim_proxy2+l,i) = dval(k,i).DValue(0);
                        }
                    }
                  else
                  
                  for (int k = 0; k < dim_proxy1; k++)
                    for (int l = 0; l < dim_proxy2; l++)
                      {
                        ud.trialfunction = proxy1;
                        ud.trial_comp = k;
                        ud.testfunction = proxy2;
                        ud.test_comp = l;
                        auto proxyrow = proxyvalues2.Row(k*dim_proxy2+l);

                        if (nonzeros(trial_cum[k1]+k, trial_cum[l1]+l))
                          {
                            // cf -> EvaluateDDeriv (mir, AFlatMatrix<>(val), AFlatMatrix<>(deriv), AFlatMatrix<> (dderiv));
                            // proxyrow = dderiv.Row(0);

                            if (k1 != l1 || k <= l)
                              {
                                cf -> Evaluate (mir, ddval);
                                for (auto i : Range(proxyrow))
                                  proxyrow(i) = ddval(i).DDValue(0);
                              }
                            if (k1 == l1 && k < l) // symmetric
                              proxyvalues2.Row(l*dim_proxy2+k) = proxyrow;
                          }
                        else
                          {
                            // dderiv = 0.0;
                            proxyrow = 0.0;
                          }
                        if (proxy1 != proxy2 || k != l)  // computed mixed second derivatives
                          {
                            proxyrow -= diags[k1].Row(k);
                            proxyrow -= diags[l1].Row(l);
                            proxyrow *= 0.5;
                          }
                      }

                  for (size_t i = 0; i < mir.Size(); i++)
                    proxyvalues2.Col(i) *= mir[i].GetWeight();
                  }
                  

                  IntRange r1 = proxy1->Evaluator()->UsedDofs(fel);
                  IntRange r2 = proxy2->Evaluator()->UsedDofs(fel);

                  FlatMatrix<SIMD<double>> bmat1(elmat.Width()*dim_proxy1, mir.Size(), lh);
                  FlatMatrix<SIMD<double>> dbmat1(elmat.Width()*dim_proxy2, mir.Size(), lh);
                  FlatMatrix<SIMD<double>> bmat2(elmat.Height()*dim_proxy2, mir.Size(), lh);

                  auto hbmat1 =  bmat1. Reshape(elmat.Width(), dim_proxy1*mir.Size());
                  auto hdbmat1 = dbmat1.Reshape(elmat.Width(), dim_proxy2*mir.Size());
                  auto hbmat2 =  bmat2. Reshape(elmat.Height(), dim_proxy2*mir.Size());
                  /*
                  FlatMatrix<SIMD<double>> hbmat1(elmat.Width(), dim_proxy1*mir.Size(), &bmat1(0,0));                  
                  FlatMatrix<SIMD<double>> hdbmat1(elmat.Width(), dim_proxy2*mir.Size(), &dbmat1(0,0));
                  FlatMatrix<SIMD<double>> hbmat2(elmat.Height(), dim_proxy2*mir.Size(), &bmat2(0,0));
                  */
                  {
                    RegionTimer reg(tbmat);
                    proxy1->Evaluator()->CalcMatrix(fel, mir, bmat1);
                  }

                  hdbmat1.Rows(r1) = 0.0; 
                  for (auto i : r1)
                    for (size_t j = 0; j < dim_proxy2; j++)
                      for (size_t k = 0; k < dim_proxy1; k++)
                        {
                          auto res = dbmat1.Row(i*dim_proxy2+j);
                          auto a = bmat1.Row(i*dim_proxy1+k);
                          auto b = proxyvalues2.Row(k*dim_proxy2+j);
                          res += pw_mult(a,b);
                        }

                  {
                    RegionTimer reg(tmult);
                    if (k1 == l1)
                      {
                        AddABt (hbmat1.Rows(r1), hdbmat1.Rows(r1), elmat.Rows(r1).Cols(r1));
                      }
                    else
                      {
                        proxy2->Evaluator()->CalcMatrix(fel, mir, bmat2);
                        AddABt (hbmat2.Rows(r2), hdbmat1.Rows(r1), elmat.Rows(r2).Cols(r1));
                        AddABt (hdbmat1.Rows(r1), hbmat2.Rows(r2), elmat.Rows(r1).Cols(r2));
                      }
                    // elmat.Rows(r1).Cols(r2) = Trans(part_elmat); // buggy if both evaluators act on same element
                  }
                }
  }

  
  double SymbolicEnergy :: Energy (const FiniteElement & fel, 
                                   const ElementTransformation & trafo, 
                                   FlatVector<double> elx, 
                                   LocalHeap & lh) const
  {
    if (simd_evaluate && (element_vb == VOL))
      {
        try
          {
            const SIMD_IntegrationRule& ir = Get_SIMD_IntegrationRule(fel, lh);            
            auto & mir = trafo(ir, lh);
        
            ProxyUserData ud(trial_proxies.Size(), lh);
            const_cast<ElementTransformation&>(trafo).userdata = &ud;
            ud.fel = &fel;
            
            for (ProxyFunction * proxy : trial_proxies)
              {
                ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);
                proxy->Evaluator()->Apply(fel, mir, elx, ud.GetAMemory(proxy));
              }
            
            FlatMatrix<SIMD<double>> values(1, mir.Size(), lh);
            cf -> Evaluate(mir, values);

            SIMD<double> sum = 0;
            for (size_t i = 0; i < mir.Size(); i++)
              sum += mir[i].GetWeight() * values(0, i);
            return HSum(sum);
          }
        catch (const ExceptionNOSIMD& e)
          {
            cout << IM(6) << e.What() << endl
                 << "switching back to standard evaluation (in SymbolicEnergy::Energy)" << endl;
            simd_evaluate = false;
            return Energy (fel, trafo, elx, lh);
          }
      }

    if (element_vb == VOL)
      {
        const IntegrationRule& ir = GetIntegrationRule(fel, lh);
        BaseMappedIntegrationRule & mir = trafo(ir, lh);
        
        ProxyUserData ud(trial_proxies.Size(), lh);
        const_cast<ElementTransformation&>(trafo).userdata = &ud;
        ud.fel = &fel;
        
        for (ProxyFunction * proxy : trial_proxies)
          {
            ud.AssignMemory (proxy, ir.Size(), proxy->Dimension(), lh);
            proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
          }
        FlatMatrix<> values(mir.Size(), 1, lh);
        cf -> Evaluate(mir, values);
        
        double sum = 0;
        for (size_t i = 0; i < mir.Size(); i++)
          sum += mir[i].GetWeight() * values(i,0);
        return sum;
      }
    else // element_boundary
      {
        double sum = 0;
        auto eltype = trafo.GetElementType();
        
        Facet2ElementTrafo transform(eltype, element_vb); 
        int nfacet = transform.GetNFacets();
        
        for (int k = 0; k < nfacet; k++)
          {
            HeapReset hr(lh);
            ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
            
            const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, 2*fel.Order()+bonus_intorder);
            IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
            BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
            // mir.ComputeNormalsAndMeasure (eltype, k);
            
            ProxyUserData ud(trial_proxies.Size(), lh);    
            const_cast<ElementTransformation&>(trafo).userdata = &ud;
            ud.fel = &fel;
            
            for (ProxyFunction * proxy : trial_proxies)
              {
                ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
                proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
              }

            FlatMatrix<> values(mir.Size(), 1, lh);
            cf -> Evaluate(mir, values);
            
            for (int i = 0; i < mir.Size(); i++)
              sum += mir[i].GetWeight() * values(i,0);
          }
        return sum;
      }
  }

  void
  SymbolicEnergy :: ApplyElementMatrix (const FiniteElement & fel, 
                                        const ElementTransformation & trafo, 
                                        const FlatVector<double> elx, 
                                        FlatVector<double> ely,
                                        void * precomputed,
                                        LocalHeap & lh) const
  {
    // static Timer t("SymbolicEnergy::ApplyElementMatrix", NoTracing); 
        
    if (simd_evaluate) //  && !element_boundary)
      {
        try
          {
            if (element_vb == VOL)
              {
                HeapReset hr(lh);
                
                ProxyUserData ud(trial_proxies.Size(), lh);        
                const_cast<ElementTransformation&>(trafo).userdata = &ud;
                ud.fel = &fel;

                const SIMD_IntegrationRule& ir = Get_SIMD_IntegrationRule(fel, lh);
                auto & mir = trafo(ir, lh);

                for (ProxyFunction * proxy : trial_proxies)
                  ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);
                for (ProxyFunction * proxy : trial_proxies)
                  proxy->Evaluator()->Apply(fel, mir, elx, ud.GetAMemory(proxy));
                
                ely = 0;
                for (auto proxy : trial_proxies)
                  {
                    HeapReset hr(lh);
                    FlatMatrix<SIMD<double>> proxyvalues(proxy->Dimension(), ir.Size(), lh);
                    FlatMatrix<AutoDiff<1,SIMD<double>>> vals(1, ir.Size(), lh);
                    for (int k = 0; k < proxy->Dimension(); k++)
                      {
                        ud.trialfunction = proxy;
                        ud.trial_comp = k;
                        cf -> Evaluate(mir, vals);
                        for (size_t j = 0; j < ir.Size(); j++)
                          proxyvalues(k,j) = vals(0,j).DValue(0);
                      }

                    for (size_t j = 0; j < ir.Size(); j++)
                      proxyvalues.Col(j) *= mir[j].GetWeight();
                    
                    proxy->Evaluator()->AddTrans(fel, mir, proxyvalues, ely);
                  }
              }
            else
              { // element_boundary

                ely = 0;
                auto eltype = trafo.GetElementType();
                
                Facet2ElementTrafo transform(eltype, element_vb); 
                int nfacet = transform.GetNFacets();
                
                for (int k = 0; k < nfacet; k++)
                  {
                    HeapReset hr(lh);
                    ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);

                    auto & ir_facet = GetSIMDIntegrationRule(etfacet, 2*fel.Order()+bonus_intorder);
                    auto & ir_facet_vol = transform(k, ir_facet, lh);
                    auto & mir = trafo(ir_facet_vol, lh);
                    // mir.ComputeNormalsAndMeasure (eltype, k);

                    ProxyUserData ud(trial_proxies.Size(), lh);    
                    const_cast<ElementTransformation&>(trafo).userdata = &ud;
                    ud.fel = &fel;

                    for (ProxyFunction * proxy : trial_proxies)
                      ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
                    for (ProxyFunction * proxy : trial_proxies)
                      proxy->Evaluator()->Apply(fel, mir, elx, ud.GetAMemory(proxy));
                
                    for (auto proxy : trial_proxies)
                      {
                        HeapReset hr(lh);
                        FlatMatrix<SIMD<double>> proxyvalues(proxy->Dimension(), mir.Size(), lh);
                        FlatMatrix<AutoDiff<1,SIMD<double>>> vals(1, mir.Size(), lh);
                        for (int k = 0; k < proxy->Dimension(); k++)
                          {
                            ud.trialfunction = proxy;
                            ud.trial_comp = k;
                            cf -> Evaluate(mir, vals);
                            for (size_t j = 0; j < mir.Size(); j++)
                              proxyvalues(k,j) = vals(0,j).DValue(0);
                          }
                        
                        for (size_t j = 0; j < mir.Size(); j++)
                          proxyvalues.Col(j) *= mir[j].GetWeight();
                        
                        proxy->Evaluator()->AddTrans(fel, mir, proxyvalues, ely);
                      }
                  }
              }
          }
        catch (const ExceptionNOSIMD& e)
          {
            cout << IM(6) << e.What() << endl
                 << "switching back to standard evaluation (in SymbolicEnergy::CalcLinearized)" << endl;              
            simd_evaluate = false;
            ApplyElementMatrix (fel, trafo, elx, ely, precomputed, lh);
          }
        return;
      }

    
    if (element_vb == VOL)
      {
        HeapReset hr(lh);
        
        ProxyUserData ud(trial_proxies.Size(), lh);        
        const_cast<ElementTransformation&>(trafo).userdata = &ud;
        ud.fel = &fel;

        const IntegrationRule& ir = GetIntegrationRule(fel, lh);
        BaseMappedIntegrationRule & mir = trafo(ir, lh);

        for (ProxyFunction * proxy : trial_proxies)
          ud.AssignMemory (proxy, ir.GetNIP(), proxy->Dimension(), lh);
        
        for (ProxyFunction * proxy : trial_proxies)
          proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
        
        ely = 0;
        FlatVector<> ely1(ely.Size(), lh);
        FlatMatrix<AutoDiff<1>> dval(mir.Size(), 1, lh);
        
        for (auto proxy : trial_proxies)
          {
            HeapReset hr(lh);
            FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
            for (int k = 0; k < proxy->Dimension(); k++)
              {
                ud.trialfunction = proxy;
                ud.trial_comp = k;
                cf -> Evaluate (mir, dval);
                for (size_t i = 0; i < mir.Size(); i++)
                  proxyvalues(i,k) = dval(i,0).DValue(0);
              }
            
            for (int i = 0; i < mir.Size(); i++)
              proxyvalues.Row(i) *= mir[i].GetWeight();
            
            proxy->Evaluator()->ApplyLinearizedTrans(fel, mir, elx, proxyvalues, ely1, lh);
            ely += ely1;
          }
      }
    else  // element_boundary
      {
        ely = 0;
        auto eltype = trafo.GetElementType();

        Facet2ElementTrafo transform(eltype, element_vb); 
        int nfacet = transform.GetNFacets();
        
        for (int k = 0; k < nfacet; k++)
          {
            HeapReset hr(lh);
            ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
            
            const IntegrationRule & ir_facet = GetIntegrationRule(etfacet, 2*fel.Order()+bonus_intorder);
            IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
            BaseMappedIntegrationRule & mir = trafo(ir_facet_vol, lh);
            // mir.ComputeNormalsAndMeasure (eltype, k);
            
            
            ProxyUserData ud(trial_proxies.Size(), lh);    
            const_cast<ElementTransformation&>(trafo).userdata = &ud;
            ud.fel = &fel;
            
            for (ProxyFunction * proxy : trial_proxies)
              {
                ud.AssignMemory (proxy, ir_facet.GetNIP(), proxy->Dimension(), lh);
                proxy->Evaluator()->Apply(fel, mir, elx, ud.GetMemory(proxy), lh);
              }
        
            FlatVector<> ely1(ely.Size(), lh);
            // FlatMatrix<> val(mir.Size(), 1,lh), deriv(mir.Size(), 1,lh);
            FlatMatrix<AutoDiff<1>> dval(mir.Size(), 1, lh);
            
            for (auto proxy : trial_proxies)
              {
                HeapReset hr(lh);
                FlatMatrix<> proxyvalues(mir.Size(), proxy->Dimension(), lh);
                for (int k = 0; k < proxy->Dimension(); k++)
                  {
                    ud.trialfunction = proxy;
                    ud.trial_comp = k;
                    // cf -> EvaluateDeriv (mir, val, deriv);
                    // proxyvalues.Col(k) = deriv.Col(0);
                    cf -> Evaluate (mir, dval);
                    for (size_t i = 0; i < mir.Size(); i++)
                      proxyvalues(i,k) = dval(i,0).DValue(0);
                  }
                
                for (int i = 0; i < mir.Size(); i++)
                  proxyvalues.Row(i) *= mir[i].GetWeight();
                
                proxy->Evaluator()->ApplyTrans(fel, mir, proxyvalues, ely1, lh);
                ely += ely1;
              }
          }
      }
  }


  shared_ptr<BilinearFormIntegrator> Integral :: MakeBilinearFormIntegrator() const
  {
    // check for DG terms
    bool has_other = false;
    cf->TraverseTree ([&has_other] (CoefficientFunction & cf)
                      {
                        if (dynamic_cast<ProxyFunction*> (&cf))
                          if (dynamic_cast<ProxyFunction&> (cf).IsOther())
                            has_other = true;
                      });
    if (has_other && (dx.element_vb != BND) && !dx.skeleton)
      throw Exception("DG-facet terms need either skeleton=True or element_boundary=True");

    shared_ptr<BilinearFormIntegrator> bfi;
    if (!has_other && !dx.skeleton)
      bfi = make_shared<SymbolicBilinearFormIntegrator> (cf, dx.vb, dx.element_vb);
    else
      bfi = make_shared<SymbolicFacetBilinearFormIntegrator> (cf, dx.vb, !dx.skeleton);
    if (dx.definedon)
      {
        if (auto definedon_bitarray = get_if<BitArray> (&*dx.definedon); definedon_bitarray)
          bfi->SetDefinedOn(*definedon_bitarray);
        /*
          // can't do that withouyt mesh
        if (auto definedon_string = get_if<string> (&*dx.definedon); definedon_string)
          {
            Region reg(self.GetFESpace()->GetMeshAccess(), dx.vb, *definedon_string);
            bfi->SetDefinedOn(reg.Mask());
          }
        */
      }
    bfi->SetDeformation(dx.deformation);               
    bfi->SetBonusIntegrationOrder(dx.bonus_intorder);
    if(dx.definedonelements)
      bfi->SetDefinedOnElements(dx.definedonelements);
    for (auto both : dx.userdefined_intrules)
      bfi->SetIntegrationRule(both.first, *both.second);

    if(linearization)
      dynamic_pointer_cast<SymbolicBilinearFormIntegrator>(bfi)->SetLinearization(linearization->MakeBilinearFormIntegrator());

    return bfi;
  }

  shared_ptr<LinearFormIntegrator> Integral :: MakeLinearFormIntegrator() const
  {
    cf -> TraverseTree ([&] (CoefficientFunction& cf) {
                          if (auto * proxy = dynamic_cast<ProxyFunction*>(&cf))
                            {
                              if (proxy->IsTrialFunction())
                                throw Exception("In MakeLinearFormIntegrator: must not have TrialFunction");
                            }
                        });
    shared_ptr<LinearFormIntegrator> lfi;
    if (!dx.skeleton)
      lfi =  make_shared<SymbolicLinearFormIntegrator> (cf, dx.vb, dx.element_vb);
    else
      lfi = make_shared<SymbolicFacetLinearFormIntegrator> (cf, dx.vb);
    if (dx.definedon)
      {
        if (auto definedon_bitarray = get_if<BitArray> (&*dx.definedon); definedon_bitarray)
          lfi->SetDefinedOn(*definedon_bitarray);
        /*
        // can't do that withouyt mesh
        if (auto definedon_string = get_if<string> (&*dx.definedon); definedon_string)
        {
          Region reg(self->GetFESpace()->GetMeshAccess(), dx.vb, *definedon_string);
          lfi->SetDefinedOn(reg.Mask());
        }
        */
      }
    lfi->SetDeformation(dx.deformation);
    lfi->SetBonusIntegrationOrder(dx.bonus_intorder);
    if(dx.definedonelements)
      lfi->SetDefinedOnElements(dx.definedonelements);
    for (auto both : dx.userdefined_intrules)
      lfi->SetIntegrationRule(both.first, *both.second);
    return lfi;
  }


}
