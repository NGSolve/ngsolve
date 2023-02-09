/*********************************************************************/
/* File:   coefficient.cpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   24. Jan. 2002                                             */
/*********************************************************************/

/* 
   Finite Element Coefficient Function
*/

#include <core/register_archive.hpp>
#include <fem.hpp>
#include <coefficient_impl.hpp>

#include <../ngstd/evalfunc.hpp>
#include <algorithm>
#ifdef NGS_PYTHON
#include <core/python_ngcore.hpp> // for shallow archive
#endif // NGS_PYTHON


namespace ngfem
{
  using ngstd::IfPos;
  
  CoefficientFunction :: ~CoefficientFunction ()
  { ; }


  void CoefficientFunction :: GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    string mycode =
      string("// GenerateCode() not overloaded for: ") + Demangle(typeid(*this).name()) + "\n"
      + R"CODE_(    typedef {scal_type} TStack{index};
    STACK_ARRAY(TStack{index}, hmem{index}, mir.Size()*{dim});
    {values_type} {values}({rows}, {cols}, reinterpret_cast<{scal_type}*>(&hmem{index}[0]));
    {
      const CoefficientFunction & cf = *reinterpret_cast<CoefficientFunction*>({this});
      {values} = {scal_type}(0.0);
      cf.Evaluate(mir, {values});
    }
    )CODE_";
    auto values = Var("values", index);
    string scal_type = code.res_type;
    string rows = ToString(Dimension());
    string cols = "mir.IR().Size()";

    std::map<string,string> variables;
    variables["scal_type"] = scal_type;
    variables["values_type"] = "FlatMatrix<"+scal_type+">";
    variables["values"] = values.S();
    variables["this"] =  code.AddPointer(this);
    variables["dim"] = ToString(Dimension());
    variables["index"] = ToString(index);
    variables["rows"] = code.is_simd ? rows : cols;
    variables["cols"] = code.is_simd ? cols : rows;
    code.header += Code::Map(mycode, variables);

    // code.Declare(code.res_type, index, Dimensions());
    code.Declare(index, Dimensions(), IsComplex());
    if(code.is_simd)
      {
        for (int i = 0; i < Dimension(); i++)
          code.body += Var(index,i,Dimensions()).Assign(values.S()+"("+ToString(i)+",i)", false);          
      }
    else
      {
        for (int i = 0; i < Dimension(); i++)
          code.body += Var(index,i,Dimensions()).Assign(values.S()+"(i,"+ToString(i)+")", false);
      }
  }

  void CoefficientFunction :: PrintReport (ostream & ost) const
  {
    // ost << "CoefficientFunction is " << typeid(*this).name() << endl;
    PrintReportRec (ost, 0);
  }
  
  void CoefficientFunction :: PrintReportRec (ostream & ost, int level) const
  {
    ost << string(2*level, ' ');
    // for (int i = 0; i < 2*level; i++)
    // ost << ' ';
    ost << "coef " << GetDescription() << ","
        << (IsComplex() ? " complex" : " real");
    if (Dimensions().Size() == 1)
      ost << ", dim=" << Dimension();
    else if (Dimensions().Size() >= 2)
      {
        ost << ", dims = " << Dimensions()[0]; //  << " x " << Dimensions()[1];
        for (int j = 1; j < Dimensions().Size(); j++)
          ost << " x " << Dimensions()[j];
      }
    ost << endl;

    Array<shared_ptr<CoefficientFunction>> input = InputCoefficientFunctions();
    for (int i = 0; i < input.Size(); i++)
      if (input[i])
        input[i] -> PrintReportRec (ost, level+1);
      else
        ost << string(2*level+2, ' ') << "none" << endl;
  }
  
  string CoefficientFunction :: GetDescription () const
  {
    if (description.length()) return description;
    return typeid(*this).name();
  }    

  shared_ptr<CoefficientFunction> CoefficientFunction :: Reshape (FlatArray<int> adims) const
  {
    return ReshapeCF (const_cast<CoefficientFunction*>(this)->shared_from_this(), adims);
  }
  
  shared_ptr<CoefficientFunction> ReshapeCF (shared_ptr<CoefficientFunction> coef,
                                             FlatArray<int> adims)
  {
    if (coef->Dimensions() == adims)
      return coef;

    if (coef->IsZeroCF())
      return ZeroCF(adims);

    Array<int> newdims (adims);

    int newdim = 1;
    int cntm1 = 0;

    for (int d : newdims) newdim *= d;
    for (int d : newdims) if (d == -1) cntm1++;

    if (cntm1 > 1)
      throw Exception("onlye one -1 allowed in reshape");

    if (cntm1 == 1)
      {
        newdim = -newdim;
        if (coef->Dimension() % newdim != 0)
          throw Exception("Reshape with -1 not a divider");

        for (int & d : newdims)
          if (d == -1)
            d = coef->Dimension() / newdim;
      }
    else                          
      {
        if (newdim != coef->Dimension())
          throw Exception("Reshape dim "+ToString(newdim)+
                          " does not fit current dimension "+ToString(coef->Dimension()));
      }
            
    if (coef -> GetDescription() == "reshape") //  && coef.use_count() == 1)
      {
        // cout << "double reshape found" << endl;
        auto inputs = coef->InputCoefficientFunctions();
        coef = nullptr;
        return ReshapeCF(inputs[0], newdims);
      }

    auto wrapper = CreateWrapperCF(coef);
    wrapper->SetDimensions(newdims);
    wrapper->SetDescription("reshape");
    return wrapper;
  }

  shared_ptr<CoefficientFunction> CoefficientFunction :: Transpose () const
  {
    auto sp = const_cast<CoefficientFunction*>(this)->shared_from_this();
    return TransposeCF (sp);
  }
  
  shared_ptr<CoefficientFunction> CoefficientFunction :: TensorTranspose (int i, int j) const
  {
    auto sp = const_cast<CoefficientFunction*>(this)->shared_from_this();
    return MakeTensorTransposeCoefficientFunction (sp, i, j);
  }
  
  shared_ptr<CoefficientFunction> CoefficientFunction ::
  Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const
  {
    throw Exception(string("Diff not implemented for CF ")+typeid(*this).name());
  }

  shared_ptr<CoefficientFunction> CoefficientFunction ::
  DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];
        
    /*
    if (var->Dimensions().Size() == 0)
      return this->Diff(var, make_shared<ConstantCoefficientFunction>(1));
    else
    */
      {
        Array<int> resultdims;
        resultdims += this->Dimensions();
        resultdims += var->Dimensions();
        
        if (this == var) // diff by me
          return IdentityCF(Dimensions());

        if (this->InputCoefficientFunctions().Size()==0)
          return ZeroCF(resultdims);          
        
        cout << IM(5) << "DiffJacobi for CoefficientFunction, type = " << typeid(*this).name() << endl;

        int dim = var->Dimension();
        Array<shared_ptr<CoefficientFunction>> ddi(dim); 
        for (int i = 0; i < dim; i++)
          {
            auto vec = UnitVectorCF(dim,i)->Reshape(var->Dimensions());
            ddi[i] = this->Diff(var, vec);
          }
        auto dvec = MakeVectorialCoefficientFunction (move(ddi));
        auto dvec1 = dvec->Reshape(var->Dimension(), this->Dimension()) -> Transpose();
        auto res = dvec1 -> Reshape(resultdims);
        cache[thisptr] = res;
        return res;
      }
  }

  
  shared_ptr<CoefficientFunction> CoefficientFunction ::
  Operator (const string & name) const
  {
    throw Exception(string("Operator ") + name + string(" not overloaded for CF ")+typeid(*this).name());
  }

  shared_ptr<CoefficientFunction> CoefficientFunction ::
  Operator (shared_ptr<DifferentialOperator> diffop) const
  {
    throw Exception(string("Operator ") + diffop->Name() + string(" not overloaded for CF ")+typeid(*this).name());
  }

  void CoefficientFunction :: SetSpaceDim (int adim)
  {
    TraverseTree ([adim](CoefficientFunction & cf) { cf.spacedim = adim; });
  }
  
  shared_ptr<CoefficientFunction> CoefficientFunctionNoDerivative ::
  Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const
  {
    if (var == this)
      return dir;
    else
      return ZeroCF(Dimensions());
  }

  shared_ptr<CoefficientFunction> CoefficientFunctionNoDerivative ::
  DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var)
      return IdentityCF(this->Dimensions());
    else
      {
        Array<int> resultdims;
        resultdims += this->Dimensions();
        resultdims += var->Dimensions();
        return ZeroCF(resultdims);
      }
  }

  bool CoefficientFunction :: IsZeroCF() const
  {
    return this->GetDescription() == "ZeroCF";
  }
  
  
  void CoefficientFunction :: TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    func(*this);
  }

  void CoefficientFunction :: 
  Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> hvalues) const
  {
    auto values = hvalues.AddSize(ir.Size(), Dimension());
    for (int i = 0; i < ir.Size(); i++)
      Evaluate (ir[i], values.Row(i)); 
  }

  void CoefficientFunction ::   
  Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
  {
    throw ExceptionNOSIMD (string("CF :: simd-Evaluate not implemented for class ") + typeid(*this).name());
  }

  
  /*
  void CoefficientFunction ::   
  Evaluate1 (const SIMD_BaseMappedIntegrationRule & ir, ABareSliceMatrix<double> values) const
  {
    static bool firsttime = true;
    if (firsttime)
    {
        cerr << "Eval1 not implemented for class " << typeid(*this).name() << endl;
        firsttime = false;
      }
    Evaluate (ir, AFlatMatrix<>(Dimension(), ir.IR().GetNIP(), &values.Get(0,0)));
  }
  */

  void CoefficientFunction ::   
  Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const
  {
    if (IsComplex())
      throw ExceptionNOSIMD (string("CF :: simd-Evaluate (complex) not implemented for class ") + typeid(*this).name());
    else
      {
        size_t nv = ir.Size();
        SliceMatrix<SIMD<double>> overlay(Dimension(), nv, 2*values.Dist(), &values(0,0).real());
        Evaluate (ir, overlay);
        for (size_t i = 0; i < Dimension(); i++)
          for (size_t j = nv; j-- > 0; )
            values(i,j) = overlay(i,j);
      }
  }

  
  void CoefficientFunction :: 
  Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const
  {
    for (int i = 0; i < ir.Size(); i++)
      Evaluate (ir[i], values.Row(i).AddSize(Dimension())); 
  }

  /*
  void CoefficientFunction ::
  EvaluateSoA (const BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
  {
    throw Exception(string ("EvaluateSoA called for ") + typeid(*this).name());
  }
    
  void CoefficientFunction ::
  EvaluateSoA (const BaseMappedIntegrationRule & ir, AFlatMatrix<Complex> values) const
  {
    throw Exception(string ("EvaluateSoAComplex called for ") + typeid(*this).name());
  }
  */

  /*
  void CoefficientFunction :: 
  NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                  FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
  {
    // cout << "CoefficientFunction::NonZeroPattern called, type = " << typeid(*this).name() << endl;
    nonzero = true;
    nonzero_deriv = false;
    nonzero_dderiv = false;
  }
  */
  void CoefficientFunction :: 
  NonZeroPattern (const class ProxyUserData & ud, FlatVector<AutoDiffDiff<1,bool>> values) const
  {
    values = AutoDiffDiff<1,bool> (true);
  }

  // shared_ptr<CoefficientFunction> shape = make_shared<ConstantCoefficientFunction>(1);




  IdentityCoefficientFunction :: ~IdentityCoefficientFunction () {  }
  
  
  ///
  ConstantCoefficientFunction ::   
  ConstantCoefficientFunction (double aval) 
    : BASE(1, false), val(aval) 
  {
    elementwise_constant = true;
  }

  ConstantCoefficientFunction ::
  ~ConstantCoefficientFunction ()
  { ; }

  void ConstantCoefficientFunction :: PrintReport (ostream & ost) const
  {
    ost << "ConstantCF, val = " << val << endl;
  }

  string ConstantCoefficientFunction :: GetDescription () const 
  {
    return ToString(val);
  }

  
  /*
  virtual string ConsantCoefficientFunction :: GetDescription() const 
  {
    return "Constant "+ToString(val);
  }
  */
  
  void ConstantCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir,
                                                BareSliceMatrix<double> values) const
  {
    values.AddSize(ir.Size(), 1) = val;
  }

  void ConstantCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir,
                                                BareSliceMatrix<Complex> values) const
  {
    values.AddSize(ir.Size(), 1) = val;
  }

  /*
  template <typename MIR, typename T, ORDERING ORD>
  void ConstantCoefficientFunction ::
  T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    size_t np = ir.Size();    
    __assume (np > 0);
    for (size_t i = 0; i < np; i++)
      values(0,i) = val;
  }
  */
  
  void ConstantCoefficientFunction :: GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    // code.body += Var(index).Declare(code.res_type);
    //code.Declare(code.res_type, index, Dimensions());
    code.Declare(index, Dimensions(), IsComplex());
    code.body += Var(index).Assign(Var(val), false);
  }

  
  shared_ptr<CoefficientFunction>
  ConstantCoefficientFunction :: DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const
  {
    return ZeroCF(var->Dimensions());
  }

  
  ///
  ConstantCoefficientFunctionC ::   
  ConstantCoefficientFunctionC (Complex aval) 
    : CoefficientFunction(1, true), val(aval) 
  { ; }

  ConstantCoefficientFunctionC ::
  ~ConstantCoefficientFunctionC ()
  { ; }
  
  double ConstantCoefficientFunctionC :: Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    throw Exception("no real evaluate for ConstantCF-Complex");
  }

  Complex ConstantCoefficientFunctionC :: EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  { 
    return val;
  }
  
  void ConstantCoefficientFunctionC :: Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const
  {
    values = val;
  }
  
  void ConstantCoefficientFunctionC :: Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const
  {
    // values = val;
    for (auto i : Range(ir.Size()))
      values(i, 0) = val;
  }
  
  void ConstantCoefficientFunctionC :: Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const
  {
    for (auto i : Range(ir.Size()))
      values(0, i) = val;
  }
  
  void ConstantCoefficientFunctionC :: PrintReport (ostream & ost) const
  {
    ost << "ConstantCFC, val = " << val << endl;
  }

  void ConstantCoefficientFunctionC :: GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    code.body += Var(index).Assign(Var(val));
  }

  shared_ptr<CoefficientFunction>
  ConstantCoefficientFunctionC :: DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const
  {
    return ZeroCF(var->Dimensions());
  }

  

  ///
  template<typename SCAL>
  ParameterCoefficientFunction<SCAL> ::
  ParameterCoefficientFunction(SCAL aval)
    : CoefficientFunctionNoDerivative(1, std::is_same_v<SCAL, Complex>), val(aval)
  { SetVariable(true); }

  template<typename SCAL>
  ParameterCoefficientFunction<SCAL> ::
  ~ParameterCoefficientFunction ()
  { ; }

  template<typename SCAL>
  void ParameterCoefficientFunction<SCAL> :: PrintReport (ostream & ost) const
  {
    ost << "ParameterCF, val = " << val << endl;
  }

  template<typename SCAL>
  void ParameterCoefficientFunction<SCAL>::
  DoArchive(Archive& ar)
  {
    CoefficientFunctionNoDerivative::DoArchive(ar);
    ar & val;
  }

  template<typename SCAL>
  double ParameterCoefficientFunction<SCAL>::
  Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    if constexpr(is_same_v<SCAL, Complex>)
      throw Exception("Using double Evaluate for complex ParameterCoefficientFunction!");
    else
      return val;
  }

  template<typename SCAL>
  void ParameterCoefficientFunction<SCAL>::
  Evaluate(const BaseMappedIntegrationRule & ir,
           BareSliceMatrix<double> values) const
  {
    if constexpr(is_same_v<SCAL, Complex>)
      throw Exception("Called double evaluate for complex ParameterCoefficientFunction!");
    else
      values.AddSize(ir.Size(), 1) = val;
  }

  template<typename SCAL>
  void ParameterCoefficientFunction<SCAL>::
  Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
            BareSliceMatrix<SIMD<double>> values) const
  {
    if constexpr(is_same_v<SCAL, Complex>)
      throw Exception("Called double evaluate for complex ParameterCoefficientFunction!");
    else
      values.AddSize(Dimension(), ir.Size()) = val;
  }

  template<typename SCAL>
  void ParameterCoefficientFunction<SCAL>::
  Evaluate(const BaseMappedIntegrationRule & ir,
           BareSliceMatrix<Complex> values) const
  {
    values.AddSize(ir.Size(), 1) = val;
  }

  template<typename SCAL>
  void ParameterCoefficientFunction<SCAL>::
  Evaluate(const SIMD_BaseMappedIntegrationRule & ir,
           BareSliceMatrix<SIMD<Complex>> values) const
  {
    values.AddSize(Dimension(), ir.Size()) = val;
  }

  template<typename SCAL>
  void ParameterCoefficientFunction<SCAL> :: GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    stringstream s;
    constexpr auto type = is_same_v<SCAL, Complex> ? "Complex" : "double";
    s << "*reinterpret_cast<" << type << "*>(" << code.AddPointer(&val) << ")";
    // code.body += Var(index).Declare(code.res_type);
    // code.Declare (code.res_type, index, this->Dimensions());
    code.Declare (index, this->Dimensions(), IsComplex());    
    code.body += Var(index).Assign(s.str(), false);
  }

  template class ParameterCoefficientFunction<double>;
  template class ParameterCoefficientFunction<Complex>;

  void PlaceholderCoefficientFunction::DoArchive(Archive& ar)
  {
    ar.Shallow(cf);
  }

  void PlaceholderCoefficientFunction::Set(shared_ptr<CoefficientFunction> _cf)
  {
    if(Dimensions().Size() != _cf->Dimensions().Size())
      throw Exception("Dimensions of function in PlaceholderCF must not change!");
    for(auto i : Range(Dimensions()))
      if(Dimensions()[i] != _cf->Dimensions()[i])
        throw Exception("Dimensions of function in PlaceholderCF must not change!");
    cf = _cf;
    is_complex = cf->IsComplex();
  }


  DomainConstantCoefficientFunction :: 
  DomainConstantCoefficientFunction (const Array<double> & aval)
    : BASE(1, false), val(aval) { ; }
  
  double DomainConstantCoefficientFunction :: Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    int elind = ip.GetTransformation().GetElementIndex();
    CheckRange (elind);
    return val[elind]; 
  }

  void DomainConstantCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir,
                                                      BareSliceMatrix<double> values) const
  {
    int elind = ir[0].GetTransformation().GetElementIndex();
    CheckRange (elind);    
    values.AddSize(ir.Size(), 1) = val[elind];
  }

  /*
  void DomainConstantCoefficientFunction :: Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
  {
    int elind = ir[0].GetTransformation().GetElementIndex();
    CheckRange (elind);        
    values.AddSize(Dimension(), ir.Size()) = val[elind];
  }
  */
  template <typename MIR, typename T, ORDERING ORD>
  void DomainConstantCoefficientFunction ::
  T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    int elind = ir[0].GetTransformation().GetElementIndex();
    CheckRange (elind);        
    // values.AddSize(Dimension(), ir.Size()) = val[elind];

    size_t np = ir.Size();    
    __assume (np > 0);
    for (size_t i = 0; i < np; i++)
      values(0,i) = val[elind];
  }
  

  void DomainConstantCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const
  {
    int elind = ir[0].GetTransformation().GetElementIndex();
    CheckRange (elind);            
    values.AddSize(ir.Size(), 1) = val[elind]; 
  }

  void DomainConstantCoefficientFunction :: GenerateCode(Code &code, FlatArray<int> inputs, int index) const
    {
      code.header += "double tmp_" + ToLiteral(index) + "["+ToLiteral(val.Size())+"] = {";
      for (auto i : Range(val))
      {
        code.header += ToLiteral(val[i]);
        if(i<val.Size()-1)
          code.header += ", ";
      }
      code.header += "};\n";
      code.header += Var(index).Assign("tmp_"+ToLiteral(index) + "[mir.GetTransformation().GetElementIndex()]");
    }


  DomainConstantCoefficientFunction :: 
  ~DomainConstantCoefficientFunction ()
  { ; }

  DomainVariableCoefficientFunction ::
  DomainVariableCoefficientFunction (const EvalFunction & afun)
    : CoefficientFunction(afun.Dimension(), afun.IsResultComplex()), fun(1)
  {
    fun[0] = make_shared<EvalFunction> (afun);
    numarg = 3;
  }

  DomainVariableCoefficientFunction ::
  DomainVariableCoefficientFunction (const EvalFunction & afun,
				     const Array<shared_ptr<CoefficientFunction>> & adepends_on)
    : CoefficientFunction(afun.Dimension(), afun.IsResultComplex()),
      fun(1), depends_on(adepends_on)
  {
    fun[0] = make_shared<EvalFunction> (afun);
    numarg = 3;
    for (int i = 0; i < depends_on.Size(); i++)
      numarg += depends_on[i]->Dimension();
  }


  DomainVariableCoefficientFunction ::
  DomainVariableCoefficientFunction (const Array<shared_ptr<EvalFunction>> & afun)
    : CoefficientFunction(1, false), fun(afun.Size())
  {
    int hdim = -1;
    for (int i = 0; i < fun.Size(); i++)
      if (afun[i])
        {
          fun[i] = afun[i];
          if (fun[i]->IsResultComplex())
            is_complex = true;
          hdim = fun[i]->Dimension();
        }
      else
        fun[i] = nullptr;
    SetDimension (hdim);
    numarg = 3;
  }

  DomainVariableCoefficientFunction ::
  DomainVariableCoefficientFunction (const Array<shared_ptr<EvalFunction>> & afun,
				     const Array<shared_ptr<CoefficientFunction>> & adepends_on)
    : CoefficientFunction(1, false), fun(afun.Size()), depends_on(adepends_on)
  {
    int hdim = -1;
    for (int i = 0; i < fun.Size(); i++)
      if (afun[i])
        {
          fun[i] = afun[i];
          if (fun[i]->IsResultComplex())
            is_complex = true;
          hdim = fun[i]->Dimension();
        }
      else
        fun[i] = nullptr;

    SetDimension (hdim);
    numarg = 3;
    for (int i = 0; i < depends_on.Size(); i++)
      numarg += depends_on[i]->Dimension();
  }


  DomainVariableCoefficientFunction ::
  ~DomainVariableCoefficientFunction ()
  {
    ;
    /*
    for (int i = 0; i < fun.Size(); i++)
      delete fun[i];
    */
  }

  double DomainVariableCoefficientFunction ::
  Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> result;
    Evaluate (ip, result);
    return result(0);
    /*
      int numarg = max2(3, depends_on.Size());
      VectorMem<10> args(numarg);
      args.Range(0,DIM) = static_cast<const DimMappedIntegrationPoint<DIM>&>(ip).GetPoint();
    
      for (int i = 3; i < depends_on.Size(); i++)
      args(i) = depends_on[i] -> Evaluate (ip);

      int elind = ip.GetTransformation().GetElementIndex();
      if (fun.Size() == 1) elind = 0;
      double val = fun[elind]->Eval (&args(0));
      return val;
    */
  }

  bool DomainVariableCoefficientFunction :: IsComplex() const 
  {
    for (int i = 0; i < fun.Size(); i++)
      if (fun[i]->IsResultComplex()) return true;
    return false;
  }
  
  int DomainVariableCoefficientFunction :: Dimension() const
  { 
    return fun[0]->Dimension(); 
  }


  Complex DomainVariableCoefficientFunction ::
  EvaluateComplex (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1, Complex> result;
    Evaluate (ip, result);
    return result(0);
    /*
      int elind = ip.GetTransformation().GetElementIndex();
      Vec<DIM, Complex> hp;
      for (int i = 0; i < DIM; i++)
      hp(i) = static_cast<const DimMappedIntegrationPoint<DIM>&>(ip).GetPoint()(i);
      return fun[elind]->Eval (&hp(0));
    */
  }
  
  void DomainVariableCoefficientFunction ::
  Evaluate(const BaseMappedIntegrationPoint & ip,
	   FlatVector<> result) const
  {
    int elind = ip.GetTransformation().GetElementIndex();
    if (fun.Size() == 1) elind = 0;
    
    if (! fun[elind] -> IsComplex ())
      {
	VectorMem<10> args(numarg);
	// args.Range(0,DIM) = static_cast<const DimMappedIntegrationPoint<DIM>&>(ip).GetPoint();
        args.Range(0,ip.DimSpace()) = ip.GetPoint();
	
	for (int i = 0, an = 3; i < depends_on.Size(); i++)
	  {
	    int dim = depends_on[i]->Dimension();
	    depends_on[i] -> Evaluate (ip, args.Range(an,an+dim));
	    an += dim;
	  }
	fun[elind]->Eval (&args(0), &result(0), result.Size());      
      }
    else
      {
	VectorMem<10, Complex> args(numarg);
	// args.Range(0,DIM) = static_cast<const DimMappedIntegrationPoint<DIM>&>(ip).GetPoint();
        args.Range(0,ip.DimSpace()) = ip.GetPoint();
	
	for (int i = 0, an = 3; i < depends_on.Size(); i++)
	  {
	    int dim = depends_on[i]->Dimension();
	    depends_on[i] -> Evaluate (ip, args.Range(an,an+dim));
	    an += dim;
	  }
	fun[elind]->Eval (&args(0), &result(0), result.Size());      
      }
  }


  void DomainVariableCoefficientFunction ::
  Evaluate(const BaseMappedIntegrationPoint & ip,
           FlatVector<Complex> result) const
  {
    VectorMem<10,Complex> args(numarg);
    args = -47;
    // args.Range(0,DIM) = static_cast<const DimMappedIntegrationPoint<DIM>&>(ip).GetPoint();
    args.Range(0,ip.DimSpace()) = ip.GetPoint();
    for (int i = 0, an = 3; i < depends_on.Size(); i++)
      {
        int dim = depends_on[i]->Dimension();
        depends_on[i] -> Evaluate (ip, args.Range(an,an+dim));
        an += dim;
      }
    
    int elind = ip.GetTransformation().GetElementIndex();
    if (fun.Size() == 1) elind = 0;
    fun[elind]->Eval (&args(0), &result(0), result.Size());
  }
  
  
void DomainVariableCoefficientFunction ::
Evaluate (const BaseMappedIntegrationRule & ir, 
	  BareSliceMatrix<double> values) const
{
  if (ir.Size() == 0) return;
  int elind = ir.GetTransformation().GetElementIndex();
  if (fun.Size() == 1) elind = 0;

  if (! fun[elind] -> IsComplex ())
    {
      ArrayMem<double,2000> mem(ir.Size()*numarg);
      FlatMatrix<> args(ir.Size(), numarg, &mem[0]);
      
      int dim = ir[0].DimSpace();
      switch (dim)
        {
        case 2:
          for (int i = 0; i < ir.Size(); i++)
            args.Row(i).Range(0,2) = ir[i].GetPoint();
          break;
        case 3:
          for (int i = 0; i < ir.Size(); i++)
            args.Row(i).Range(0,3) = ir[i].GetPoint();
          break;
        default:
          for (int i = 0; i < ir.Size(); i++)
            args.Row(i).Range(0,dim) = ir[i].GetPoint();
        }
      

      /*
	args.Row(i).Range(0,DIM) = 
	  static_cast<const DimMappedIntegrationPoint<DIM> & > (ir[i]).GetPoint();
      */
      for (int i = 0, an = 3; i < depends_on.Size(); i++)
	{
	  int dim = depends_on[i]->Dimension();
	  Matrix<> hmat(ir.Size(), dim);
	  depends_on[i] -> Evaluate (ir, hmat);
	  args.Cols(an,an+dim) = hmat;
	  an += dim;
	}
      for (int i = 0; i < ir.Size(); i++)
	fun[elind]->Eval (&args(i,0), &values(i,0), values.Dist());
    }
  else
    {
      Matrix<Complex> args(ir.Size(), numarg);
      for (int i = 0; i < ir.Size(); i++)
	args.Row(i).Range(0,ir[i].DimSpace()) = ir[i].GetPoint();
      
      for (int i = 0, an = 3; i < depends_on.Size(); i++)
	{
	  int dim = depends_on[i]->Dimension();
	  Matrix<Complex> hmat(ir.Size(), dim);
	  depends_on[i] -> Evaluate (ir, hmat);
	  args.Cols(an,an+dim) = hmat;
	  an += dim;
	}
    
      for (int i = 0; i < ir.Size(); i++)
	fun[elind]->Eval (&args(i,0), &values(i,0), values.Dist());
    }
}

  
void DomainVariableCoefficientFunction :: PrintReport (ostream & ost) const
{
  *testout << "DomainVariableCoefficientFunction, functions are: " << endl;
  for (int i = 0; i < fun.Size(); i++)
    fun[i] -> Print(ost);
}

void DomainVariableCoefficientFunction :: GenerateCode(Code &code, FlatArray<int> inputs, int index) const
{
  code.body += "// DomainVariableCoefficientFunction: not implemented";
}

/*
  template class DomainVariableCoefficientFunction<1>;
  template class DomainVariableCoefficientFunction<2>;
  template class DomainVariableCoefficientFunction<3>;
*/

PolynomialCoefficientFunction::
PolynomialCoefficientFunction(const Array < Array< Array<double>* >* > & polycoeffs_in,
                              const Array < Array<double>* > & polybounds_in)
  : CoefficientFunction(1, false), polycoeffs(polycoeffs_in), polybounds(polybounds_in)
{ ; }

PolynomialCoefficientFunction::
PolynomialCoefficientFunction(const Array < Array<double>* > & polycoeffs_in)
  : CoefficientFunction(1, false)
{
  polycoeffs.SetSize(polycoeffs_in.Size());
  polybounds.SetSize(polycoeffs_in.Size());
  
  for(int i=0; i<polycoeffs_in.Size(); i++)
    {
      polycoeffs[i] = new Array< Array<double>* >(1);
      (*polycoeffs[i])[0] = polycoeffs_in[i];
      polybounds[i] = new Array<double>(0);
    } 
}


PolynomialCoefficientFunction::~PolynomialCoefficientFunction()
{
  for(int i=0; i<polycoeffs.Size(); i++)
    {
      delete polybounds[i];
      for(int j=0; j<polycoeffs[i]->Size(); j++)
	{
	  delete (*polycoeffs[i])[j];
	}
      delete polycoeffs[i];
    }
  polycoeffs.DeleteAll();
  polybounds.DeleteAll();
}
    
  
  
double PolynomialCoefficientFunction::Evaluate (const BaseMappedIntegrationPoint & ip) const
{
  return Evaluate(ip,0);
}



double PolynomialCoefficientFunction::EvalPoly(const double t, const Array<double> & coeffs) const
{
  const int last = coeffs.Size()-1;
    
  double retval = coeffs[last];
  for(int i=last-1; i>=0; i--)
    {
      retval *= t;
      retval += coeffs[i];
    }

  return retval;    
}


double PolynomialCoefficientFunction::EvalPolyDeri(const double t, const Array<double> & coeffs) const
{
  const int last = coeffs.Size()-1;

  double retval = last*coeffs[last];
  for(int i=last-1; i>=1; i--)
    {
      retval *= t;
      retval += i*coeffs[i];
    }  

  return retval;    
}


double PolynomialCoefficientFunction::Evaluate (const BaseMappedIntegrationPoint & ip, const double & t) const
{
  const int elind = ip.GetTransformation().GetElementIndex();
    
  if (elind < 0 || elind >= polycoeffs.Size())
    {
      ostringstream ost;
      ost << "PolynomialCoefficientFunction: Element index "
	  << elind << " out of range 0 - " << polycoeffs.Size()-1 << endl;
      throw Exception (ost.str());
    }
 
  int pos;
  for(pos=0; pos < polybounds[elind]->Size() && t > (*polybounds[elind])[pos]; pos++){}
   
  return EvalPoly(t,*((*(polycoeffs[elind]))[pos]));

    
}


 
double PolynomialCoefficientFunction::EvaluateDeri (const BaseMappedIntegrationPoint & ip, const double & t) const
{
  const int elind = ip.GetTransformation().GetElementIndex();
    
  if (elind < 0 || elind >= polycoeffs.Size())
    {
      ostringstream ost;
      ost << "PolynomialCoefficientFunction: Element index "
	  << elind << " out of range 0 - " << polycoeffs.Size()-1 << endl;
      throw Exception (ost.str());
    }

  int pos;
  for(pos=0; pos < polybounds[elind]->Size() && t > (*polybounds[elind])[pos]; pos++){}

  return EvalPolyDeri(t,*((*(polycoeffs[elind]))[pos]));
}


double PolynomialCoefficientFunction::EvaluateConst () const
{
  return (*(*polycoeffs[0])[0])[0];
}



//////////////////

FileCoefficientFunction :: FileCoefficientFunction ()
  : CoefficientFunction(1, false)
{
  writeips = false;
}

  
FileCoefficientFunction :: FileCoefficientFunction (const string & filename)
  : CoefficientFunction(1, false)  
{
  StartWriteIps(filename);
}

FileCoefficientFunction :: FileCoefficientFunction (const string & aipfilename,
						    const string & ainfofilename,
						    const string & avaluesfilename,
						    const bool loadvalues)
  : CoefficientFunction(1, false)  
{
  ipfilename = aipfilename;
  infofilename = ainfofilename;
  valuesfilename = avaluesfilename;

  if(loadvalues)
    {
      writeips = false;
      LoadValues();
    }
  else
    StartWriteIps();
}
    

  
void FileCoefficientFunction :: EmptyValues(void)
{
  for(int i=0; i<ValuesAtIps.Size(); i++)
    delete ValuesAtIps[i];

  ValuesAtIps.SetSize(0);
}

void FileCoefficientFunction :: Reset(void)
{
  EmptyValues();
}

FileCoefficientFunction :: ~FileCoefficientFunction()
{
  if(writeips)
    StopWriteIps(); 

  EmptyValues();
}


void FileCoefficientFunction :: LoadValues(const string & filename)
{
  cout << "Loading values for coefficient function ..."; cout.flush();

  if(writeips) cerr << "WARNING: CoefficientFunction still writing points to \"" 
		    << ipfilename << "\"" << endl;

  ifstream infile(filename.c_str());
    
  int numels,numips,numentries,eln,ipn;
  double val;

  infile >> numels;
  infile >> numips;
  infile >> numentries;
    
  EmptyValues();
    
  ValuesAtIps.SetSize(numels);
    
  for(int i=0; i<numels; i++)
    {
      ValuesAtIps[i] = new Array<double>(numips);
      *(ValuesAtIps[i]) = 0.;
    }

  for(int i=0; i<numentries; i++)
    {
      infile >> eln;
      infile >> ipn;
      infile >> val;
      (*(ValuesAtIps[eln]))[ipn] = val;
    }

  infile.close();
  cout << "done\n";
}



double FileCoefficientFunction :: Evaluate (const BaseMappedIntegrationPoint & ip) const
{
  const ElementTransformation & eltrans = ip.GetTransformation();
  const int elnum = eltrans.GetElementNr();
  const int ipnum = ip.GetIPNr();

  if(writeips)
    {
      if(elnum > maxelnum) const_cast<int&> (maxelnum) = elnum;
      if(ipnum > maxipnum) const_cast<int&> (maxipnum) = ipnum;
      const_cast<int&> (totalipnum)++;

      Vec<3> point;
      eltrans.CalcPoint(ip.IP(),point);

      const_cast<ofstream&> (outfile) << elnum << " " << ipnum << " " << point << "\n";
    }

  if(elnum < ValuesAtIps.Size())
    {
      return (*(ValuesAtIps[elnum]))[ipnum];
    }

  return 0.;
}

void FileCoefficientFunction :: StartWriteIps(const string & filename)
{
  writeips = true;
  maxelnum = 0;
  maxipnum = 0;
  totalipnum = 0;

  outfile.open(filename.c_str());
  outfile.precision(12);
    
}

void FileCoefficientFunction :: StopWriteIps(const string & infofilename)
{
  writeips = false;

  outfile.close();

    
  cout << "Stopped writing to " << ipfilename << endl;
  cout << "Writing info file to " << infofilename << endl;

  ofstream info(infofilename.c_str());

  info << "numelts " << maxelnum+1 << endl
       << "maxnumips " << maxipnum+1 << endl
       << "totalipnum " << totalipnum << endl;

  info.close();

}


class UnitVectorCoefficientFunction : public T_CoefficientFunction<UnitVectorCoefficientFunction>
{
  using BASE = T_CoefficientFunction<UnitVectorCoefficientFunction>;
  int coord;
public:
  UnitVectorCoefficientFunction () : T_CoefficientFunction<UnitVectorCoefficientFunction>(1, false)
  {
    SetDimension(1);
  }
  
  UnitVectorCoefficientFunction (int dim, int _coord) : T_CoefficientFunction<UnitVectorCoefficientFunction>(dim, false), coord(_coord)
  {
    if (coord >= dim)
      throw Exception("In UnitVectorCoefficientFunction: coord >= dim");
  }

  int GetCoordinate() const { return coord; }


  void DoArchive (Archive & archive) override
  {
    BASE::DoArchive(archive);
    archive & coord;
  }

  virtual void PrintReport (ostream & ost) const override
  {
    ost << "UnitVectorCoefficientFunction";
  }

  virtual string GetDescription() const override
  {
    return "UnitVectorCF " + ToString(coord);
  }


  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    func(*this);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    for (int i : Range(Dimension()))
        {
          /*
          if (i == coord)
            code.body += Var(index,i, Dimensions()).Assign(string("1.0"));
          else
            code.body += Var(index,i, Dimensions()).Assign(string("0.0"));
          */
          code.body += Var(index,i, Dimensions()).Assign(string( (i==coord) ? "1.0" : "0.0"));          
        }
  }

  using T_CoefficientFunction<UnitVectorCoefficientFunction>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    if (Dimension() > 1)
      throw Exception ("UnitVectorCF:: scalar evaluate for non scalar called");
    return 1.0;
  }
 
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    result = 0.0;
    result(coord) = 1.0;
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   BareSliceMatrix<T,ORD> values) const
  {
    values.AddSize(Dimension(), ir.Size()) = T(0.0);
    values.Row(coord).Range(ir.Size()) = T(1.0);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    values.AddSize(Dimension(), ir.Size()) = T(0.0);
    values.Row(coord).Range(ir.Size()) = T(1.0);
  }



  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    values = AutoDiffDiff<1,bool>(false);
    values(coord) = AutoDiffDiff<1,bool>(true);
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    values = AutoDiffDiff<1,bool>(false);
    values(coord) = AutoDiffDiff<1,bool>(true);
  }

  using CoefficientFunction::Operator;
  shared_ptr<CoefficientFunction> Operator (const string & name) const override
  {
    if (spacedim == -1)
      throw Exception("cannot differentiate constant since we don't know the space dimension, use 'coef.spacedim=dim'");
    if (name != "grad")
      throw Exception ("cannot apply operator "+name+" for constant");
    return ZeroCF ( Array( { spacedim } ) );
  }
  
  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    return ZeroCF(Dimensions());
  }
  
};

class ZeroCoefficientFunction : public T_CoefficientFunction<ZeroCoefficientFunction>
{
  using BASE = T_CoefficientFunction<ZeroCoefficientFunction>;
public:
  ZeroCoefficientFunction () : T_CoefficientFunction<ZeroCoefficientFunction>(1, false)
  {
    SetDimension(1);
  }
  
  ZeroCoefficientFunction (int dim) : T_CoefficientFunction<ZeroCoefficientFunction>(1, false)
  {
    SetDimensions(Array({dim}));
  }

  ZeroCoefficientFunction (int dim1, int dim2) : T_CoefficientFunction<ZeroCoefficientFunction>(1, false)
  {
    SetDimensions(Array( {dim1,dim2} ));
  }

  ZeroCoefficientFunction (Array<int> dims) : T_CoefficientFunction<ZeroCoefficientFunction>(1, false)
  {
    SetDimensions(move(dims));
  }

  

  void DoArchive (Archive & archive) override
  {
    BASE::DoArchive(archive);
  }

  virtual void PrintReport (ostream & ost) const override
  {
    ost << "ZeroCoefficientFunction";
  }

  virtual string GetDescription() const override
  {
    return "ZeroCF";
  }

  virtual bool IsZeroCF() const override { return true; }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    func(*this);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    // code.Declare (code.res_type, index, this->Dimensions());
    code.Declare (index, this->Dimensions(), IsComplex());
    for (int i = 0; i < Dimension(); i++)
      code.body += Var(index,i,this->Dimensions()).Assign(string("0.0"), false);      
  }

  using T_CoefficientFunction<ZeroCoefficientFunction>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    if (Dimension() > 1)
      throw Exception ("ZeroCF:: scalar evaluate for non scalar called");
    return 0.0;
  }
 
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    result = 0.0;
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   BareSliceMatrix<T,ORD> values) const
  {
    values.AddSize(Dimension(), ir.Size()) = T(0.0);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    values.AddSize(Dimension(), ir.Size()) = T(0.0);
  }



  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    values = AutoDiffDiff<1,bool>(false);
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    values = AutoDiffDiff<1,bool>(false);
  }

  using CoefficientFunction::Operator;
  shared_ptr<CoefficientFunction> Operator (const string & name) const override
  {
    if (spacedim == -1)
      throw Exception("cannot differentiate constant since we don't know the space dimension, use 'coef.spacedim=dim'");
    if (name != "grad")
      throw Exception ("cannot apply operator "+name+" for constant");
    return ZeroCF ( Array( { spacedim } ) );
  }
  
  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    return const_cast<ZeroCoefficientFunction*>(this)->shared_from_this();
  }
  
};


  
class ScaleCoefficientFunction : public T_CoefficientFunction<ScaleCoefficientFunction>
{
  double scal;
  shared_ptr<CoefficientFunction> c1;
  typedef T_CoefficientFunction<ScaleCoefficientFunction> BASE;
public:
  ScaleCoefficientFunction() = default;
  ScaleCoefficientFunction (double ascal, 
                            shared_ptr<CoefficientFunction> ac1)
    : BASE(ac1->Dimension(), ac1->IsComplex()),
      scal(ascal), c1(ac1)
  {
    SetDimensions(c1->Dimensions());
    elementwise_constant = c1->ElementwiseConstant();
  }
  
  void DoArchive (Archive & archive) override
  {
    BASE::DoArchive(archive);
    archive.Shallow(c1) & scal;
  }

  virtual void PrintReport (ostream & ost) const override
  {
    ost << scal << "*(";
    c1->PrintReport(ost);
    ost << ")";
  }

  virtual string GetDescription() const override
  {
    return "scale "+ToString(scal);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    // code.Declare (code.res_type, index, this->Dimensions());
    code.Declare (index, this->Dimensions(), IsComplex());  
    for (int i = 0; i < Dimension(); i++)
      code.body += Var(index,i,this->Dimensions()).Assign(Var(scal) * Var(inputs[0],i,c1->Dimensions()), false);      
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 }); }

  virtual bool DefinedOn (const ElementTransformation & trafo) override
  { return c1->DefinedOn(trafo); }
    
  using BASE::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    return scal * c1->Evaluate(ip);
  }
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const override
  {
    return scal * c1->EvaluateComplex(ip);
  }
  virtual double EvaluateConst () const override
  {
    return scal * c1->EvaluateConst();
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
  
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         BareSliceMatrix<double> values) const override
  {
    c1->Evaluate (ir, values);
    values.AddSize(ir.Size(), Dimension()) *= scal;
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   BareSliceMatrix<T,ORD> values) const
  {
    c1->Evaluate (ir, values);
    values.AddSize(Dimension(), ir.Size()) *= scal;
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto in0 = input[0];
    // values.AddSize(Dimension(), ir.Size()) = scal * in0;
    // SliceMatrix<T,ORD> sval(Dimension(), ir.Size(), values.Dist(), values.Data());
    // sval = scal * in0;   // failing on WIN-AVX

    /*
    // working on WIN-AVX
    for (int i = 0; i < Dimension(); i++)
      for (int j = 0; j < ir.Size(); j++)
        sval(i,j) = scal * in0(i,j);
    */

    /*
      working:
      SliceMatrix<T,ORD> sval(Dimension(), ir.Size(), values.Dist(), values.Data());
    auto prod = scal*in0;
    for (int i = 0; i < Dimension(); i++)
      for (int j = 0; j < ir.Size(); j++)
        sval(i,j) = prod(i,j);
    */

    /*
      // working
    auto sval = values.AddSize(Dimension(), ir.Size());
    auto prod = scal*in0;
    for (int i = 0; i < Dimension(); i++)
      for (int j = 0; j < ir.Size(); j++)
        sval(i,j) = prod(i,j);
    */

    /*
    // failing
    auto sval = values.AddSize(Dimension(), ir.Size());
    auto prod = scal*in0;
    for (int i = 0; i < sval.Height(); i++)
      for (int j = 0; j < sval.Width(); j++)
        sval(i,j) = prod(i,j);
    */

    // failing, but works with check
    auto sval = values.AddSize(Dimension(), ir.Size());
    auto prod = scal*in0;
    sval = prod;
    // assert(sval.Height() == Dimension());
    // assert(sval.Width() == ir.Size());
    
    // if (sval.Height() != Dimension()) throw Exception("wrong height");
    // if (sval.Width() != ir.Size()) throw Exception("wrong width");
    
    // values.AddSize(Dimension(), ir.Size()) = scal*in0;
    
    /*
      // working on WIN-AVX
    for (int i = 0; i < Dimension(); i++)
      for (int j = 0; j < ir.Size(); j++)
        values(i,j) = scal * in0(i,j);
    */
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         BareSliceMatrix<Complex> values) const override
  {
    c1->Evaluate (ir, values);
    values.AddSize(ir.Size(), Dimension()) *= scal;
  }

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    c1->NonZeroPattern (ud, nonzero, nonzero_deriv, nonzero_dderiv);
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    c1->NonZeroPattern (ud, values);
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    values = input[0];
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return scal * c1->Diff(var, dir);
  }
  
  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];
             
    if (this == var)
      {
        if (this -> Dimensions().Size() == 0)
            return make_shared<ConstantCoefficientFunction>(1);
        return IdentityCF(this->Dimensions());
      }
    auto res = scal * c1->DiffJacobi(var, cache);
    cache[thisptr] = res;
    return res;
  }
};


class ScaleCoefficientFunctionC : public CoefficientFunction
{
  Complex scal;
  shared_ptr<CoefficientFunction> c1;
public:
  ScaleCoefficientFunctionC() = default;
  ScaleCoefficientFunctionC (Complex ascal, 
                            shared_ptr<CoefficientFunction> ac1)
    : CoefficientFunction(ac1->Dimension(), true), scal(ascal), c1(ac1)
  {
    SetDimensions (c1->Dimensions());
  }

  void DoArchive(Archive& ar) override
  {
    CoefficientFunction::DoArchive(ar);
    ar.Shallow(c1) & scal;
  }
  
  // virtual bool IsComplex() const { return true; }
  // virtual int Dimension() const { return c1->Dimension(); }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 }); }
  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    for (int i = 0; i < c1->Dimension(); i++)
      code.body += Var(index,i,Dimensions()).Assign( Var(scal) * Var(inputs[0],i,c1->Dimensions()) );
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const  override
  {
    throw Exception ("real Evaluate called for complex ScaleCF");
  }
  
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const override
  {
    return scal * c1->EvaluateComplex(ip);    
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        BareSliceMatrix<Complex> result) const override
  {
    c1->Evaluate (ir, result);
    result.AddSize(ir.Size(), Dimension()) *= scal;
  }

  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         BareSliceMatrix<SIMD<Complex>> values) const override
  {
    c1->Evaluate (ir, values);
    values.AddSize(Dimension(), ir.Size()) *= scal;
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
                         BareSliceMatrix<AutoDiffDiff<1,double>> values) const override
  {
    throw Exception ("can't diff complex CF (ScaleCoefficientFunctionC)");
  }
  

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    c1->NonZeroPattern (ud, nonzero, nonzero_deriv, nonzero_dderiv);
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    c1->NonZeroPattern (ud, values);
  }
  
  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return scal * c1->Diff(var, dir);
  }

  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];
             
    if (this == var)
      {
        if (this -> Dimensions().Size() == 0)
            return make_shared<ConstantCoefficientFunction>(1);
        return IdentityCF(this->Dimensions());
      }
    auto res = scal * c1->DiffJacobi(var, cache);
    cache[thisptr] = res;
    return res;
  }
  
};

// ***********************************************************************************

class MultScalVecCoefficientFunction : public T_CoefficientFunction<MultScalVecCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;  // scalar
  shared_ptr<CoefficientFunction> c2;  // vector
  typedef T_CoefficientFunction<MultScalVecCoefficientFunction> BASE;
public:
  MultScalVecCoefficientFunction() = default;
  MultScalVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                  shared_ptr<CoefficientFunction> ac2)
    : BASE(ac2->Dimension(), ac1->IsComplex() || ac2->IsComplex()),
      c1(ac1), c2(ac2)
  {
    SetDimensions (c2->Dimensions());
  }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1).Shallow(c2);
  }
  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return { c1, c2 }; }

  virtual string GetDescription () const override
  {
    switch (Dimensions().Size())
      {
      case 1: return "scalar-vector multiply";
      case 2: return "scalar-matrix multiply";    
      default: return "scalar-tensor multiply";
      }
  }
  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    // code.Declare (code.res_type, index, Dimensions());
    code.Declare (index, Dimensions(), IsComplex());
    if (code_uses_tensors)
      {
        code.body += "for (size_t i = 0; i < "+ToString(this->Dimension())+"; i++)\n";
        code.body += "var_"+ToString(index)+"[i] = var_"+ToString(inputs[0])+"[0]*var_"+ToString(inputs[1])+"[i];\n";
      }
    else
      for (int i = 0; i < Dimension(); i++)
        code.body += Var(index,i,Dimensions()).Assign( Var(inputs[0]) * Var(inputs[1],i,Dimensions()), false );      
  }

  using BASE::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    throw Exception ("double MultScalVecCF::Evaluate called");
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    Vec<1> v1;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, result);
    result *= v1(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    Vec<1,Complex> v1;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, result);
    result *= v1(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        BareSliceMatrix<Complex> result) const override
  {
    STACK_ARRAY(double, hmem1, 2*ir.Size());
    FlatMatrix<Complex> temp1(ir.Size(), 1, reinterpret_cast<Complex*> (&hmem1[0]));
    
    c1->Evaluate(ir, temp1);
    c2->Evaluate(ir, result);
    for (int i = 0; i < ir.Size(); i++)
      result.Row(i).AddSize(Dimension()) *= temp1(i,0);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    size_t w = ir.Size();
    __assume (w > 0);
    STACK_ARRAY(T, hmem1, w);
    FlatMatrix<T,ORD> temp1(1, w, &hmem1[0]);
    
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, values);

    for (size_t j = 0; j < Dimension(); j++)
      for (size_t i = 0; i < w; i++)
        values(j,i) *= temp1(0,i);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto in0 = input[0];
    auto in1 = input[1];
    size_t dim = Dimension();
    size_t np = ir.Size();
    
    for (size_t j = 0; j < dim; j++)
      for (size_t i = 0; i < np; i++)
        values(j,i) = in0(0,i) * in1(j,i);
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return c1->Diff(var,dir)*c2 + c1 * c2->Diff(var,dir);
  }

  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];
    int dimvar = var->Dimension();
    int mydim = this->Dimension();    
      
    Array<int> dimres { this -> Dimensions() + var->Dimensions() };

    if (this == var)
      return IdentityCF(this->Dimensions());

    auto diffc1 = c1->DiffJacobi(var, cache);
    auto diffc2 = c2->DiffJacobi(var, cache);
    
    auto prod1_ = c2 -> Reshape(mydim, 1) * diffc1 -> Reshape(1, dimvar);
    auto prod1 = prod1_ -> Reshape(dimres);

    auto prod2 = c1 * diffc2;
    auto res = prod1 + prod2;
    cache[thisptr] = res;
    return res;
  }

  
  
  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    int dim = Dimension();
    Vector<bool> v1(1), d1(1), dd1(1);
    Vector<bool> v2(dim), d2(dim), dd2(dim);
    c1->NonZeroPattern (ud, v1, d1, dd1);
    c2->NonZeroPattern (ud, v2, d2, dd2);
    for (auto i : Range(dim))
      {
        nonzero(i) = v1(0) && v2(i);
        nonzero_deriv(i) = (v1(0) && d2(i)) || (d1(0) && v2(i));
        nonzero_dderiv(i) = (v1(0) && dd2(i)) || (d1(0) && d2(i)) || (dd1(0) && v2(i));
      }
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    int dim = Dimension();
    Vector<AutoDiffDiff<1,bool>> v1(1), v2(dim);
    c1->NonZeroPattern (ud, v1);
    c2->NonZeroPattern (ud, v2);
    for (size_t j = 0; j < dim; j++)
      values(j) = v1(0) * v2(j);
  }


  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto in0 = input[0];
    auto in1 = input[1];
    size_t dim = Dimension();
    
    for (size_t j = 0; j < dim; j++)
      values(j) = in0(0) * in1(j);
  }
};


class MultVecVecCoefficientFunction : public T_CoefficientFunction<MultVecVecCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
  int dim1;
  using BASE = T_CoefficientFunction<MultVecVecCoefficientFunction>;
public:
  MultVecVecCoefficientFunction() = default;
  MultVecVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                 shared_ptr<CoefficientFunction> ac2)
    : T_CoefficientFunction<MultVecVecCoefficientFunction>(1, ac1->IsComplex() || ac2->IsComplex()), c1(ac1), c2(ac2)
  {
    elementwise_constant = c1->ElementwiseConstant() && c2->ElementwiseConstant();
    dim1 = c1->Dimension();
    if (dim1 != c2->Dimension())
      throw Exception("MultVecVec : dimensions don't fit");
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1).Shallow(c2) & dim1;
  }
  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    CodeExpr result;
    for (int i = 0; i < c1->Dimension(); i++)
      result += Var(inputs[0], i, c1->Dimensions()) * Var(inputs[1], i, c2->Dimensions());
    // code.Declare (code.res_type, index, Dimensions());
    code.Declare (index, Dimensions(), IsComplex());  
    code.body += Var(index).Assign(result.S(), false);
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1, c2 }); }  
  
    using T_CoefficientFunction<MultVecVecCoefficientFunction>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    STACK_ARRAY(double, hmem1, dim1);
    FlatVector<> v1(dim1, hmem1);
    STACK_ARRAY(double, hmem2, dim1);
    FlatVector<> v2(dim1, hmem2);

    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    Vector<Complex> v1(dim1), v2(dim1);
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    size_t w = ir.Size();
    __assume (w > 0);

    size_t dim = dim1; // Dimension();
    STACK_ARRAY(T, hmem, 2*dim*w);
    FlatMatrix<T,ORD> temp1(dim, w, &hmem[0]);
    FlatMatrix<T,ORD> temp2(dim, w, &hmem[dim*w]);
    
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, temp2);

    for (size_t i = 0; i < w; i++)
      {
        T sum{0.0};
        for (size_t j = 0; j < dim; j++)
          sum += temp1(j,i) * temp2(j,i);
        values(0,i) = sum; 
      }
  }


  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto in0 = input[0];
    auto in1 = input[1];
    size_t dim = dim1;
    size_t np = ir.Size();

    for (size_t i = 0; i < np; i++)
      {
        T sum{0.0};
        for (size_t j = 0; j < dim; j++)
          sum += in0(j,i) * in1(j,i);
        values(0,i) = sum; 
      }    
  }  

  /*
  virtual bool ElementwiseConstant () const override
  { return c1->ElementwiseConstant() && c2->ElementwiseConstant(); }
  */
  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    Vector<bool> v1(dim1), v2(dim1), d1(dim1), d2(dim1), dd1(dim1), dd2(dim1);
    c1->NonZeroPattern (ud, v1, d1, dd1);
    c2->NonZeroPattern (ud, v2, d2, dd2);
    bool nz = false, nzd = false, nzdd = false;
    for (int i = 0; i < dim1; i++)
      {
        if (v1(i) && v2(i)) nz = true;
        if ((v1(i) && d2(i)) || (d1(i) && v2(i))) nzd = true;
        if ((v1(i) && dd2(i)) || (d1(i) && d2(i)) || (dd1(i) && v2(i))) nzdd = true;
      }
    nonzero = nz;
    nonzero_deriv = nzd;
    nonzero_dderiv = nzdd;
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> v1(dim1), v2(dim1);
    c1->NonZeroPattern (ud, v1);
    c2->NonZeroPattern (ud, v2);
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < dim1; i++)
      sum += v1(i)*v2(i);
    values(0) = sum;
  }

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto v1 = input[0];
    auto v2 = input[1];
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < dim1; i++)
      sum += v1(i)*v2(i);
    values(0) = sum;
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
      if (this == var) return dir;
      return InnerProduct(c1->Diff(var,dir),c2) + InnerProduct(c1,c2->Diff(var,dir));
  }

  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    if (this == var)
      return make_shared<ConstantCoefficientFunction> (1);

    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    shared_ptr<CoefficientFunction> dv1v2, dv2v1;
    int dimip = c1->Dimension();
    int dimvar = var->Dimension();

    auto vc1 = c1->Reshape( dimip );
    auto vc2 = c2->Reshape( dimip );

    if (c1.get() == var)
      dv1v2 = c2;
    else
      {
        auto dvc1 = c1->DiffJacobi (var, cache);
        dv1v2 = dvc1 -> Reshape(dimip, dimvar) -> Transpose() * vc2;
        dv1v2 = dv1v2 -> Reshape (var->Dimensions());
      }

    if (c2.get() == var)
      dv2v1 = c1;
    else
      {
        auto dvc2 = c2->DiffJacobi (var, cache);
        dv2v1 = dvc2 -> Reshape(dimip, dimvar) -> Transpose() * vc1;
        dv2v1 = dv2v1 -> Reshape (var->Dimensions());
      }
    auto res = dv1v2 + dv2v1;
    cache[thisptr] = res;
    return res;
  }
  
};




template <int DIM>
class T_MultVecVecCoefficientFunction : public T_CoefficientFunction<T_MultVecVecCoefficientFunction<DIM>>
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
  using BASE = T_CoefficientFunction<T_MultVecVecCoefficientFunction<DIM>>;
  using typename BASE::T_DJC;
public:
  T_MultVecVecCoefficientFunction() = default;
  T_MultVecVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                   shared_ptr<CoefficientFunction> ac2)
    : T_CoefficientFunction<T_MultVecVecCoefficientFunction<DIM>>(1, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
  {
    this->elementwise_constant = c1->ElementwiseConstant() && c2->ElementwiseConstant();
    if (DIM != c1->Dimension() || DIM != c2->Dimension())
      throw Exception("T_MultVecVec : dimensions don't fit");
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1).Shallow(c2);
  }

  virtual string GetDescription () const override
  { return "innerproduct, fix size = "+ToString(DIM); }

  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    // code.Declare (code.res_type, index, this->Dimensions());
    code.Declare (index, this->Dimensions(), this->IsComplex());
    CodeExpr result;
    for (int i = 0; i < c1->Dimension(); i++)
      result += Var(inputs[0], i, c1->Dimensions()) * Var(inputs[1], i, c2->Dimensions());
    code.body += Var(index).Assign(result.S(), false);
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1, c2 }); }  

  using T_CoefficientFunction<T_MultVecVecCoefficientFunction<DIM>>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    Vec<DIM> v1, v2;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    Vec<DIM,Complex> v1, v2;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    size_t w = ir.Size();
    __assume (w > 0);
    
    STACK_ARRAY(T, hmem, 2*DIM*w);
    FlatMatrix<T,ORD> temp1(DIM, w, &hmem[0]);
    FlatMatrix<T,ORD> temp2(DIM, w, &hmem[DIM*w]);
    
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, temp2);

    for (size_t i = 0; i < w; i++)
      {
        T sum{0.0};
        for (size_t j = 0; j < DIM; j++)
          sum += temp1(j,i) * temp2(j,i);
        values(0,i) = sum; 
      }
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto in0 = input[0];
    auto in1 = input[1];
    size_t np = ir.Size();

    for (size_t i = 0; i < np; i++)
      {
        T sum{0.0};
        for (size_t j = 0; j < DIM; j++)
          sum += in0(j,i) * in1(j,i);
        values(0,i) = sum; 
      }    
  }  
  
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        BareSliceMatrix<Complex> result) const override
  {
    STACK_ARRAY(double, hmem1, 2*ir.Size()*DIM);
    FlatMatrix<Complex> temp1(ir.Size(), DIM, (Complex*)hmem1);
    STACK_ARRAY(double, hmem2, 2*ir.Size()*DIM);
    FlatMatrix<Complex> temp2(ir.Size(), DIM, (Complex*)hmem2);

    c1->Evaluate(ir, temp1);
    c2->Evaluate(ir, temp2);
    for (int i = 0; i < ir.Size(); i++)
      result(i,0) = InnerProduct(temp1.Row(i), temp2.Row(i));
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return InnerProduct(c1->Diff(var,dir),c2) + InnerProduct(c1,c2->Diff(var,dir));
  }
  
  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    if (this == var)
      return make_shared<ConstantCoefficientFunction> (1);

    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    shared_ptr<CoefficientFunction> dv1v2, dv2v1;
    int dimip = c1->Dimension();
    int dimvar = var->Dimension();
    
    auto vc1 = c1->Reshape( dimip );
    auto vc2 = c2->Reshape( dimip );
       
    if (c1.get() == var)
      dv1v2 = c2;
    else
      {
        auto dvc1 = c1->DiffJacobi (var, cache);
        dv1v2 = dvc1 -> Reshape(dimip, dimvar) -> Transpose() * vc2;
        dv1v2 = dv1v2 -> Reshape (var->Dimensions());
      }

    if (c2.get() == var)
      dv2v1 = c1;
    else
      {
        auto dvc2 = c2->DiffJacobi (var, cache);
        dv2v1 = dvc2 -> Reshape(dimip, dimvar) -> Transpose() * vc1;
        dv2v1 = dv2v1 -> Reshape (var->Dimensions());        
      }
    auto res = dv1v2 + dv2v1;
    cache[thisptr] = res;
    return res;
  }
  
  

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> v1(DIM), v2(DIM);
    c1->NonZeroPattern (ud, v1);
    c2->NonZeroPattern (ud, v2);
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < DIM; i++)
      sum += v1(i)*v2(i);
    values(0) = sum;
  }

  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto v1 = input[0];
    auto v2 = input[1];
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < DIM; i++)
      sum += v1(i)*v2(i);
    values(0) = sum;
  }
};







template <int DIM>
class T_MultVecVecSameCoefficientFunction : public T_CoefficientFunction<T_MultVecVecSameCoefficientFunction<DIM>>
{
  shared_ptr<CoefficientFunction> c1;
  using BASE = T_CoefficientFunction<T_MultVecVecSameCoefficientFunction<DIM>>;
  using typename BASE::T_DJC;
public:
  T_MultVecVecSameCoefficientFunction() = default;
  T_MultVecVecSameCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : T_CoefficientFunction<T_MultVecVecSameCoefficientFunction<DIM>>(1, ac1->IsComplex()), c1(ac1)
  {
    this->elementwise_constant = c1->ElementwiseConstant();
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1);
  }

  virtual string GetDescription () const override
  { return "innerproduct, same vectors, fix size = "+ToString(DIM); }

  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    // code.Declare (code.res_type, index, this->Dimensions());
    code.Declare (index, this->Dimensions(), this->IsComplex()); 
    
    CodeExpr result;
    for (int i = 0; i < c1->Dimension(); i++)
      result += Var(inputs[0], i, c1->Dimensions()) * Var(inputs[0], i, c1->Dimensions());
    code.body += Var(index).Assign(result.S(), false);
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 }); }  

  using T_CoefficientFunction<T_MultVecVecSameCoefficientFunction<DIM>>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    Vec<DIM> v1, v2;
    c1->Evaluate (ip, v1);
    result(0) = InnerProduct (v1, v1);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    Vec<DIM,Complex> v1;
    c1->Evaluate (ip, v1);
    result(0) = InnerProduct (v1, v1);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    size_t w = ir.Size();
    __assume (w > 0);
    
    STACK_ARRAY(T, hmem, DIM*w);
    FlatMatrix<T,ORD> temp1(DIM, w, &hmem[0]);
    
    c1->Evaluate (ir, temp1);

    for (size_t i = 0; i < w; i++)
      {
        T sum{0.0};
        for (size_t j = 0; j < DIM; j++)
          sum += temp1(j,i) * temp1(j,i);
        values(0,i) = sum; 
      }
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto in0 = input[0];
    size_t np = ir.Size();

    for (size_t i = 0; i < np; i++)
      {
        T sum{0.0};
        for (size_t j = 0; j < DIM; j++)
          sum += in0(j,i) * in0(j,i);
        values(0,i) = sum; 
      }    
  }  
  
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        BareSliceMatrix<Complex> result) const override
  {
    STACK_ARRAY(double, hmem1, 2*ir.Size()*DIM);
    FlatMatrix<Complex> temp1(ir.Size(), DIM, (Complex*)hmem1);

    c1->Evaluate(ir, temp1);
    for (int i = 0; i < ir.Size(); i++)
      result(i,0) = InnerProduct(temp1.Row(i), temp1.Row(i));
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return 2*InnerProduct(c1->Diff(var,dir),c1);
  }

  
  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];
    
    if (this == var)
      return make_shared<ConstantCoefficientFunction> (1);

    shared_ptr<CoefficientFunction> dv1v1;
    int dimip = c1->Dimension();
    int dimvar = var->Dimension();
    
    auto vc1 = c1->Reshape( dimip );
       
    if (c1.get() == var)
      dv1v1 = c1;
    else
      {
        auto dvc1 = c1->DiffJacobi (var, cache);
        dv1v1 = dvc1 -> Reshape(dimip, dimvar) -> Transpose() * vc1;
        dv1v1 = dv1v1 -> Reshape (var->Dimensions());
      }
    auto res = 2*dv1v1;
    cache[thisptr] = res;
    return res;
  }
  

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> v1(DIM);
    c1->NonZeroPattern (ud, v1);
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < DIM; i++)
      sum += v1(i)*v1(i);
    values(0) = sum;
  }

  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto v1 = input[0];
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < DIM; i++)
      sum += v1(i)*v1(i);
    values(0) = sum;
  }
};










class EigCoefficientFunction : public CoefficientFunctionNoDerivative
{
  shared_ptr<CoefficientFunction> cfmat;
  int dim1;
  int vecdim;
  
public:
  EigCoefficientFunction() = default;
  EigCoefficientFunction (shared_ptr<CoefficientFunction> ac1) : CoefficientFunctionNoDerivative(ac1->Dimension() + ac1->Dimensions()[0],false), cfmat(ac1)
  {
    vecdim = cfmat->Dimensions()[0];
    dim1 = cfmat->Dimension();
  }

  void DoArchive(Archive& ar) override
  {
    CoefficientFunctionNoDerivative::DoArchive(ar);
    ar.Shallow(cfmat) & dim1 & vecdim;
  }
  
  using CoefficientFunctionNoDerivative::Evaluate;
  double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    return 0;
  }
  void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> res) const override
  {
    STACK_ARRAY(double,mem, dim1);
    FlatVector<double> vec(dim1, &mem[0]);
    
    cfmat->Evaluate (ip, vec);

    FlatMatrix<double> mat(vecdim,vecdim, &mem[0]);
    FlatVector<double> lami(vecdim, &res[dim1]);
    FlatMatrix<double> eigenvecs(vecdim,vecdim,&res[0]);
    
    CalcEigenSystem(mat,lami,eigenvecs);
  }
};



class NormCoefficientFunction : public T_CoefficientFunction<NormCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  int dim1;
  typedef double TIN;
  using BASE = T_CoefficientFunction<NormCoefficientFunction>;
public:
  NormCoefficientFunction() = default;
  NormCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : T_CoefficientFunction<NormCoefficientFunction> (1, false), c1(ac1)
  {
    dim1 = c1->Dimension();
    elementwise_constant = c1->ElementwiseConstant();
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1) & dim1;
  }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 }); }  
  
    using T_CoefficientFunction<NormCoefficientFunction>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    VectorMem<10,TIN> v1(dim1);
    c1->Evaluate (ip, v1);
    result(0) = L2Norm(v1);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    Vec<1> res;
    Evaluate (ip, res);
    result(0) = res(0);
  }


  /*
  virtual bool ElementwiseConstant () const override
  { return c1->ElementwiseConstant(); }
  */

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    size_t np = ir.Size();
    size_t dim1 = c1->Dimension();
    STACK_ARRAY(T,mem, np*dim1);
    FlatMatrix<T,ORD> m1(dim1, np, &mem[0]);
    c1->Evaluate (ir, m1);
    
    for (size_t i = 0; i < np; i++)
      {
        T sum{0.0};
        for (size_t j = 0; j < dim1; j++)
          sum += sqr(m1(j,i));
        values(0,i) = sqrt(sum);
      }
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    size_t np = ir.Size();
    auto in = input[0];
    for (size_t i = 0; i < np; i++)
      {
        T sum{0.0};
        for (size_t j = 0; j < dim1; j++)
          sum += sqr(in(j,i));
        values(0,i) = sqrt(sum);
      }
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    auto res = CodeExpr();
    for (int i = 0; i < c1->Dimension(); i++)
      res += Var(inputs[0],i,c1->Dimensions()).Func("L2Norm2");

    // code.Declare (code.res_type, index, this->Dimensions());
    code.Declare (index, this->Dimensions(), this->IsComplex());
    code.body += Var(index).Assign( res.Func("sqrt"), false);
  }

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    Vector<bool> v1(dim1), d1(dim1), dd1(dim1);
    c1->NonZeroPattern (ud, v1, d1, dd1);
    bool nz = false, nzd = false, nzdd = false;
    for (int i = 0; i < dim1; i++)
      {
        if (v1(i)) nz = true;
        if (d1(i)) nzd = true;
        if (dd1(i)) nzdd = true;
      }
    nonzero = nz;
    nonzero_deriv = nzd;
    nonzero_dderiv = nzd || nzdd;
  }
  */

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> v1(dim1);
    c1->NonZeroPattern (ud, v1);

    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < dim1; i++)
      sum += v1(i);
    values(0).Value() = sum.Value();
    values(0).DValue(0) = sum.DValue(0);
    values(0).DDValue(0) = sum.DValue(0) || sum.DDValue(0);
  }

  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto v1 = input[0];
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < dim1; i++)
      sum += v1(i);
    values(0).Value() = sum.Value();
    values(0).DValue(0) = sum.DValue(0);
    values(0).DDValue(0) = sum.DValue(0) || sum.DDValue(0);
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    if (var == this) return dir;
    return make_shared<ConstantCoefficientFunction>(1.0)/NormCF(c1) * InnerProduct(c1,c1->Diff(var,dir));
  }


  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    if (this == var)
      return make_shared<ConstantCoefficientFunction> (1);

    shared_ptr<CoefficientFunction> dv1v1;
    int dimip = c1->Dimension();
    int dimvar = var->Dimension();
    
    auto vc1 = c1->Reshape( dimip );
       
    if (c1.get() == var)
      dv1v1 = c1;
    else
      {
        auto dvc1 = c1->DiffJacobi (var, cache);
        dv1v1 = dvc1 -> Reshape(dimip, dimvar) -> Transpose() * vc1;
        dv1v1 = dv1v1 -> Reshape (var->Dimensions());
      }
    auto res = (1.0 / const_cast<CoefficientFunction*>((CoefficientFunction*)this)->shared_from_this()) * dv1v1;
    cache[thisptr] = res;
    return res;
  }
  

  
};




class NormCoefficientFunctionC : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  int dim1;
  typedef Complex TIN;
public:
  NormCoefficientFunctionC() = default;
  NormCoefficientFunctionC (shared_ptr<CoefficientFunction> ac1)
    : CoefficientFunction (1, false), c1(ac1)
  {
    dim1 = c1->Dimension();
    elementwise_constant = c1->ElementwiseConstant(); 
  }

  void DoArchive(Archive& ar) override
  {
    CoefficientFunction::DoArchive(ar);
    ar.Shallow(c1) & dim1;
  }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 }); }  
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    VectorMem<10,TIN> v1(dim1);
    c1->Evaluate (ip, v1);
    result(0) = L2Norm(v1);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    Vec<1> res;
    Evaluate (ip, res);
    result(0) = res(0);
  }


  /*
  virtual bool ElementwiseConstant () const override
  { return c1->ElementwiseConstant(); }
  */
  
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        BareSliceMatrix<> result) const override
  {
    STACK_ARRAY(double,hmem,ir.Size()*dim1*sizeof(TIN)/sizeof(double));
    FlatMatrix<TIN> inval(ir.IR().GetNIP(), dim1, reinterpret_cast<TIN*>(&hmem[0]));
    c1->Evaluate (ir, inval);
    for (size_t i = 0; i < ir.Size(); i++)
      result(i,0) = L2Norm(inval.Row(i));
  }


  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override
  {
    STACK_ARRAY(SIMD<Complex>,hmem,ir.Size()*dim1);
    FlatMatrix<SIMD<Complex>> inval(dim1, ir.Size(), &hmem[0]);
    c1->Evaluate (ir, inval);
    for (size_t i = 0; i < ir.Size(); i++)
      {
        SIMD<double> sum = 0;
        for (size_t j = 0; j < dim1; j++)
          sum += sqr(inval(j,i).real())+sqr(inval(j,i).imag());
        values(0,i) = sqrt(sum);
      }
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    auto res = CodeExpr();
    for (int i = 0; i < c1->Dimension(); i++)
      res += Var(inputs[0],i,c1->Dimensions()).Func("L2Norm2");
    code.body += Var(index).Assign( res.Func("sqrt"));
  }

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    Vector<bool> v1(dim1), d1(dim1), dd1(dim1);
    c1->NonZeroPattern (ud, v1, d1, dd1);
    bool nz = false, nzd = false, nzdd = false;
    for (int i = 0; i < dim1; i++)
      {
        if (v1(i)) nz = true;
        if (d1(i)) nzd = true;
        if (dd1(i)) nzdd = true;
      }
    nonzero = nz;
    nonzero_deriv = nzd;
    nonzero_dderiv = nzd || nzdd;
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> v1(dim1);
    c1->NonZeroPattern (ud, v1);
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < dim1; i++)
      sum = sum + v1(i);
    values(0) = sum;
  }
};

  

class MultMatMatCoefficientFunction : public T_CoefficientFunction<MultMatMatCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
  int inner_dim;
  using BASE = T_CoefficientFunction<MultMatMatCoefficientFunction>;
public:
  MultMatMatCoefficientFunction() = default;
  MultMatMatCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                 shared_ptr<CoefficientFunction> ac2)
    : T_CoefficientFunction<MultMatMatCoefficientFunction>(1, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
  {
    auto dims_c1 = c1 -> Dimensions();
    auto dims_c2 = c2 -> Dimensions();
    if (dims_c1.Size() != 2 || dims_c2.Size() != 2)
      throw Exception("Mult of non-matrices called");
    if (dims_c1[1] != dims_c2[0])
      throw Exception(string("Matrix dimensions don't fit: m1 is ") +
                      ToLiteral(dims_c1[0]) + " x " + ToLiteral(dims_c1[1]) +
                      ", m2 is " + ToLiteral(dims_c2[0]) + " x " + ToLiteral(dims_c2[1]) );
    SetDimensions( ngstd::INT<2> (dims_c1[0], dims_c2[1]) );
    inner_dim = dims_c1[1];
  }

  virtual string GetDescription () const override
  { return "matrix-matrix multiply"; }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1).Shallow(c2) & inner_dim;
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
    FlatArray<int> hdims = Dimensions();
    // code.Declare (code.res_type, index, hdims);
    code.Declare (index, hdims, IsComplex());

    if (code_uses_tensors)
      {
        code.body += "for (size_t i = 0; i < "+ToString(hdims[0])+"; i++)\n";
        code.body += "for (size_t j = 0; j < "+ToString(hdims[1])+"; j++) { \n";
        code.body += "auto sum = var_" + ToString(inputs[0]) + "(i,0) * var_" + ToString(inputs[1]) + "(0,j); \n";
        code.body += "for (size_t k = 1; k < "+ToString(inner_dim)+"; k++) \n";
        code.body += "sum += var_" + ToString(inputs[0]) + "(i,k) * var_" + ToString(inputs[1]) + "(k,j); \n";
        code.body += "var_" + ToString(index) + "(i,j) = sum; } \n";
      }
    else
      {
        for (int i : Range(hdims[0]))
          for (int j : Range(hdims[1])) {
            CodeExpr s;
            for (int k : Range(inner_dim))
              s += Var(inputs[0], i, k) * Var(inputs[1], k, j);
            code.body += Var(index, i, j).Assign(s, false);
          }
      }
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1, c2 }); }  


  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    FlatArray<int> hdims = Dimensions();
    Vector<bool> v1(hdims[0]*inner_dim), v2(hdims[1]*inner_dim);
    Vector<bool> d1(hdims[0]*inner_dim), d2(hdims[1]*inner_dim);
    Vector<bool> dd1(hdims[0]*inner_dim), dd2(hdims[1]*inner_dim);
    c1->NonZeroPattern (ud, v1, d1, dd1);
    c2->NonZeroPattern (ud, v2, d2, dd2);
    nonzero = false;
    nonzero_deriv = false;
    nonzero_dderiv = false;
    FlatMatrix<bool> m1(hdims[0], inner_dim, &v1(0));
    FlatMatrix<bool> m2(inner_dim, hdims[1], &v2(0));
    FlatMatrix<bool> md1(hdims[0], inner_dim, &d1(0));
    FlatMatrix<bool> md2(inner_dim, hdims[1], &d2(0));
    FlatMatrix<bool> mdd1(hdims[0], inner_dim, &dd1(0));
    FlatMatrix<bool> mdd2(inner_dim, hdims[1], &dd2(0));

    for (int i = 0; i < hdims[0]; i++)
      for (int j = 0; j < hdims[1]; j++)
        for (int k = 0; k < inner_dim; k++)
          {
            nonzero(i*hdims[1]+j) |= m1(i,k) && m2(k,j);
            nonzero_deriv(i*hdims[1]+j) |= (m1(i,k) && md2(k,j)) || (md1(i,k) && m2(k,j));
            nonzero_dderiv(i*hdims[1]+j) |= (m1(i,k) && mdd2(k,j)) || (md1(i,k) && md2(k,j)) || (mdd1(i,k) && m2(k,j));
          }
  }
  */

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    FlatArray<int> hdims = Dimensions();
    Vector<AutoDiffDiff<1,bool>> va(hdims[0]*inner_dim), vb(hdims[1]*inner_dim);
    c1->NonZeroPattern (ud, va);
    c2->NonZeroPattern (ud, vb);
    
    size_t d1 = hdims[1];

    values = false;
    
    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < hdims[1]; k++)
        for (size_t l = 0; l < inner_dim; l++)
          values(j*d1+k) += va(j*inner_dim+l) * vb(l*d1+k);
  }

  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto va = input[0];
    auto vb = input[1];

    FlatArray<int> hdims = Dimensions();    
    size_t d1 = hdims[1];

    values = false;
    
    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < hdims[1]; k++)
        for (size_t l = 0; l < inner_dim; l++)
          values(j*d1+k) += va(j*inner_dim+l) * vb(l*d1+k);
  }

    using T_CoefficientFunction<MultMatMatCoefficientFunction>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    throw Exception ("MultMatMatCF:: scalar evaluate for matrix called");
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const override
  {
    FlatArray<int> hdims = Dimensions();
    Vector<> va(hdims[0]*inner_dim);
    Vector<> vb(hdims[1]*inner_dim);
    FlatMatrix<> a(hdims[0], inner_dim, &va[0]);
    FlatMatrix<> b(inner_dim, hdims[1], &vb[0]);
    
    c1->Evaluate (ip, va);
    c2->Evaluate (ip, vb);

    FlatMatrix<> c(hdims[0], hdims[1], &result(0));
    c = a*b;
  }  

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<Complex> result) const override
  {
    FlatArray<int> hdims = Dimensions();
    STACK_ARRAY(double,mema,2*hdims[0]*inner_dim);
    STACK_ARRAY(double,memb,2*hdims[1]*inner_dim);
    FlatVector<Complex> va(hdims[0]*inner_dim, reinterpret_cast<Complex*>(&mema[0]));
    FlatVector<Complex> vb(inner_dim*hdims[1], reinterpret_cast<Complex*>(&memb[0]));
    
    c1->Evaluate (ip, va);
    c2->Evaluate (ip, vb);

    FlatMatrix<Complex> a(hdims[0], inner_dim, &va(0));
    FlatMatrix<Complex> b(inner_dim, hdims[1], &vb(0));

    FlatMatrix<Complex> c(hdims[0], hdims[1], &result(0));
    c = a*b;
    //cout << "MultMatMat: complex not implemented" << endl;
  }  

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & mir, BareSliceMatrix<T,ORD> values) const 
  {
    FlatArray<int> hdims = Dimensions();    
    STACK_ARRAY(T, hmem1, mir.Size()*hdims[0]*inner_dim);
    STACK_ARRAY(T, hmem2, mir.Size()*hdims[1]*inner_dim);
    FlatMatrix<T,ORD> va(hdims[0]*inner_dim, mir.Size(), &hmem1[0]);
    FlatMatrix<T,ORD> vb(hdims[1]*inner_dim, mir.Size(), &hmem2[0]);

    c1->Evaluate (mir, va);
    c2->Evaluate (mir, vb);

    values.AddSize(Dimension(),mir.Size()) = T(0.0);

    size_t d1 = hdims[1];
    size_t mir_size = mir.Size();
    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < hdims[1]; k++)
        for (size_t l = 0; l < inner_dim; l++)
          {
            auto row_a = va.Row(j*inner_dim+l);
            auto row_b = vb.Row(l*d1+k);
            auto row_c = values.Row(j*d1+k);
            for (size_t i = 0; i < mir_size; i++)
              row_c(i) += row_a(i) * row_b(i);
            // row_c = pw_mult (row_a, row_b);
          }
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto va = input[0];
    auto vb = input[1];

    FlatArray<int> hdims = Dimensions();    
    size_t d1 = hdims[1];
    size_t np = ir.Size();

    values.AddSize(Dimension(),np) = T(0.0);
    
    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < hdims[1]; k++)
        for (size_t l = 0; l < inner_dim; l++)
          {
            auto row_a = va.Row(j*inner_dim+l);
            auto row_b = vb.Row(l*d1+k);
            auto row_c = values.Row(j*d1+k);
            for (size_t i = 0; i < np; i++)
              row_c(i) += row_a(i) * row_b(i);
            // row_c = pw_mult (row_a, row_b);
          }
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    if (var == this) return dir;
    return c1->Diff(var,dir)*c2 + c1 * c2->Diff(var,dir);
  }
  
  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {    
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    int h = Dimensions()[0];
    int w = Dimensions()[1];
    int dimvar = var->Dimension();

    if (this == var)
      return IdentityCF(this->Dimensions());

    Array<int> dimres{h,w};
    dimres += var->Dimensions();

    auto diffc1 = c1->DiffJacobi (var, cache);
    auto diffc2 = c2->DiffJacobi (var, cache);
    
    auto diffc1_trans = diffc1 -> TensorTranspose( 0, 1 ) -> Reshape( inner_dim, h*dimvar );
    auto prod1 = c2->Transpose() * diffc1_trans;
    Array<int> dimtmp{w, h};
    dimtmp += var->Dimensions();
    auto prod1trans = prod1->Reshape(dimtmp) -> TensorTranspose( 0, 1 );
    
    auto diffc2_trans = diffc2 -> Reshape( Array<int>{inner_dim,w*dimvar} );
    auto prod2 = c1 * diffc2_trans;
    auto prod2trans = prod2->Reshape(dimres);

    auto res = prod1trans + prod2trans;
    cache[thisptr] = res;
    return res;
  }
};







class MultMatVecCoefficientFunction : public T_CoefficientFunction<MultMatVecCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
  // Array<int> dims;
  int inner_dim;
  using BASE = T_CoefficientFunction<MultMatVecCoefficientFunction>;
public:
  MultMatVecCoefficientFunction() = default;
  MultMatVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                 shared_ptr<CoefficientFunction> ac2)
    : T_CoefficientFunction(1, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
  {
    auto dims_c1 = c1 -> Dimensions();
    auto dims_c2 = c2 -> Dimensions();
    if (dims_c1.Size() < 2 || dims_c2.Size() != 1)
      throw Exception("Not a mat-vec multiplication");
    if (dims_c1.Last() != dims_c2[0])
      throw Exception(string ("Matrix dimensions don't fit: mat is ") +
                      ToLiteral(dims_c1[0]) + " x " + ToLiteral(dims_c1[1]) + ", vec is " + ToLiteral(dims_c2[0]));
    // SetDimensions (ngstd::INT<1>(dims_c1[0]));
    // SetDimensions (dims_c1.Range(0, dims_c1.Size()-1)); // ngstd::INT<1>(dims_c1[0]));
    SetDimensions (dims_c1.Range(0, END-1));
    inner_dim = dims_c1.Last(); // [1];
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1).Shallow(c2) & inner_dim;
  }

  virtual string GetDescription () const override
  { return "matrix-vector multiply"; }

  
  // virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  // virtual int Dimension() const { return dims[0]; }
  // virtual Array<int> Dimensions() const { return Array<int> (dims); } 

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1, c2 }); }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
    // code.Declare (code.res_type, index, Dimensions());
    code.Declare (index, Dimensions(), IsComplex()); 
      auto dims = c1->Dimensions();
      for (int i : Range(dims[0])) {
        CodeExpr s;
        for (int j : Range(dims[1]))
            s += Var(inputs[0], i, j) * Var(inputs[1], j);
	code.body += Var(index, i).Assign(s, false);
      }
  }

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    FlatArray<int> hdims = Dimensions();
    Vector<bool> v1(hdims[0]*inner_dim), v2(inner_dim);
    Vector<bool> d1(hdims[0]*inner_dim), d2(inner_dim);
    Vector<bool> dd1(hdims[0]*inner_dim), dd2(inner_dim);
    c1->NonZeroPattern (ud, v1, d1, dd1);
    c2->NonZeroPattern (ud, v2, d2, dd2);
    nonzero = false;
    nonzero_deriv = false;
    nonzero_dderiv = false;
    FlatMatrix<bool> m1(hdims[0], inner_dim, &v1(0));
    FlatMatrix<bool> md1(hdims[0], inner_dim, &d1(0));
    FlatMatrix<bool> mdd1(hdims[0], inner_dim, &dd1(0));
    for (int i = 0; i < hdims[0]; i++)
      for (int j = 0; j < inner_dim; j++)
        {
          nonzero(i) |= (m1(i,j) && v2(j));
          nonzero_deriv(i) |= ((m1(i,j) && d2(j)) || (md1(i,j) && v2(j)));
          nonzero_dderiv(i) |= ((m1(i,j) && dd2(j)) || (md1(i,j) && d2(j)) || (mdd1(i,j) && v2(j)));
        }
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    // FlatArray<int> hdims = Dimensions();
    Vector<AutoDiffDiff<1,bool>> va(Dimension()*inner_dim), vb(inner_dim);
    c1->NonZeroPattern (ud, va);
    c2->NonZeroPattern (ud, vb);
    
    values = false;
    
    for (size_t i = 0; i < Dimension(); i++)
      for (size_t j = 0; j < inner_dim; j++)
        values(i) += va(i*inner_dim+j) * vb(j);
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto va = input[0];
    auto vb = input[1];
    
    // FlatArray<int> hdims = Dimensions();    
    values = false;
    
    for (size_t i = 0; i < Dimension(); i++)
      for (size_t j = 0; j < inner_dim; j++)
        values(i) += va(i*inner_dim+j) * vb(j);
  }
    using T_CoefficientFunction<MultMatVecCoefficientFunction>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    throw Exception ("MultMatVecCF:: scalar evaluate for matrix called");
  }

  /*
  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const override
  {
    FlatArray<int> hdims = Dimensions();
    VectorMem<20> va(hdims[0]*inner_dim);
    VectorMem<20> vb(inner_dim);
    FlatMatrix<> a(hdims[0], inner_dim, &va[0]);

    c1->Evaluate (ip, va);
    c2->Evaluate (ip, vb);

    result = a * vb;
  }  

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<Complex> result) const override
  {
    FlatArray<int> hdims = Dimensions();
    STACK_ARRAY(double,mema,2*hdims[0]*inner_dim);
    STACK_ARRAY(double,memb,2*inner_dim);
    FlatVector<Complex> va(hdims[0]*inner_dim,reinterpret_cast<Complex*>(&mema[0]));
    FlatVector<Complex> vb(inner_dim,reinterpret_cast<Complex*>(&memb[0]));
    FlatMatrix<Complex> a(hdims[0], inner_dim, &va(0));

    c1->Evaluate (ip, va);
    c2->Evaluate (ip, vb);

    result = a * vb;
    //cout << "MultMatMat: complex not implemented" << endl;
  }  
  */

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    // FlatArray<int> hdims = Dimensions();
    int dim = Dimension();
    STACK_ARRAY(T, hmem1, ir.Size()*dim*inner_dim);
    STACK_ARRAY(T, hmem2, ir.Size()*inner_dim);
    FlatMatrix<T,ORD> temp1(dim*inner_dim, ir.Size(), &hmem1[0]);
    FlatMatrix<T,ORD> temp2(inner_dim, ir.Size(), &hmem2[0]);
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, temp2);
    values.AddSize(Dimension(),ir.Size()) = T(0.0);
    for (size_t i = 0; i < dim; i++)
      for (size_t j = 0; j < inner_dim; j++)
        for (size_t k = 0; k < ir.Size(); k++)
          values(i,k) += temp1(i*inner_dim+j, k) * temp2(j,k);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto va = input[0];
    auto vb = input[1];
    
    // FlatArray<int> hdims = Dimensions();
    int dim = Dimension();
    values.AddSize(Dimension(),ir.Size()) = T(0.0);
    
    for (size_t i = 0; i < dim; i++)
      for (size_t j = 0; j < inner_dim; j++)
        for (size_t k = 0; k < ir.Size(); k++)
          values(i,k) += va(i*inner_dim+j, k) * vb(j,k);
  }


  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return c1->Diff(var,dir)*c2 + c1 * c2->Diff(var,dir);
  }
  
  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    int h = Dimensions()[0];
    int dimvar = var->Dimension();

    if (this == var)
      return IdentityCF(h);

    Array<int> dimres{h};
    dimres += var->Dimensions();

    auto diffc1 = c1->DiffJacobi (var, cache);
    auto diffc2 = c2->DiffJacobi (var, cache);
    
    auto diffc1_trans = diffc1 -> TensorTranspose(0,1) -> Reshape(inner_dim, h*dimvar) -> Transpose();
    auto prod1 = (diffc1_trans * c2) -> Reshape(h,dimvar);
    auto prod1trans = prod1 -> Reshape(dimres);

    auto diffc2_trans = diffc2 -> Reshape(inner_dim, dimvar);
    auto prod2 = c1 * diffc2_trans;
    auto prod2trans = prod2 -> Reshape(dimres);

    auto res = prod1trans + prod2trans;
    cache[thisptr] = res;
    return res;
  }
  
};



class CrossProductCoefficientFunction : public T_CoefficientFunction<CrossProductCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
  using BASE = T_CoefficientFunction<CrossProductCoefficientFunction>;
public:
  CrossProductCoefficientFunction() = default;
  CrossProductCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                   shared_ptr<CoefficientFunction> ac2)
    : T_CoefficientFunction(3, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
  {
    if (c1->Dimension() != 3)
      throw Exception("first factor of cross product does not have dim=3");
    if (c2->Dimension() != 3)
      throw Exception("second factor of cross product does not have dim=3");
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1).Shallow(c2);
  }

  virtual string GetDescription () const override
  { return "cross-product"; }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1, c2 }); }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
    // code.Declare (code.res_type, index, this->Dimensions());
    code.Declare (index, this->Dimensions(), this->IsComplex());
    code.body += Var(index, 0).Assign (Var(inputs[0],1)*Var(inputs[1],2)-Var(inputs[0],2)*Var(inputs[1],1), false);
    code.body += Var(index, 1).Assign (Var(inputs[0],2)*Var(inputs[1],0)-Var(inputs[0],0)*Var(inputs[1],2), false);
    code.body += Var(index, 2).Assign (Var(inputs[0],0)*Var(inputs[1],1)-Var(inputs[0],1)*Var(inputs[1],0), false);
  }

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> va(3), vb(3);
    c1->NonZeroPattern (ud, va);
    c2->NonZeroPattern (ud, vb);
    
    values(0) = va(1)*vb(2)+va(2)*vb(1);
    values(1) = va(2)*vb(0)+va(0)*vb(2);
    values(2) = va(0)*vb(1)+va(1)*vb(0);    
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto va = input[0];
    auto vb = input[1];
    
    values(0) = va(1)*vb(2)+va(2)*vb(1);
    values(1) = va(2)*vb(0)+va(0)*vb(2);
    values(2) = va(0)*vb(1)+va(1)*vb(0);    
  }
  
  using T_CoefficientFunction<CrossProductCoefficientFunction>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    throw Exception ("MultMatVecCF:: scalar evaluate for vector called");
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const override
  {
    Vec<3> va, vb;
    c1->Evaluate (ip, va);
    c2->Evaluate (ip, vb);
    result = Cross(va, vb);
  }  

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<Complex> result) const override
  {
    Vec<3,Complex> va, vb;
    c1->Evaluate (ip, va);
    c2->Evaluate (ip, vb);
    result = Cross(va, vb);
  }  


  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    STACK_ARRAY(T, hmem1, ir.Size()*3);
    STACK_ARRAY(T, hmem2, ir.Size()*3);
    FlatMatrix<T,ORD> temp1(3, ir.Size(), &hmem1[0]);
    FlatMatrix<T,ORD> temp2(3, ir.Size(), &hmem2[0]);
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, temp2);
    for (size_t k = 0; k < ir.Size(); k++)
      {
        Vec<3,T> ai = temp1.Col(k);
        Vec<3,T> bi = temp2.Col(k);
        values.Col(k).Range(3) = Cross(ai, bi);
      }
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto va = input[0];
    auto vb = input[1];

    for (size_t k = 0; k < ir.Size(); k++)
      {
        Vec<3,T> ai = va.Col(k);
        Vec<3,T> bi = vb.Col(k);
        values.Col(k).Range(3) = Cross(ai, bi);
      }
  }


  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return CrossProduct(c1->Diff(var,dir),c2) + CrossProduct(c1, c2->Diff(var,dir));
  }


  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, typename BASE::T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];
    
    int h = c1->Dimension();
    int dimvar = var->Dimension();
    
    if (this == var)
      return IdentityCF(h);

    Array<int> dimres{h};
    dimres += var->Dimensions();

    shared_ptr<CoefficientFunction> res;
    if (false) //first version, a little bit slower
      {
        auto diffc1 = c1 -> DiffJacobi (var, cache);
        auto diffc2 = c2 -> DiffJacobi (var, cache);
        
        Array<shared_ptr<CoefficientFunction>> cross1(dimvar);
        Array<shared_ptr<CoefficientFunction>> cross2(dimvar);
        for (auto i : Range(dimvar))
          {
            cross1[i] = CrossProduct(MakeSubTensorCoefficientFunction(diffc1,i,Array<int>({3}),Array<int>({dimvar})),c2);
            cross2[i] = CrossProduct(c1,MakeSubTensorCoefficientFunction(diffc2,i,Array<int>({3}),Array<int>({dimvar})));
          }
        
        auto dcross1 = MakeVectorialCoefficientFunction (move(cross1))->Reshape(dimvar,h) -> Transpose() -> Reshape(dimres);
        auto dcross2 = MakeVectorialCoefficientFunction (move(cross2))->Reshape(dimvar,h) -> Transpose() -> Reshape(dimres);
    
        res = dcross1 + dcross2;
      }
    else //second version, a little bit faster
      {
        auto c11 = MakeComponentCoefficientFunction(c1,0);
        auto c12 = MakeComponentCoefficientFunction(c1,1);
        auto c13 = MakeComponentCoefficientFunction(c1,2);
        auto c21 = MakeComponentCoefficientFunction(c2,0);
        auto c22 = MakeComponentCoefficientFunction(c2,1);
        auto c23 = MakeComponentCoefficientFunction(c2,2);

        Array<shared_ptr<CoefficientFunction>> cross(3);
        cross[0] = c12*c23 - c13*c22;
        cross[1] = c13*c21 - c11*c23;
        cross[2] = c11*c22 - c12*c21;

        res = MakeVectorialCoefficientFunction(move(cross)) -> DiffJacobi (var, cache);
        
      }
    cache[thisptr] = res;
    return res;
  }
  
};



class TransposeCoefficientFunction : public T_CoefficientFunction<TransposeCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  using BASE = T_CoefficientFunction<TransposeCoefficientFunction>;
public:
  TransposeCoefficientFunction() = default;
  TransposeCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : T_CoefficientFunction<TransposeCoefficientFunction>(1, ac1->IsComplex()), c1(ac1)
  {
    auto dims_c1 = c1 -> Dimensions();
    if (dims_c1.Size() != 2)
      throw Exception("Transpose of non-matrix called");

    SetDimensions (ngstd::INT<2> (dims_c1[1], dims_c1[0]) );
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1);
  }

  virtual string GetDescription () const override
  { return "Matrix transpose"; }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      FlatArray<int> hdims = Dimensions();
      // code.Declare (code.res_type, index, this->Dimensions());
      code.Declare (index, this->Dimensions(), this->IsComplex());
      
      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1]))
          code.body += Var(index,i,j).Assign( Var(inputs[0],j,i), false );
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    FlatArray<int> hdims = Dimensions();    
    Vector<bool> v1(hdims[0]*hdims[1]), d1(hdims[0]*hdims[1]), dd1(hdims[0]*hdims[1]);
    c1->NonZeroPattern (ud, v1, d1, dd1);
    {
      FlatMatrix<bool> m1(hdims[1], hdims[0], &v1(0));
      FlatMatrix<bool> m2(hdims[0], hdims[1], &nonzero(0));
      m2 = Trans(m1);
    }
    {
      FlatMatrix<bool> m1(hdims[1], hdims[0], &d1(0));
      FlatMatrix<bool> m2(hdims[0], hdims[1], &nonzero_deriv(0));
      m2 = Trans(m1);
    }
    {
      FlatMatrix<bool> m1(hdims[1], hdims[0], &dd1(0));
      FlatMatrix<bool> m2(hdims[0], hdims[1], &nonzero_dderiv(0));
      m2 = Trans(m1);
    }
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    FlatArray<int> hdims = Dimensions();    
    Vector<AutoDiffDiff<1,bool>> v1(hdims[0]*hdims[1]);
    c1->NonZeroPattern (ud, v1);

    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < hdims[1]; k++)
        values(j*hdims[1]+k) = v1(k*hdims[0]+j);
  }

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    FlatArray<int> hdims = Dimensions();    
    auto in0 = input[0];
    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < hdims[1]; k++)
        values(j*hdims[1]+k) = in0(k*hdims[0]+j);
  }
    using T_CoefficientFunction<TransposeCoefficientFunction>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    throw Exception ("TransposeCF:: scalar evaluate for matrix called");
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const override
  {
    FlatArray<int> hdims = Dimensions();        
    VectorMem<20> input(result.Size());
    c1->Evaluate (ip, input);    
    FlatMatrix<> reshape1(hdims[1], hdims[0], &input(0));  // source matrix format
    FlatMatrix<> reshape2(hdims[0], hdims[1], &result(0));  // range matrix format
    reshape2 = Trans(reshape1);
    
    /*
    c1->Evaluate (ip, result);
    static Timer t("Transpose - evaluate");
    RegionTimer reg(t);
    FlatMatrix<> reshape(dims[1], dims[0], &result(0));  // source matrix format
    Matrix<> tmp = Trans(reshape);
    FlatMatrix<> reshape2(dims[0], dims[1], &result(0));  // range matrix format
    reshape2 = tmp;
    */
  }  

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<Complex> result) const override
  {
    FlatArray<int> hdims = Dimensions();        
    STACK_ARRAY(double,meminput,2*hdims[0]*hdims[1]);
    FlatVector<Complex> input(hdims[0]*hdims[1],reinterpret_cast<Complex*>(&meminput[0]));
    c1->Evaluate (ip, input);    
    FlatMatrix<Complex> reshape1(hdims[1], hdims[0], &input(0));  // source matrix format
    FlatMatrix<Complex> reshape2(hdims[0], hdims[1], &result(0));  // range matrix format
    reshape2 = Trans(reshape1);
    //cout << "Transpose: complex not implemented" << endl;
  }  

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & mir,
                   BareSliceMatrix<T,ORD> result) const
  {
    FlatArray<int> hdims = Dimensions();    
    c1->Evaluate (mir, result);
    STACK_ARRAY(T, hmem, hdims[0]*hdims[1]);
    FlatMatrix<T,ORD> tmp (hdims[0], hdims[1], &hmem[0]);

    for (size_t i = 0; i < mir.Size(); i++)
      {
        for (int j = 0; j < hdims[0]; j++)
          for (int k = 0; k < hdims[1]; k++)
            tmp(j,k) = result(k*hdims[0]+j, i);
        for (int j = 0; j < hdims[0]; j++)
          for (int k = 0; k < hdims[1]; k++)
            result(j*hdims[1]+k, i) = tmp(j,k);
      }
  }  

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    FlatArray<int> hdims = Dimensions();
    size_t np = ir.Size();
    
    auto in0 = input[0];
    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < hdims[1]; k++)
        for (size_t i = 0; i < np; i++)
          values(j*hdims[1]+k, i) = in0(k*hdims[0]+j, i);
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return TransposeCF (c1->Diff(var, dir));
  }

  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    if (this == var)
      return IdentityCF(this->Dimensions());
  
    auto diffc1 = c1->DiffJacobi (var, cache);
    auto res = diffc1 -> TensorTranspose( 0, 1 );
    cache[thisptr] = res;
    return res;
  }
};



template <int D>
class InverseCoefficientFunction : public T_CoefficientFunction<InverseCoefficientFunction<D>>
{
  shared_ptr<CoefficientFunction> c1;
  using BASE = T_CoefficientFunction<InverseCoefficientFunction<D>>;
  using typename BASE::T_DJC;
public:
  InverseCoefficientFunction() = default;
  InverseCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : T_CoefficientFunction<InverseCoefficientFunction>(D*D, ac1->IsComplex()), c1(ac1)
  {
    this->SetDimensions (ngstd::INT<2> (D,D));
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1);
  }

  virtual string GetDescription () const override
  { return "inverse"; }

  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
    // auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.res_type+">";
    auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.GetType(this->IsComplex())+">";
    auto mat_var = Var("mat", index);
    auto inv_var = Var("inv", index);
    code.body += mat_var.Declare(mat_type);
    code.body += inv_var.Declare(mat_type);
    for (int j = 0; j < D; j++)
      for (int k = 0; k < D; k++)
        code.body += mat_var(j,k).Assign(Var(inputs[0], j, k), false);

    code.body += inv_var.Assign(mat_var.Func("Inv"), false);

    for (int j = 0; j < D; j++)
      for (int k = 0; k < D; k++)
        code.body += Var(index, j, k).Assign(inv_var(j,k));
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    nonzero = true;
    nonzero_deriv = true;
    nonzero_dderiv = true;
  }
  */

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> v1(c1->Dimension());
    c1->NonZeroPattern (ud, v1);
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < v1.Size(); i++)
      sum += v1(i);
    values = sum;
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto v1 = input[0];
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < v1.Size(); i++)
      sum += v1(i);
    values = sum;
    /*
    AutoDiffDiff<1,bool> add(true);
    add.DValue(0) = true;
    add.DDValue(0,0) = true;
    values = add;
    */
    /*
    FlatArray<int> hdims = Dimensions();    
    auto in0 = input[0];
    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < hdims[1]; k++)
        values(j*hdims[1]+k) = in0(k*hdims[0]+j);
    */
  }
  using T_CoefficientFunction<InverseCoefficientFunction<D>>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    throw Exception ("InverseCF:: scalar evaluate for matrix called");
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & mir,
                   BareSliceMatrix<T,ORD> result) const
  {
    c1->Evaluate (mir, result);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        Mat<D,D,T> hm;
        for (int j = 0; j < D; j++)
          for (int k = 0; k < D; k++)
            hm(j,k) = result(j*D+k, i);
        hm = Inv(hm);
        for (int j = 0; j < D; j++)
          for (int k = 0; k < D; k++)
            result(j*D+k, i) = hm(j,k);
      }
  }  

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    size_t np = ir.Size();
    auto in0 = input[0];

    for (size_t i = 0; i < np; i++)
      {
        Mat<D,D,T> hm;
        for (int j = 0; j < D; j++)
          for (int k = 0; k < D; k++)
            hm(j,k) = in0(j*D+k, i);
        hm = Inv(hm);
        for (int j = 0; j < D; j++)
          for (int k = 0; k < D; k++)
            values(j*D+k, i) = hm(j,k);
      }
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;

    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    return (-1) * thisptr * c1->Diff(var,dir) * thisptr;
  }

  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    if (this == var)
      return IdentityCF(this->Dimensions());

    auto diffc1 = c1->DiffJacobi (var, cache);
    auto inv1 = thisptr;

    Array<int> dimres { D, D };
    dimres += var->Dimensions();
    
    auto prod1 = -inv1 * diffc1->Reshape( D, -1 );
    auto prod1r = prod1 -> Reshape( Array<int> (dimres) );
    auto trans = prod1r -> TensorTranspose( 0, 1 );
    auto prod2 = inv1->Transpose() * trans -> Reshape( D, -1 );
    auto res = prod2 -> Reshape( Array<int> (dimres) ) -> TensorTranspose( 0, 1);
    cache[thisptr] = res;
    return res;
  }
};


class InverseCoefficientFunctionAnyDim : public CoefficientFunction
{
  // This is implemented in a separate class because of limitations of msvc in the context of "constexpr if".
  // Otherwise, it could be easily integrated in InverseCoefficientFunction<D>, eg. for D == -1.

  shared_ptr<CoefficientFunction> c1;
  using T_DJC = CoefficientFunction::T_DJC;

public:
  InverseCoefficientFunctionAnyDim() = default;
  InverseCoefficientFunctionAnyDim(shared_ptr<CoefficientFunction> ac1)
          : CoefficientFunction(ac1->Dimension() * ac1->Dimension(), ac1->IsComplex()), c1(ac1)
  {
    this->SetDimensions(ac1->Dimensions());
  }

  void DoArchive(Archive &ar) override
  {
    CoefficientFunction::DoArchive(ar);
    ar.Shallow(c1);
  }

  virtual void TraverseTree(const function<void(CoefficientFunction &)> &func) override
  {
    c1->TraverseTree(func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>>
  InputCoefficientFunctions() const override
  {
    return Array<shared_ptr<CoefficientFunction>>({c1});
  }

  virtual void NonZeroPattern(const class ProxyUserData &ud,
                              FlatVector<AutoDiffDiff<1, bool>> values) const override
  {
    Vector<AutoDiffDiff<1, bool>> v1(c1->Dimension());
    c1->NonZeroPattern(ud, v1);
    AutoDiffDiff<1, bool> sum(false);
    for (auto i : Range(v1))
      sum += v1(i);
    values = sum;
  }

  virtual void NonZeroPattern(const class ProxyUserData &ud,
                              FlatArray<FlatVector<AutoDiffDiff<1, bool>>> input,
                              FlatVector<AutoDiffDiff<1, bool>> values) const override
  {
    auto v1 = input[0];
    AutoDiffDiff<1, bool> sum(false);
    for (auto i : Range(v1))
      sum += v1(i);
    values = sum;
  }

  virtual double Evaluate(const BaseMappedIntegrationPoint &ip) const override
  {
    throw Exception("InverseAnyDimCF:: scalar evaluate for matrix called");
  }

  void Evaluate(const BaseMappedIntegrationPoint & ip,
                FlatVector<> result) const override
  {
    FlatMatrix<double,ColMajor> mat(Dimension(), 1, &result(0));
    ip.IntegrationRuleFromPoint([this,mat] (const BaseMappedIntegrationRule & ir)
                                { this -> T_Evaluate (ir, BareSliceMatrix<double,ColMajor>(mat)); });
  }

  void Evaluate(const BaseMappedIntegrationPoint & ip,
                FlatVector<Complex> result) const override
  {
    FlatMatrix<Complex,ColMajor> mat(Dimension(), 1, &result(0));
    ip.IntegrationRuleFromPoint([this,mat] (const BaseMappedIntegrationRule & ir)
                                { this -> T_Evaluate (ir, BareSliceMatrix<Complex,ColMajor>(mat)); });
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const override
  { this -> /* template */ T_Evaluate (ir, Trans(values)); }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         FlatArray<BareSliceMatrix<double,ColMajor>> input,
                         BareSliceMatrix<double,ColMajor> values) const override
  { this -> T_Evaluate (ir, input, values); }


  virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const override
  {
    if (!IsComplex())
      {
        BareSliceMatrix<double> realvalues(2*values.Dist(), (double*)values.Data(), DummySize(values.Height(), values.Width()));
        Evaluate (ir, realvalues);
        for (size_t i = 0; i < ir.Size(); i++)
          for (size_t j = Dimension(); j-- > 0; )
            values(i,j) = realvalues(i,j);
        return;
      }
    this -> /* template */ T_Evaluate (ir, Trans(values));
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate(const MIR &mir, BareSliceMatrix<T, ORD> result) const
  {
    c1->Evaluate(mir, result);

    const int D = c1->Dimensions()[0];
    ArrayMem<T, 1000> mem(D * D);
    FlatMatrix<T> hm(D, D, mem.Data());

    for (auto i : Range(mir))
      {
        for (auto j : Range(D))
          for (auto k : Range(D))
            hm(j, k) = result(j * D + k, i);
        hm = Inv(hm);
        for (auto j : Range(D))
          for (auto k : Range(D))
            result(j * D + k, i) = hm(j, k);
      }
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate(const MIR &ir, FlatArray<BareSliceMatrix<T, ORD>> input,
                  BareSliceMatrix<T, ORD> values) const
  {
    auto in0 = input[0];

    const int D = c1->Dimensions()[0];
    ArrayMem<T, 1000> mem(D * D);
    FlatMatrix<T> hm(D, D, mem.Data());

    for (auto i : Range(ir))
      {
        for (auto j : Range(D))
          for (auto k : Range(D))
            hm(j, k) = in0(j * D + k, i);
        hm = Inv(hm);
        for (auto j : Range(D))
          for (auto k : Range(D))
            values(j * D + k, i) = hm(j, k);
      }
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    if (code.is_simd)
        return CoefficientFunction::GenerateCode(code, inputs, index);

    const auto h = Dimensions()[0];
    const auto h_str = ToString(h);

    const auto mem_var = Var("mem", index);
    stringstream _code;
    _code << "ArrayMem<" << code.res_type << ", 400> " << mem_var.code << "("
          << Dimension() * 2 << ");"
          << "\n";
    code.body += _code.str();

    auto mat_var = Var("mat", index);
    auto inv_var = Var("inv", index);
    stringstream _code2;
    _code2 << "FlatMatrix<" + code.res_type + "> " << mat_var.code << "(" << h_str
           << ", " << mem_var.code << ".Data()"
           << ");"
           << "\n";
    _code2 << "FlatMatrix<" + code.res_type + "> " << inv_var.code << "(" << h_str
           << ", " << mem_var.code << ".Data() + " << Dimension() << ");"
           << "\n";

    code.body += _code2.str();
    for (auto j : Range(h))
      for (auto k : Range(h))
        code.body += mat_var(j, k).Assign(Var(inputs[0], j, k), false);

    code.body += inv_var.Assign(mat_var.Func("Inv"), false);

    for (auto j : Range(h))
      for (auto k : Range(h))
        code.body += Var(index, j, k).Assign(inv_var(j, k));
  }

  shared_ptr<CoefficientFunction>
  Diff(const CoefficientFunction *var,
       shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var)
      return c1->Diff(c1.get(), dir);

    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    return (-1) * thisptr * c1->Diff(var, dir) * thisptr;
  }

  shared_ptr<CoefficientFunction> DiffJacobi(const CoefficientFunction *var,
                                             T_DJC &cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    if (this == var)
      return IdentityCF(this->Dimensions());

    auto diffc1 = c1->DiffJacobi (var, cache);
    auto inv1 = thisptr;

    const int D = c1->Dimensions()[0];
    Array<int> dimres { D, D };
    dimres += var->Dimensions();

    auto prod1 = -inv1 * diffc1->Reshape( D, -1 );
    auto prod1r = prod1 -> Reshape( Array<int> (dimres) );
    auto trans = prod1r -> TensorTranspose( 0, 1 );
    auto prod2 = inv1->Transpose() * trans -> Reshape( D, -1 );
    auto res = prod2 -> Reshape( Array<int> (dimres) ) -> TensorTranspose( 0, 1);
    cache[thisptr] = res;
    return res;
  }
};



template <int D>
class DeterminantCoefficientFunction : public T_CoefficientFunction<DeterminantCoefficientFunction<D>>
{
  shared_ptr<CoefficientFunction> c1;
  using BASE = T_CoefficientFunction<DeterminantCoefficientFunction<D>>;
public:
  DeterminantCoefficientFunction() = default;
  DeterminantCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : T_CoefficientFunction<DeterminantCoefficientFunction>(1, ac1->IsComplex()), c1(ac1)
  {
    auto dims_c1 = c1 -> Dimensions();
    if (dims_c1.Size() != 2)
      throw Exception("Determinant of non-matrix called");
    if (dims_c1[0] != dims_c1[1])
      throw Exception("Determinant of non-symmetric matrix called");
  }

  virtual string GetDescription () const override
  { return "Determinant"; }

  
  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1);
  }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
    // auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.res_type+">";
    auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.GetType(this->IsComplex())+">";
    auto mat_var = Var("mat", index);
    code.body += mat_var.Declare(mat_type);
    for (int j = 0; j < D; j++)
      for (int k = 0; k < D; k++)
        code.body += mat_var(j,k).Assign(Var(inputs[0], j, k), false);

    // code.Declare (code.res_type, index, this->Dimensions());
    code.Declare (index, this->Dimensions(), this->IsComplex());
    
    code.body += Var(index).Assign(mat_var.Func("Det"), false);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    nonzero = true;
    nonzero_deriv = true;
    nonzero_dderiv = true;
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    /*
    AutoDiffDiff<1,bool> add(true);
    add.DValue(0) = true;
    add.DDValue(0,0) = true;
    values = add;
    */
    Vector<AutoDiffDiff<1,bool>> in(D*D);
    c1->NonZeroPattern (ud, in);
    Array<FlatVector<AutoDiffDiff<1,bool>>> input{1UL};
    input[0].AssignMemory(D*D, &in(0));
    NonZeroPattern (ud, input, values);
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto sm = input[0];
    AutoDiffDiff<1,bool> res;
    switch (D)
      {
      case 1: 
        res = sm(0);
        break;
      case 2:
        res = sm(0)*sm(3)+sm(1)*sm(2);
        break;
      case 3:
        res = 
          sm(0) * (sm(4) * sm(8) + sm(5) * sm(7)) +
          sm(1) * (sm(5) * sm(6) + sm(3) * sm(8)) +
          sm(2) * (sm(3) * sm(7) + sm(4) * sm(6));
        break;
      default:
	{
	  cerr << "general det not implemented" << endl;
	}
      }

    /*
    Mat<D,D,AutoDiffDiff<1,bool>> hm;
    for (int j = 0; j < D; j++)
      for (int k = 0; k < D; k++)
        hm(j,k) = in0(j*D+k);
    cout << "Nonzero mat:" << endl << hm << endl;
    values(0) = Det(hm);
    cout << "nonzero det: " << values(0) << endl;
    */
    values(0) = res;
  }
  
  using T_CoefficientFunction<DeterminantCoefficientFunction<D>>::Evaluate;
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & mir,
                   BareSliceMatrix<T,ORD> result) const
  {
    STACK_ARRAY(T, hmem, mir.Size()*D*D);
    FlatMatrix<T,ORD> hv(D*D, mir.Size(), &hmem[0]);
    c1->Evaluate (mir, hv);
    
    for (size_t i = 0; i < mir.Size(); i++)
      {
        Mat<D,D,T> hm;
        for (int j = 0; j < D; j++)
          for (int k = 0; k < D; k++)
            hm(j,k) = hv(j*D+k, i);
        result(0,i) = Det(hm);
      }
  }  

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    size_t np = ir.Size();
    auto in0 = input[0];

    for (size_t i = 0; i < np; i++)
      {
        Mat<D,D,T> hm;
        for (int j = 0; j < D; j++)
          for (int k = 0; k < D; k++)
            hm(j,k) = in0(j*D+k, i);
        values(0,i) = Det(hm);        
      }
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    // return DeterminantCF(c1) * InnerProduct( TransposeCF(InverseCF(c1)), c1->Diff(var,dir) );
    return InnerProduct( CofactorCF(c1), c1->Diff(var,dir) );
  }

  shared_ptr<CoefficientFunction>
  DiffJacobi (const CoefficientFunction * var, typename BASE::T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    if (c1.get() == var) return CofactorCF (c1);
    auto input = c1->InputCoefficientFunctions();
    if (input.Size() == 0) return ZeroCF(var->Dimensions());

    auto cof = CofactorCF(c1) -> Reshape( 1, D*D );
    auto diffc1 = c1->DiffJacobi (var, cache) -> Reshape( D*D, var->Dimension() );
    auto prod = cof * diffc1;
    auto res = prod->Reshape(var->Dimensions());

    /*
    auto cof = CofactorCF(c1) -> Reshape( D*D );
    auto diffc1 = c1->DiffJacobi (var, cache) -> Reshape( D*D, var->Dimension() );
    auto res = diffc1->Transpose() * cof;
    */
    
    cache[thisptr] = res;
    return res;
  }

  
};








template <int D>
class CofactorCoefficientFunction : public T_CoefficientFunction<CofactorCoefficientFunction<D>>
{
  shared_ptr<CoefficientFunction> c1;
  using BASE = T_CoefficientFunction<CofactorCoefficientFunction<D>>;
public:
  CofactorCoefficientFunction() = default;
  CofactorCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : T_CoefficientFunction<CofactorCoefficientFunction>(D*D, ac1->IsComplex()), c1(ac1)
  {
    this->SetDimensions (ngstd::INT<2> (D,D));
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1);
  }

  virtual string GetDescription () const override
  { return "cofactor"; }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
    // auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.res_type+">";
    auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.GetType(this->IsComplex())+">";
    auto mat_var = Var("mat", index);
    auto cof_var = Var("cof", index);
    code.body += mat_var.Declare(mat_type);
    code.body += cof_var.Declare(mat_type);
    for (int j = 0; j < D; j++)
      for (int k = 0; k < D; k++)
        code.body += mat_var(j,k).Assign(Var(inputs[0], j, k), false);

    code.body += cof_var.Assign(mat_var.Func("Cof"), false);
    
    // code.Declare (code.res_type, index, this->Dimensions());
    code.Declare (index, this->Dimensions(), this->IsComplex()); 
    for (int j = 0; j < D; j++)
      for (int k = 0; k < D; k++)
        code.body += Var(index, j, k).Assign(cof_var(j,k), false);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    nonzero = true;
    nonzero_deriv = true;
    nonzero_dderiv = true;
  }
  */

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> v1(c1->Dimension());
    c1->NonZeroPattern (ud, v1);
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < v1.Size(); i++)
      sum += v1(i);
    values = sum;
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto v1 = input[0];
    AutoDiffDiff<1,bool> sum(false);
    for (int i = 0; i < v1.Size(); i++)
      sum += v1(i);
    values = sum;
    /*
    AutoDiffDiff<1,bool> add(true);
    add.DValue(0) = true;
    add.DDValue(0,0) = true;
    values = add;
    */
    /*
    FlatArray<int> hdims = Dimensions();    
    auto in0 = input[0];
    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < hdims[1]; k++)
        values(j*hdims[1]+k) = in0(k*hdims[0]+j);
    */
  }
  using T_CoefficientFunction<CofactorCoefficientFunction<D>>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    throw Exception ("CofactorCF:: scalar evaluate for matrix called");
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & mir,
                   BareSliceMatrix<T,ORD> result) const
  {
    c1->Evaluate (mir, result);
    for (size_t i = 0; i < mir.Size(); i++)
      {
        Mat<D,D,T> hm;
        for (int j = 0; j < D; j++)
          for (int k = 0; k < D; k++)
            hm(j,k) = result(j*D+k, i);
        hm = Cof(hm);
        for (int j = 0; j < D; j++)
          for (int k = 0; k < D; k++)
            result(j*D+k, i) = hm(j,k);
      }
  }  

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    size_t np = ir.Size();
    auto in0 = input[0];

    for (size_t i = 0; i < np; i++)
      {
        Mat<D,D,T> hm;
        for (int j = 0; j < D; j++)
          for (int k = 0; k < D; k++)
            hm(j,k) = in0(j*D+k, i);
        hm = Cof(hm);
        for (int j = 0; j < D; j++)
          for (int k = 0; k < D; k++)
            values(j*D+k, i) = hm(j,k);
      }
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    if (this->Dimensions()[0] <= 2)
      {
        //Cofactor Matrix linear in 2d (in 1d Cofactor Matrix = 0)
        return CofactorCF(c1->Diff(var,dir));
      }
    else if (this->Dimensions()[0] == 3) //3d
      {
        //formula follows from Cayley–Hamilton
        //Cof(A) = 0.5*(tr(A)**2 - tr(A**2))I - tr(A)A^T +(AA)^T

        //return (0.5*(TraceCF(c1)*TraceCF(c1) - TraceCF(c1*c1))*IdentityCF(3) - TraceCF(c1)*TransposeCF(c1) + TransposeCF(c1*c1))->Diff(var,dir);
        return  0.5*(2*TraceCF(c1)*TraceCF(c1->Diff(var,dir)) - TraceCF(c1->Diff(var,dir)*c1 + c1 * c1->Diff(var,dir)))*IdentityCF(3)- TraceCF(c1->Diff(var,dir))*TransposeCF(c1) - TraceCF(c1)*TransposeCF(c1->Diff(var,dir)) + TransposeCF(c1->Diff(var,dir)*c1 + c1 * c1->Diff(var,dir));
      }
    else
      throw Exception("CofactorCF Diff only implemented for dim <=3");
  }

  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, typename BASE::T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    if (this == var)
      return IdentityCF(this->Dimensions());

    shared_ptr<CoefficientFunction> res;

    if (this->Dimensions()[0] == 2)
      res = (TraceCF(c1)*IdentityCF(2)-TransposeCF(c1)) -> DiffJacobi (var, cache);
    else if (this->Dimensions()[0] == 3)
      {
        auto trcf = TraceCF(c1);
        auto prodcf = c1*c1;
        res = (0.5*(trcf*trcf - TraceCF(prodcf))*IdentityCF(3) - trcf*TransposeCF(c1) + TransposeCF(prodcf)) -> DiffJacobi (var, cache);
      }
    else
      res = (DeterminantCF(c1) * InverseCF(c1)->Transpose() ) -> DiffJacobi (var, cache);
    cache[thisptr] = res;
    return res;
  }

  
};








class SymmetricCoefficientFunction : public T_CoefficientFunction<SymmetricCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  using BASE = T_CoefficientFunction<SymmetricCoefficientFunction>;
public:
  SymmetricCoefficientFunction() = default;
  SymmetricCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : T_CoefficientFunction<SymmetricCoefficientFunction>(1, ac1->IsComplex()), c1(ac1)
  {
    auto dims_c1 = c1 -> Dimensions();
    if (dims_c1.Size() != 2)
      throw Exception("Sym of non-matrix called");
    if (dims_c1[0] != dims_c1[1])
      throw Exception("Sym of non-square matrix called");
    
    SetDimensions (ngstd::INT<2> (dims_c1[0], dims_c1[0]) );
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1);
  }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual string GetDescription () const override
  { return "symmetric"; }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      FlatArray<int> hdims = Dimensions();        
      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1]))
          code.body += Var(index,i,j).Assign("0.5*("+Var(inputs[0],i,j).S()+"+"+Var(inputs[0],j,i).S()+")");
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    int hd = Dimensions()[0];
    c1->NonZeroPattern (ud, values);

    for (int i = 0; i < hd; i++)
      for (int j = 0; j < hd; j++)
        {
          int ii = i*hd+j;
          int jj = j*hd+i;
          values(ii) = values(ii)+values(jj);  // logical or
        }
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    int hd = Dimensions()[0];    
    auto in0 = input[0];
    for (int i = 0; i < hd; i++)
      for (int j = 0; j < hd; j++)
        {
          int ii = i*hd+j;
          int jj = j*hd+i;
          values(ii) = in0(ii)+in0(jj);   // logical or 
        }
  }
  using T_CoefficientFunction<SymmetricCoefficientFunction>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    throw Exception ("SymCF:: scalar evaluate for matrix called");
  }
  
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & mir,
                   BareSliceMatrix<T,ORD> result) const
  {
    int hd = Dimensions()[0];
    c1->Evaluate (mir, result);
    STACK_ARRAY(T, hmem, hd*hd);
    FlatMatrix<T,ORD> tmp (hd, hd, &hmem[0]);

    for (size_t i = 0; i < mir.Size(); i++)
      {
        for (int j = 0; j < hd; j++)
          for (int k = 0; k < hd; k++)
            tmp(j,k) = result(k*hd+j, i);
        for (int j = 0; j < hd; j++)
          for (int k = 0; k < hd; k++)
            result(j*hd+k, i) = 0.5*(tmp(j,k)+tmp(k,j));
      }
  }  

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    int hd = Dimensions()[0];
    size_t np = ir.Size();
    
    auto in0 = input[0];
    for (size_t j = 0; j < hd; j++)
      for (size_t k = 0; k < hd; k++)
        for (size_t i = 0; i < np; i++)
          values(j*hd+k, i) = 0.5 * (in0(k*hd+j, i)+in0(j*hd+k, i));
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return SymmetricCF(c1->Diff(var, dir));
  }

  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    if (this == var)
      return IdentityCF(this->Dimensions());
  
    auto diffc1 = c1->DiffJacobi (var, cache);
    auto res = 0.5*(diffc1 + diffc1 -> TensorTranspose( 0, 1 ));
    cache[thisptr] = res;
    return res;
  }
};


class SkewCoefficientFunction : public T_CoefficientFunction<SkewCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  using BASE = T_CoefficientFunction<SkewCoefficientFunction>;
public:
  SkewCoefficientFunction() = default;
  SkewCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : T_CoefficientFunction<SkewCoefficientFunction>(1, ac1->IsComplex()), c1(ac1)
  {
    auto dims_c1 = c1 -> Dimensions();
    if (dims_c1.Size() != 2)
      throw Exception("Skew of non-matrix called");
    if (dims_c1[0] != dims_c1[1])
      throw Exception("Skew of non-square matrix called");
    
    SetDimensions (ngstd::INT<2> (dims_c1[0], dims_c1[0]) );
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1);
  }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual string GetDescription () const override
  { return "skew"; }
  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      FlatArray<int> hdims = Dimensions();        
      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1]))
          code.body += Var(index,i,j).Assign("0.5*("+Var(inputs[0],i,j).S()+"-"+Var(inputs[0],j,i).S()+")");
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    int hd = Dimensions()[0];    
    c1->NonZeroPattern (ud, nonzero, nonzero_deriv, nonzero_dderiv);
    for (int i = 0; i < hd; i++)
      for (int j = 0; j < hd; j++)
        {
          int ii = i*hd+j;
          int jj = j*hd+i;
          nonzero(ii) |= nonzero(jj);
          nonzero_deriv(ii) |= nonzero_deriv(jj);
          nonzero_dderiv(ii) |= nonzero_dderiv(jj);
        }
  }
  */

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    int hd = Dimensions()[0];
    c1->NonZeroPattern (ud, values);
    for (int i = 0; i < hd; i++)
      for (int j = 0; j < hd; j++)
        {
          int ii = i*hd+j;
          int jj = j*hd+i;
          values(ii) = values(ii)+values(jj);   // logical or 
        }
  }

  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    int hd = Dimensions()[0];    
    auto in0 = input[0];
    for (int i = 0; i < hd; i++)
      for (int j = 0; j < hd; j++)
        {
          int ii = i*hd+j;
          int jj = j*hd+i;
          values(ii) = in0(ii)+in0(jj);   // logical or 
        }
  }
  using T_CoefficientFunction<SkewCoefficientFunction>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    throw Exception ("SkewCF:: scalar evaluate for matrix called");
  }

  
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & mir,
                   BareSliceMatrix<T,ORD> result) const
  {
    int hd = Dimensions()[0];
    c1->Evaluate (mir, result);
    STACK_ARRAY(T, hmem, hd*hd);
    FlatMatrix<T,ORD> tmp (hd, hd, &hmem[0]);

    for (size_t i = 0; i < mir.Size(); i++)
      {
        for (int j = 0; j < hd; j++)
          for (int k = 0; k < hd; k++)
            tmp(j,k) = result(j*hd+k, i);
        for (int j = 0; j < hd; j++)
          for (int k = 0; k < hd; k++)
            result(j*hd+k, i) = 0.5*(tmp(j,k)-tmp(k,j));
      }
  }  

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    int hd = Dimensions()[0];
    size_t np = ir.Size();
    
    auto in0 = input[0];
    for (size_t j = 0; j < hd; j++)
      for (size_t k = 0; k < hd; k++)
        for (size_t i = 0; i < np; i++)
          values(j*hd+k, i) = 0.5 * (in0(j*hd+k, i)-in0(k*hd+j, i));
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return SkewCF(c1->Diff(var, dir));
  }

  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    if (this == var)
      return IdentityCF(this->Dimensions());
  
    auto diffc1 = c1->DiffJacobi (var, cache);
    auto res = 0.5*(diffc1 - diffc1 -> TensorTranspose( 0, 1 ));
    cache[thisptr] = res;
    return res;
  }
};





class TraceCoefficientFunction : public T_CoefficientFunction<TraceCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  using BASE = T_CoefficientFunction<TraceCoefficientFunction>;
public:
  TraceCoefficientFunction() = default;
  TraceCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : T_CoefficientFunction<TraceCoefficientFunction>(1, ac1->IsComplex()), c1(ac1)
  {
    auto dims_c1 = c1 -> Dimensions();
    if (dims_c1.Size() != 2)
      throw Exception("Trace of non-matrix called");
    if (dims_c1[0] != dims_c1[1])
      throw Exception("Trace of non-square matrix called");
  }

  virtual string GetDescription () const override
  { return "trace"; }
  
  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1);
  }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
    CodeExpr result;
    int dim1 = c1->Dimensions()[0];
    code.Declare (index, Array<int>(), IsComplex()); 
    for (int i = 0; i < dim1; i++)
      result += Var(inputs[0],i,i);
    code.body += Var(index).Assign(result.S(), false);
  }
  
  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    int dim1 = c1->Dimension();
    Vector<bool> v1(dim1), d1(dim1), dd1(dim1);
    c1->NonZeroPattern (ud, v1, d1, dd1);

    bool v = false, d = false, dd = false;
    for (int i = 0; i < dim1; i++)
      {
        v |= v1(i);
        d |= d1(i);
        dd |= dd1(i);
      }
    nonzero = v;
    nonzero_deriv = d;
    nonzero_dderiv = dd;
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    int dim1 = c1->Dimension();
    int d = c1->Dimensions()[0];
    Vector<AutoDiffDiff<1,bool>> v1(dim1);
    c1->NonZeroPattern (ud, v1);
    values(0) = false;
    /*
    for (int i = 0; i < dim1; i++)
      values(0) = values(0)+v1(i);   // logical or
    */
    for (int i = 0; i < d; i++)
      values(0) = values(0)+v1(i*(d+1));   // logical or
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    // int dim1 = c1->Dimension();
    int d = c1->Dimensions()[0];
    auto in0 = input[0];
    values(0) = false;
    /*
    for (int i = 0; i < dim1; i++)
      values(0) = values(0)+in0(i);   // logical or
    */
    for (int i = 0; i < d; i++)
      values(0) = values(0)+in0(i*(d+1));   // logical or
  }

  using T_CoefficientFunction<TraceCoefficientFunction>::Evaluate;
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & mir,
                   BareSliceMatrix<T,ORD> result) const
  {
    int hd = c1->Dimensions()[0];
    STACK_ARRAY(T, hmem, hd*hd*mir.Size());
    FlatMatrix<T,ORD> m1(hd*hd, mir.Size(), &hmem[0]);
    c1->Evaluate (mir, m1);
    
    for (size_t i = 0; i < mir.Size(); i++)
      {
        T sum{0.0};
        for (int j = 0; j < hd; j++)
          sum += m1(j*(hd+1), i);
        result(0, i) = sum;
      }
  }  

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    int hd = c1->Dimensions()[0];
    size_t np = ir.Size();
    
    auto in0 = input[0];
    for (size_t i = 0; i < np; i++)
      {
        T sum{0.0};
        for (size_t j = 0; j < hd; j++)
          sum += in0(j*(hd+1), i);
        values(0,i) = sum;
      }
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return TraceCF(c1->Diff(var, dir));
  }

  shared_ptr<CoefficientFunction>
  DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    if (c1.get() == var) return IdentityCF (c1->Dimensions()[0]);
    auto input = c1->InputCoefficientFunctions();
    if (input.Size() == 0) return ZeroCF(c1->Dimensions());

    shared_ptr<CoefficientFunction> res;
    if (c1->GetDescription()=="binary operation '-'")
      res = TraceCF(input[0])->DiffJacobi (var, cache) - TraceCF(input[1])->DiffJacobi (var, cache);
    else if (dynamic_pointer_cast<MultMatMatCoefficientFunction>(c1) && !c1->IsComplex())
      {
        auto AB = c1->InputCoefficientFunctions();
        res = InnerProduct(AB[0], AB[1]->Transpose() )->DiffJacobi (var, cache);
      }
    else
      res = MakeTensorTraceCoefficientFunction (c1->DiffJacobi (var, cache), 0, 1);
    cache[thisptr] = res;
    return res;
  }
  
};




  

  
  // ///////////////////////////// operators  /////////////////////////

  struct GenericPlus {
    template <typename T> T operator() (T x, T y) const { return x+y; }
    void DoArchive(Archive& ar){}
  };
  struct GenericMinus {
    template <typename T> T operator() (T x, T y) const { return x-y; }
    void DoArchive(Archive& ar){}
  };
  struct GenericMult {
    template <typename T> T operator() (T x, T y) const { return x*y; }
    void DoArchive(Archive& ar){}
  };
  struct GenericDiv {
    template <typename T> T operator() (T x, T y) const { return x/y; }
    void DoArchive(Archive& ar){}
  };
  GenericPlus gen_plus;
  GenericMinus gen_minus;
  GenericMult gen_mult;
  GenericDiv gen_div;

template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericPlus>::Diff(const CoefficientFunction * var,
                                   shared_ptr<CoefficientFunction> dir) const
{
  if (var == this) return dir;
  return c1->Diff(var,dir) + c2->Diff(var,dir);
}

template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericPlus>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
{
  auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
  if (cache.find(thisptr) != cache.end())
    return cache[thisptr];

  if (var == this) return IdentityCF(Dimensions());
  auto res = c1->DiffJacobi (var, cache) + c2->DiffJacobi (var, cache);
  cache[thisptr] = res;
  return res;
}

template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericPlus>::Operator(const string & name) const
{
  return c1->Operator(name) + c2->Operator(name);
}

shared_ptr<CoefficientFunction> operator+ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2)
{
  if (c1->IsZeroCF())
    {
      if (c2->IsZeroCF())
        return c1;
      else
        return c2;
    }
  else if (c2->IsZeroCF())
    return c1;
  
  return BinaryOpCF (c1, c2, gen_plus, "+");
}


template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericMinus>::Diff(const CoefficientFunction * var,
                                    shared_ptr<CoefficientFunction> dir) const
{
  if (var == this) return dir;      
  return c1->Diff(var,dir) - c2->Diff(var,dir);
}

template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericMinus>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
{
  auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
  if (cache.find(thisptr) != cache.end())
    return cache[thisptr];

  if (var == this) return IdentityCF(Dimensions());
  auto res = c1->DiffJacobi (var, cache) - c2->DiffJacobi (var, cache);
  cache[thisptr] = res;
  return res;
}


template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericMinus>::Operator(const string & name) const
{
  return c1->Operator(name) - c2->Operator(name);
}


shared_ptr<CoefficientFunction> operator- (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2)
{
  if (c1->IsZeroCF())
    {
      if (c2->IsZeroCF())
        return c1;
      else
        return -c2;
    }
  else if (c2->IsZeroCF())
    return c1;

  return BinaryOpCF (c1, c2, gen_minus, "-");
}

template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericMult>::Operator(const string & name) const
{
  if (c1->Dimension() != 1 || c2->Dimension() != 1)
    throw Exception ("can only differentiate scalar multiplication");
  if (name != "grad")
    throw Exception ("can only from 'grad' operator of product");
  return c1->Operator(name) * c2  + c1 * c2->Operator(name);
}


shared_ptr<CoefficientFunction> CWMult (shared_ptr<CoefficientFunction> cf1,
                                        shared_ptr<CoefficientFunction> cf2)
{
  if (cf1->IsZeroCF() || cf2->IsZeroCF())
    return ZeroCF( cf1->Dimensions() );
  return BinaryOpCF (cf1, cf2, gen_mult, "*");
}


template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericMult>::Diff(const CoefficientFunction * var,
                                   shared_ptr<CoefficientFunction> dir) const
{
  if (var == this) return dir;    
  //return c1->Diff(var,dir)*c2 + c1*c2->Diff(var,dir); // would replace point-wise mult by InnerProduct
  return CWMult (c1->Diff(var,dir), c2) + 
    CWMult(c1, c2->Diff(var,dir));
}

template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericMult>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
{
  auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
  if (cache.find(thisptr) != cache.end())
    return cache[thisptr];

  if (var == this) return make_shared<ConstantCoefficientFunction>(1);
  shared_ptr<CoefficientFunction> res;
  if (c1 == c2)
    res = (2*c1) * c1->DiffJacobi (var, cache);
  else
    res = c1 * c2->DiffJacobi (var, cache) + c2 * c1->DiffJacobi (var, cache);
  cache[thisptr] = res;
  return res;
}




template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericDiv>::Diff(const CoefficientFunction * var,
                                shared_ptr<CoefficientFunction> dir) const 
{
  if (var == this) return dir;
  // return (c1->Diff(var,dir)*c2 - c1*c2->Diff(var,dir)) / (c2*c2);
  /*return (BinaryOpCF (c1->Diff(var,dir), c2, gen_mult, "*") -
    BinaryOpCF (c1, c2->Diff(var,dir), gen_mult, "*")) /
    BinaryOpCF (c2, c2, gen_mult, "*");
    BinaryOpCF (c2, c2, gen_mult, "*");*/
  return (CWMult (c1->Diff(var,dir), c2) - CWMult (c1, c2->Diff(var,dir))) / CWMult (c2, c2);
}

template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericDiv>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const 
{
  auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
  if (cache.find(thisptr) != cache.end())
    return cache[thisptr];

  if (Dimensions().Size() > 0)
    return CoefficientFunction::DiffJacobi (var, cache);
  
  if (var == this)
    return make_shared<ConstantCoefficientFunction>(1);
  
  auto res = (c2*c1->DiffJacobi (var, cache) - c1*c2->DiffJacobi (var, cache)) / (c2*c2);
  cache[thisptr] = res;
  return res;
}




shared_ptr<CoefficientFunction> operator* (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2)
  {
    if (c1->IsZeroCF() || c2->IsZeroCF())
      {
        if (c1->Dimensions().Size() == 2 && c2->Dimensions().Size() == 2)
          return ZeroCF(Array( {c1->Dimensions()[0],c2->Dimensions()[1]} ));
        if (c1->Dimensions().Size() == 2 && c2->Dimensions().Size() == 1)
          return ZeroCF(Array( {c1->Dimensions()[0]} ));
        if (c1->Dimension() > 1 && c2->Dimension() > 1)
          return ZeroCF(Array<int>( {} ));
        if ( (c1->Dimension() == 1 && c2->Dimension() > 1) || (c1->Dimension() > 1 && c2->Dimension() == 1))
          return ZeroCF(Array( {int(c1->Dimension()*c2->Dimension())} ));

        return ZeroCF(Array<int>( {} ));
          
      }
    if (c1->Dimensions().Size() == 2 && c2->Dimensions().Size() == 2)
      {
        if (dynamic_pointer_cast<IdentityCoefficientFunction>(c1) && !c1->IsVariable())
          return c2;
        if (dynamic_pointer_cast<IdentityCoefficientFunction>(c2) && !c2->IsVariable())
          return c1;
        return make_shared<MultMatMatCoefficientFunction> (c1, c2);
      }
    if (c1->Dimensions().Size() >= 2 && c2->Dimensions().Size() == 1)
      {
        if (dynamic_pointer_cast<IdentityCoefficientFunction>(c1) && !c1->IsVariable())
          return c2;
        return make_shared<MultMatVecCoefficientFunction> (c1, c2);
      }

    if (c1->Dimensions().Size()==1 && c2->Dimensions().Size()==1)
      if (c1->Dimension() > 1 && c1->Dimension()==c2->Dimension())
        {
          switch (c1->Dimension())
          {
          case 2:
            return make_shared<T_MultVecVecCoefficientFunction<2>> (c1, c2);
          case 3:
            return make_shared<T_MultVecVecCoefficientFunction<3>> (c1, c2);
          case 4:
            return make_shared<T_MultVecVecCoefficientFunction<4>> (c1, c2);
          case 5:
            return make_shared<T_MultVecVecCoefficientFunction<5>> (c1, c2);
          default:
            return make_shared<MultVecVecCoefficientFunction> (c1, c2);
          }
      }
    if (c1->Dimension() == 1 && c2->Dimension() > 1)
      {
        if (c1->Dimensions().Size())
          return make_shared<MultScalVecCoefficientFunction> (MakeComponentCoefficientFunction(c1,0), c2);
        else
          return make_shared<MultScalVecCoefficientFunction> (c1, c2);
      }
    if (c1->Dimension() > 1 && c2->Dimension() == 1)
      {
        if (c2->Dimensions().Size())
          return make_shared<MultScalVecCoefficientFunction> (MakeComponentCoefficientFunction(c2,0), c1);
        else
          return make_shared<MultScalVecCoefficientFunction> (c2, c1);
      }
    return BinaryOpCF (c1, c2, gen_mult,"*");
  }

  shared_ptr<CoefficientFunction> operator* (double v1, shared_ptr<CoefficientFunction> c2)
  {
    if (c2->IsZeroCF())
      return c2;
    else if (v1 == double(0))
      return ZeroCF(c2->Dimensions());

    return make_shared<ScaleCoefficientFunction> (v1, c2); 
  }
  
  shared_ptr<CoefficientFunction> operator* (Complex v1, shared_ptr<CoefficientFunction> c2)
  {
    if (c2->IsZeroCF())
      return c2;
    else if (v1 == Complex(0))
      return ZeroCF(c2->Dimensions());

    return make_shared<ScaleCoefficientFunctionC> (v1, c2); 
  }


  struct GenericConj {
    template <typename T> T operator() (T x) const { return Conj(x); } // from bla
    static string Name() { return "conj"; }
    SIMD<double> operator() (SIMD<double> x) const { return x; }
    template<typename T>
    AutoDiff<1,T> operator() (AutoDiff<1,T> x) const { throw Exception ("Conj(..) is not complex differentiable"); }
    template<typename T>
    AutoDiffDiff<1,T> operator() (AutoDiffDiff<1,T> x) const { throw Exception ("Conj(..) is not complex differentiable"); }
    void DoArchive(Archive& ar) {}
  };

  template <> 
  shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericConj>::Diff(const CoefficientFunction * var,
                                   shared_ptr<CoefficientFunction> dir) const
  {
    if (var == this) return dir;
    cout << "Warning: differentiate conjugate by taking conjugate of derivative" << endl;
    return ConjCF(c1->Diff(var, dir));
  }


  shared_ptr<CoefficientFunction> ConjCF (shared_ptr<CoefficientFunction> c1)
  {
    if (c1->IsZeroCF())
      return c1;
    return UnaryOpCF(c1, GenericConj(), GenericConj::Name());
  }


template <>
shared_ptr<CoefficientFunction>
cl_UnaryOpCF<GenericIdentity>::Diff(const CoefficientFunction * var,
                                      shared_ptr<CoefficientFunction> dir) const
{
  if (var == this) return dir;
  auto hcf = c1->Diff(var, dir);
  if (! (this->Dimensions() == hcf->Dimensions()) )
    { // reshaping requires wrapper, e.g. for code generation
      hcf = UnaryOpCF (hcf, GenericIdentity{}, " ");
      hcf->SetDimensions(Dimensions());
    }
  return hcf;
}

template <>
shared_ptr<CoefficientFunction>
cl_UnaryOpCF<GenericIdentity>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
{
  auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
  if (cache.find(thisptr) != cache.end())
    return cache[thisptr];

  if (var == this)
    {
      if (this->Dimensions().Size()==0)
        return make_shared<ConstantCoefficientFunction>(1);
      return IdentityCF(this->Dimensions());
    }

  Array<int> ndims(this->Dimensions());
  ndims += var->Dimensions();
  auto res = c1->DiffJacobi (var, cache) -> Reshape(ndims);
  cache[thisptr] = res;
  return res;
}





template <>
shared_ptr<CoefficientFunction>
cl_UnaryOpCF<GenericIdentity>::Operator(const string & name) const
{
  return c1->Operator(name);
}

  shared_ptr<CoefficientFunction> CreateWrapperCF (shared_ptr<CoefficientFunction> cf)
  {
    return UnaryOpCF (cf, GenericIdentity{}, " ");
  }


  shared_ptr<CoefficientFunction> InnerProduct (shared_ptr<CoefficientFunction> c1,
                                                shared_ptr<CoefficientFunction> c2)
  {
    if (c1->IsZeroCF() || c2->IsZeroCF())
      return ZeroCF( Array<int>() );
        
    if (c2->IsComplex())
      {
        auto conj = ConjCF(c2);
        if (conj->GetDescription() == c2->GetDescription())
          cout << "Info: InnerProduct has been changed and takes now conjugate" << endl
               << "since c2 is already a Conjugate operation, we don't take conjugate" << endl
               << "is you don't want conjugate, use a*b" << endl;
        else
          c2 = conj;
      }

    // if (c1->GetDescription() == "UnitCF")
    // if (typeid(*c1) == typeid(UnitVectorCoefficientFunctionX))
    if (auto uv = dynamic_pointer_cast<UnitVectorCoefficientFunction>(c1))
      return MakeComponentCoefficientFunction(c2, uv->GetCoordinate());
    
    // if (c2->GetDescription() == "UnitCF")
    // return MakeComponentCoefficientFunction(move(c1),static_cast<UnitVectorCoefficientFunction*>(c2.get())->GetCoordinate());
    if (auto uv = dynamic_pointer_cast<UnitVectorCoefficientFunction>(c2))
      return MakeComponentCoefficientFunction(c1, uv->GetCoordinate());

    if (auto transc1 = dynamic_pointer_cast<TransposeCoefficientFunction>(c1))
      if (auto transc2 = dynamic_pointer_cast<TransposeCoefficientFunction>(c2))
        {
          cout << IM(6) << "simplify double transpose" << endl;
          auto sub1 = c1->InputCoefficientFunctions();
          auto sub2 = c2->InputCoefficientFunctions();
          return InnerProduct(sub1[0], sub2[0]);
        }
      
    
    if (c1 == c2)
      {
        switch (c1->Dimension())
          {
          case 1:
            return make_shared<T_MultVecVecSameCoefficientFunction<1>> (c1);
          case 2:
            return make_shared<T_MultVecVecSameCoefficientFunction<2>> (c1);
          case 3:
            return make_shared<T_MultVecVecSameCoefficientFunction<3>> (c1);
          case 4:
            return make_shared<T_MultVecVecSameCoefficientFunction<4>> (c1);
          case 5:
            return make_shared<T_MultVecVecSameCoefficientFunction<5>> (c1);
          case 6:
            return make_shared<T_MultVecVecSameCoefficientFunction<6>> (c1);
          case 8:
            return make_shared<T_MultVecVecSameCoefficientFunction<8>> (c1);
          case 9:
            return make_shared<T_MultVecVecSameCoefficientFunction<9>> (c1);
          default:
            ;
          }
      }
    
    switch (c1->Dimension())
      {
      case 1:
        return make_shared<T_MultVecVecCoefficientFunction<1>> (c1, c2);
      case 2:
        return make_shared<T_MultVecVecCoefficientFunction<2>> (c1, c2);
      case 3:
        return make_shared<T_MultVecVecCoefficientFunction<3>> (c1, c2);
      case 4:
        return make_shared<T_MultVecVecCoefficientFunction<4>> (c1, c2);
      case 5:
        return make_shared<T_MultVecVecCoefficientFunction<5>> (c1, c2);
      case 6:
        return make_shared<T_MultVecVecCoefficientFunction<6>> (c1, c2);
      case 8:
        return make_shared<T_MultVecVecCoefficientFunction<8>> (c1, c2);
      case 9:
        return make_shared<T_MultVecVecCoefficientFunction<9>> (c1, c2);
      default:
        return make_shared<MultVecVecCoefficientFunction> (c1, c2);
      }
    
    // return make_shared<MultVecVecCoefficientFunction> (c1, c2);
  }


  shared_ptr<CoefficientFunction> CrossProduct (shared_ptr<CoefficientFunction> c1,
                                                shared_ptr<CoefficientFunction> c2)
  {
    if (c1->IsZeroCF() || c2->IsZeroCF())
      return ZeroCF( c1->Dimensions() );
    
    return make_shared<CrossProductCoefficientFunction> (c1, c2);
  }

  shared_ptr<CoefficientFunction> IdentityCF (int dim)
  {
     return make_shared<IdentityCoefficientFunction> (dim);
  }

  shared_ptr<CoefficientFunction> IdentityCF (FlatArray<int> dims)
  {

    if (dims.Size() == 0)
      return ConstantCF(1);

    int dim = 1;
    for (auto d : dims)
      dim *= d;
    
    Array<int> tensor_dims;
    tensor_dims.Append(dims);
    tensor_dims.Append(dims);
    return make_shared<IdentityCoefficientFunction>(dim)->Reshape(tensor_dims);
  }

  shared_ptr<CoefficientFunction> UnitVectorCF (int dim, int coord)
  {
    return make_shared<UnitVectorCoefficientFunction> (dim, coord);
  }


  shared_ptr<CoefficientFunction> ZeroCF (FlatArray<int> dims)
  {
    if (dims.Size() == 2)
      return make_shared<ZeroCoefficientFunction> (dims[0], dims[1]);
    else if (dims.Size() == 1)
      return make_shared<ZeroCoefficientFunction> (dims[0]);
    else if (dims.Size() == 0)
      return make_shared<ZeroCoefficientFunction> ();
    else
      return make_shared<ZeroCoefficientFunction> (Array<int>(dims));      
  }

  shared_ptr<CoefficientFunction> ConstantCF (double val)
  {
    return make_shared<ConstantCoefficientFunction>(val);
  }

  shared_ptr<CoefficientFunction> TransposeCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsZeroCF())
      {
        auto dims = coef->Dimensions();
        coef->SetDimensions( Array( {dims[1], dims[0]}) );
        return coef;
      }

    // if (coef.use_count() == 1)
    if (dynamic_pointer_cast<IdentityCoefficientFunction> (coef) && !coef->IsVariable())
      return coef;

    return make_shared<TransposeCoefficientFunction> (coef);
  }

  shared_ptr<CoefficientFunction> InverseCF (shared_ptr<CoefficientFunction> coef)
  {
    auto dims = coef->Dimensions();
    if (dims.Size() != 2) throw Exception("Inverse of non-matrix");
    if (dims[0] != dims[1]) throw Exception("Inverse of non-quadratic matrix");
    switch (dims[0])
      {
      case 1: return make_shared<InverseCoefficientFunction<1>> (coef);
      case 2: return make_shared<InverseCoefficientFunction<2>> (coef);
      case 3: return make_shared<InverseCoefficientFunction<3>> (coef);
      default:
        return make_shared<InverseCoefficientFunctionAnyDim>(coef);
      }
  }

  shared_ptr<CoefficientFunction> DeterminantCF (shared_ptr<CoefficientFunction> coef)
  {
    auto dims = coef->Dimensions();
    if (dims.Size() != 2) throw Exception("Inverse of non-matrix");
    if (dims[0] != dims[1]) throw Exception("Inverse of non-quadratic matrix");

    if (coef->IsZeroCF())
      return ZeroCF(Array<int>());

    if (dynamic_pointer_cast<IdentityCoefficientFunction> (coef) && !coef->IsVariable())
      return make_shared<ConstantCoefficientFunction>(1);


    // common pattern : Det (F^T F) = Det(F)**2, for F square
    if (!coef->IsVariable())
      if (auto mmm = dynamic_pointer_cast<MultMatMatCoefficientFunction> (coef))
        {
          auto AB = mmm->InputCoefficientFunctions();
          if (!AB[0]->IsVariable())
            if (auto trans = dynamic_pointer_cast<TransposeCoefficientFunction> (AB[0]))
              {
                auto At = trans->InputCoefficientFunctions()[0];
                if (At->Dimensions()[0] == At->Dimensions()[1])
                  {
                    auto detF = DeterminantCF (At);
                    return detF*detF;
                  }
              }
        }
    
    switch (dims[0])
      {
      case 1: return make_shared<DeterminantCoefficientFunction<1>> (coef);
      case 2: return make_shared<DeterminantCoefficientFunction<2>> (coef);
      case 3: return make_shared<DeterminantCoefficientFunction<3>> (coef);
      default:
        throw Exception("Determinant of matrix of size "+ToString(dims[0]) + " not available");
      }
  }

  shared_ptr<CoefficientFunction> CofactorCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsZeroCF())
      return coef;
    
    auto dims = coef->Dimensions();
    if (dims.Size() != 2) throw Exception("Cofactor of non-matrix");
    if (dims[0] != dims[1]) throw Exception("Cofactor of non-quadratic matrix");
    switch (dims[0])
      {
      case 1: return make_shared<CofactorCoefficientFunction<1>> (coef);
      case 2: return make_shared<CofactorCoefficientFunction<2>> (coef);
      case 3: return make_shared<CofactorCoefficientFunction<3>> (coef);
      case 4: return make_shared<CofactorCoefficientFunction<4>> (coef);
      default:
        throw Exception("Cofactor of matrix of size "+ToString(dims[0]) + " not available");
      }
  }

  shared_ptr<CoefficientFunction> SymmetricCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsZeroCF())
      return coef;

    return make_shared<SymmetricCoefficientFunction> (coef);
  }

  shared_ptr<CoefficientFunction> SkewCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsZeroCF())
      return coef;

    return make_shared<SkewCoefficientFunction> (coef);
  }

  shared_ptr<CoefficientFunction> TraceCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsZeroCF())
      return ZeroCF(Array<int>());
    
    return make_shared<TraceCoefficientFunction> (coef);
  }

  shared_ptr<CoefficientFunction> NormCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsZeroCF())
      return ZeroCF(Array<int>());
    
    if (coef->IsComplex())
      return make_shared<NormCoefficientFunctionC> (coef);
    else
      return make_shared<NormCoefficientFunction> (coef);
  }

  shared_ptr<CoefficientFunction> EigCF (shared_ptr<CoefficientFunction> coef)
  {
    return make_shared<EigCoefficientFunction> (coef);
  }
  
  shared_ptr<CoefficientFunction> operator/ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2)
  {
    if (c1->IsZeroCF())
      return c1;
    if (c2->Dimensions().Size() == 0 && c1->Dimensions().Size() > 0)
      return (make_shared<ConstantCoefficientFunction>(1.0)/c2)*c1;
    return BinaryOpCF (c1, c2, gen_div, "/");
  }

  shared_ptr<CoefficientFunction> operator/ (double val, shared_ptr<CoefficientFunction> c2)
  {
    return ConstantCF(val) / c2;
  }
  
class ComponentCoefficientFunction : public T_CoefficientFunction<ComponentCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  int dim1;
  int comp;
  typedef T_CoefficientFunction<ComponentCoefficientFunction> BASE;
public:
  ComponentCoefficientFunction() = default;
  ComponentCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                int acomp)
    : BASE(1, ac1->IsComplex()), c1(ac1), comp(acomp)
  {
    dim1 = c1->Dimension();
    elementwise_constant = c1->ElementwiseConstant();
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1) & dim1 & comp;
  }

  virtual string GetDescription () const override
  { return "ComponentCoefficientFunction " + ToString(comp); }

  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    /*
    auto dims = c1->Dimensions();
    int i,j;
    GetIndex(dims, comp, i, j);
    code.body += Var(index).Assign( Var(inputs[0], i, j ));
    */
    code.Declare (index, Dimensions(), IsComplex()); 
    code.body += Var(index).Assign( Var(inputs[0], comp, c1->Dimensions() ), false);    
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 }); }

  
  using BASE::Evaluate;
  /*
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    VectorMem<20> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    return v1(comp);
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const override
  {
    VectorMem<20> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    result(0) = v1(comp);
  }  
  
  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<Complex> result) const override
  {
    Vector<Complex> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    result(0) = v1(comp);
  }
  */

  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         BareSliceMatrix<Complex> result) const override
  {
    // int dim1 = c1->Dimension();
    STACK_ARRAY(double, hmem, 2*ir.Size()*dim1);
    FlatMatrix<Complex> temp(ir.Size(), dim1, (Complex*)hmem);
    c1->Evaluate (ir, temp);
    result.Col(0).Range(0,ir.Size()) = temp.Col(comp);
  }  

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    STACK_ARRAY(T, hmem, ir.Size()*dim1);
    FlatMatrix<T,ORD> temp(dim1, ir.Size(), &hmem[0]);
    
    c1->Evaluate (ir, temp);
    size_t nv = ir.Size();
    __assume(nv > 0);
    for (size_t i = 0; i < nv; i++)
      values(0,i) = temp(comp, i);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto in0 = input[0];    
    values.Row(0).Range(ir.Size()) = in0.Row(comp);
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return MakeComponentCoefficientFunction (c1->Diff(var, dir), comp);
  }  

  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];
    
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    auto diffc1 = c1->DiffJacobi (var, cache);
    Array<int> vardims = Array<int> (var->Dimensions());
    Array<int> dist(vardims.Size());
    int prod = 1;
    for (int i = dist.Size()-1; i >= 0; i--)
      {
        dist[i] = prod;
        prod *= vardims[i];
      }
    auto res = MakeSubTensorCoefficientFunction(diffc1, comp*var->Dimension(), move(vardims), move(dist));
    cache[thisptr] = res;
    return res;
  }  
  
  
  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    Vector<bool> v1(c1->Dimension()), d1(c1->Dimension()), dd1(c1->Dimension());
    c1->NonZeroPattern (ud, v1, d1, dd1);
    nonzero(0) = v1(comp);
    nonzero_deriv(0) = d1(comp);
    nonzero_dderiv(0) = dd1(comp);
  }  
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> v1(c1->Dimension());
    c1->NonZeroPattern (ud, v1);
    values(0) = v1(comp);
  }

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    values(0) = input[0](comp);
  }
};

shared_ptr<CoefficientFunction>
MakeComponentCoefficientFunction (shared_ptr<CoefficientFunction> c1, int comp)
{
  if (c1->IsZeroCF())
    return ZeroCF( Array<int>({}) );

  if (c1->GetDescription() == "unary operation ' '"  && !c1->IsVariable())
    return MakeComponentCoefficientFunction(c1->InputCoefficientFunctions()[0], comp);

  if (c1->GetDescription() == "VectorialCoefficientFunction")
  {
    int ci = 0;
    int dim = 0;
    auto coefs = c1->InputCoefficientFunctions();
    while(dim<=comp)
    {
      int dim_last = dim;
      if(dim == comp && coefs[ci]->Dimension() == 1)
        return coefs[ci];
      dim += coefs[ci]->Dimension();
      if(dim>comp)
        return MakeComponentCoefficientFunction( coefs[ci], comp-dim_last );
      ci++;
    }
  }
  return make_shared<ComponentCoefficientFunction> (c1, comp);
}


// ********************** SubTensorCoefficientFunction **********************

  
class SubTensorCoefficientFunction : public T_CoefficientFunction<SubTensorCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  int dim1;
  int first;
  Array<int> num, dist;
  typedef T_CoefficientFunction<SubTensorCoefficientFunction> BASE;
  Array<int> mapping; // output index -> input index
public:
  SubTensorCoefficientFunction() = default;
  SubTensorCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                int afirst, Array<int> anum, Array<int> adist)
    : BASE(1, ac1->IsComplex()), c1(ac1), first(afirst), num(anum), dist(adist)
  {
    SetDimensions(anum);
    dim1 = c1->Dimension();    
    elementwise_constant = c1->ElementwiseConstant();

    stringstream descr{};
    auto append_array_str = [&](const auto& array) {
      bool append{false};
      auto array_str = ToString(array);
      for (auto c : array_str)
        if (c == ':')
          append = true;
        else if (c == '\n')
          {
            descr << ',';
            append = false;
          }
        else if (append)
          descr << c;
        else
          continue;
    };

    descr << "subtensor [ first: " << first << ", num: (";
    append_array_str(num);
    descr << "), dist: (";
    append_array_str(dist);
    descr << ") ]";
    SetDescription(descr.str());
    
    for (int i = 0; i < Dimension(); i++)
      {
        int index = i;
        int inputindex = first;
        for (int j = num.Size()-1; j >= 0; j--)
          {
            int indexj = index % num[j];
            inputindex += indexj*dist[j];
            index /= num[j];
          }
        mapping.Append (inputindex);
      }
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1) & dim1 & first & num & dist & mapping;
  }

  /*
  virtual string GetDescription () const override
  { return "subtensor"; }
  */

  int First() const { return first; }
  FlatArray<int> Num() const { return num; }
  FlatArray<int> Dist() const { return dist; }
  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    // code.Declare (code.res_type, index, Dimensions());
    code.Declare (index, Dimensions(), IsComplex());
    auto dims1 = c1->Dimensions();
    for (auto i : Range(mapping))
      code.body += Var(index, i, num).Assign( Var(inputs[0], mapping[i], dims1), false);
  }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 }); }

  
  using BASE::Evaluate;

  /*
  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         BareSliceMatrix<Complex> result) const override
  {
    // int dim1 = c1->Dimension();
    STACK_ARRAY(double, hmem, 2*ir.Size()*dim1);
    FlatMatrix<Complex> temp(ir.Size(), dim1, (Complex*)hmem);
    c1->Evaluate (ir, temp);
    result.Col(0).Range(0,ir.Size()) = temp.Col(comp);
  }  
  */

  
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    STACK_ARRAY(T, hmem, ir.Size()*dim1);
    FlatMatrix<T,ORD> temp(dim1, ir.Size(), &hmem[0]);
    
    c1->Evaluate (ir, temp);
    for (size_t i = 0; i < mapping.Size(); i++)
      values.Row(i).Range(ir.Size()) = temp.Row(mapping[i]);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto in0 = input[0];
    for (auto i : Range(mapping))
      values.Row(i).Range(ir.Size()) = in0.Row(mapping[i]);
  }
  
  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return MakeSubTensorCoefficientFunction (c1->Diff(var, dir), first, Array<int> (num), Array<int> (dist));
  }  


  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];
    
    if (this == var)
      return IdentityCF(Dimensions());
    
    Array<int> dimres { this->Dimensions() };
    dimres += var->Dimensions();
    
    auto diffc1 = c1->DiffJacobi (var, cache);

    int vardim = var->Dimension();
    Array<int> vardims = Array<int> (var->Dimensions());
    Array<int> dist(vardims.Size());
    int prod = 1;
    for (int i = dist.Size()-1; i >= 0; i--)
      {
        dist[i] = prod;
        prod *= vardims[i];
      }
    Array<int> firstdist { this -> dist };
    for (auto & fd : firstdist)
      fd *= vardim;
      
    Array<int> alldist { firstdist + dist };
    
    auto res = MakeSubTensorCoefficientFunction(diffc1, first*vardim, move(dimres), move(alldist));
    cache[thisptr] = res;
    return res;
  }  
  

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> v1(c1->Dimension());
    c1->NonZeroPattern (ud, v1);
    for (auto i : Range(mapping))
      values(i) = v1(mapping[i]);
  }


  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    c1->NonZeroPattern (ud, input[0]);
    for (auto i : Range(mapping))
      values(i) = input[0](mapping[i]);
  }
};

shared_ptr<CoefficientFunction>
MakeSubTensorCoefficientFunction (shared_ptr<CoefficientFunction> c1, int first, Array<int> num, Array<int> dist)
{
  //further optimization possible by checking if only subtensor of c1 is zero..
  if (c1->IsZeroCF())
    return ZeroCF(num);

  // trivial sub-tensor ?
  const auto dims = c1->Dimensions();
  bool trivial = (first == 0) && (num.Size() == dims.Size()) && (num == dims);
  for (int i = 0; i+1 < dist.Size(); i++)
    if (dist[i] != num[i]*dist[i+1]) trivial = false;
  if (dist.Size() >= 1)
    if (dist.Last() != 1) trivial = false;

  if (trivial)
    {
      cout << IM(2) << "optimizing out trivial sub-tensor" << endl;
      return c1;
    }
  
  return make_shared<SubTensorCoefficientFunction> (c1, first, move(num), move(dist));
}


shared_ptr<CoefficientFunction>
MakeTensorTransposeCoefficientFunction (shared_ptr<CoefficientFunction> c1, Array<int> ordering)
{
  auto dims1 = c1->Dimensions();
  if (dims1.Size() != ordering.Size())
    throw Exception("TensorTranspose - tensor dimensions don't match");

  if (!c1->IsVariable())
    if (auto subcf = dynamic_pointer_cast<SubTensorCoefficientFunction>(c1))
      {
        cout << IM(2) << "Optimization: Tensor-Transpose of subtensor is a subtensor" << endl;
        Array<int> dims(dims1.Size()), dist(dims1.Size());
        for (int i = 0; i < dims.Size(); i++)
          {
            if (ordering[i] < 0 || ordering[i] >= dims1.Size())
              throw Exception ("ordering out of range");
            dims[i] = dims1[ordering[i]];
            dist[i] = subcf->Dist()[ordering[i]];
          }
        
        auto input = c1->InputCoefficientFunctions();
        auto res = MakeSubTensorCoefficientFunction (input[0], subcf->First(),
                                                     std::move(dims), std::move(dist));
        return res;
      }
  
  Array<int> dist1(dims1.Size());
  
  int disti = 1;
  for (int i = dims1.Size()-1; i >= 0; i--)
    {
      dist1[i] = disti;
      disti *= dims1[i];
        }
  
  Array<int> dims(dims1.Size()), dist(dims1.Size());
  for (int i = 0; i < dims.Size(); i++)
    {
      if (ordering[i] < 0 || ordering[i] >= dims1.Size())
        throw Exception ("ordering out of range");
      dims[i] = dims1[ordering[i]];
      dist[i] = dist1[ordering[i]];
    }
      
  auto res = MakeSubTensorCoefficientFunction (c1, 0, std::move(dims), std::move(dist));
  stringstream descr;
  descr << "tensor-transpose [";
  for (auto i : Range(ordering.Size() - 1))
    descr << " " << ordering[i] << ",";
  descr << " " << ordering.Last() << " ]";
  res -> SetDescription(descr.str());
  return res;
}

shared_ptr<CoefficientFunction>
MakeTensorTransposeCoefficientFunction (shared_ptr<CoefficientFunction> c1, int i1, int i2)
{
  Array<int> ia(c1->Dimensions().Size());
  for (int i = 0; i < ia.Size(); i++)
    ia[i] = i;
  
  if (i1 < 0 || i1 >= ia.Size())
    throw Exception ("TensorTranspose, i1 out of range");
  if (i2 < 0 || i2 >= ia.Size())
    throw Exception ("TensorTranspose, i2 out of range");
  
  ia[i1] = i2;
  ia[i2] = i1;

  return MakeTensorTransposeCoefficientFunction (c1, move(ia));
}


shared_ptr<CoefficientFunction>
MakeTensorTraceCoefficientFunction (shared_ptr<CoefficientFunction> c1, int i1, int i2)
{
  if (i1 < 0 || i1 >= c1->Dimensions().Size())
    throw Exception ("TensorTrace, i1 out of range");
  if (i2 < 0 || i2 >= c1->Dimensions().Size())
      throw Exception ("TensorTrace, i2 out of range");
  
  // a simple workaround for the beginning:
  if ( (i1 == 0 && i2 == 1) || (i1 == 1 && i2 == 0) )
    {
      Array<int> dims1 { c1->Dimensions() };
      auto dimsres = dims1.Range(2, END);

      // auto shapematT = TransposeCF (ReshapeCF (std::move(c1), sqr(dims1[0]), -1) );
      // auto prod = shapematT * ReshapeCF ( IdentityCF(dims1[0]), sqr (dims1[0]) );

      auto shapematT = c1 -> Reshape (sqr(dims1[0]), -1) -> Transpose(); 
      auto prod = shapematT * IdentityCF(dims1[0]) -> Reshape (sqr (dims1[0]));

      return prod -> Reshape (Array<int> (dimsres));
    }

  throw Exception("MakeTensorTraceCF not implemented for general case");
}



// ********************** ExtendDimensionCoefficientFunction ***************************

  
class ExtendDimensionCoefficientFunction : public T_CoefficientFunction<ExtendDimensionCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  typedef T_CoefficientFunction<ExtendDimensionCoefficientFunction> BASE;
  Array<int> mapping; // input index -> output index
  Array<int> dims, pos, stride;
  int dim1;
public:
  ExtendDimensionCoefficientFunction() = default;
  ExtendDimensionCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                  Array<int> adims, Array<int> apos, Array<int> astride)
    : BASE(1, ac1->IsComplex()), c1(ac1), dims(adims), pos(apos), stride(astride)
  {
    SetDimensions(dims);
    elementwise_constant = c1->ElementwiseConstant();

    auto dims1 = c1->Dimensions();
    dim1 = c1->Dimension();
      
    if (dims1.Size() != dims.Size())
      throw Exception("ExtendDimension needs same tensor dimension");

    for (int i = pos.Size(); i < dims.Size(); i++)
      pos.Append(0);

    if (stride.Size() > 0 && stride.Size() != dims.Size())
      throw Exception("stride must be either of size zero or the same size as dims");

    if (stride.Size() == 0)
      {
        stride.SetSize(dims.Size());
        stride = 1;
        for (int i = dims.Size() - 1; i >= 0; i--)
          for (int j = 0; j < i; ++j)
            stride[j] *= dims[i];
      }

    stringstream descr;
    descr << "extend-dimension [";
    descr << " input dims: ";
    for (auto i : Range(dims1.Size()-1))
      descr << dims1[i] << ", ";
    descr << dims1.Last() << " | ";
    descr << " pos: ";
    for (auto i : Range(pos.Size()-1))
      descr << pos[i] << ", ";
    descr << pos.Last() << " | ";
    descr << " stride: ";
    for (auto i : Range(stride.Size()-1))
      descr << stride[i] << ", ";
    descr << stride.Last() << " ]";
    SetDescription(descr.str());  

    int firstoutput = 0;
    for (int i = 0; i < dims.Size(); i++)
      firstoutput += pos[i] * stride[i];

    for (int i = 0; i < c1->Dimension(); i++)
      {
        int index = i;
        int outputindex = firstoutput;
        for (int j = dims1.Size()-1; j >= 0; j--)
          {
            int indexj = index % dims1[j];
            outputindex += indexj * stride[j];
            index /= dims1[j];
          }
        if (outputindex > Dimension())
          throw Exception("illegal output index "+ToString(outputindex));
        mapping.Append (outputindex);
      }
//    cout << "output indices = " << mapping << endl;
  }

  

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1) & mapping & dim1;
  }


  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    auto dims1 = c1->Dimensions();

    Array<int> imapping(Dimension());
    imapping = -1;
    for (auto i : mapping.Range())
      imapping[mapping[i]] = i;

    for (int i = 0; i < imapping.Size(); i++)
      if (imapping[i] == -1)
        code.body += Var(index, i, Dimensions()).Assign( string("0.0"));
      else
        code.body += Var(index, i, Dimensions()).Assign( Var(inputs[0], imapping[i], dims1));
  }

  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

//  virtual string GetDescription () const override
//  { return "extend-dimension"; }
  
  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 }); }

  
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    STACK_ARRAY(T, hmem, ir.Size()*dim1);
    FlatMatrix<T,ORD> temp(dim1, ir.Size(), &hmem[0]);
    
    c1->Evaluate (ir, temp);

    values.AddSize(Dimension(), ir.Size()) = T(0.0);
    for (size_t i = 0; i < mapping.Size(); i++)
      values.Row(mapping[i]).Range(ir.Size()) = temp.Row(i);
  }


  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto in0 = input[0];
    values.AddSize(Dimension(), ir.Size()) = T(0.0);
    
    for (size_t i = 0; i < mapping.Size(); i++)
      values.Row(mapping[i]).Range(ir.Size()) = in0.Row(i);
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return MakeExtendDimensionCoefficientFunction (c1->Diff(var, dir), Array<int> (dims), Array<int> (pos),
                                               Array<int>(stride));
  }

  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];
    
    if (this == var) return IdentityCF(this->Dimensions());

    Array<int> resdims;
    resdims += Dimensions();
    resdims += var->Dimensions();

    Array<int> resstride(resdims.Size());
    resstride.Range(0, dims.Size()) = stride;
    resstride.Range(dims.Size(), END) = 1;
    for (int i = resdims.Size() - 1; i >= dims.Size(); i--)
      for (int j = 0; j < i; ++j)
        resstride[j] *= resdims[i];

//    cout << "new stride: " << resstride << endl;
    auto diffc1 = c1->DiffJacobi (var, cache);
    auto res = MakeExtendDimensionCoefficientFunction (diffc1, move(resdims), Array<int>(pos), move(resstride));
    cache[thisptr] = res;
    return res;
  }

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> v1(c1->Dimension());
    c1->NonZeroPattern (ud, v1);
    values = false;
    for (int i = 0; i < mapping.Size(); i++)
      values[mapping[i]] = v1[i];
  }

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto in0 = input[0];
    values = false;
    for (int i = 0; i < mapping.Size(); i++)
      values[mapping[i]] = in0[i];
  }
  
};

shared_ptr<CoefficientFunction>
MakeExtendDimensionCoefficientFunction (shared_ptr<CoefficientFunction> c1,
                                    Array<int> dims, Array<int> pos, Array<int> stride)
{
  if (c1->IsZeroCF())
    return ZeroCF(dims);
  return make_shared<ExtendDimensionCoefficientFunction> (c1, move(dims), move(pos), move(stride));
}


// ********************** VectorContractionCoefficientFunction **********************

  
class VectorContractionCoefficientFunction : public T_CoefficientFunction<VectorContractionCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  typedef T_CoefficientFunction<VectorContractionCoefficientFunction> BASE;
  Array<shared_ptr<CoefficientFunction>> vectors;
public:
  VectorContractionCoefficientFunction() = default;
  VectorContractionCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                        Array<shared_ptr<CoefficientFunction>> avectors)
    : BASE(1, ac1->IsComplex()), c1(ac1), vectors(move(avectors))
  {
    elementwise_constant = c1->ElementwiseConstant();
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1) & vectors;
  }

  /*
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    auto dims1 = c1->Dimensions();

    if(num.Size()==1)
      {
        for (int i = 0; i < num[0]; i++)
         {
            int i1,k1;
            auto comp = first+i*dist[0];
            GetIndex(dims1, comp, i1, k1);
            code.body += Var(index, i).Assign( Var(inputs[0], i1, k1 ));
          }
      }

    if(num.Size()==2)
      {
        for (int i = 0; i < num[0]; i++)
          for (int j = 0; j < num[1]; j++)
             {
                int i1,j1;
                auto comp = first+i*dist[0]+j*dist[1];
                GetIndex(dims1, comp, i1, j1);
                code.body += Var(index, i, j).Assign( Var(inputs[0], i1, j1 ));
              }
      }
  }
  */
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    for (auto v : vectors)
      v->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  {
    
    Array<shared_ptr<CoefficientFunction>> inputs;
    inputs.Append(c1);
    inputs += vectors;
    return inputs;
  }
  
  using BASE::Evaluate;

  
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    STACK_ARRAY(T, hmem, ir.Size()*c1->Dimension());
    FlatMatrix<T,ORD> temp(c1->Dimension(), ir.Size(), &hmem[0]);
    STACK_ARRAY(T, hmem2, ir.Size()*c1->Dimension());  // 
    
    c1->Evaluate (ir, temp);

    size_t nv = ir.Size();
    size_t actdim = c1->Dimension();

    for (int dir = 0; dir < vectors.Size(); dir++)
      {
        FlatMatrix<T,ORD> vi(vectors[dir]->Dimension(), ir.Size(), &hmem2[0]);
        vectors[dir]->Evaluate(ir, vi);

        size_t newdim = actdim / vi.Height();
        for (size_t i = 0; i < newdim; i++)
          temp.Row(i) = pw_mult(temp.Row(i), vi.Row(0));
        for (size_t j = 1; j < vi.Height(); j++)
          for (size_t i = 0; i < newdim; i++)
            temp.Row(i) += pw_mult(temp.Row(j*newdim+i), vi.Row(j));

        /*
        for (size_t k = 0; k < nv; k++)
          for (size_t i = 0; i < newdim; i++)
            {
              T sum{0.0};
              for (size_t j = 0; j < vi.Height(); j++)
                sum += temp(i+j*newdim, k) * vi(j,k);
              temp(i,k) = sum;
            }
        */
        actdim = newdim;
      }
    for (size_t k = 0; k < nv; k++)
      values(0,k) = temp(0,k);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    STACK_ARRAY(T, hmem, ir.Size()*c1->Dimension());
    FlatMatrix<T,ORD> temp(c1->Dimension(), ir.Size(), &hmem[0]);
    
    // c1->Evaluate (ir, temp);
    temp = input[0];

    size_t nv = ir.Size();
    size_t actdim = c1->Dimension();

    for (int dir = 0; dir < vectors.Size(); dir++)
      {
        // FlatMatrix<T,ORD> vi(vectors[dir]->Dimension(), ir.Size(), &hmem2[0]);
        // vectors[dir]->Evaluate(ir, vi);
        auto vi = input[dir+1].AddSize(vectors[dir]->Dimension(), ir.Size());
          
        size_t newdim = actdim / vi.Height();
        for (size_t i = 0; i < newdim; i++)
          temp.Row(i) = pw_mult(temp.Row(i), vi.Row(0));
        for (size_t j = 1; j < vi.Height(); j++)
          for (size_t i = 0; i < newdim; i++)
            temp.Row(i) += pw_mult(temp.Row(j*newdim+i), vi.Row(j));

        actdim = newdim;
      }
    for (size_t k = 0; k < nv; k++)
      values(0,k) = temp(0,k);
  }

  /*
  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return MakeVectorContractionCoefficientFunction (c1->Diff(var, dir), first, Array<int> (num), Array<int> (dist));
  }  

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    Vector<AutoDiffDiff<1,bool>> v1(c1->Dimension());
    c1->NonZeroPattern (ud, v1);
    switch (num.Size())
      {
      case 1:
        for (int i = 0; i < num[0]; i++)
          values(i) = v1(first+i*dist[0]);
        break;
      case 2:
        for (int i = 0, ii = 0; i < num[0]; i++)
          for (int j = 0; j < num[1]; j++, ii++)
            values(ii) = v1(first+i*dist[0]+j*dist[1]);
        break;

      default:
        throw Exception("subtensor of order "+ToString(num.Size())+" not supported");
      }
  }

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    c1->NonZeroPattern (ud, values);
    switch (num.Size())
      {
      case 1:
        for (int i = 0; i < num[0]; i++)
          values(i) = values(first+i*dist[0]);
        break;
      case 2:
        for (int i = 0, ii = 0; i < num[0]; i++)
          for (int j = 0; j < num[1]; j++, ii++)
            values(ii) = values(first+i*dist[0]+j*dist[1]);
        break;
      default:
        throw Exception("subtensor of order "+ToString(num.Size())+" not supported");
      }
  }
  */
};

shared_ptr<CoefficientFunction>
MakeVectorContractionCoefficientFunction (shared_ptr<CoefficientFunction> c1,
                                          Array<shared_ptr<CoefficientFunction>> vectors)
{
  return make_shared<VectorContractionCoefficientFunction> (c1, move(vectors));
}





// ********************** SingleContractionCoefficientFunction **********************

  
class SingleContractionCoefficientFunction : public T_CoefficientFunction<SingleContractionCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  typedef T_CoefficientFunction<SingleContractionCoefficientFunction> BASE;
  shared_ptr<CoefficientFunction> vec;
  int index, dimbefore, dimafter;
public:
  SingleContractionCoefficientFunction() = default;
  SingleContractionCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                        shared_ptr<CoefficientFunction> avec,
                                        int aindex)
    : BASE(ac1->Dimension()/avec->Dimension(), ac1->IsComplex()),
      c1(ac1), vec(avec), index(aindex)
  {
    elementwise_constant = c1->ElementwiseConstant() && vec->ElementwiseConstant();
    dimbefore = 1;
    dimafter = 1;
    Array<int> dims;
    dims.SetSize(c1->Dimensions().Size()-1);
    for (int j = 0; j < index; j++)
      {
        dims[j] = c1->Dimensions()[j];
        dimbefore *= dims[j];
      }
    for (int j = index; j < dims.Size(); j++)
      {
        dims[j] = c1->Dimensions()[j+1];
        dimafter *= dims[j];
      }
    SetDimensions(dims);
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1) & vec;
  }

  /*
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  */
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    vec->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1, vec }); }
  
  using BASE::Evaluate;

  
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> hvalues) const
  {
    STACK_ARRAY(T, hmem, ir.Size()*c1->Dimension());
    FlatMatrix<T,ORD> temp(c1->Dimension(), ir.Size(), &hmem[0]);
    c1->Evaluate (ir, temp);

    STACK_ARRAY(T, hmem2, ir.Size()*vec->Dimension());
    FlatMatrix<T,ORD> vi(vec->Dimension(), ir.Size(), &hmem2[0]);
    vec->Evaluate(ir, vi);

    size_t nv = ir.Size();

    auto values = hvalues.AddSize(Dimension(), nv);
    values = T(0);

    for (int i = 0, ii = 0; i < dimbefore; i++)
      for (int j = 0; j < vec->Dimension(); j++)
        for (int k = 0; k < dimafter; k++, ii++)
          values.Row(i*dimafter+k) += pw_mult(vi.Row(j), temp.Row(ii));
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> hvalues) const
  {
    auto temp = input[0];
    auto vi = input[1];
    
    size_t nv = ir.Size();

    auto values = hvalues.AddSize(Dimension(), nv);
    values = T(0);

    for (int i = 0, ii = 0; i < dimbefore; i++)
      for (int j = 0; j < vec->Dimension(); j++)
        for (int k = 0; k < dimafter; k++, ii++)
          values.Row(i*dimafter+k) += pw_mult(vi.Row(j), temp.Row(ii));
  }

  /*
  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return MakeVectorContractionCoefficientFunction (c1->Diff(var, dir), first, Array<int> (num), Array<int> (dist));
  }  

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
  }

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
  }
  */
};

shared_ptr<CoefficientFunction>
MakeSingleContractionCoefficientFunction (shared_ptr<CoefficientFunction> c1,
                                          shared_ptr<CoefficientFunction> vec,
                                          int index)
{
  return make_shared<SingleContractionCoefficientFunction> (c1, vec, index);
}










// ************************ DomainWiseCoefficientFunction *************************************

class DomainWiseCoefficientFunction : public T_CoefficientFunction<DomainWiseCoefficientFunction>
{
  Array<shared_ptr<CoefficientFunction>> ci;
  typedef T_CoefficientFunction<DomainWiseCoefficientFunction> BASE;
  using BASE::Evaluate;
public:
  DomainWiseCoefficientFunction() = default;
  DomainWiseCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
    : BASE(1, false), ci(aci) 
  { 
    for (auto & cf : ci)
      if (cf && cf->IsComplex()) is_complex = true;
    for (auto & cf : ci)
      if (cf)
        SetDimensions(cf->Dimensions());

    elementwise_constant = true;
    for (auto cf : ci)
      if (cf && !cf->ElementwiseConstant())
        elementwise_constant = false;
  }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    auto size = ci.Size();
    ar & size;
    ci.SetSize(size);
    for(auto& cf : ci)
      ar.Shallow(cf);
  }

  virtual bool DefinedOn (const ElementTransformation & trafo) override
  {
    int matindex = trafo.GetElementIndex();
    return (matindex < ci.Size() && ci[matindex]);
  }

  /*
  bool ElementwiseConstant() const override
  {
    for(auto cf : ci)
      if(cf && !cf->ElementwiseConstant())
        return false;
    return true;
  }
  */
  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    code.body += "// DomainWiseCoefficientFunction:\n";
    string type = "decltype(0.0";
    for(int in : inputs)
      type += "+decltype("+Var(in,0,Dimensions()).S()+")()";
    type += ")";
    for (int i = 0; i < Dimension(); i++)
      code.body += Var(index,i,this->Dimensions()).Declare(type);
    code.body += "switch(domain_index) {\n";
    for(int domain : Range(inputs))
    {
        code.body += "case " + ToLiteral(domain) + ": \n";
        for (int i = 0; i < Dimension(); i++)
          code.body += "  "+Var(index, i, Dimensions()).Assign(Var(inputs[domain], i, Dimensions()), false);          
        code.body += "  break;\n";
    }
    code.body += "default: \n";
    for (int i = 0; i < Dimension(); i++)
      code.body += "  "+Var(index, i, Dimensions()).Assign(string("0.0"), false);      
    code.body += "  break;\n";
    code.body += "}\n";
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    for (auto & cf : ci)
      if (cf)
        cf->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  {
    Array<shared_ptr<CoefficientFunction>> cfa;
    for (auto cf : ci)
      cfa.Append (cf);
    return Array<shared_ptr<CoefficientFunction>>(cfa);
  } 

  shared_ptr<CoefficientFunction> Operator (const string & name) const override
  {
    Array<shared_ptr<CoefficientFunction>> cfop;
    
    for (auto & cf : ci)
      if (cf)
        cfop.Append (cf->Operator(name));
      else
        cfop.Append (nullptr);
    return MakeDomainWiseCoefficientFunction(move (cfop));
  }
  
  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    Array<shared_ptr<CoefficientFunction>> ci_deriv;
    for (auto & cf : ci)
      if (cf)
        ci_deriv.Append (cf->Diff(var, dir));
      else
        ci_deriv.Append (nullptr);
    return MakeDomainWiseCoefficientFunction(move (ci_deriv));
  }  
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    result = 0;
    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ip, result);
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const override
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ir, values);
    else
      values.AddSize(ir.Size(), Dimension()) = 0.0;
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ir, values);
    else
      values.AddSize(Dimension(), ir.Size()) = T(0.0);
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      values.AddSize(Dimension(), ir.Size()) = input[matindex];
    else
      values.AddSize(Dimension(), ir.Size()) = T(0.0);
  }  

  
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    result = 0;
    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ip, result);
  }
  
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const override
  {
    Vec<1,Complex> res;
    Evaluate (ip, res);
    return res(0);
  }

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv,
                               FlatVector<bool> nonzero_dderiv) const override
  {
    size_t dim = Dimension();
    STACK_ARRAY(bool, mem, 3*dim);
    FlatVector<bool> nzi(dim, &mem[0]);
    FlatVector<bool> nzdi(dim, &mem[dim]);
    FlatVector<bool> nzddi(dim, &mem[2*dim]);
    nonzero = false;
    nonzero_deriv = false;
    nonzero_dderiv = false;
    for (auto & aci : ci)
      if (aci)
        {
          aci -> NonZeroPattern(ud, nzi, nzdi, nzddi);
          for (size_t i = 0; i < nonzero.Size(); i++)
            {
              nonzero(i) |= nzi(i);
              nonzero_deriv(i) |= nzdi(i);
              nonzero_dderiv(i) |= nzddi(i);
            }
        }
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override 
  {
    size_t dim = Dimension();
    Vector<AutoDiffDiff<1,bool>> nzi(dim);
    values = AutoDiffDiff<1,bool> (false);
    for (auto & aci : ci)
      if (aci)
        {
          aci -> NonZeroPattern(ud, nzi);
          for (size_t i = 0; i < values.Size(); i++)
            values(i) += nzi(i);
        }
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override 
  {
    values = AutoDiffDiff<1,bool> (false);
    for (auto j : Range(input))
      if (ci[j])
        for (size_t i = 0; i < values.Size(); i++)
          values(i) += input[j](i);
  }
};

  shared_ptr<CoefficientFunction>
  MakeDomainWiseCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
  {
    for(auto cf : aci)
      if (cf && !cf->IsZeroCF())
        return make_shared<DomainWiseCoefficientFunction> (move (aci));
    for(auto cf : aci)
      if(cf)
        return ZeroCF(cf->Dimensions());
    return nullptr;
  }






// ************************ OtherCoefficientFunction *************************************

class OtherCoefficientFunction : public T_CoefficientFunction<OtherCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  typedef T_CoefficientFunction<OtherCoefficientFunction> BASE;
  using BASE::Evaluate;
public:
  OtherCoefficientFunction() = default;
  OtherCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : BASE(ac1->Dimension(), ac1->IsComplex()), c1(ac1)
  { SetDimensions(ac1->Dimensions()); }

  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    ar.Shallow(c1);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    throw Exception ("OtherCF::GenerateCode not available");
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  {
    Array<shared_ptr<CoefficientFunction>> cfa;
    cfa.Append (c1);
    return Array<shared_ptr<CoefficientFunction>>(cfa);
  } 
  
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    throw Exception ("OtherCF::Evaluated (mip) not available");    
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    throw Exception ("OtherCF::Evaluated (mip) not available");        
  }

  /*
  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
    if (!ir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");
    c1->Evaluate (*ir.GetOtherMIR(), values);
  }
  */
  
  virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const override
  {
    if (!ir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");    
    c1->Evaluate (*ir.GetOtherMIR(), values);    
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    if (!ir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");    
    c1->Evaluate (*ir.GetOtherMIR(), values);    
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    if (!ir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");    
    c1->Evaluate (*ir.GetOtherMIR(), values);    
  }

  /*
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const 
  {
    // compile not available
    if (!ir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");    
    c1->Evaluate (*ir.GetOtherMIR(), values);        
  }
  */
  
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    throw Exception ("OtherCF::Evaluated (mip) not available");        
  }
  
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const override
  {
    throw Exception ("OtherCF::Evaluated (mip) not available");            
  }

  /*
  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    if (!mir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");    
    c1->EvaluateDeriv (*mir.GetOtherMIR(), result, deriv);            
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    if (!mir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");    
    c1->EvaluateDDeriv (*mir.GetOtherMIR(), result, deriv, dderiv);                
  }
  */ 
};

shared_ptr<CoefficientFunction>
MakeOtherCoefficientFunction (shared_ptr<CoefficientFunction> me)
{
  me->TraverseTree
    ( [&] (CoefficientFunction & nodecf)
      {
        if (dynamic_cast<const ProxyFunction*> (&nodecf))
          throw Exception ("Other() can be applied either to a proxy, or to an expression without any proxy\n  ---> use the Other()-operator on sub-trees");
      }
      );
  return make_shared<OtherCoefficientFunction> (me);
}

bool IsOtherCoefficientFunction (CoefficientFunction & coef)
{
  return (dynamic_cast<OtherCoefficientFunction*> (&coef) != nullptr);
}






  

  // ///////////////////////////// IfPos   ////////////////////////////////  

  
class IfPosCoefficientFunction : public T_CoefficientFunction<IfPosCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> cf_if;
    shared_ptr<CoefficientFunction> cf_then;
    shared_ptr<CoefficientFunction> cf_else;
    typedef T_CoefficientFunction<IfPosCoefficientFunction> BASE;
  public:
    IfPosCoefficientFunction() = default;
    IfPosCoefficientFunction (shared_ptr<CoefficientFunction> acf_if,
                              shared_ptr<CoefficientFunction> acf_then,
                              shared_ptr<CoefficientFunction> acf_else)
      : BASE(acf_then->Dimension(),
             acf_then->IsComplex() || acf_else->IsComplex()),
        cf_if(acf_if), cf_then(acf_then), cf_else(acf_else)
    {
      if (acf_then->Dimension() != acf_else->Dimension())
        throw Exception(string("In IfPosCoefficientFunction: dim(cf_then) == ") + ToLiteral(acf_then->Dimension()) + string(" != dim(cf_else) == ") + ToLiteral(acf_else->Dimension()));
      SetDimensions(cf_then->Dimensions());
    }

    void DoArchive(Archive& ar) override
    {
      BASE::DoArchive(ar);
      ar.Shallow(cf_if).Shallow(cf_then).Shallow(cf_else);
    }

    virtual ~IfPosCoefficientFunction () { ; }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      Vec<1> val;
      cf_if->Evaluate(ip, val);
      if (val(0) > 0)
        return cf_then->Evaluate(ip);
      else
        return cf_else->Evaluate(ip);      
    }

    virtual void Evaluate (const BaseMappedIntegrationPoint& ip, FlatVector<double> values) const override
    {
      Vec<1> val;
      cf_if->Evaluate(ip, val);
      if(val(0) > 0)
        cf_then->Evaluate(ip,values);
      else
        cf_else->Evaluate(ip,values);
    }

    template <typename MIR, typename T, ORDERING ORD>    
    void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
    {
      size_t np = ir.Size();
      size_t dim = Dimension();
      
      STACK_ARRAY(T, hmem1, np);
      FlatMatrix<T,ORD> if_values(1, np, hmem1);
      STACK_ARRAY(T, hmem2, np*dim);
      FlatMatrix<T,ORD> then_values(dim, np, hmem2);
      STACK_ARRAY(T, hmem3, np*dim);
      FlatMatrix<T,ORD> else_values(dim, np, hmem3);
      
      cf_if->Evaluate (ir, if_values);
      cf_then->Evaluate (ir, then_values);
      cf_else->Evaluate (ir, else_values);
      
      for (size_t i = 0; i < np; i++)
        for (size_t j = 0; j < dim; j++)
          values(j,i) = IfPos(if_values(0,i), then_values(j,i), else_values(j,i));
    }
    
    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    {
      size_t np = ir.Size();
      size_t dim = Dimension();

      auto if_values = input[0];
      auto then_values = input[1];
      auto else_values = input[2];
      
      for (size_t i = 0; i < np; i++)
        for (size_t j = 0; j < dim; j++)
          values(j,i) = IfPos(if_values(0,i), then_values(j,i), else_values(j,i));
    }
    
    /*
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> hvalues) const override
    {
      auto values = hvalues.AddSize(ir.Size(), Dimension());
      
      STACK_ARRAY(double, hmem1, ir.Size());
      FlatMatrix<> if_values(ir.Size(), 1, hmem1);
      STACK_ARRAY(double, hmem2, ir.Size()*values.Width());
      FlatMatrix<> then_values(ir.Size(), values.Width(), hmem2);
      STACK_ARRAY(double, hmem3, ir.Size()*values.Width());
      FlatMatrix<> else_values(ir.Size(), values.Width(), hmem3);
      
      cf_if->Evaluate (ir, if_values);
      cf_then->Evaluate (ir, then_values);
      cf_else->Evaluate (ir, else_values);
      
      for (int i = 0; i < ir.Size(); i++)
        if (if_values(i) > 0)
          values.Row(i) = then_values.Row(i);
        else
          values.Row(i) = else_values.Row(i);

      // for (int i = 0; i < ir.Size(); i++)
      //   values(i) = (if_values(i) > 0) ? then_values(i) : else_values(i);
    }
    */
    
      using T_CoefficientFunction<IfPosCoefficientFunction>::Evaluate;
    virtual void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<Complex> values) const override
    {
      if(cf_if->Evaluate(ip)>0)
        cf_then->Evaluate(ip,values);
      else
        cf_else->Evaluate(ip,values);
    }

    /*
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override
    {
      size_t nv = ir.Size(), dim = Dimension();
      STACK_ARRAY(SIMD<double>, hmem1, nv);
      ABareMatrix<double> if_values(&hmem1[0], nv);
      STACK_ARRAY(SIMD<double>, hmem2, nv*dim);
      ABareMatrix<double> then_values(&hmem2[0], nv);
      STACK_ARRAY(SIMD<double>, hmem3, nv*dim);
      ABareMatrix<double> else_values(&hmem3[0], nv);
      
      cf_if->Evaluate (ir, if_values);
      cf_then->Evaluate (ir, then_values);
      cf_else->Evaluate (ir, else_values);
      for (size_t k = 0; k < dim; k++)
        for (size_t i = 0; i < nv; i++)
          values(k,i) = ngstd::IfPos (if_values.Get(i),
                                      then_values.Get(k,i),
                                      else_values.Get(k,i));
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const override
    {
      size_t nv = ir.Size(), dim = Dimension();      
      auto if_values = *input[0];
      auto then_values = *input[1];
      auto else_values = *input[2];
      
      for (size_t k = 0; k < dim; k++)
        for (size_t i = 0; i < nv; i++)
          values.Get(k,i) = ngstd::IfPos (if_values.Get(i),
                                          then_values.Get(k,i),
                                          else_values.Get(k,i)); 
    }
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatArray<FlatMatrix<>*> input,
                           FlatMatrix<double> values) const 
    {
      FlatMatrix<> if_values = *input[0];
      FlatMatrix<> then_values = *input[1];
      FlatMatrix<> else_values = *input[2];
      for (int i = 0; i < if_values.Height(); i++)
        values.Row(i) = (if_values(i) > 0) ? then_values.Row(i) : else_values.Row(i);
    }
    */

    // virtual bool IsComplex() const { return cf_then->IsComplex() | cf_else->IsComplex(); }
    // virtual int Dimension() const { return cf_then->Dimension(); }

    void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
    {
      /*
      auto cast_value = [&] (int input, int i, int j) {
          return code.res_type + "(" + Var(inputs[input], i, j).S() + ")";
      };
      */
      auto cast_value = [&] (int input, int i, FlatArray<int> dims) {
        // return code.res_type + "(" + Var(inputs[input], i, dims).S() + ")";
        return code.GetType(this->IsComplex()) + "(" + Var(inputs[input], i, dims).S() + ")";
      };

      auto var_if = Var(inputs[0]);
      // for (int i = 0; i < cf_then->Dimension(); i++)
      // code.body += Var(index,i,cf_then->Dimensions()).Declare(code.res_type);        
      // code.Declare (code.res_type, index, Dimensions());
      code.Declare (index, Dimensions(), IsComplex());      
      
      if(code.is_simd) {
        for (int i = 0; i < cf_then->Dimension(); i++)
            // cast all input parameters of IfPos to enforce the right overload (f.i. intermediate results could be double instead of AutoDiff<>)
          code.body += Var(index,i,cf_then->Dimensions()).Assign("IfPos("+cast_value(0, 0, cf_if->Dimensions()) + ',' + cast_value(1, i, cf_then->Dimensions()) +
                                                                 ',' + cast_value(2, i, cf_else->Dimensions())+')', false);
      } else {
        code.body += "if (" + var_if.S() + ">0.0) {\n";
        for (int i = 0; i < cf_then->Dimension(); i++)
          code.body += Var(index,i,cf_then->Dimensions()).Assign( Var(inputs[1],i,cf_then->Dimensions()), false );          
        code.body += "} else {\n";
        for (int i = 0; i < cf_then->Dimension(); i++)
          code.body += Var(index,i,cf_then->Dimensions()).Assign( Var(inputs[2],i,cf_then->Dimensions()), false );          
        
        code.body += "}\n";
      }
    }

    /*
    virtual Array<int> Dimensions() const
    {
      return cf_then->Dimensions();
    }
    */

    /*
    [[deprecated]]
    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                FlatMatrix<> values,
                                FlatMatrix<> deriv) const override
    {
      STACK_ARRAY(double, hmem1, ir.Size());
      FlatMatrix<> if_values(ir.Size(), 1, hmem1);
      STACK_ARRAY(double, hmem2, ir.Size()*values.Width());
      FlatMatrix<> then_values(ir.Size(), values.Width(), hmem2);
      STACK_ARRAY(double, hmem3, ir.Size()*values.Width());
      FlatMatrix<> else_values(ir.Size(), values.Width(), hmem3);
      STACK_ARRAY(double, hmem4, ir.Size()*values.Width());
      FlatMatrix<> then_deriv(ir.Size(), values.Width(), hmem4);
      STACK_ARRAY(double, hmem5, ir.Size()*values.Width());
      FlatMatrix<> else_deriv(ir.Size(), values.Width(), hmem5);

      
      cf_if->Evaluate (ir, if_values);
      cf_then->EvaluateDeriv (ir, then_values, then_deriv);
      cf_else->EvaluateDeriv (ir, else_values, else_deriv);
      
      for (int i = 0; i < ir.Size(); i++)
        if (if_values(i) > 0)
          {
            values.Row(i) = then_values.Row(i);
            deriv.Row(i) = then_deriv.Row(i);
          }
        else
          {
            values.Row(i) = else_values.Row(i);
            deriv.Row(i) = else_deriv.Row(i);
          }
    }
    */
    
    /*
    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                FlatMatrix<Complex> result,
                                FlatMatrix<Complex> deriv) const
    {
      Evaluate (ir, result);
      deriv = 0;
    }

    virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatMatrix<> result,
                                 FlatMatrix<> deriv,
                                 FlatMatrix<> dderiv) const
    {
      EvaluateDeriv (ir, result, deriv);
      dderiv = 0;
    }

    virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatMatrix<Complex> result,
                                 FlatMatrix<Complex> deriv,
                                 FlatMatrix<Complex> dderiv) const
    {
      EvaluateDeriv (ir, result, deriv);
      dderiv = 0;
    }

    
    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatArray<FlatMatrix<>*> input,
                                 FlatArray<FlatMatrix<>*> dinput,
                                 FlatMatrix<> result,
                                 FlatMatrix<> deriv) const
    {
      EvaluateDeriv (ir, result, deriv);
    }

    virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatArray<FlatMatrix<>*> input,
                                 FlatArray<FlatMatrix<>*> dinput,
                                 FlatArray<FlatMatrix<>*> ddinput,
                                 FlatMatrix<> result,
                                 FlatMatrix<> deriv,
                                 FlatMatrix<> dderiv) const
    {
      EvaluateDDeriv (ir, result, deriv, dderiv);
    }
    */

    // virtual bool ElementwiseConstant () const { return false; }
    
    // virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const;

    /*
    virtual void PrintReport (ostream & ost) const;
    virtual void PrintReportRec (ostream & ost, int level) const;
    virtual string GetName () const;
    */

    /*
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                AFlatMatrix<> values,
                                AFlatMatrix<> deriv) const
    {
      int dim = Dimension();
      STACK_ARRAY(SIMD<double>, hmem1, ir.Size());
      AFlatMatrix<> if_values(1, ir.IR().GetNIP(), hmem1);
      STACK_ARRAY(SIMD<double>, hmem2, ir.Size()*dim);
      AFlatMatrix<> then_values(dim, ir.IR().GetNIP(), hmem2);
      STACK_ARRAY(SIMD<double>, hmem3, ir.Size()*dim);
      AFlatMatrix<> else_values(dim, ir.IR().GetNIP(), hmem3);
      STACK_ARRAY(SIMD<double>, hmem4, ir.Size()*dim);
      AFlatMatrix<> then_deriv(dim, ir.IR().GetNIP(), hmem4);
      STACK_ARRAY(SIMD<double>, hmem5, ir.Size()*dim);
      AFlatMatrix<> else_deriv(dim, ir.IR().GetNIP(), hmem5);

      cf_if->Evaluate (ir, if_values);
      cf_then->EvaluateDeriv (ir, then_values, then_deriv);
      cf_else->EvaluateDeriv (ir, else_values, else_deriv);
      
      for (int i = 0; i < ir.Size(); i++)
        for (int j = 0; j < dim; j++)
          {
            values.Get(j,i) = IfPos(if_values.Get(0,i), then_values.Get(j,i), else_values.Get(j,i));
            deriv.Get(j,i) = IfPos(if_values.Get(0,i), then_deriv.Get(j,i), else_deriv.Get(j,i));
          }
    }

    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                FlatArray<AFlatMatrix<>*> input,
                                FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result,
                                AFlatMatrix<> deriv) const
    {
      size_t nv = ir.Size(), dim = Dimension();      
      auto if_values = *input[0];
      auto then_values = *input[1];
      auto else_values = *input[2];
      auto then_deriv = *dinput[1];
      auto else_deriv = *dinput[2];
      
      for (size_t k = 0; k < dim; k++)
        for (size_t i = 0; i < nv; i++)
          {
            result.Get(k,i) = ngstd::IfPos (if_values.Get(i),
                                            then_values.Get(k,i),
                                            else_values.Get(k,i));
            deriv.Get(k,i) = ngstd::IfPos (if_values.Get(i),
                                           then_deriv.Get(k,i),
                                           else_deriv.Get(k,i));
          }
    }
    */
    
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      cf_if->TraverseTree (func);
      cf_then->TraverseTree (func);
      cf_else->TraverseTree (func);
      func(*this);
    }
    
    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    {
      return Array<shared_ptr<CoefficientFunction>>( { cf_if, cf_then, cf_else } );
    }

    /*
    virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                                 FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
    {
      int dim = Dimension();
      Vector<bool> v1(dim), d1(dim), dd1(dim);
      Vector<bool> v2(dim), d2(dim), dd2(dim);
      cf_then->NonZeroPattern (ud, v1, d1, dd1);
      cf_else->NonZeroPattern (ud, v2, d2, dd2);
      for (int i = 0; i < dim; i++)
        {
          nonzero(i) = v1(i) || v2(i);
          nonzero_deriv(i) = d1(i) || d2(i);
          nonzero_dderiv(i) = dd1(i) || dd2(i);
        }
    }  
    */

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<AutoDiffDiff<1,bool>> values) const override
    {
      int dim = Dimension();
      Vector<AutoDiffDiff<1,bool>> v1(dim), v2(dim);
      cf_then->NonZeroPattern (ud, v1);
      cf_else->NonZeroPattern (ud, v2);
      values = v1+v2;
    }
    
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                                 FlatVector<AutoDiffDiff<1,bool>> values) const override
    {
      auto v1 = input[1];
      auto v2 = input[2];
      values = v1+v2;
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (this == var) return dir;
      return IfPos (cf_if, cf_then->Diff(var, dir), cf_else->Diff(var, dir));
    }

    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];
      
      if (this == var)
        return IdentityCF(cf_then->Dimensions());
      
      auto res = IfPos (cf_if, cf_then->DiffJacobi (var, cache), cf_else->DiffJacobi (var, cache));
      cache[thisptr] = res;
      return res;
    }
  };
  
  extern
  shared_ptr<CoefficientFunction> IfPos (shared_ptr<CoefficientFunction> cf_if,
                                         shared_ptr<CoefficientFunction> cf_then,
                                         shared_ptr<CoefficientFunction> cf_else)
  {
    if(cf_if->Dimension() != 1)
      throw Exception("Dimension of first component in IfPos must be 1!");
    if (cf_then->IsZeroCF() && cf_else->IsZeroCF())
      return cf_then;
    return make_shared<IfPosCoefficientFunction> (cf_if, cf_then, cf_else);
  }


class VectorialCoefficientFunction : public T_CoefficientFunction<VectorialCoefficientFunction>
{
  Array<shared_ptr<CoefficientFunction>> ci;
  Array<size_t> dimi;  // dimensions of components
  typedef T_CoefficientFunction<VectorialCoefficientFunction> BASE;
public:
  VectorialCoefficientFunction() = default;
  VectorialCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
    : BASE(0, false), ci(aci), dimi(aci.Size())
  {
    int hdim = 0;
    for (int i : Range(ci))
      {
        dimi[i] = ci[i]->Dimension();
        hdim += dimi[i];
      }
    
    for (auto cf : ci)
      if (cf && cf->IsComplex())
        is_complex = true;

    SetDimension(hdim);

    elementwise_constant = true;
    for (auto cf : ci)
      if (!cf->ElementwiseConstant())
        elementwise_constant = false;
    // dims = Array<int> ( { dimension } ); 
  }
  
  void DoArchive(Archive& ar) override
  {
    BASE::DoArchive(ar);
    auto size = ci.Size();
    ar & size;
    ci.SetSize(size);
    for(auto& cf : ci)
      ar.Shallow(cf);
    ar & dimi;
  }

  virtual string GetDescription () const override
  { return "VectorialCoefficientFunction"; }
  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override;

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    for (auto cf : ci)
      cf->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  {
    Array<shared_ptr<CoefficientFunction>> cfa;
    for (auto cf : ci)
      cfa.Append (cf);
    return Array<shared_ptr<CoefficientFunction>>(cfa);
  } 


  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        cf->NonZeroPattern(ud,
                           nonzero.Range(base,base+dimi),
                           nonzero_deriv.Range(base,base+dimi),
                           nonzero_dderiv.Range(base,base+dimi));
        base += dimi;
      }
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        cf->NonZeroPattern(ud, values.Range(base,base+dimi));
        base += dimi;
      }
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    size_t base = 0;
    for (size_t i : Range(ci))
      {
        values.Range(base,base+dimi[i]) = input[i];
        base += dimi[i];
      }    
  }
  
  
  virtual bool DefinedOn (const ElementTransformation & trafo) override
  {
    for (auto & cf : ci)
      if (!cf->DefinedOn(trafo)) return false;
    return true;
  }
  

  using BASE::Evaluate;  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    int base = 0;
    for (auto & cf : ci)
      {
        int dimi = cf->Dimension();
        cf->Evaluate(ip, result.Range(base,base+dimi));
        base += dimi;
      }
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        cf->Evaluate(ip, result.Range(base,base+dimi));
        base += dimi;
      }

    // for (int i : Range(ci))
    // ci[i]->Evaluate(ip, result.Range(i,i+1));
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    size_t base = 0;
    for (auto i : Range(ci.Size()))
      {
        ci[i]->Evaluate(ir, values.Rows(base, base + dimi[i]));
        base += dimi[i];
      }
  }
 
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    size_t base = 0;
    size_t np = ir.Size();
    for (size_t i : Range(ci))
      {
        values.Rows(base,base+dimi[i]).AddSize(dimi[i], np) = input[i];
        base += dimi[i];
      }
  }

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        BareSliceMatrix<Complex> result) const override
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        STACK_ARRAY(double, hmem, 2*ir.Size()*dimi);
        FlatMatrix<Complex> temp(ir.Size(), dimi, (Complex*)hmem);
        cf->Evaluate(ir, temp);
        result.Cols(base,base+dimi).AddSize(ir.Size(), dimi) = temp;
        base += dimi;
      }
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    Array<shared_ptr<CoefficientFunction>> diff_ci;
    for (auto & cf : ci)
      if (cf)
        diff_ci.Append (cf->Diff(var, dir));
      else
        diff_ci.Append (nullptr);
    auto veccf = MakeVectorialCoefficientFunction (move(diff_ci));
    veccf->SetDimensions(Dimensions());
    return veccf;
  }

  shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, typename BASE::T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    if (this == var)
      return IdentityCF(this->Dimensions());

    int dimvar = var->Dimension();
    Array<int> dimres;
    dimres += this->Dimensions();
    dimres += var->Dimensions();

    Array<shared_ptr<CoefficientFunction>> diff_ci;
    for (auto & cf : ci)
      if (cf)
        diff_ci.Append (cf->DiffJacobi (var, cache) -> Reshape(cf->Dimension()*dimvar));
      else
        diff_ci.Append (nullptr);
    auto res = MakeVectorialCoefficientFunction (move(diff_ci)) -> Reshape(dimvar,this->Dimension()) -> Reshape(dimres);

    cache[thisptr] = res;
    return res;

  }

};

  void VectorialCoefficientFunction::GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    int input = 0;
    int input_index = 0;
    // code.Declare (code.res_type, index, Dimensions());
    code.Declare (index, Dimensions(), IsComplex());
    for (int i = 0; i < Dimension(); i++)
      {
        auto cfi = ci[input];
        code.body += Var(index, i, this->Dimensions())
          .Assign(Var(inputs[input], input_index, cfi->Dimensions()), false);
        input_index++;
        if (input_index == cfi->Dimension() )
          {
            input++;
            input_index = 0;
          }
      }

  }


  shared_ptr<CoefficientFunction>
  MakeVectorialCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
  {
    int totdim = 0;
    for (auto cf : aci)
      if (!cf->IsZeroCF())
        return make_shared<VectorialCoefficientFunction> (move (aci));
      else
        totdim += cf->Dimension();
    
    return ZeroCF( Array<int>( {totdim} ));
  }


// ////////////////////////// Coordinate CF ////////////////////////

  class CoordCoefficientFunction
    : public T_CoefficientFunction<CoordCoefficientFunction, CoefficientFunctionNoDerivative>
  {
    int dir;
    typedef T_CoefficientFunction<CoordCoefficientFunction, CoefficientFunctionNoDerivative> BASE;
  public:
    CoordCoefficientFunction() = default;
    CoordCoefficientFunction (int adir) : BASE(1, false), dir(adir) { SetVariable(true); }

    void DoArchive(Archive& ar) override
    {
      BASE::DoArchive(ar);
      ar & dir;
    }

    virtual string GetDescription () const override
    {
      string dirname;
      switch (dir)
        {
        case 0: dirname = "x"; break;
        case 1: dirname = "y"; break;
        case 2: dirname = "z"; break;
        default: dirname = ToLiteral(dir);
        }
      return string("coordinate ")+dirname;
    }

    using BASE::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      if (!ip.IsComplex())
        return ip.GetPoint()(dir);
      else
        return ip.GetPointComplex()(dir).real();
    }
    /*
    virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                          FlatMatrix<> result) const
    {
      if (!ir.IsComplex())
        result.Col(0) = ir.GetPoints().Col(dir);
      else
      {
        auto pnts = ir.GetPointsComplex().Col(dir);
        for (auto i : Range(ir.Size()))
          result(i,0) = pnts(i).real();
      }
    }
    virtual void Evaluate(const BaseMappedIntegrationRule & ir,
			  FlatMatrix<Complex> result) const override
    {
      result.Col(0) = ir.GetPoints().Col(dir);
    }
    */

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
        auto v = Var(index);
        // code.body += v.Assign(CodeExpr(string("mir.GetPoints()(i,")+ToLiteral(dir)+")"));

        // code.Declare(code.res_type, index, Dimensions());
        code.Declare(index, Dimensions(), IsComplex());
        code.body += v.Assign(CodeExpr(string("points(i,")+ToLiteral(dir)+")"), false);
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
    {
      size_t nv = ir.Size();
      __assume (nv > 0);

      if(dir>=ir.DimSpace())
      {
        for (size_t i = 0; i < nv; i++)
          values(0,i) = 0.0;
        return;
      }

      if(!ir.IsComplex())
        {
          auto points = ir.GetPoints();
          for (size_t i = 0; i < nv; i++)
            values(0,i) = points(i, dir);
        }
      else
        {
          auto cpoints = ir.GetPointsComplex();
          __assume (nv > 0);
          for(auto i : Range(nv))
            values(0,i) = cpoints(i,dir).real();
        }
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    { T_Evaluate (ir, values); }

    /*
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    {
      Evaluate (ir, values);
    }
    */

    using CoefficientFunction::Operator;
    shared_ptr<CoefficientFunction> Operator (const string & name) const override
    {
      if (spacedim == -1)
        throw Exception("cannot differentiate coordinate since we don't know the space dimension, use 'coef.spacedim=dim'");
      if (name != "grad")
        throw Exception ("cannot apply operator "+name+" for coordinate");
      
      Array<shared_ptr<CoefficientFunction>> funcs(spacedim);
      funcs = ZeroCF(Array<int>());
      funcs[dir] = make_shared<ConstantCoefficientFunction> (1);
      return MakeVectorialCoefficientFunction (move(funcs));
    }
    
    
    
    shared_ptr<CoefficientFunction>
    Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dirdiff) const override
    {
      // if (var == shape.get())
      if (dynamic_cast<const DiffShapeCF*>(var))
        return MakeComponentCoefficientFunction (dirdiff, dir);
      // return BASE::Diff (var, dirdiff);
      
      if (auto coordcf = dynamic_cast<const CoordCoefficientFunction*>(var))
        if (coordcf->dir == this->dir)
          return dirdiff;
      
      return ZeroCF(Dimensions());
    }
    
    
  };


shared_ptr<CoefficientFunction> MakeCoordinateCoefficientFunction (int comp)
{
  return make_shared<CoordCoefficientFunction> (comp);
}



class NGS_DLL_HEADER FrozenCoefficientFunction
  : public T_CoefficientFunction<FrozenCoefficientFunction, CoefficientFunctionNoDerivative>
{
  shared_ptr<CoefficientFunction> cf;
 public:
  FrozenCoefficientFunction (shared_ptr<CoefficientFunction> acf) : cf(acf)
    {
      SetDimensions (cf->Dimensions());
    }
  
  ~FrozenCoefficientFunction() { ; }
  
  void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    cf->TraverseTree (func);
    func(*this);
  }

  string GetDescription () const override
  {
    return "frozen";
  }    

  
  Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  {
    return { cf };
  }
  

  
  double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    return cf -> Evaluate(ip); 
  }

  
  template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    cf -> Evaluate (ir, values);
  }
  
  template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
  {
    cf -> Evaluate (ir, values);
  }

  shared_ptr<CoefficientFunction>
    Operator (const string & name) const override
  {
    return Freeze (cf->Operator(name));
  }
  
};



shared_ptr<CoefficientFunction> Freeze (shared_ptr<CoefficientFunction> cf)
{
  return make_shared<FrozenCoefficientFunction> (cf);
}




  // ///////////////////////////// Compiled CF /////////////////////////
class CompiledCoefficientFunction : public CompiledCoefficientFunctionInterface //, public std::enable_shared_from_this<CompiledCoefficientFunction>
  {
    typedef void (*lib_function)(const ngfem::BaseMappedIntegrationRule &, ngbla::BareSliceMatrix<double>);
    typedef void (*lib_function_simd)(const ngfem::SIMD_BaseMappedIntegrationRule &, BareSliceMatrix<SIMD<double>>);
    typedef void (*lib_function_deriv)(const ngfem::BaseMappedIntegrationRule &, ngbla::BareSliceMatrix<AutoDiff<1,double>>);
    typedef void (*lib_function_simd_deriv)(const ngfem::SIMD_BaseMappedIntegrationRule &, BareSliceMatrix<AutoDiff<1,SIMD<double>>>);
    typedef void (*lib_function_dderiv)(const ngfem::BaseMappedIntegrationRule &, ngbla::BareSliceMatrix<AutoDiffDiff<1,double>>);
    typedef void (*lib_function_simd_dderiv)(const ngfem::SIMD_BaseMappedIntegrationRule &, BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>>);

    typedef void (*lib_function_complex)(const ngfem::BaseMappedIntegrationRule &, ngbla::BareSliceMatrix<Complex>);
    typedef void (*lib_function_simd_complex)(const ngfem::SIMD_BaseMappedIntegrationRule &, BareSliceMatrix<SIMD<Complex>>);

    shared_ptr<CoefficientFunction> cf;
    DynamicTable<int> inputs;
    size_t max_inputsize;
    Array<int> dim;
    int totdim;
    Array<bool> is_complex;
    // Array<Timer*> timers;
    unique_ptr<SharedLibrary> library;
    lib_function compiled_function = nullptr;
    lib_function_simd compiled_function_simd = nullptr;
    lib_function_deriv compiled_function_deriv = nullptr;
    lib_function_simd_deriv compiled_function_simd_deriv = nullptr;
    lib_function_dderiv compiled_function_dderiv = nullptr;
    lib_function_simd_dderiv compiled_function_simd_dderiv = nullptr;

    lib_function_complex compiled_function_complex = nullptr;
    lib_function_simd_complex compiled_function_simd_complex = nullptr;

    bool _real_compile = false;
    int _maxderiv = 2;
    bool _wait = false;
    bool _keep_files = false;
    
  public:
    CompiledCoefficientFunction() = default;
    CompiledCoefficientFunction (shared_ptr<CoefficientFunction> acf)
      : CompiledCoefficientFunctionInterface(acf->Dimension(), acf->IsComplex()), cf(acf) // , compiled_function(nullptr), compiled_function_simd(nullptr)
    {
      SetDimensions (cf->Dimensions());
      cf -> TraverseTree
        ([&] (CoefficientFunction & stepcf)
         {
           if (!steps.Contains(&stepcf))
             {
               steps.Append (&stepcf);
               // timers.Append (new Timer(string("CompiledCF")+typeid(stepcf).name()));
               dim.Append (stepcf.Dimension());
               is_complex.Append (stepcf.IsComplex());
             }
         });
      
      totdim = 0;
      for (int d : dim) totdim += d;
      
      cout << IM(3) << "Compiled CF:" << endl;
      for (auto cf : steps)
        cout << IM(3) << typeid(*cf).name() << endl;
      
      inputs = DynamicTable<int> (steps.Size());
      max_inputsize = 0;
      
      cf -> TraverseTree
        ([&] (CoefficientFunction & stepcf)
         {
           int mypos = steps.Pos (&stepcf);
           if (!inputs[mypos].Size())
             {
               Array<shared_ptr<CoefficientFunction>> in = stepcf.InputCoefficientFunctions();
               max_inputsize = max2(in.Size(), max_inputsize);
               for (auto incf : in)
                 inputs.Add (mypos, steps.Pos(incf.get()));
             }
         });
      cout << IM(3) << "inputs = " << endl << inputs << endl;

    }


  void PrintReport (ostream & ost) const override
  {
    ost << "Compiled CF:" << endl;
    for (int i : Range(steps))
      {
        auto & cf = steps[i];
        ost << "Step " << i << ": " << cf->GetDescription();
        if (cf->Dimensions().Size() == 1)
          ost << ", dim=" << cf->Dimension();
        else if (cf->Dimensions().Size() >= 2)
          {
            // ost << ", dims = " << cf->Dimensions()[0] << " x " << cf->Dimensions()[1];
            ost << ", dims = " << cf->Dimensions()[0];
            for (int j = 1; j < cf->Dimensions().Size(); j++)
              ost << " x " << cf->Dimensions()[j];
          }
        ost << endl;
        if (inputs[i].Size() > 0)
          {
            ost << "     input: ";
            for (auto innr : inputs[i])
              ost << innr << " ";
            ost << endl;
          }
      }
    /*
    for (auto cf : steps)
      ost << cf -> GetDescription() << endl;
    ost << "inputs = " << endl << inputs << endl;
    */
    
      /*
    for (int i = 0; i < 2*level; i++)
      ost << ' ';
    ost << "coef " << GetDescription() << ","
        << (IsComplex() ? " complex" : " real");
    if (Dimensions().Size() == 1)
      ost << ", dim=" << Dimension();
    else if (Dimensions().Size() == 2)
      ost << ", dims = " << Dimensions()[0] << " x " << Dimensions()[1];
    ost << endl;

    Array<shared_ptr<CoefficientFunction>> input = InputCoefficientFunctions();
    for (int i = 0; i < input.Size(); i++)
      input[i] -> PrintReportRec (ost, level+1);
    */
  }

    
    void DoArchive(Archive& ar) override
    {
      CoefficientFunction::DoArchive(ar);
      ar.Shallow(cf);
      if(ar.Input())
        {
          cf -> TraverseTree
            ([&] (CoefficientFunction & stepcf)
             {
               if (!steps.Contains(&stepcf))
                 {
                   steps.Append (&stepcf);
                   // timers.Append (new Timer(string("CompiledCF")+typeid(stepcf).name()));
                   dim.Append (stepcf.Dimension());
                   is_complex.Append (stepcf.IsComplex());
                 }
             });
          totdim = 0;
          for (int d : dim) totdim += d;
          inputs = DynamicTable<int> (steps.Size());
          max_inputsize = 0;

          cf -> TraverseTree
            ([&] (CoefficientFunction & stepcf)
             {
               int mypos = steps.Pos (&stepcf);
               if (!inputs[mypos].Size())
                 {
                   Array<shared_ptr<CoefficientFunction>> in = stepcf.InputCoefficientFunctions();
                   max_inputsize = max2(in.Size(), max_inputsize);
                   for (auto incf : in)
                     inputs.Add (mypos, steps.Pos(incf.get()));
                 }
             });
        }
    }


    virtual shared_ptr<CoefficientFunction>
    Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override
    {
      if (var == this) return Compile (dir, _real_compile, _maxderiv, _wait, _keep_files);
      auto diff_cf = cf->Diff(var, dir);
      // return Compile (diff_cf, false, 0, 0);
      return Compile (diff_cf, _real_compile, _maxderiv, _wait, _keep_files);
    }

    virtual shared_ptr<CoefficientFunction>
    DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];

      if (var == this) return Compile (IdentityCF(Dimensions()), _real_compile, _maxderiv, _wait);
      auto diff_cf = cf->DiffJacobi (var, cache);
      // return Compile (diff_cf, false, 0, 0);
      auto res = Compile (diff_cf, _real_compile, _maxderiv, _wait, _keep_files);
      cache[thisptr] = res;
      return res;
    }


    Code GenerateProgram (int deriv, bool simd) const override
    {
      Code code;
      code.is_simd = simd;
      code.deriv = deriv;

      string res_type = cf->IsComplex() ? "Complex" : "double";
      if(simd) res_type = "SIMD<" + res_type + ">";
      if(deriv==1) res_type = "AutoDiff<1," + res_type + ">";
      if(deriv==2) res_type = "AutoDiffDiff<1," + res_type + ">";
      code.res_type = res_type;
      
      for (auto i : Range(steps))
        {
          auto& step = *steps[i];
          if (!simd && deriv == 0)              
            cout << IM(3) << "step " << i << ": " << typeid(step).name() << endl;
          code.body += string("// step ")+ToString(i)+": "+step.GetDescription()+"\n";
          step.GenerateCode(code, inputs[i],i);
        }
      
      return code;
    }
    
    void RealCompile(int maxderiv, bool wait, bool keep_files)
    {
      _real_compile = true;
      _maxderiv = maxderiv;
      _wait = wait;
      _keep_files = keep_files;
      
        std::vector<string> link_flags;
        if(cf->IsComplex())
            maxderiv = 0;
        stringstream s;
        string pointer_code;
        string top_code = ""
          "#include<fem.hpp>\n"
          "#if defined(__GNUC__) || defined(__clang__)\n"
          "#pragma GCC diagnostic ignored \"-Wunused-but-set-variable\"\n"
          "#endif\n"
             "using namespace ngfem;\n"
             "extern \"C\" {\n"
             ;

        string parameters[3] = {"results", "deriv", "dderiv"};

        for (int deriv : Range(maxderiv+1))
          for (auto simd : {false,true}) {
            if (!simd && deriv == 0)
              cout << IM(3) << "Compiled CF:" << endl;

            Code code = GenerateProgram(deriv, simd);
            
            /*
            Code code;
            code.is_simd = simd;
            code.deriv = deriv;

            string res_type = cf->IsComplex() ? "Complex" : "double";
            if(simd) res_type = "SIMD<" + res_type + ">";
            if(deriv==1) res_type = "AutoDiff<1," + res_type + ">";
            if(deriv==2) res_type = "AutoDiffDiff<1," + res_type + ">";
            code.res_type = res_type;
            */
            
            string res_type = code.res_type;

            /*
            for (auto i : Range(steps)) {
              auto& step = *steps[i];
              if (!simd && deriv == 0)              
                cout << IM(3) << "step " << i << ": " << typeid(step).name() << endl;
              code.body += string("// step ")+ToString(i)+": "+step.GetDescription()+"\n";
              step.GenerateCode(code, inputs[i],i);
            }
            */
            
            pointer_code += code.pointer;
            top_code += code.top;

            // set results
            string scal_type = cf->IsComplex() ? "Complex" : "double";
            int ii = 0;

            // code.Declare (code.res_type, steps.Size(), cf->Dimensions());
            code.Declare (steps.Size(), cf->Dimensions(), cf->IsComplex());
            for (int ind = 0; ind < cf->Dimension(); ind++) {
              // code.body += Var(steps.Size(),ind,cf->Dimensions()).Declare(res_type);
              code.body += Var(steps.Size(),ind,cf->Dimensions()).Assign(Var(steps.Size()-1,ind, cf->Dimensions()),false);
              string sget = "(i," + ToLiteral(ii) + ") =";
              if(simd) sget = "(" + ToLiteral(ii) + ",i) =";
              
              // for (auto ideriv : Range(simd ? 1 : deriv+1))
              for (auto ideriv : Range(1))
                {
                  code.body += parameters[ideriv] + sget + Var(steps.Size(),ind,cf->Dimensions()).code;
                  code.body += ";\n";
                }
              ii++;
            };
            

            if(code.header.find("gridfunction_local_heap") != std::string::npos)
            {
                code.header.insert(0, "LocalHeapMem<100000> gridfunction_local_heap(\"compiled_cf_gfheap\");\n");
                code.header.insert(0, "ArrayMem<int, 100> gridfunction_dnums;\n");
                code.header.insert(0, "ArrayMem<double, 100> gridfunction_elu;\n");
            }

            // Function name
#ifdef WIN32
            s << "__declspec(dllexport) ";
#endif
            s << "void CompiledEvaluate";
            if(deriv==2) s << "D";
            if(deriv>=1) s << "Deriv";
            if(simd) s << "SIMD";

            // Function parameters
            if (simd)
              s << "(SIMD_BaseMappedIntegrationRule & mir, BareSliceMatrix<" << res_type << "> results";
            else
              s << "(BaseMappedIntegrationRule & mir, BareSliceMatrix<" << res_type << "> results";
            s << " ) {" << endl;
            s << code.header << endl;
            s << "[[maybe_unused]] auto points = mir.GetPoints();" << endl;
            s << "[[maybe_unused]] auto domain_index = mir.GetTransformation().GetElementIndex();" << endl;
            s << "for ( auto i : Range(mir)) {" << endl;
            s << "[[maybe_unused]] auto & ip = mir[i];" << endl;
            s << code.body << endl;
            s << "}\n}" << endl << endl;

            for(const auto &lib : code.link_flags)
                if(std::find(std::begin(link_flags), std::end(link_flags), lib) == std::end(link_flags))
                    link_flags.push_back(lib);

        }
        s << "}" << endl;
        string file_code = top_code + s.str();
        std::vector<std::variant<filesystem::path, string>> codes;
        codes.push_back(file_code);
        if(pointer_code.size()) {
          pointer_code = "extern \"C\" {\n" + pointer_code;
          pointer_code += "}\n";
          codes.push_back(pointer_code);
        }

        auto self = dynamic_pointer_cast<CompiledCoefficientFunction>(shared_from_this());
        auto compile_func = [self, codes, link_flags, maxderiv, keep_files] () {
              self->library = CompileCode( codes, link_flags, keep_files );
              if(self->cf->IsComplex())
              {
                  self->compiled_function_simd_complex = self->library->GetFunction<lib_function_simd_complex>("CompiledEvaluateSIMD");
                  self->compiled_function_complex = self->library->GetFunction<lib_function_complex>("CompiledEvaluate");
              }
              else
              {
                  self->compiled_function_simd = self->library->GetFunction<lib_function_simd>("CompiledEvaluateSIMD");
                  self->compiled_function = self->library->GetFunction<lib_function>("CompiledEvaluate");
                  if(maxderiv>0)
                  {
                      self->compiled_function_simd_deriv = self->library->GetFunction<lib_function_simd_deriv>("CompiledEvaluateDerivSIMD");
                      self->compiled_function_deriv = self->library->GetFunction<lib_function_deriv>("CompiledEvaluateDeriv");
                  }
                  if(maxderiv>1)
                  {
                      self->compiled_function_simd_dderiv = self->library->GetFunction<lib_function_simd_dderiv>("CompiledEvaluateDDerivSIMD");
                      self->compiled_function_dderiv = self->library->GetFunction<lib_function_dderiv>("CompiledEvaluateDDeriv");
                  }
              }
              cout << IM(7) << "Compilation done" << endl;
        };
        if(wait)
            compile_func();
        else
        {
          std::thread(
            [=](){
              try {compile_func();}
              catch (const std::exception &e) {
                cout << "error: compilation to cpp failed" << endl;
                cerr << IM(3) << "Compilation of CoefficientFunction failed: " << e.what() << endl;
              }
            }
          ).detach();
        }
    }

    void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      cf -> TraverseTree (func);
      // func(*cf);
      func(*this);      
    }

    // virtual bool IsComplex() const { return cf->IsComplex(); }
    // virtual int Dimension() const { return cf->Dimension(); }
    // virtual Array<int> Dimensions() const  { return cf->Dimensions(); } 
    
    
    // bool ElementwiseConstant () const override { return false; }
    /*
    virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                                 FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
    { cf->NonZeroPattern (ud, nonzero, nonzero_deriv, nonzero_dderiv); }
    */
    
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      return cf -> Evaluate(ip);
    }

    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<> result) const override
    {
      cf->Evaluate (ip, result);      
    }

    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<Complex> result) const override
    {
      cf->Evaluate (ip, result);
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     BareSliceMatrix<T,ORD> values) const
    {
      ArrayMem<T, 1000> hmem(ir.Size()*totdim);
      size_t mem_ptr = 0;
      ArrayMem<BareSliceMatrix<T,ORD>,100> temp(steps.Size());
      ArrayMem<BareSliceMatrix<T,ORD>, 100> in(max_inputsize);

      // used if nullptr are in tree
      BareSliceMatrix<T,ORD> zeromat(FlatMatrix<T, ORD>(0, ir.Size(), hmem.Data()));
      for (size_t i = 0; i < steps.Size()-1; i++)
        {
          new (&temp[i]) BareSliceMatrix<T,ORD> (FlatMatrix<T,ORD> (dim[i], ir.Size(), &hmem[mem_ptr]));
          mem_ptr += ir.Size()*dim[i];
        }
      
      new (&temp.Last()) BareSliceMatrix<T,ORD>(values);

      for (size_t i = 0; i < steps.Size(); i++)
        {
          auto inputi = inputs[i];
          for (int nr : Range(inputi))
          {
            if(inputi[nr] == -1)
              new (&in[nr]) BareSliceMatrix<T,ORD> (zeromat);
            else
              new (&in[nr]) BareSliceMatrix<T,ORD> (temp[inputi[nr]]);
          }
          steps[i] -> Evaluate (ir, in.Range(0, inputi.Size()), temp[i]);
        }
    }

    /*
    void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                         FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
    {
      typedef AutoDiffDiff<1,bool> T;
      ArrayMem<T, 1000> hmem(totdim);
      size_t mem_ptr = 0;
      ArrayMem<FlatVector<T>,100> temp(steps.Size());
      ArrayMem<FlatVector<T>,100> in(max_inputsize);
      for (size_t i = 0; i < steps.Size(); i++)
        {
          new (&temp[i]) FlatVector<T> (dim[i], &hmem[mem_ptr]);
          mem_ptr += dim[i];
        }
      
      for (size_t i = 0; i < steps.Size(); i++)
        {
          auto inputi = inputs[i];
          for (int nr : Range(inputi))
            new (&in[nr]) FlatVector<T> (temp[inputi[nr]]);
          steps[i] -> NonZeroPattern (ud, in.Range(0, inputi.Size()), temp[i]);
        }
      auto & last = temp.Last();
      for (size_t i = 0; i < nonzero.Size(); i++)
        {
          nonzero(i) = last(i).Value();
          nonzero_deriv(i) = last(i).DValue(0);
          nonzero_dderiv(i) = last(i).DDValue(0);
        }
    }
    */
    void NonZeroPattern (const class ProxyUserData & ud, FlatVector<AutoDiffDiff<1,bool>> nonzero) const override
    {
      typedef AutoDiffDiff<1,bool> T;
      ArrayMem<T, 1000> hmem(totdim);
      size_t mem_ptr = 0;
      ArrayMem<FlatVector<T>,100> temp(steps.Size());
      ArrayMem<FlatVector<T>,100> in(max_inputsize);

      // used if nullptr are in tree
      FlatVector<T> zerovec(0, hmem.Data());

      for (size_t i = 0; i < steps.Size(); i++)
        {
          new (&temp[i]) FlatVector<T> (dim[i], &hmem[mem_ptr]);
          mem_ptr += dim[i];
        }
      
      for (size_t i = 0; i < steps.Size(); i++)
        {
          auto inputi = inputs[i];
          for (int nr : Range(inputi))
          {
            if(inputi[nr] == -1)
              new (&in[nr]) FlatVector<T> (zerovec);
            else
              new (&in[nr]) FlatVector<T> (temp[inputi[nr]]);
          }
          steps[i] -> NonZeroPattern (ud, in.Range(0, inputi.Size()), temp[i]);
        }
      auto & last = temp.Last();
      for (size_t i = 0; i < nonzero.Size(); i++)
        nonzero(i) = last(i);
    }
    
    void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const override
    {
      if(compiled_function)
      {
        compiled_function(ir,values);
        return;
      }

      T_Evaluate (ir, Trans(values));
    }



    void Evaluate (const BaseMappedIntegrationRule & ir,
                   BareSliceMatrix<AutoDiff<1,double>> values) const override
    {
      if(compiled_function_deriv)
        {
          compiled_function_deriv(ir, values);
          return;
        }

      T_Evaluate (ir, Trans(values));
    }



    void Evaluate (const BaseMappedIntegrationRule & ir,
                   BareSliceMatrix<AutoDiffDiff<1,double>> values) const override
    {
      if(compiled_function_dderiv)
      {
        compiled_function_dderiv(ir, values);
        return;
      }

      T_Evaluate (ir, Trans(values));
    }


    
    void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const override
    {
      if(compiled_function_simd_deriv)
        {
          compiled_function_simd_deriv(ir, values);
          return;
        }

      T_Evaluate (ir, values);
    }


    
    void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const override
    {
      if(compiled_function_simd_dderiv)
      {
        compiled_function_simd_dderiv(ir, values);
        return;
      }
      
      T_Evaluate (ir, values);
    }

    
    void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   BareSliceMatrix<SIMD<double>> values) const override
    {
      if(compiled_function_simd)
      {
        compiled_function_simd(ir, values);
        return;
      }

      T_Evaluate (ir, values);
    }

    void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const override
    {
      if(compiled_function_complex)
      {
          compiled_function_complex(ir,values);
          return;
      }
      else
      {
          cf->Evaluate (ir, values);
      }
    }

    void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   BareSliceMatrix<SIMD<Complex>> values) const override
    {
      if(compiled_function_simd_complex)
      {
        compiled_function_simd_complex(ir,values);
        return;
      }
      else
      {
          cf->Evaluate (ir, values);
      }
    }

    void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
    {
      for (int i = 0; i < cf->Dimension(); i++)
        code.body += Var(index,i,cf->Dimensions()).Assign( Var(inputs[0], i, cf->Dimensions()) );
    }

    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    {
      return Array<shared_ptr<CoefficientFunction>>({ cf });
    }
  };

class RealCF : public CoefficientFunctionNoDerivative
  {
    shared_ptr<CoefficientFunction> cf;
    bool cf_is_complex;
  public:
    RealCF() = default;
    RealCF(shared_ptr<CoefficientFunction> _cf)
      : CoefficientFunctionNoDerivative(_cf->Dimension(),false), cf(_cf)
    {
      cf_is_complex = cf->IsComplex();
      SetDimensions(cf->Dimensions());
    }

    void DoArchive(Archive& ar) override
    {
      CoefficientFunctionNoDerivative::DoArchive(ar);
      ar.Shallow(cf) & cf_is_complex;
    }

    virtual string GetDescription() const override
    {
      return "RealCF";
    }

    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      cf->TraverseTree (func);
      func(*this);
    }

    void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
    {
      for (int i = 0; i < Dimension(); i++)
        code.body += Var(index,i,Dimensions()).Assign( string(Var(inputs[0],i,cf->Dimensions())) + ".real()");
    }
    
    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ cf }); }
    
      using CoefficientFunctionNoDerivative::Evaluate;
    virtual double Evaluate(const BaseMappedIntegrationPoint& ip) const override
    {
      if(cf->IsComplex())
        {
          Vec<1,Complex> val;
          cf->Evaluate(ip,val);
          return val(0).real();
        }
      return cf->Evaluate(ip);
    }

    virtual void Evaluate(const BaseMappedIntegrationPoint& ip, FlatVector<> vec) const override
    {
      if(cf->IsComplex())
        {
          VectorMem<10,Complex> complex_vec(vec.Size());
          cf->Evaluate(ip,complex_vec);
          vec = Real(complex_vec);
        }
      else
          cf->Evaluate(ip,vec);
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const override
    {
      if (!cf_is_complex)
        {
          cf->Evaluate(ir, values);
          return;
        }

      STACK_ARRAY(Complex, mem, ir.Size()*Dimension());
      FlatMatrix<Complex> cvalues(ir.Size(), Dimension(), &mem[0]);
      cf->Evaluate (ir, cvalues);
      values.AddSize(ir.Size(), Dimension()) = Real(cvalues);
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           BareSliceMatrix<SIMD<double>> values) const override
    {
      if (!cf_is_complex)
        {
          cf->Evaluate(ir, values);
          return;
        }

      STACK_ARRAY(SIMD<Complex>, mem, ir.Size()*Dimension());
      FlatMatrix<SIMD<Complex>> cvalues(Dimension(), ir.Size(), &mem[0]);
      cf->Evaluate (ir, cvalues);
      values.AddSize(Dimension(), ir.Size()) = Real(cvalues);
    }
  };

  class ImagCF : public CoefficientFunctionNoDerivative
  {
    shared_ptr<CoefficientFunction> cf;
  public:
    ImagCF() = default;
    ImagCF(shared_ptr<CoefficientFunction> _cf) : CoefficientFunctionNoDerivative(_cf->Dimension(),false), cf(_cf)
    {
      SetDimensions(cf->Dimensions());
    }

    void DoArchive(Archive& ar) override
    {
      CoefficientFunctionNoDerivative::DoArchive(ar);
      ar.Shallow(cf);
    }

    virtual string GetDescription() const override
    {
      return "ImagCF";
    }

    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      cf->TraverseTree (func);
      func(*this);
    }
    
    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ cf }); }
    
    using CoefficientFunctionNoDerivative::Evaluate;
    virtual double Evaluate(const BaseMappedIntegrationPoint& ip) const override
    {
      if(!cf->IsComplex())
          throw Exception("real cf has no imag part!");

      VectorMem<10,Complex> val(cf->Dimension());
      cf->Evaluate(ip,val);
      return val(0).imag();
    }

    virtual void Evaluate(const BaseMappedIntegrationPoint& ip, FlatVector<> vec) const override
    {
      if(cf->IsComplex())
        {
          VectorMem<10,Complex> complex_vec(vec.Size());
          cf->Evaluate(ip,complex_vec);
          vec = Imag(complex_vec);
        }
      else
          cf->Evaluate(ip,vec);
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const override
    {
      if (cf->IsComplex())
        {
          STACK_ARRAY(Complex, mem, ir.Size()*Dimension());
          FlatMatrix<Complex> cvalues(ir.Size(), Dimension(), &mem[0]);
          cf->Evaluate (ir, cvalues);
          values.AddSize(ir.Size(),Dimension()) = Imag(cvalues);
        }
      else
        values.AddSize(ir.Size(),Dimension()) = 0.;
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           BareSliceMatrix<SIMD<double>> values) const override
    {
      if (!cf->IsComplex())
          throw Exception("real cf has no imag part!");

      STACK_ARRAY(SIMD<Complex>, mem, ir.Size()*Dimension());
      FlatMatrix<SIMD<Complex>> cvalues(Dimension(), ir.Size(), &mem[0]);
      cf->Evaluate (ir, cvalues);
      values.AddSize(Dimension(), ir.Size()) = Imag(cvalues);
    }
  };

  shared_ptr<CoefficientFunction> Real(shared_ptr<CoefficientFunction> cf)
  {
    return make_shared<RealCF>(cf);
  }
  shared_ptr<CoefficientFunction> Imag(shared_ptr<CoefficientFunction> cf)
  {
    return make_shared<ImagCF>(cf);
  }

  shared_ptr<CompiledCoefficientFunctionInterface> Compile (shared_ptr<CoefficientFunction> c, bool realcompile, int maxderiv, bool wait, bool keep_files)
  {
    auto compiledcf = dynamic_pointer_cast<CompiledCoefficientFunction>(c);
    auto cf = compiledcf ? compiledcf : make_shared<CompiledCoefficientFunction> (c);
    if(realcompile)
      cf->RealCompile(maxderiv, wait, keep_files);
    return cf;
  }

class LoggingCoefficientFunction : public T_CoefficientFunction<LoggingCoefficientFunction>
{
protected:
  shared_ptr<CoefficientFunction> func;
  unique_ptr<ostream> out;
public:
  LoggingCoefficientFunction (shared_ptr<CoefficientFunction> f, string logfile ) :
      T_CoefficientFunction<LoggingCoefficientFunction>(f->Dimension(), f->IsComplex()),
      func(f)
    {
      this->SetDimensions (func->Dimensions());
      this->elementwise_constant = func->ElementwiseConstant();

      if(logfile=="stdout")
          out = make_unique<ostream>(cout.rdbuf());
      else if(logfile=="stderr")
          out = make_unique<ostream>(cerr.rdbuf());
      else
          out = make_unique<ofstream>(logfile);
    }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
    {
      *out << "======== Evaluate(" << Demangle(typeid(ir).name()) << ", " << Demangle(typeid(values).name()) << ")\n";
      *out << ir;
      func->Evaluate(ir, values);
      *out << "result = \n" << values.AddSize(Dimension(), ir.Size()) << '\n';
    }

  template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,
                     BareSliceMatrix<T,ORD> values) const
  {
    *out << "======== Evaluate(" << Demangle(typeid(ir).name()) << ", " << Demangle(typeid(input).name()) << ", " << Demangle(typeid(values).name()) << ")\n";
    *out << ir;
    *out << "input = \n" << input;
    func->Evaluate(ir, input, values);
    *out << "result = \n" << values.AddSize(Dimension(), ir.Size()) << '\n';
  }

  Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ func }); }

  void PrintReport (ostream & ost) const override
    {
      ost << "LoggingCF(";
      func->PrintReport(ost);
      ost << ")";
    }

  string GetDescription() const override
    {
      return "LoggingCF";
    }

  void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,bool>> nonzero) const override
    {
      func->NonZeroPattern(ud, nonzero);
    }

  void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
    {
      func->NonZeroPattern(ud, input, values);
    }

  void TraverseTree (const function<void(CoefficientFunction&)> & func_) override
    {
      func->TraverseTree(func_);
      func_(*this);
    }
};

shared_ptr<CoefficientFunction> LoggingCF(shared_ptr<CoefficientFunction> func, string logfile)
{
  return make_shared<LoggingCoefficientFunction>(func, logfile);
}









class CacheCoefficientFunction : public T_CoefficientFunction<CacheCoefficientFunction>
{
protected:
  shared_ptr<CoefficientFunction> func;
public:
  CacheCoefficientFunction (shared_ptr<CoefficientFunction> f) 
    : T_CoefficientFunction<CacheCoefficientFunction>(f->Dimension(), f->IsComplex()),
    func(f)
  {
    func->TraverseTree([&](CoefficientFunction &nodecf) {
      if (dynamic_cast<ProxyFunction *>(&nodecf))
        throw Exception("CacheCoefficientFunction: func to be cache must not contain proxy functions");
    });

    this->SetDimensions (func->Dimensions());
    this->elementwise_constant = func->ElementwiseConstant();
  }
  shared_ptr<CoefficientFunction> & Func() { return func; }
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    // should have caches ...
    auto & trafo = ir.GetTransformation();

    if (ProxyUserData * ud = static_cast<ProxyUserData*> (trafo.userdata))
      for (auto [cf,data] : ud->caches)
        {
          if (cf == this)
            {
              if constexpr (std::is_same<T,double>::value)
                {
                  // cout << "d" << endl;
                  auto * mat = static_cast<FlatMatrix<double>*> (data);
                  values.AddSize(mat->Width(), mat->Height()) = 1.0*Trans(*mat);
                  return;
                }
              if constexpr (std::is_same<T,Complex>::value)
                {
                  // static Timer t("CacheCF::Eval c");
                  // RegionTracer reg(TaskManager::GetThreadId(), t);
                  
                  // cout << "c" << endl;
                  if (IsComplex())
                    {
                      auto * mat = static_cast<FlatMatrix<Complex>*> (data);
                      values.AddSize(mat->Width(), mat->Height()) = Trans(*mat);
                    }
                  else
                    {
                      auto * mat = static_cast<FlatMatrix<double>*> (data);
                      values.AddSize(mat->Width(), mat->Height()) = Trans(*mat);
                    }
                  return;
                }
              if constexpr (std::is_same<T,SIMD<double>>::value)
                {
                  // cout << "simd<d>" << endl;
                  auto * mat = static_cast<FlatMatrix<SIMD<double>>*> (data);
                  values.AddSize(mat->Height(), mat->Width()) = *mat;
                  return;
                }
              if constexpr (std::is_same<T,SIMD<Complex>>::value)
                {
                  if (!IsComplex())
                    cout << "simd complex eval of real" << endl;
                  // cout << "simd<c>" << endl;
                  auto * mat = static_cast<FlatMatrix<SIMD<Complex>>*> (data);
                  values.AddSize(mat->Height(), mat->Width()) = *mat;
                  return;
                }
            }
        }
    func->Evaluate(ir, values);
  }
  
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,
                   BareSliceMatrix<T,ORD> values) const
  {
    func->Evaluate(ir, input, values);
  }

  Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ func }); }
  
  void PrintReport (ostream & ost) const override
  {
    ost << "CacheCF(";
    func->PrintReport(ost);
    ost << ")";
  }

  string GetDescription() const override
  {
    return "CacheCF";
  }
  
  void NonZeroPattern (const class ProxyUserData & ud,
                       FlatVector<AutoDiffDiff<1,bool>> nonzero) const override
  {
    func->NonZeroPattern(ud, nonzero);
  }
  
  void NonZeroPattern (const class ProxyUserData & ud,
                       FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                       FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    func->NonZeroPattern(ud, input, values);
  }
  
  void TraverseTree (const function<void(CoefficientFunction&)> & func_) override
  {
    func->TraverseTree(func_);
    func_(*this);
  }
};

shared_ptr<CoefficientFunction> CacheCF(shared_ptr<CoefficientFunction> func)
{
  return make_shared<CacheCoefficientFunction>(func);
}


Array<CoefficientFunction*> FindCacheCF (CoefficientFunction & func)
{
  Array<CoefficientFunction*> cachecf;
  func.TraverseTree
    ( [&] (CoefficientFunction & nodecf)
      {
        if (dynamic_cast<CacheCoefficientFunction*> (&nodecf))
          {
            if (cachecf.Contains(&nodecf)) return;
            cachecf.Append (&nodecf);
          }
      });
  return cachecf;
}


void PrecomputeCacheCF (const Array<CoefficientFunction*> & cachecfs, BaseMappedIntegrationRule & mir,
                        LocalHeap & lh)
{
  // static Timer t("Precompute CacheCF");
  // RegionTracer reg(TaskManager::GetThreadId(), t);
  
  auto & trafo = mir.GetTransformation();
  ProxyUserData * ud = static_cast<ProxyUserData*> (trafo.userdata);
  if (!ud) throw Exception ("don't have a userdata");
  
  new (&ud->caches) FlatArray<pair<const CoefficientFunction*, void*>> (cachecfs.Size(), lh);

  for (int i = 0; i < cachecfs.Size(); i++)
    {
      ud->caches[i].first = cachecfs[i];
      // or complex ...
      if (cachecfs[i]->IsComplex())
        {
          auto * mat = new (lh) FlatMatrix<Complex>(mir.Size(), cachecfs[i]->Dimension(), lh);
          static_cast<CacheCoefficientFunction*>(cachecfs[i])->Func() -> Evaluate (mir, *mat);
          ud->caches[i].second = mat;
        }
      else
        {
          auto * mat = new (lh) FlatMatrix<double>(mir.Size(), cachecfs[i]->Dimension(), lh);
          static_cast<CacheCoefficientFunction*>(cachecfs[i])->Func() -> Evaluate (mir, *mat);
          ud->caches[i].second = mat;
        }
        
    }

}



void PrecomputeCacheCF (const Array<CoefficientFunction*> & cachecfs, SIMD_BaseMappedIntegrationRule & mir,
                        LocalHeap & lh)
{
  auto & trafo = mir.GetTransformation();
  ProxyUserData * ud = static_cast<ProxyUserData*> (trafo.userdata);
  if (!ud) throw Exception ("don't have a userdata");
  
  new (&ud->caches) FlatArray<pair<const CoefficientFunction*, void*>> (cachecfs.Size(), lh);

  for (int i = 0; i < cachecfs.Size(); i++)
    {
      ud->caches[i].first = cachecfs[i];
      // or complex ...
      if (cachecfs[i]->IsComplex())
        {
          auto * mat = new (lh) FlatMatrix<SIMD<Complex>>(cachecfs[i]->Dimension(), mir.Size(), lh);
          static_cast<CacheCoefficientFunction*>(cachecfs[i])->Func() -> Evaluate (mir, *mat);
          ud->caches[i].second = mat;
        }
      else
        {
          auto * mat = new (lh) FlatMatrix<SIMD<double>>(cachecfs[i]->Dimension(), mir.Size(), lh);
          static_cast<CacheCoefficientFunction*>(cachecfs[i])->Func() -> Evaluate (mir, *mat);
          ud->caches[i].second = mat;
        }
        
    }

}



void PrecomputeCacheCF (CoefficientFunction & func, SIMD_BaseMappedIntegrationRule & mir,
                        LocalHeap & lh)
{
  // cout << "precompute cachecf" << endl;
  // first we cnt number of Caches:
  ArrayMem<CacheCoefficientFunction*,10> cachecf;
  func.TraverseTree
    ( [&] (CoefficientFunction & nodecf)
      {
        if (auto ccf = dynamic_cast<CacheCoefficientFunction*> (&nodecf))
          {
            if (cachecf.Contains(ccf)) return;
            cachecf.Append (ccf);
          }
      });

  auto & trafo = mir.GetTransformation();
  ProxyUserData & ud = *static_cast<ProxyUserData*> (trafo.userdata);

  new (&ud.caches) FlatArray<pair<const CoefficientFunction*, void*>> (cachecf.Size(), lh);
  // new FlatArray<pair<const CoefficientFunction*, void*>> (cachecf.Size(), lh);
  for (int i = 0; i < cachecf.Size(); i++)
    {
      ud.caches[i].first = cachecf[i];
      // or complex ...
      if (cachecf[i]->IsComplex())
        {
          auto * mat = new (lh) FlatMatrix<SIMD<Complex>>(cachecf[i]->Dimension(), mir.Size(), lh);
          cachecf[i]->Func() -> Evaluate (mir, *mat);
          ud.caches[i].second = mat;
        }
      else
        {
          auto * mat = new (lh) FlatMatrix<SIMD<double>>(cachecf[i]->Dimension(), mir.Size(), lh);
          cachecf[i]->Func() -> Evaluate (mir, *mat);
          ud.caches[i].second = mat;
        }
        
    }
}







static RegisterClassForArchive<CoefficientFunction> regcf;
static RegisterClassForArchive<ConstantCoefficientFunction, CoefficientFunction> regccf;
static RegisterClassForArchive<ConstantCoefficientFunctionC, CoefficientFunction> regccfc;
static RegisterClassForArchive<ParameterCoefficientFunction<double>, CoefficientFunction> regpard;
static RegisterClassForArchive<ParameterCoefficientFunction<Complex>, CoefficientFunction> regparc;
static RegisterClassForArchive<PlaceholderCoefficientFunction, CoefficientFunction> regplacholdercf;
static RegisterClassForArchive<DomainConstantCoefficientFunction, CoefficientFunction> regdccf;
static RegisterClassForArchive<DomainVariableCoefficientFunction, CoefficientFunction> regdvcf;
static RegisterClassForArchive<IntegrationPointCoefficientFunction, CoefficientFunction> regipcf;
static RegisterClassForArchive<PolynomialCoefficientFunction, CoefficientFunction> regpolcf;
static RegisterClassForArchive<FileCoefficientFunction, CoefficientFunction> regfilecf;
static RegisterClassForArchive<CoordCoefficientFunction, CoefficientFunction> regcoocf;
static RegisterClassForArchive<DomainWiseCoefficientFunction, CoefficientFunction> regdwcf;
static RegisterClassForArchive<VectorialCoefficientFunction, CoefficientFunction> regveccf;
static RegisterClassForArchive<ComponentCoefficientFunction, CoefficientFunction> regcompcf;
static RegisterClassForArchive<ScaleCoefficientFunction, CoefficientFunction> regscale;
static RegisterClassForArchive<ScaleCoefficientFunctionC, CoefficientFunction> regscalec;
static RegisterClassForArchive<MultScalVecCoefficientFunction, CoefficientFunction> regscalvec;
static RegisterClassForArchive<MultVecVecCoefficientFunction, CoefficientFunction> regmultvecvec;
static RegisterClassForArchive<T_MultVecVecCoefficientFunction<1>, CoefficientFunction> regtmultvecvec1;
static RegisterClassForArchive<T_MultVecVecCoefficientFunction<2>, CoefficientFunction> regtmultvecvec2;
static RegisterClassForArchive<T_MultVecVecCoefficientFunction<3>, CoefficientFunction> regtmultvecvec3;
static RegisterClassForArchive<EigCoefficientFunction, CoefficientFunction> regeigcf;
static RegisterClassForArchive<NormCoefficientFunction, CoefficientFunction> regnormcf;
static RegisterClassForArchive<NormCoefficientFunctionC, CoefficientFunction> regnormcfc;
static RegisterClassForArchive<MultMatMatCoefficientFunction, CoefficientFunction> regmultmatmatcf;
static RegisterClassForArchive<MultMatVecCoefficientFunction, CoefficientFunction> regmultmatveccf;
static RegisterClassForArchive<ZeroCoefficientFunction, CoefficientFunction> regzerocf;
static RegisterClassForArchive<UnitVectorCoefficientFunction, CoefficientFunction> regunitcf;
static RegisterClassForArchive<IdentityCoefficientFunction, CoefficientFunction> regidentitycf;
static RegisterClassForArchive<TransposeCoefficientFunction, CoefficientFunction> regtransposecf;
static RegisterClassForArchive<SymmetricCoefficientFunction, CoefficientFunction> regsymmetriccf;
static RegisterClassForArchive<SkewCoefficientFunction, CoefficientFunction> regskewcf;
static RegisterClassForArchive<TraceCoefficientFunction, CoefficientFunction> regtracecf;
static RegisterClassForArchive<InverseCoefficientFunction<1>, CoefficientFunction> reginversecf1;
static RegisterClassForArchive<InverseCoefficientFunction<2>, CoefficientFunction> reginversecf2;
static RegisterClassForArchive<InverseCoefficientFunction<3>, CoefficientFunction> reginversecf3;
static RegisterClassForArchive<InverseCoefficientFunctionAnyDim, CoefficientFunction> reginverseanydimcf;
static RegisterClassForArchive<DeterminantCoefficientFunction<1>, CoefficientFunction> regdetcf1;
static RegisterClassForArchive<DeterminantCoefficientFunction<2>, CoefficientFunction> regdetcf2;
static RegisterClassForArchive<DeterminantCoefficientFunction<3>, CoefficientFunction> regdetcf3;
static RegisterClassForArchive<CofactorCoefficientFunction<1>, CoefficientFunction> regcof1;
static RegisterClassForArchive<CofactorCoefficientFunction<2>, CoefficientFunction> regcof2;
static RegisterClassForArchive<CofactorCoefficientFunction<3>, CoefficientFunction> regcof3;
static RegisterClassForArchive<cl_BinaryOpCF<GenericPlus>, CoefficientFunction> regcfplus;
static RegisterClassForArchive<cl_BinaryOpCF<GenericMinus>, CoefficientFunction> regcfminus;
static RegisterClassForArchive<cl_BinaryOpCF<GenericMult>, CoefficientFunction> regcfmult;
static RegisterClassForArchive<cl_BinaryOpCF<GenericDiv>, CoefficientFunction> regcfdiv;
static RegisterClassForArchive<cl_UnaryOpCF<GenericIdentity>, CoefficientFunction> regcfid;
static RegisterClassForArchive<IfPosCoefficientFunction, CoefficientFunction> regfifpos;
static RegisterClassForArchive<RealCF, CoefficientFunction> regrealcf;
static RegisterClassForArchive<ImagCF, CoefficientFunction> regimagcf;
static RegisterClassForArchive<CompiledCoefficientFunction, CoefficientFunction> regcompiledcf;
static RegisterClassForArchive<OtherCoefficientFunction, CoefficientFunction> regothercf;
static RegisterClassForArchive<LoggingCoefficientFunction, CoefficientFunction> regloggingcf;
static RegisterClassForArchive<CacheCoefficientFunction, CoefficientFunction> regcachingcf;
static RegisterClassForArchive<SubTensorCoefficientFunction, CoefficientFunction> regsubtensorcf;
}

