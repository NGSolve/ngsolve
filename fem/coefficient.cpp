/*********************************************************************/
/* File:   coefficient.cpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   24. Jan. 2002                                             */
/*********************************************************************/

/* 
   Finite Element Coefficient Function
*/

#include <fem.hpp>
#include <../ngstd/evalfunc.hpp>
#include <algorithm>



namespace ngfem
{
  
  CoefficientFunction :: ~CoefficientFunction ()
  { ; }


  void CoefficientFunction :: GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    string mycode =
      string("// GenerateCode() not overloaded for: ") + Demangle(typeid(*this).name()) + "\n"
      + R"CODE_(    STACK_ARRAY({stack_type}, {hmem}, {stack_size});
    {values_type} {values}({rows}, {cols}, reinterpret_cast<{vscal_type}*>(&{hmem}[0]));
    {
      const CoefficientFunction & cf = *reinterpret_cast<CoefficientFunction*>({this});
      {values} = 0.0;
      cf.Evaluate(mir, {values});
    }
    )CODE_";
    auto values = Var("values", index);
    string scal_type = IsComplex() ? "Complex" : "double";
    string vscal_type = code.is_simd ? "SIMD<"+scal_type+">" : scal_type;
    string rows = ToString(Dimension());
    string cols = "mir.IR().Size()";

    std::map<string,string> variables;
    variables["scal_type"] = scal_type;
    variables["vscal_type"] = vscal_type;
    variables["values_type"] = "FlatMatrix<"+vscal_type+">";
    variables["values"] = values.S();
    variables["this"] =  code.AddPointer(this);
    variables["stack_type"] = code.is_simd ? "SIMD<double>" : "double";
    variables["stack_size"] = "mir.Size()*sizeof("+scal_type+")/sizeof(double)*"+ToString(Dimension());
    variables["rows"] = code.is_simd ? rows : cols;
    variables["cols"] = code.is_simd ? cols : rows;
    variables["hmem"] = Var("hmem", index).S();
    code.header += Code::Map(mycode, variables);
    if(code.is_simd)
      {
        TraverseDimensions( Dimensions(), [&](int ind, int i, int j) {
            code.body += Var(index,i,j).Assign(values.S()+"("+ToString(ind)+",i)");
          });
      }
    else
      {
        TraverseDimensions( Dimensions(), [&](int ind, int i, int j) {
            code.body += Var(index,i,j).Assign(values.S()+"(i,"+ToString(ind)+")");
          });
      }
  }

  void CoefficientFunction :: PrintReport (ostream & ost) const
  {
    // ost << "CoefficientFunction is " << typeid(*this).name() << endl;
    PrintReportRec (ost, 0);
  }
  
  void CoefficientFunction :: PrintReportRec (ostream & ost, int level) const
  {
    for (int i = 0; i < 2*level; i++)
      ost << ' ';
    ost << "coef " << GetDescription() << ","
        << (IsComplex() ? " complex" : " real");
    if (Dimensions().Size() == 1)
      ost << ", dim=" << Dimension();
    else if (Dimensions().Size() == 2)
      ost << ", dims = " << Dimensions()[0] << " x " << Dimensions()[1];
    ost << endl;

    Array<CoefficientFunction*> input = InputCoefficientFunctions();
    for (int i = 0; i < input.Size(); i++)
      input[i] -> PrintReportRec (ost, level+1);
  }
  
  string CoefficientFunction :: GetDescription () const
  {
    return typeid(*this).name();
  }    

  
  void CoefficientFunction :: TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    func(*this);
  }
  
  void CoefficientFunction :: 
  Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
    // cout << "switching from rule to point, cf = " << typeid(*this).name() << endl;
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
    throw ExceptionNOSIMD (string("CF :: simd-Evaluate (complex) not implemented for class ") + typeid(*this).name());
  }

  
  void CoefficientFunction :: 
  Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const
  {
    for (int i = 0; i < ir.Size(); i++)
      Evaluate (ir[i], values.Row(i)); 
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

  void CoefficientFunction :: 
  NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                  FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
  {
    // cout << "CoefficientFunction::NonZeroPattern called, type = " << typeid(*this).name() << endl;
    nonzero = true;
    nonzero_deriv = false;
    nonzero_dderiv = false;
  }
  

  ///
  ConstantCoefficientFunction ::   
  ConstantCoefficientFunction (double aval) 
    : BASE(1, false), val(aval) 
  { ; }

  ConstantCoefficientFunction ::
  ~ConstantCoefficientFunction ()
  { ; }

  void ConstantCoefficientFunction :: PrintReport (ostream & ost) const
  {
    ost << "ConstantCF, val = " << val << endl;
  }

  void ConstantCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir,
                                                FlatMatrix<double> values) const
  {
    values = val;
  }

  void ConstantCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir,
                                                FlatMatrix<Complex> values) const
  {
    values = val;
  }

  template <typename T>
  void ConstantCoefficientFunction ::
  T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<T> values) const
  {
    size_t nv = ir.Size();    
    __assume (nv > 0);
    for (size_t i = 0; i < nv; i++)
      values(0,i) = val;
  }
  
  void ConstantCoefficientFunction :: GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    string type = "double";
    if(code.is_simd) type = "SIMD<double>";
    if(code.deriv==1) type = "AutoDiff<1,"+type+">";
    if(code.deriv==2) type = "AutoDiffDiff<1,"+type+">";
    code.body += Var(index).Declare(type);
    code.body += Var(index).Assign(Var(val), false);
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
  
  void ConstantCoefficientFunctionC :: Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const
  {
    values = val;
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


  ///
  ParameterCoefficientFunction ::   
  ParameterCoefficientFunction (double aval) 
    : CoefficientFunctionNoDerivative(1, false), val(aval)
  { ; }

  ParameterCoefficientFunction ::
  ~ParameterCoefficientFunction ()
  { ; }

  void ParameterCoefficientFunction :: PrintReport (ostream & ost) const
  {
    ost << "ParameterCF, val = " << val << endl;
  }

  void ParameterCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir,
                                                FlatMatrix<double> values) const
  {
    values = val;
  }

  void ParameterCoefficientFunction :: GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    string type = "double";
    if(code.is_simd) type = "SIMD<double>";
    if(code.deriv==1) type = "AutoDiff<1,"+type+">";
    if(code.deriv==2) type = "AutoDiffDiff<1,"+type+">";
    stringstream s;
    s << "*reinterpret_cast<double*>(" << code.AddPointer(&val) << ")";
    code.body += Var(index).Declare(type);
    code.body += Var(index).Assign(s.str(), false);
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

  void DomainConstantCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
    int elind = ir[0].GetTransformation().GetElementIndex();
    CheckRange (elind);    
    values = val[elind];
  }

  /*
  void DomainConstantCoefficientFunction :: Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
  {
    int elind = ir[0].GetTransformation().GetElementIndex();
    CheckRange (elind);        
    values.AddSize(Dimension(), ir.Size()) = val[elind];
  }
  */
  template <typename T>
  void DomainConstantCoefficientFunction ::
  T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<T> values) const
  {
    int elind = ir[0].GetTransformation().GetElementIndex();
    CheckRange (elind);        
    // values.AddSize(Dimension(), ir.Size()) = val[elind];

    size_t nv = ir.Size();    
    __assume (nv > 0);
    for (size_t i = 0; i < nv; i++)
      values(0,i) = val[elind];
  }
  

  void DomainConstantCoefficientFunction :: Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const
  {
    int elind = ir[0].GetTransformation().GetElementIndex();
    CheckRange (elind);            
    values = val[elind]; 
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
        args.Range(0,ip.Dim()) = ip.GetPoint();
	
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
        args.Range(0,ip.Dim()) = ip.GetPoint();
	
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
    args.Range(0,ip.Dim()) = ip.GetPoint();
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
	  FlatMatrix<double> values) const
{
  if (ir.Size() == 0) return;
  int elind = ir.GetTransformation().GetElementIndex();
  if (fun.Size() == 1) elind = 0;

  if (! fun[elind] -> IsComplex ())
    {
      ArrayMem<double,2000> mem(ir.Size()*numarg);
      FlatMatrix<> args(ir.Size(), numarg, &mem[0]);
      
      int dim = ir[0].Dim();
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
	fun[elind]->Eval (&args(i,0), &values(i,0), values.Width());
    }
  else
    {
      Matrix<Complex> args(ir.Size(), numarg);
      for (int i = 0; i < ir.Size(); i++)
	args.Row(i).Range(0,ir[i].Dim()) = ir[i].GetPoint();
      
      for (int i = 0, an = 3; i < depends_on.Size(); i++)
	{
	  int dim = depends_on[i]->Dimension();
	  Matrix<Complex> hmat(ir.Size(), dim);
	  depends_on[i] -> Evaluate (ir, hmat);
	  args.Cols(an,an+dim) = hmat;
	  an += dim;
	}
    
      for (int i = 0; i < ir.Size(); i++)
	fun[elind]->Eval (&args(i,0), &values(i,0), values.Width());
    }
}

  
void DomainVariableCoefficientFunction :: PrintReport (ostream & ost) const
{
  *testout << "DomainVariableCoefficientFunction, functios are: " << endl;
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




  
class ScaleCoefficientFunction : public T_CoefficientFunction<ScaleCoefficientFunction>
{
  double scal;
  shared_ptr<CoefficientFunction> c1;
  typedef T_CoefficientFunction<ScaleCoefficientFunction> BASE;
public:
  ScaleCoefficientFunction (double ascal, 
                            shared_ptr<CoefficientFunction> ac1)
    : BASE(ac1->Dimension(), ac1->IsComplex()),
      scal(ascal), c1(ac1)
  {
    SetDimensions(c1->Dimensions());
  }
  
  virtual void PrintReport (ostream & ost) const
  {
    ost << scal << "*(";
    c1->PrintReport(ost);
    ost << ")";
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    TraverseDimensions( c1->Dimensions(), [&](int ind, int i, int j) {
        code.body += Var(index,i,j).Assign(Var(scal) * Var(inputs[0],i,j));
    });
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get() }); }

  virtual bool DefinedOn (const ElementTransformation & trafo)
  { return c1->DefinedOn(trafo); }
    
  using BASE::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->Evaluate(ip);
  }
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->EvaluateComplex(ip);
  }
  virtual double EvaluateConst () const
  {
    return scal * c1->EvaluateConst();
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
  
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         FlatMatrix<double> values) const
  {
    c1->Evaluate (ir, values);
    values *= scal;
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   BareSliceMatrix<T> values) const
  {
    c1->Evaluate (ir, values);
    values.AddSize(Dimension(), ir.Size()) *= scal;
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   FlatArray<BareSliceMatrix<T>> input,                       
                   BareSliceMatrix<T> values) const
  {
    auto in0 = input[0];
    values.AddSize(Dimension(), ir.Size()) = scal * in0;
  }
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    auto in0 = *input[0];
    values = scal * in0;
  }
  
  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         FlatMatrix<Complex> values) const
  {
    c1->Evaluate (ir, values);
    values *= scal;
  }
  
  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                              FlatMatrix<> result, FlatMatrix<> deriv) const
  {
    c1->EvaluateDeriv (ir, result, deriv);
    result *= scal;
    deriv *= scal;
  }
  
  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                               FlatMatrix<> result, FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const
  {
    c1->EvaluateDDeriv (ir, result, deriv, dderiv);
    result *= scal;
    deriv *= scal;
    dderiv *= scal;
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & mir,
                         FlatArray<FlatMatrix<>*> input,
                         FlatMatrix<> result) const
  {
    FlatMatrix<> v1 = *input[0];
    result = scal * v1;
  }

  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const
  {
    FlatMatrix<> v1 = *input[0];
    FlatMatrix<> dv1 = *dinput[0];

    result = scal * v1;
    deriv = scal * dv1;
  }

  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                               FlatArray<FlatMatrix<>*> input,
                               FlatArray<FlatMatrix<>*> dinput,
                               FlatArray<FlatMatrix<>*> ddinput,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const
  {
    FlatMatrix<> v1 = *input[0];
    FlatMatrix<> dv1 = *dinput[0];
    FlatMatrix<> ddv1 = *ddinput[0];

    result = scal * v1;
    deriv = scal * dv1;
    dderiv = scal * ddv1;
  }

  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & mir, 
                              AFlatMatrix<double> values, AFlatMatrix<double> deriv) const
  {
    c1 -> EvaluateDeriv (mir, values, deriv);
    values *= scal;
    deriv *= scal;
  }  
  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & mir, 
                               AFlatMatrix<double> values, AFlatMatrix<double> deriv,
                               AFlatMatrix<double> dderiv) const
  {
    c1 -> EvaluateDDeriv (mir, values, deriv, dderiv);
    values *= scal;
    deriv *= scal;
    dderiv *= scal;
  }
  
  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                              FlatArray<AFlatMatrix<>*> input,
                              FlatArray<AFlatMatrix<>*> dinput,
                              AFlatMatrix<> result,
                              AFlatMatrix<> deriv) const
  {
    result = scal * (*input[0]);
    deriv = scal * (*dinput[0]);
  }
  
  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                               FlatArray<AFlatMatrix<>*> input,
                               FlatArray<AFlatMatrix<>*> dinput,
                               FlatArray<AFlatMatrix<>*> ddinput,
                               AFlatMatrix<> result,
                               AFlatMatrix<> deriv,
                               AFlatMatrix<> dderiv) const
  {
    result = scal * (*input[0]);
    deriv = scal * (*dinput[0]);
    dderiv = scal * (*ddinput[0]);
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
  {
    c1->NonZeroPattern (ud, nonzero, nonzero_deriv, nonzero_dderiv);
  }  
};


class ScaleCoefficientFunctionC : public CoefficientFunction
{
  Complex scal;
  shared_ptr<CoefficientFunction> c1;
public:
  ScaleCoefficientFunctionC (Complex ascal, 
                            shared_ptr<CoefficientFunction> ac1)
    : CoefficientFunction(ac1->Dimension(), true), scal(ascal), c1(ac1)
  {
    SetDimensions (c1->Dimensions());
  }
  
  // virtual bool IsComplex() const { return true; }
  // virtual int Dimension() const { return c1->Dimension(); }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get() }); }
  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    TraverseDimensions( c1->Dimensions(), [&](int ind, int i, int j) {
        code.body += Var(index,i,j).Assign(Var(scal) * Var(inputs[0],i,j));
    });
  }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    throw Exception ("real Evaluate called for complex ScaleCF");
  }
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->EvaluateComplex(ip);    
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<Complex> result) const
  {
    c1->Evaluate (ir, result);
    result *= scal;
  }

  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                         BareSliceMatrix<SIMD<Complex>> values) const
  {
    c1->Evaluate (ir, values);
    values.AddSize(Dimension(), ir.Size()) *= scal;
  }
  
  
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
  {
    c1->NonZeroPattern (ud, nonzero, nonzero_deriv, nonzero_dderiv);
  }  
    
};

// ***********************************************************************************

class MultScalVecCoefficientFunction : public T_CoefficientFunction<MultScalVecCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;  // scalar
  shared_ptr<CoefficientFunction> c2;  // vector
  typedef T_CoefficientFunction<MultScalVecCoefficientFunction> BASE;
public:
  MultScalVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                  shared_ptr<CoefficientFunction> ac2)
    : BASE(ac2->Dimension(), ac1->IsComplex() || ac2->IsComplex()),
      c1(ac1), c2(ac2)
  {
    SetDimensions (c2->Dimensions());
  }
  
  // virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  // virtual int Dimension() const { return c2->Dimension(); }
  // virtual Array<int> Dimensions() const { return c2->Dimensions(); }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get(), c2.get() }); }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    TraverseDimensions( c2->Dimensions(), [&](int ind, int i, int j) {
      code.body += Var(index,i,j).Assign( Var(inputs[0]) * Var(inputs[1],i,j) );
    });
  }

  using BASE::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    throw Exception ("double MultScalVecCF::Evaluate called");
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    Vec<1> v1;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, result);
    result *= v1(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    Vec<1,Complex> v1;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, result);
    result *= v1(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
    STACK_ARRAY(double, hmem1, ir.Size());
    FlatMatrix<> temp1(ir.Size(), 1, hmem1);
    
    c1->Evaluate(ir, temp1);
    c2->Evaluate(ir, result);
    for (int i = 0; i < ir.Size(); i++)
      result.Row(i) *= temp1(i,0);
  }

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<Complex> result) const
  {
    STACK_ARRAY(double, hmem1, 2*ir.Size());
    FlatMatrix<Complex> temp1(ir.Size(), 1, reinterpret_cast<Complex*> (&hmem1[0]));
    
    c1->Evaluate(ir, temp1);
    c2->Evaluate(ir, result);
    for (int i = 0; i < ir.Size(); i++)
      result.Row(i) *= temp1(i,0);
  }

  /*
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
  {
    STACK_ARRAY(SIMD<double>, hmem1, values.Width());
    AFlatMatrix<double> temp1(1, values.Width(), &hmem1[0]);
    
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, values);

    for (size_t j = 0; j < values.Height(); j++)
      for (size_t i = 0; i < values.VWidth(); i++)
        values.Get(j,i) *= temp1.Get(0,i);
  }
  */
  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<T> values) const
  {
    size_t w = ir.Size();
    __assume (w > 0);
    STACK_ARRAY(T, hmem1, w);
    FlatMatrix<T> temp1(1, w, &hmem1[0]);
    
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, values);

    for (size_t j = 0; j < Dimension(); j++)
      for (size_t i = 0; i < w; i++)
        values(j,i) *= temp1(0,i);
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   FlatArray<BareSliceMatrix<T>> input,                       
                   BareSliceMatrix<T> values) const
  {
    auto in0 = input[0];
    auto in1 = input[1];
    size_t dim = Dimension();
    size_t np = ir.Size();
    
    for (size_t j = 0; j < dim; j++)
      for (size_t i = 0; i < np; i++)
        values(j,i) = in0(0,i) * in1(j,i);
  }
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    auto in0 = *input[0];
    auto in1 = *input[1];

    for (size_t j = 0; j < values.Height(); j++)
      for (size_t i = 0; i < values.VWidth(); i++)
        values.Get(j,i) = in0.Get(0,i) * in1.Get(j,i);
  }



  
  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                              FlatMatrix<> result, FlatMatrix<> deriv) const
  {
    STACK_ARRAY(double, hmem1, ir.Size());
    FlatMatrix<> temp1(ir.Size(), 1, hmem1);
    STACK_ARRAY(double, hmem2, ir.Size());
    FlatMatrix<> deriv1(ir.Size(), 1, hmem2);
    c1->EvaluateDeriv(ir, temp1, deriv1);
    c2->EvaluateDeriv(ir, result, deriv);
    for (int i = 0; i < ir.Size(); i++)
      {
        deriv.Row(i) *= temp1(i,0);
        deriv.Row(i) += deriv1(i,0) * result.Row(i);
        result.Row(i) *= temp1(i,0);
      }
  }


  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                               FlatMatrix<> result, FlatMatrix<> deriv, FlatMatrix<> dderiv) const
  {
    STACK_ARRAY(double, hmem1, ir.Size());
    FlatMatrix<> temp1(ir.Size(), 1, hmem1);
    STACK_ARRAY(double, hmem2, ir.Size());
    FlatMatrix<> deriv1(ir.Size(), 1, hmem2);
    STACK_ARRAY(double, hmem3, ir.Size());
    FlatMatrix<> dderiv1(ir.Size(), 1, hmem3);

    c1->EvaluateDDeriv(ir, temp1, deriv1, dderiv1);
    c2->EvaluateDDeriv(ir, result, deriv, dderiv);
    for (int i = 0; i < ir.Size(); i++)
      {
        dderiv.Row(i) *= temp1(i,0);
        dderiv.Row(i) += 2*deriv1(i,0) * deriv.Row(i);
        dderiv.Row(i) += dderiv1(i,0) * result.Row(i);
        deriv.Row(i) *= temp1(i,0);
        deriv.Row(i) += deriv1(i,0) * result.Row(i);
        result.Row(i) *= temp1(i,0);
      }
  }


  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatArray<FlatMatrix<>*> input,
                         FlatMatrix<double> result) const
  {
    FlatMatrix<> temp1 = *input[0];
    FlatMatrix<> temp2 = *input[1];
    for (int i = 0; i < ir.Size(); i++)
      result.Row(i) = temp1(i,0) * temp2.Row(i);
  }

  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const
  {
    FlatMatrix<> v1 = *input[0], v2 = *input[1];
    FlatMatrix<> dv1 = *dinput[0], dv2 = *dinput[1];
    
    for (int k = 0; k < mir.Size(); k++)
      {
        result.Row(k) = v1(k,0)*v2.Row(k);
        deriv.Row(k) = v1(k,0)*dv2.Row(k)+dv1(k,0)*v2.Row(k);
      }
  }

  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & mir,
                              AFlatMatrix<> result,
                              AFlatMatrix<> deriv) const
  {
    STACK_ARRAY(SIMD<double>, hmem1, mir.Size());
    AFlatMatrix<> temp1(1, mir.IR().GetNIP(), hmem1);
    STACK_ARRAY(SIMD<double>, hmem2, mir.Size());
    AFlatMatrix<> deriv1(1, mir.IR().GetNIP(), hmem2);
    c1->EvaluateDeriv(mir, temp1, deriv1);
    c2->EvaluateDeriv(mir, result, deriv);

    for (int i = 0; i < result.Height(); i++)
      for (int k = 0; k < mir.Size(); k++)
        {
          deriv.Get(i,k) = deriv.Get(i,k)*temp1.Get(0,k) + result.Get(i,k) * deriv1.Get(0,k);
          result.Get(i,k) *= temp1.Get(0,k);
        }
  }
  
  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & mir,
                              FlatArray<AFlatMatrix<>*> input,
                              FlatArray<AFlatMatrix<>*> dinput,
                              AFlatMatrix<> result,
                              AFlatMatrix<> deriv) const
  {
    AFlatMatrix<> v1 = *input[0], v2 = *input[1];
    AFlatMatrix<> dv1 = *dinput[0], dv2 = *dinput[1];

    for (int i = 0; i < result.Height(); i++)
      for (int k = 0; k < mir.Size(); k++)
        {
          result.Get(i,k) = v1.Get(0,k)*v2.Get(i,k);
          deriv.Get(i,k) = v1.Get(0,k)*dv2.Get(i,k)+dv1.Get(0,k)*v2.Get(i,k);
        }
  }
  

  
  
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
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
};


class MultVecVecCoefficientFunction : public T_CoefficientFunction<MultVecVecCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
  int dim1;
public:
  MultVecVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                 shared_ptr<CoefficientFunction> ac2)
    : T_CoefficientFunction<MultVecVecCoefficientFunction>(1, ac1->IsComplex() || ac2->IsComplex()), c1(ac1), c2(ac2)
  {
    dim1 = c1->Dimension();
    if (dim1 != c2->Dimension())
      throw Exception("MultVecVec : dimensions don't fit");
  }
  
  // virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  // virtual int Dimension() const { return 1; }
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    CodeExpr result;
    TraverseDimensions( c1->Dimensions(), [&](int ind, int i, int j) {
        int i2, j2;
        GetIndex( c2->Dimensions(), ind, i2, j2 );
        result += Var(inputs[0],i,j) * Var(inputs[1],i2,j2);
    });
    code.body += Var(index).Assign(result.S());
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get(), c2.get() }); }  
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
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
                        FlatVector<Complex> result) const
  {
    Vector<Complex> v1(dim1), v2(dim1);
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
    STACK_ARRAY(double, hmem1, ir.Size()*dim1);
    FlatMatrix<> temp1(ir.Size(), dim1, hmem1);
    STACK_ARRAY(double, hmem2, ir.Size()*dim1);
    FlatMatrix<> temp2(ir.Size(), dim1, hmem2);

    c1->Evaluate(ir, temp1);
    c2->Evaluate(ir, temp2);
    for (int i = 0; i < ir.Size(); i++)
      result(i,0) = InnerProduct(temp1.Row(i), temp2.Row(i));
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatArray<FlatMatrix<>*> input,
                         FlatMatrix<double> result) const
  {
    FlatMatrix<> temp1 = *input[0];
    FlatMatrix<> temp2 = *input[1];
    for (int i = 0; i < ir.Size(); i++)
      result(i,0) = InnerProduct(temp1.Row(i), temp2.Row(i));
  }


  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<T> values) const
  {
    size_t w = ir.Size();
    __assume (w > 0);

    size_t dim = Dimension();
    STACK_ARRAY(T, hmem, 2*dim*w);
    FlatMatrix<T> temp1(dim, w, &hmem[0]);
    FlatMatrix<T> temp2(dim, w, &hmem[dim*w]);
    
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


  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   FlatArray<BareSliceMatrix<T>> input,                       
                   BareSliceMatrix<T> values) const
  {
    auto in0 = input[0];
    auto in1 = input[1];
    size_t dim = Dimension();
    size_t np = ir.Size();

    for (size_t i = 0; i < np; i++)
      {
        T sum{0.0};
        for (size_t j = 0; j < dim; j++)
          sum += in0(j,i) * in1(j,i);
        values(0,i) = sum; 
      }    
  }  
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    auto in0 = *input[0];
    auto in1 = *input[1];
    
    for (size_t i = 0; i < values.VWidth(); i++)
      {
        SIMD<double> sum = 0.0;
        for (size_t j = 0; j < dim1; j++)
          sum += in0.Get(j,i) * in1.Get(j,i);
        values.Get(i) = sum;
      }
  }

  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    Matrix<> v1(mir.Size(), dim1), v2(mir.Size(),dim1);
    Matrix<> dv1(mir.Size(), dim1), dv2(mir.Size(), dim1);
    c1->EvaluateDeriv (mir, v1, dv1);
    c2->EvaluateDeriv (mir, v2, dv2);
    for (int k = 0; k < mir.Size(); k++)
      {
        result(k,0) = InnerProduct (v1.Row(k), v2.Row(k));
        deriv(k,0) = InnerProduct (v1.Row(k), dv2.Row(k))+InnerProduct(v2.Row(k),dv1.Row(k));
      }
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    Matrix<> v1(mir.Size(), dim1), v2(mir.Size(), dim1);
    Matrix<> dv1(mir.Size(), dim1), dv2(mir.Size(), dim1);
    Matrix<> ddv1(mir.Size(), dim1), ddv2(mir.Size(), dim1);
    c1->EvaluateDDeriv (mir, v1, dv1, ddv1);
    c2->EvaluateDDeriv (mir, v2, dv2, ddv2);

    for (int k = 0; k < mir.Size(); k++)
      {
        result(k,0) = InnerProduct (v1.Row(k), v2.Row(k));
        deriv(k,0) = InnerProduct (v1.Row(k), dv2.Row(k))+InnerProduct(v2.Row(k),dv1.Row(k));
        dderiv(k,0) = InnerProduct (v1.Row(k), ddv2.Row(k))+
          2*InnerProduct(dv1.Row(k),dv2.Row(k))+InnerProduct(ddv1.Row(k),v2.Row(k));
      }

  }


  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const
  {
    FlatMatrix<> v1 = *input[0], v2 = *input[1];
    FlatMatrix<> dv1 = *dinput[0], dv2 = *dinput[1];
    
    for (int k = 0; k < mir.Size(); k++)
      {
        result(k,0) = InnerProduct (v1.Row(k), v2.Row(k));
        deriv(k,0) = InnerProduct (v1.Row(k), dv2.Row(k))+InnerProduct(v2.Row(k),dv1.Row(k));
      }
  }

  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                               FlatArray<FlatMatrix<>*> input,
                               FlatArray<FlatMatrix<>*> dinput,
                               FlatArray<FlatMatrix<>*> ddinput,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const
  {
    FlatMatrix<> v1 = *input[0], v2 = *input[1];
    FlatMatrix<> dv1 = *dinput[0], dv2 = *dinput[1];
    FlatMatrix<> ddv1 = *ddinput[0], ddv2 = *ddinput[1];
    
    for (int k = 0; k < mir.Size(); k++)
      {
        result(k,0) = InnerProduct (v1.Row(k), v2.Row(k));
        deriv(k,0) = InnerProduct (v1.Row(k), dv2.Row(k))+InnerProduct(v2.Row(k),dv1.Row(k));
        dderiv(k,0) = InnerProduct (v1.Row(k), ddv2.Row(k))+
          2*InnerProduct(dv1.Row(k),dv2.Row(k))+InnerProduct(ddv1.Row(k),v2.Row(k));
      }
  }

  virtual bool ElementwiseConstant () const
  { return c1->ElementwiseConstant() && c2->ElementwiseConstant(); }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
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

};

template <int DIM>
class T_MultVecVecCoefficientFunction : public T_CoefficientFunction<T_MultVecVecCoefficientFunction<DIM>>
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
public:
  T_MultVecVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                   shared_ptr<CoefficientFunction> ac2)
    : T_CoefficientFunction<T_MultVecVecCoefficientFunction<DIM>>(1, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
  {
    if (DIM != c1->Dimension() || DIM != c2->Dimension())
      throw Exception("T_MultVecVec : dimensions don't fit");
  }
  
  // virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  // virtual int Dimension() const { return 1; }
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    CodeExpr result;
    TraverseDimensions( c1->Dimensions(), [&](int ind, int i, int j) {
        int i2, j2;
        GetIndex( c2->Dimensions(), ind, i2, j2 );
        result += Var(inputs[0],i,j) * Var(inputs[1],i2,j2);
    });
    code.body += Var(index).Assign(result.S());
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get(), c2.get() }); }  

  using T_CoefficientFunction<T_MultVecVecCoefficientFunction<DIM>>::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    Vec<DIM> v1, v2;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    Vec<DIM,Complex> v1, v2;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
    STACK_ARRAY(double, hmem1, ir.Size()*DIM);
    FlatMatrixFixWidth<DIM> temp1(ir.Size(), hmem1);
    STACK_ARRAY(double, hmem2, ir.Size()*DIM);
    FlatMatrixFixWidth<DIM> temp2(ir.Size(), hmem2);

    c1->Evaluate(ir, temp1);
    c2->Evaluate(ir, temp2);

    for (size_t i = 0; i < ir.Size(); i++)
      result(i,0) = InnerProduct(temp1.Row(i), temp2.Row(i));
  }
  /*
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
  {
    STACK_ARRAY(SIMD<double>, hmem1, DIM*values.Width());
    STACK_ARRAY(SIMD<double>, hmem2, DIM*values.Width());
    AFlatMatrix<double> temp1(DIM, values.Width(), &hmem1[0]);
    AFlatMatrix<double> temp2(DIM, values.Width(), &hmem2[0]);
    
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, temp2);

    for (size_t i = 0; i < values.VWidth(); i++)
      {
        SIMD<double> sum = 0.0;
        for (size_t j = 0; j < DIM; j++)
          sum += temp1.Get(j,i) * temp2.Get(j,i);
        values.Get(i) = sum; 
      }
  }
  */



  /*
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, ABareSliceMatrix<> values) const
  {
    size_t w = ir.Size();
    __assume (w > 0);
    
    STACK_ARRAY(SIMD<double>, hmem, 2*DIM*w);
    AFlatMatrix<double> temp1(DIM, w*SIMD<double>::Size(), &hmem[0]);
    AFlatMatrix<double> temp2(DIM, w*SIMD<double>::Size(), &hmem[DIM*w]);
    
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, temp2);

    for (size_t i = 0; i < w; i++)
      {
        SIMD<double> sum = 0.0;
        for (size_t j = 0; j < DIM; j++)
          sum += temp1.Get(j,i) * temp2.Get(j,i);
        values.Get(0,i) = sum; 
      }
  }

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, ABareSliceMatrix<Complex> values) const
  {
    size_t w = ir.Size();
    __assume (w > 0);
    
    STACK_ARRAY(SIMD<Complex>, hmem, 2*DIM*w);
    AFlatMatrix<Complex> temp1(DIM, w*SIMD<double>::Size(), &hmem[0]);
    AFlatMatrix<Complex> temp2(DIM, w*SIMD<double>::Size(), &hmem[DIM*w]);
    
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, temp2);

    for (size_t i = 0; i < w; i++)
      {
        SIMD<Complex> sum = Complex(0.0);
        for (size_t j = 0; j < DIM; j++)
          sum += temp1.Get(j,i) * temp2.Get(j,i);
        values.Get(0,i) = sum; 
      }
  }
  */

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<T> values) const
  {
    size_t w = ir.Size();
    __assume (w > 0);
    
    STACK_ARRAY(T, hmem, 2*DIM*w);
    FlatMatrix<T> temp1(DIM, w, &hmem[0]);
    FlatMatrix<T> temp2(DIM, w, &hmem[DIM*w]);
    
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

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   FlatArray<BareSliceMatrix<T>> input,                       
                   BareSliceMatrix<T> values) const
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
  
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    auto in0 = *input[0];
    auto in1 = *input[1];
    
    for (size_t i = 0; i < values.VWidth(); i++)
      {
        SIMD<double> sum = 0.0;
        for (size_t j = 0; j < DIM; j++)
          sum += in0.Get(j,i) * in1.Get(j,i);
        values.Get(i) = sum; 
      }
  }

  
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<Complex> result) const
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



  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatArray<FlatMatrix<>*> input,
                         FlatMatrix<double> result) const
  {
    FlatMatrix<> temp1 = *input[0];
    FlatMatrix<> temp2 = *input[1];
    for (int i = 0; i < ir.Size(); i++)
      result(i,0) = InnerProduct(temp1.Row(i), temp2.Row(i));
  }

  

  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    size_t si = DIM*mir.Size();
    STACK_ARRAY(double, hmem, 4*si);
    FlatMatrix<> v1(mir.Size(), DIM, &hmem[0]), v2(mir.Size(),DIM, &hmem[si]);
    FlatMatrix<> dv1(mir.Size(), DIM, &hmem[2*si]), dv2(mir.Size(), DIM, &hmem[3*si]);
    
    c1->EvaluateDeriv (mir, v1, dv1);
    c2->EvaluateDeriv (mir, v2, dv2);
    for (int k = 0; k < mir.Size(); k++)
      {
        result(k,0) = InnerProduct (v1.Row(k), v2.Row(k));
        deriv(k,0) = InnerProduct (v1.Row(k), dv2.Row(k))+InnerProduct(v2.Row(k),dv1.Row(k));
      }
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    size_t si = DIM*mir.Size();
    STACK_ARRAY(double, hmem, 6*si);
    FlatMatrix<> v1(mir.Size(), DIM, &hmem[0]), v2(mir.Size(),DIM, &hmem[si]);
    FlatMatrix<> dv1(mir.Size(), DIM, &hmem[2*si]), dv2(mir.Size(), DIM, &hmem[3*si]);
    FlatMatrix<> ddv1(mir.Size(), DIM, &hmem[4*si]), ddv2(mir.Size(), DIM, &hmem[5*si]);

    c1->EvaluateDDeriv (mir, v1, dv1, ddv1);
    c2->EvaluateDDeriv (mir, v2, dv2, ddv2);

    /*
    cout << "eval dderiv:" << endl
         << "v1 = " << v1 << ", dv1 = " << dv1 << "; ddv1 = " << ddv1 << endl
         << "v2 = " << v2 << ", dv2 = " << dv2 << "; ddv2 = " << ddv2 << endl;
    */
    for (int k = 0; k < mir.Size(); k++)
      {
        result(k,0) = InnerProduct (v1.Row(k), v2.Row(k));
        deriv(k,0) = InnerProduct (v1.Row(k), dv2.Row(k))+InnerProduct(v2.Row(k),dv1.Row(k));
        dderiv(k,0) = InnerProduct (v1.Row(k), ddv2.Row(k))+
          2*InnerProduct(dv1.Row(k),dv2.Row(k))+InnerProduct(ddv1.Row(k),v2.Row(k));
      }
    // cout << "res = " << result << ", deriv = " << deriv << ", dderiv = " << dderiv << endl;
  }



  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const
  {
    FlatMatrix<> v1 = *input[0], v2 = *input[1];
    FlatMatrix<> dv1 = *dinput[0], dv2 = *dinput[1];
    
    for (int k = 0; k < mir.Size(); k++)
      {
        result(k,0) = InnerProduct (v1.Row(k), v2.Row(k));
        deriv(k,0) = InnerProduct (v1.Row(k), dv2.Row(k))+InnerProduct(v2.Row(k),dv1.Row(k));
      }
  }

  
  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                               FlatArray<FlatMatrix<>*> input,
                               FlatArray<FlatMatrix<>*> dinput,
                               FlatArray<FlatMatrix<>*> ddinput,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const
  {
    FlatMatrix<> v1 = *input[0], v2 = *input[1];
    FlatMatrix<> dv1 = *dinput[0], dv2 = *dinput[1];
    FlatMatrix<> ddv1 = *ddinput[0], ddv2 = *ddinput[1];
    
    for (int k = 0; k < mir.Size(); k++)
      {
        result(k,0) = InnerProduct (v1.Row(k), v2.Row(k));
        deriv(k,0) = InnerProduct (v1.Row(k), dv2.Row(k))+InnerProduct(v2.Row(k),dv1.Row(k));
        dderiv(k,0) = InnerProduct (v1.Row(k), ddv2.Row(k))+
          2*InnerProduct(dv1.Row(k),dv2.Row(k))+InnerProduct(ddv1.Row(k),v2.Row(k));
      }
  }



  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & mir,
                              AFlatMatrix<> result,
                              AFlatMatrix<> deriv) const
  {
    size_t si = DIM*mir.Size();
    STACK_ARRAY(SIMD<double>, hmem, 4*si);
    AFlatMatrix<> v1(DIM, mir.IR().GetNIP(), &hmem[0]), v2(DIM, mir.IR().GetNIP(), &hmem[si]);
    AFlatMatrix<> dv1(DIM, mir.IR().GetNIP(), &hmem[2*si]), dv2(DIM, mir.IR().GetNIP(), &hmem[3*si]);
    
    c1->EvaluateDeriv (mir, v1, dv1);
    c2->EvaluateDeriv (mir, v2, dv2);
    
    for (size_t k = 0; k < mir.Size(); k++)
      {
        SIMD<double> sum = 0.0, dsum = 0.0;
        for (size_t i = 0; i < DIM; i++)
          {
            sum += v1.Get(i,k) * v2.Get(i,k);
            dsum += v1.Get(i,k) * dv2.Get(i,k) + dv1.Get(i,k) * v2.Get(i,k);
          }
        result.Get(0,k) = sum;
        deriv.Get(0,k) = dsum;
      }
  }
  



  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & mir,
                              FlatArray<AFlatMatrix<>*> input,
                              FlatArray<AFlatMatrix<>*> dinput,
                              AFlatMatrix<> result,
                              AFlatMatrix<> deriv) const
  {
    AFlatMatrix<> v1 = *input[0], v2 = *input[1];
    AFlatMatrix<> dv1 = *dinput[0], dv2 = *dinput[1];
    
    for (size_t k = 0; k < mir.Size(); k++)
      {
        SIMD<double> sum = 0.0, dsum = 0.0;
        for (size_t i = 0; i < DIM; i++)
          {
            sum += v1.Get(i,k) * v2.Get(i,k);
            dsum += v1.Get(i,k) * dv2.Get(i,k) + dv1.Get(i,k) * v2.Get(i,k);
          }
        result.Get(0,k) = sum;
        deriv.Get(0,k) = dsum;
      }
  }
  
  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & mir,
                               FlatArray<AFlatMatrix<>*> input,
                               FlatArray<AFlatMatrix<>*> dinput,
                               FlatArray<AFlatMatrix<>*> ddinput,
                               AFlatMatrix<> result,
                               AFlatMatrix<> deriv,
                               AFlatMatrix<> dderiv) const
  {
    AFlatMatrix<> v1 = *input[0], v2 = *input[1];
    AFlatMatrix<> dv1 = *dinput[0], dv2 = *dinput[1];
    AFlatMatrix<> ddv1 = *ddinput[0], ddv2 = *ddinput[1];
    
    for (size_t k = 0; k < mir.Size(); k++)
      {
        SIMD<double> sum = 0.0, dsum = 0.0, ddsum = 0.0;
        for (size_t i = 0; i < DIM; i++)
          {
            sum += v1.Get(i,k) * v2.Get(i,k);
            dsum += v1.Get(i,k) * dv2.Get(i,k) + dv1.Get(i,k) * v2.Get(i,k);
            ddsum += v1.Get(i,k) * ddv2.Get(i,k) + 2*dv1.Get(i,k) * dv2.Get(i,k)
              + ddv1.Get(i,k) * v2.Get(i,k);
          }
        result.Get(0,k) = sum;
        deriv.Get(0,k) = dsum;
        dderiv.Get(0,k) = ddsum;
      }
  }


  
  virtual bool ElementwiseConstant () const
  { return c1->ElementwiseConstant() && c2->ElementwiseConstant(); }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
  {
    Vector<bool> v1(DIM), v2(DIM), d1(DIM), dd1(DIM), d2(DIM), dd2(DIM);
    c1->NonZeroPattern (ud, v1, d1, dd1);
    c2->NonZeroPattern (ud, v2, d2, dd2);
    // cout << "nonzero, v1 = " << v1 << ", d1 = " << d1 << ", dd1 = " << dd1 << endl;
    // cout << "nonzero, v2 = " << v2 << ", d2 = " << d2 << ", dd2 = " << dd2 << endl;
    bool nz = false, nzd = false, nzdd = false;
    for (int i = 0; i < DIM; i++)
      {
        if (v1(i) && v2(i)) nz = true;
        if ((v1(i) && d2(i)) || (d1(i) && v2(i))) nzd = true;
        if ((v1(i) && dd2(i)) || (d1(i) && d2(i)) || (dd1(i) && v2(i))) nzdd = true;
      }
    // cout << "nz = " << nz << ", nzd = " << nzd << ", nzdd = " << nzdd << endl;
    nonzero = nz;
    nonzero_deriv = nzd;
    nonzero_dderiv = nzdd;
  }

};






template <typename TIN>
class NormCoefficientFunction : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  int dim1;
public:
  NormCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : CoefficientFunction(1, false), c1(ac1)
  {
    dim1 = c1->Dimension();
  }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get() }); }  
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    VectorMem<10,TIN> v1(dim1);
    c1->Evaluate (ip, v1);
    result(0) = L2Norm(v1);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    result(0) = res(0);
  }

  
  virtual bool ElementwiseConstant () const
  { return c1->ElementwiseConstant(); }

  
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix <> result) const
  {
    STACK_ARRAY(double,hmem,ir.Size()*dim1*sizeof(TIN)/sizeof(double));
    FlatMatrix<TIN> inval(ir.IR().GetNIP(), dim1, reinterpret_cast<TIN*>(&hmem[0]));
    c1->Evaluate (ir, inval);
    for (size_t i = 0; i < result.Height(); i++)
      result(i,0) = L2Norm(inval.Row(i));
  }


  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
  {
    if (typeid(TIN)==typeid(Complex)) throw ExceptionNOSIMD("CF Norm of complex cannot use simds");
    STACK_ARRAY(SIMD<double>,hmem,ir.Size()*dim1);
    AFlatMatrix<> inval(dim1, ir.IR().GetNIP(), &hmem[0]);
    c1->Evaluate (ir, inval);
    for (size_t i = 0; i < ir.Size(); i++)
      {
        SIMD<double> sum = 0;
        for (size_t j = 0; j < dim1; j++)
          sum += inval.Get(j,i)*inval.Get(j,i);
        values(0,i) = sqrt(sum);
      }
  }


  
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
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

};



  

class MultMatMatCoefficientFunction : public T_CoefficientFunction<MultMatMatCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
  // Array<int> dims;
  int inner_dim;
public:
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
    // dims = { dims_c1[0], dims_c2[1] };
    SetDimensions( Array<int> ({ dims_c1[0], dims_c2[1] }));
    inner_dim = dims_c1[1];
  }
  
  // virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  // virtual int Dimension() const { return dims[0]*dims[1]; }
  // virtual Array<int> Dimensions() const { return Array<int> (dims); } 

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
    FlatArray<int> hdims = Dimensions();
      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1])) {
          CodeExpr s;
          for (int k : Range(inner_dim))
            s += Var(inputs[0], i, k) * Var(inputs[1], k, j);
          code.body += Var(index, i, j).Assign(s);
        }
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get(), c2.get() }); }  


  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
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

  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    throw Exception ("MultMatMatCF:: scalar evaluate for matrix called");
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const
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
                         FlatVector<Complex> result) const
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

  virtual void Evaluate (const BaseMappedIntegrationRule & mir,
                         FlatMatrix<> result) const
  {
    // Matrix<> va(mir.Size(), dims[0]*inner_dim);
    // Matrix<> vb(mir.Size(), dims[1]*inner_dim);
    FlatArray<int> hdims = Dimensions();
    STACK_ARRAY(double, mema, mir.Size()*hdims[0]*inner_dim);
    STACK_ARRAY(double, memb, mir.Size()*hdims[1]*inner_dim);
    FlatMatrix<> va(mir.Size(), hdims[0]*inner_dim, mema);
    FlatMatrix<> vb(mir.Size(), hdims[1]*inner_dim, memb);

    c1->Evaluate (mir, va);
    c2->Evaluate (mir, vb);

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        FlatMatrix<> b(inner_dim, hdims[1], &vb(i,0));
        FlatMatrix<> c(hdims[0], hdims[1], &result(i,0));
        c = a*b;
      }
  }

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & mir, BareSliceMatrix<SIMD<double>> values) const
  {
    FlatArray<int> hdims = Dimensions();    
    STACK_ARRAY(SIMD<double>, hmem1, mir.IR().Size()*hdims[0]*inner_dim);
    STACK_ARRAY(SIMD<double>, hmem2, mir.IR().Size()*hdims[1]*inner_dim);
    AFlatMatrix<double> va(hdims[0]*inner_dim, mir.IR().GetNIP(), &hmem1[0]);
    AFlatMatrix<double> vb(hdims[1]*inner_dim, mir.IR().GetNIP(), &hmem2[0]);
    c1->Evaluate (mir, va);
    c2->Evaluate (mir, vb);
    values.AddSize(Dimension(),mir.Size()) = 0.0;

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
              row_c(i) += row_a.Get(i) * row_b.Get(i);
            // row_c = pw_mult (row_a, row_b);
          }
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & mir, BareSliceMatrix<T> values) const
  {
    FlatArray<int> hdims = Dimensions();    
    STACK_ARRAY(T, hmem1, mir.Size()*hdims[0]*inner_dim);
    STACK_ARRAY(T, hmem2, mir.Size()*hdims[1]*inner_dim);
    FlatMatrix<T> va(hdims[0]*inner_dim, mir.Size(), &hmem1[0]);
    FlatMatrix<T> vb(hdims[1]*inner_dim, mir.Size(), &hmem2[0]);

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

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   FlatArray<BareSliceMatrix<T>> input,                       
                   BareSliceMatrix<T> values) const
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

  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & mir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    FlatArray<int> hdims = Dimensions();    
    auto va = *input[0];
    auto vb = *input[1];
    values = 0.0;

    size_t dims0 = hdims[0];
    size_t dims1 = hdims[1];
    size_t idim = inner_dim;
    for (size_t j = 0; j < dims0; j++)
      for (size_t k = 0; k < dims1; k++)
        for (size_t l = 0; l < idim; l++)
          {
            auto row_a = va.Row(j*inner_dim+l);
            auto row_b = vb.Row(l*dims1+k);
            auto row_c = values.Row(j*dims1+k);
            for (size_t i = 0; i < mir.Size(); i++)
              row_c.Get(i) += row_a.Get(i) * row_b.Get(i);
          }    
  }
  
  
  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    FlatArray<int> hdims = Dimensions();    
    Matrix<> va(mir.Size(), hdims[0]*inner_dim);
    Matrix<> vb(mir.Size(), hdims[1]*inner_dim);
    Matrix<> vda(mir.Size(), hdims[0]*inner_dim);
    Matrix<> vdb(mir.Size(), hdims[1]*inner_dim);
    c1->EvaluateDeriv (mir, va, vda);
    c2->EvaluateDeriv (mir, vb, vdb);

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        FlatMatrix<> b(inner_dim, hdims[1], &vb(i,0));
        FlatMatrix<> da(hdims[0], inner_dim, &vda(i,0));
        FlatMatrix<> db(inner_dim, hdims[1], &vdb(i,0));
        FlatMatrix<> c(hdims[0], hdims[1], &result(i,0));
        FlatMatrix<> dc(hdims[0], hdims[1], &deriv(i,0));
        c = a*b;
        dc = a*db+da*b;
      }
  }
  
  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    FlatArray<int> hdims = Dimensions();        
    Matrix<> va(mir.Size(), hdims[0]*inner_dim);
    Matrix<> vb(mir.Size(), hdims[1]*inner_dim);
    Matrix<> vda(mir.Size(), hdims[0]*inner_dim);
    Matrix<> vdb(mir.Size(), hdims[1]*inner_dim);
    Matrix<> vdda(mir.Size(), hdims[0]*inner_dim);
    Matrix<> vddb(mir.Size(), hdims[1]*inner_dim);
    c1->EvaluateDDeriv (mir, va, vda, vdda);
    c2->EvaluateDDeriv (mir, vb, vdb, vddb);

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        FlatMatrix<> b(inner_dim, hdims[1], &vb(i,0));
        FlatMatrix<> da(hdims[0], inner_dim, &vda(i,0));
        FlatMatrix<> db(inner_dim, hdims[1], &vdb(i,0));
        FlatMatrix<> dda(hdims[0], inner_dim, &vdda(i,0));
        FlatMatrix<> ddb(inner_dim, hdims[1], &vddb(i,0));
        FlatMatrix<> c(hdims[0], hdims[1], &result(i,0));
        FlatMatrix<> dc(hdims[0], hdims[1], &deriv(i,0));
        FlatMatrix<> ddc(hdims[0], hdims[1], &dderiv(i,0));
        c = a*b;
        dc = a*db+da*b;
        ddc = a*ddb+2*da*db+dda*b;
      }
  }


  



  virtual void Evaluate(const BaseMappedIntegrationRule & mir,
                        FlatArray<FlatMatrix<>*> input,
                        FlatMatrix<> result) const
  {
    FlatArray<int> hdims = Dimensions();        
    FlatMatrix<> va = *input[0], vb = *input[1];

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        FlatMatrix<> b(inner_dim, hdims[1], &vb(i,0));
        FlatMatrix<> c(hdims[0], hdims[1], &result(i,0));
        c = a*b;
      }
  }



  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatArray<FlatMatrix<>*> input,
                             FlatArray<FlatMatrix<>*> dinput,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    FlatArray<int> hdims = Dimensions();        
    FlatMatrix<> va = *input[0], vb = *input[1];
    FlatMatrix<> vda = *dinput[0], vdb = *dinput[1];

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        FlatMatrix<> b(inner_dim, hdims[1], &vb(i,0));
        FlatMatrix<> da(hdims[0], inner_dim, &vda(i,0));
        FlatMatrix<> db(inner_dim, hdims[1], &vdb(i,0));
        FlatMatrix<> c(hdims[0], hdims[1], &result(i,0));
        FlatMatrix<> dc(hdims[0], hdims[1], &deriv(i,0));
        c = a*b;
        dc = a*db+da*b;
      }
  }


  
  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatArray<FlatMatrix<>*> ddinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    FlatArray<int> hdims = Dimensions();        
    FlatMatrix<> va = *input[0], vb = *input[1];
    FlatMatrix<> vda = *dinput[0], vdb = *dinput[1];
    FlatMatrix<> vdda = *ddinput[0], vddb = *ddinput[1];

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        FlatMatrix<> b(inner_dim, hdims[1], &vb(i,0));
        FlatMatrix<> da(hdims[0], inner_dim, &vda(i,0));
        FlatMatrix<> db(inner_dim, hdims[1], &vdb(i,0));
        FlatMatrix<> dda(hdims[0], inner_dim, &vdda(i,0));
        FlatMatrix<> ddb(inner_dim, hdims[1], &vddb(i,0));
        FlatMatrix<> c(hdims[0], hdims[1], &result(i,0));
        FlatMatrix<> dc(hdims[0], hdims[1], &deriv(i,0));
        FlatMatrix<> ddc(hdims[0], hdims[1], &dderiv(i,0));
        c = a*b;
        dc = a*db+da*b;
        ddc = a*ddb+2*da*db+dda*b;
      }
  }


  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & mir,
                              FlatArray<AFlatMatrix<double>*> input,
                              FlatArray<AFlatMatrix<double>*> dinput,
                              AFlatMatrix<double> values,
                              AFlatMatrix<double> deriv) const
  {
    FlatArray<int> hdims = Dimensions();        
    auto va = *input[0];
    auto dva = *dinput[0];
    auto vb = *input[1];
    auto dvb = *dinput[1];
    values = 0.0;
    deriv = 0.0;
    
    for (int j = 0; j < hdims[0]; j++)
      for (int k = 0; k < hdims[1]; k++)
        for (int l = 0; l < inner_dim; l++)
          {
            auto row_a = va.Row(j*inner_dim+l);
            auto row_da = dva.Row(j*inner_dim+l);
            auto row_b = vb.Row(l*hdims[1]+k);
            auto row_db = dvb.Row(l*hdims[1]+k);
            auto row_c = values.Row(j*hdims[1]+k);
            auto row_dc = deriv.Row(j*hdims[1]+k);

            for (int i = 0; i < mir.Size(); i++)
              {
                row_c.Get(i) += row_a.Get(i) * row_b.Get(i);
                row_dc.Get(i) += row_a.Get(i) * row_db.Get(i) + row_da.Get(i) * row_b.Get(i);
              }    
          }
  }


  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & mir,
                               FlatArray<AFlatMatrix<double>*> input,
                               FlatArray<AFlatMatrix<double>*> dinput,
                               FlatArray<AFlatMatrix<double>*> ddinput,
                               AFlatMatrix<double> values,
                               AFlatMatrix<double> deriv,
                               AFlatMatrix<double> dderiv) const
  {
    FlatArray<int> hdims = Dimensions();        
    auto va = *input[0];
    auto dva = *dinput[0];
    auto ddva = *ddinput[0];
    auto vb = *input[1];
    auto dvb = *dinput[1];
    auto ddvb = *ddinput[1];
    values = 0.0;
    deriv = 0.0;
    dderiv = 0.0;
    
    for (int j = 0; j < hdims[0]; j++)
      for (int k = 0; k < hdims[1]; k++)
        for (int l = 0; l < inner_dim; l++)
          {
            auto row_a = va.Row(j*inner_dim+l);
            auto row_da = dva.Row(j*inner_dim+l);
            auto row_dda = ddva.Row(j*inner_dim+l);
            auto row_b = vb.Row(l*hdims[1]+k);
            auto row_db = dvb.Row(l*hdims[1]+k);
            auto row_ddb = ddvb.Row(l*hdims[1]+k);
            auto row_c = values.Row(j*hdims[1]+k);
            auto row_dc = deriv.Row(j*hdims[1]+k);
            auto row_ddc = dderiv.Row(j*hdims[1]+k);
            
            for (int i = 0; i < mir.Size(); i++)
              {
                row_c.Get(i) += row_a.Get(i) * row_b.Get(i);
                row_dc.Get(i) += row_a.Get(i) * row_db.Get(i) + row_da.Get(i) * row_b.Get(i);
                row_ddc.Get(i) += row_a.Get(i) * row_ddb.Get(i) +
                  2 * row_da.Get(i) * row_db.Get(i) +
                  row_dda.Get(i) * row_b.Get(i);
              }
          }
  }

    
};







class MultMatVecCoefficientFunction : public T_CoefficientFunction<MultMatVecCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
  // Array<int> dims;
  int inner_dim;
public:
  MultMatVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                 shared_ptr<CoefficientFunction> ac2)
    : T_CoefficientFunction(1, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
  {
    auto dims_c1 = c1 -> Dimensions();
    auto dims_c2 = c2 -> Dimensions();
    if (dims_c1.Size() != 2 || dims_c2.Size() != 1)
      throw Exception("Not a mat-vec multiplication");
    if (dims_c1[1] != dims_c2[0])
      throw Exception(string ("Matrix dimensions don't fit: mat is ") +
                      ToLiteral(dims_c1[0]) + " x " + ToLiteral(dims_c1[1]) + ", vec is " + ToLiteral(dims_c2[0]));
    // dims = Array<int> ({ dims_c1[0] });
    SetDimensions (Array<int> ({ dims_c1[0] }));
    inner_dim = dims_c1[1];
  }
  
  // virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  // virtual int Dimension() const { return dims[0]; }
  // virtual Array<int> Dimensions() const { return Array<int> (dims); } 

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get(), c2.get() }); }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
      auto dims = c1->Dimensions();
      for (int i : Range(dims[0])) {
        CodeExpr s;
        for (int j : Range(dims[1]))
            s += Var(inputs[0], i, j) * Var(inputs[1], j);
	code.body += Var(index, i).Assign(s);
      }
  }

  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
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

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    throw Exception ("MultMatVecCF:: scalar evaluate for matrix called");
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const
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
                         FlatVector<Complex> result) const
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

  virtual void Evaluate (const BaseMappedIntegrationRule & mir,
                         FlatMatrix<> result) const
  {
    FlatArray<int> hdims = Dimensions();    
    Matrix<> va(mir.Size(), hdims[0]*inner_dim);
    Matrix<> vb(mir.Size(), inner_dim);
    c1->Evaluate (mir, va);
    c2->Evaluate (mir, vb);

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        result.Row(i) = a * vb.Row(i);
      }
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<T> values) const
  {
    FlatArray<int> hdims = Dimensions();    
    STACK_ARRAY(T, hmem1, ir.Size()*hdims[0]*inner_dim);
    STACK_ARRAY(T, hmem2, ir.Size()*inner_dim);
    FlatMatrix<T> temp1(hdims[0]*inner_dim, ir.Size(), &hmem1[0]);
    FlatMatrix<T> temp2(inner_dim, ir.Size(), &hmem2[0]);
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, temp2);
    values.AddSize(Dimension(),ir.Size()) = T(0.0);
    for (size_t i = 0; i < hdims[0]; i++)
      for (size_t j = 0; j < inner_dim; j++)
        for (size_t k = 0; k < ir.Size(); k++)
          values(i,k) += temp1(i*inner_dim+j, k) * temp2(j,k);
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   FlatArray<BareSliceMatrix<T>> input,                       
                   BareSliceMatrix<T> values) const
  {
    auto va = input[0];
    auto vb = input[1];
    
    FlatArray<int> hdims = Dimensions();    
    values.AddSize(Dimension(),ir.Size()) = T(0.0);
    
    for (size_t i = 0; i < hdims[0]; i++)
      for (size_t j = 0; j < inner_dim; j++)
        for (size_t k = 0; k < ir.Size(); k++)
          values(i,k) += va(i*inner_dim+j, k) * vb(j,k);
  }
    
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
  {
    FlatArray<int> hdims = Dimensions();    
    STACK_ARRAY(SIMD<double>, hmem1, (ir.IR().GetNIP()+8)*hdims[0]*inner_dim);
    STACK_ARRAY(SIMD<double>, hmem2, (ir.IR().GetNIP()+8)*inner_dim);
    AFlatMatrix<double> temp1(hdims[0]*inner_dim, ir.IR().GetNIP(), &hmem1[0]);
    AFlatMatrix<double> temp2(inner_dim, ir.IR().GetNIP(), &hmem2[0]);
    c1->Evaluate (ir, temp1);
    c2->Evaluate (ir, temp2);
    values.AddSize(Dimension(),ir.Size()) = 0.0;
    for (size_t i = 0; i < hdims[0]; i++)
      for (size_t j = 0; j < inner_dim; j++)
        for (size_t k = 0; k < ir.Size(); k++)
          values(i,k) += temp1.Get(i*inner_dim+j, k) * temp2.Get(j,k);
  }
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & mir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    FlatArray<int> hdims = Dimensions();    
    auto in0 = *input[0];
    auto in1 = *input[1];
    values = 0.0;
    auto _inner_dim = inner_dim;
    if (_inner_dim <= 0) return;
    int ii = 0;
    for (auto i : Range(hdims[0]))
      for (auto j : Range(_inner_dim))
        // values.Row(i) += pw_mult (in0.Row(i*_inner_dim+j), in1.Row(j));
        values.Row(i) += pw_mult (in0.Row(ii++), in1.Row(j));
  }

  
  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    FlatArray<int> hdims = Dimensions();
    Matrix<> va(mir.Size(), hdims[0]*inner_dim);
    Matrix<> vb(mir.Size(), inner_dim);
    Matrix<> vda(mir.Size(), hdims[0]*inner_dim);
    Matrix<> vdb(mir.Size(), inner_dim);
    c1->EvaluateDeriv (mir, va, vda);
    c2->EvaluateDeriv (mir, vb, vdb);

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        FlatMatrix<> da(hdims[0], inner_dim, &vda(i,0));

        result.Row(i) = a*vb.Row(i);
        deriv.Row(i) = a*vdb.Row(i) + da*vb.Row(i);
      }
  }
  
  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    FlatArray<int> hdims = Dimensions();    
    Matrix<> va(mir.Size(), hdims[0]*inner_dim);
    Matrix<> vb(mir.Size(), inner_dim);
    Matrix<> vda(mir.Size(), hdims[0]*inner_dim);
    Matrix<> vdb(mir.Size(), inner_dim);
    Matrix<> vdda(mir.Size(), hdims[0]*inner_dim);
    Matrix<> vddb(mir.Size(), inner_dim);
    c1->EvaluateDDeriv (mir, va, vda, vdda);
    c2->EvaluateDDeriv (mir, vb, vdb, vddb);

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        FlatMatrix<> da(hdims[0], inner_dim, &vda(i,0));
        FlatMatrix<> dda(hdims[0], inner_dim, &vdda(i,0));

        result.Row(i) = a*vb.Row(i);
        deriv.Row(i) = a*vdb.Row(i) + da*vb.Row(i);
        dderiv.Row(i) = a*vddb.Row(i) + 2*da*vdb.Row(i) + dda*vb.Row(i);
      }
  }


  



  virtual void Evaluate(const BaseMappedIntegrationRule & mir,
                        FlatArray<FlatMatrix<>*> input,
                        FlatMatrix<> result) const
  {
    FlatArray<int> hdims = Dimensions();    
    FlatMatrix<> va = *input[0], vb = *input[1];
    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        result.Row(i) = a * vb.Row(i);
      }
  }

  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatArray<FlatMatrix<>*> input,
                             FlatArray<FlatMatrix<>*> dinput,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    FlatArray<int> hdims = Dimensions();    
    FlatMatrix<> va = *input[0], vb = *input[1];
    FlatMatrix<> vda = *dinput[0], vdb = *dinput[1];

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        FlatMatrix<> da(hdims[0], inner_dim, &vda(i,0));

        // FlatMatrix<> c(dims[0], dims[1], &result(i,0));
        // FlatMatrix<> dc(dims[0], dims[1], &deriv(i,0));
        // c = a*b;
        // dc = a*db+da*b;
        result.Row(i) = a * vb.Row(i);
        deriv.Row(i) = da * vb.Row(i) + a*vdb.Row(i);        
      }
  }

  virtual void EvaluateDeriv(const SIMD_BaseMappedIntegrationRule & mir,
                             AFlatMatrix<> result,
                             AFlatMatrix<> deriv) const
  {
    FlatArray<int> hdims = Dimensions();

    STACK_ARRAY(SIMD<double>, hmema, mir.Size()*hdims[0]*inner_dim);
    STACK_ARRAY(SIMD<double>, hmemb, mir.Size()*inner_dim);
    AFlatMatrix<double> va(hdims[0]*inner_dim, mir.IR().GetNIP(), &hmema[0]);
    AFlatMatrix<double> vb(inner_dim, mir.IR().GetNIP(), &hmemb[0]);
    STACK_ARRAY(SIMD<double>, hmemda, mir.Size()*hdims[0]*inner_dim);
    STACK_ARRAY(SIMD<double>, hmemdb, mir.Size()*inner_dim);
    AFlatMatrix<double> dva(hdims[0]*inner_dim, mir.IR().GetNIP(), &hmemda[0]);
    AFlatMatrix<double> dvb(inner_dim, mir.IR().GetNIP(), &hmemdb[0]);
    c1->EvaluateDeriv (mir, va, dva);
    c2->EvaluateDeriv (mir, vb, dvb);
    
    // AFlatMatrix<> va = *input[0], vb = *input[1];
    // AFlatMatrix<> dva = *dinput[0], dvb = *dinput[1];

    result = 0.0;
    deriv = 0.0;
    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < inner_dim; k++)
        {
          size_t row = j*inner_dim+k;
          for (size_t i = 0; i < mir.Size(); i++)
            result.Get(j,i) += va.Get(row,i)*vb.Get(k,i);
          for (size_t i = 0; i < mir.Size(); i++)
            deriv.Get(j,i) += dva.Get(row,i)*vb.Get(k,i) + va.Get(row,i)*dvb.Get(k,i);
        }
  }

  

  virtual void EvaluateDeriv(const SIMD_BaseMappedIntegrationRule & mir,
                             FlatArray<AFlatMatrix<>*> input,
                             FlatArray<AFlatMatrix<>*> dinput,
                             AFlatMatrix<> result,
                             AFlatMatrix<> deriv) const
  {
    FlatArray<int> hdims = Dimensions();    
    AFlatMatrix<> va = *input[0], vb = *input[1];
    AFlatMatrix<> dva = *dinput[0], dvb = *dinput[1];

    result = 0.0;
    deriv = 0.0;
    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < inner_dim; k++)
        {
          size_t row = j*inner_dim+k;
          for (size_t i = 0; i < mir.Size(); i++)
            result.Get(j,i) += va.Get(row,i)*vb.Get(k,i);
          for (size_t i = 0; i < mir.Size(); i++)
            deriv.Get(j,i) += dva.Get(row,i)*vb.Get(k,i) + va.Get(row,i)*dvb.Get(k,i);
        }
  }


  
  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatArray<FlatMatrix<>*> ddinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    FlatArray<int> hdims = Dimensions();    
    FlatMatrix<> va = *input[0], vb = *input[1];
    FlatMatrix<> vda = *dinput[0], vdb = *dinput[1];
    FlatMatrix<> vdda = *ddinput[0], vddb = *ddinput[1];

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> a(hdims[0], inner_dim, &va(i,0));
        FlatMatrix<> da(hdims[0], inner_dim, &vda(i,0));
        FlatMatrix<> dda(hdims[0], inner_dim, &vdda(i,0));

        result.Row(i) = a*vb.Row(i);
        deriv.Row(i) = a*vdb.Row(i) + da*vb.Row(i);
        dderiv.Row(i) = a*vddb.Row(i) + 2*da*vdb.Row(i) + dda*vb.Row(i);
      }
  }


  
};



  
class TransposeCoefficientFunction : public T_CoefficientFunction<TransposeCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  // Array<int> dims;
public:
  TransposeCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : T_CoefficientFunction<TransposeCoefficientFunction>(1, ac1->IsComplex()), c1(ac1)
  {
    auto dims_c1 = c1 -> Dimensions();
    if (dims_c1.Size() != 2)
      throw Exception("Transpose of non-matrix called");
    // dims = { dims_c1[1], dims_c1[0] };
    SetDimensions (Array<int> ({ dims_c1[1], dims_c1[0] }));
  }
  
  // virtual bool IsComplex() const { return c1->IsComplex(); }
  // virtual int Dimension() const { return c1->Dimension(); }
  // virtual Array<int> Dimensions() const { return Array<int> (dims); } 

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
      FlatArray<int> hdims = Dimensions();        
      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1]))
          code.body += Var(index,i,j).Assign( Var(inputs[0],j,i) );
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get() } ); }  

  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
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

  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    throw Exception ("TransposeCF:: scalar evaluate for matrix called");
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const
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
                         FlatVector<Complex> result) const
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

  virtual void Evaluate (const BaseMappedIntegrationRule & mir,
                         FlatMatrix<> result) const
  {
    FlatArray<int> hdims = Dimensions();        
    c1->Evaluate (mir, result);
    STACK_ARRAY(double, hmem, hdims[0]*hdims[1]);
    FlatMatrix<> tmp (hdims[0], hdims[1], hmem);

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(hdims[1], hdims[0], &result(i,0));  // source matrix format
        tmp = Trans(reshape);
        FlatMatrix<> reshape2(hdims[0], hdims[1], &result(i,0));  // range matrix format
        reshape2 = tmp;
      }
  }  

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
                         BareSliceMatrix<SIMD<double>> result) const
  {
    FlatArray<int> hdims = Dimensions();    
    c1->Evaluate (mir, result);
    STACK_ARRAY(SIMD<double>, hmem, hdims[0]*hdims[1]);
    AFlatMatrix<double> tmp (hdims[0], hdims[1]*SIMD<double>::Size(), &hmem[0]);

    for (int i = 0; i < mir.Size(); i++)
      {
        for (int j = 0; j < hdims[0]; j++)
          for (int k = 0; k < hdims[1]; k++)
            tmp.Get(j,k) = result(k*hdims[0]+j, i);
        for (int j = 0; j < hdims[0]; j++)
          for (int k = 0; k < hdims[1]; k++)
            result(j*hdims[1]+k, i) = tmp.Get(j,k);
      }
  }  

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & mir,
                   BareSliceMatrix<T> result) const
  {
    FlatArray<int> hdims = Dimensions();    
    c1->Evaluate (mir, result);
    STACK_ARRAY(T, hmem, hdims[0]*hdims[1]);
    FlatMatrix<T> tmp (hdims[0], hdims[1], &hmem[0]);

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

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   FlatArray<BareSliceMatrix<T>> input,                       
                   BareSliceMatrix<T> values) const
  {
    FlatArray<int> hdims = Dimensions();
    size_t np = ir.Size();
    
    auto in0 = input[0];
    for (size_t j = 0; j < hdims[0]; j++)
      for (size_t k = 0; k < hdims[1]; k++)
        for (size_t i = 0; i < np; i++)
          values(j*hdims[1]+k, i) = in0(k*hdims[0]+j, i);
  }
  

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & mir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    FlatArray<int> hdims = Dimensions();        
    auto in0 = *input[0];
    for (int i = 0; i < mir.Size(); i++)
      {
        for (int j = 0; j < hdims[0]; j++)
          for (int k = 0; k < hdims[1]; k++)
            values.Get(j*hdims[1]+k, i) = in0.Get(k*hdims[0]+j, i);
      }
  }

  
  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    FlatArray<int> dims = Dimensions();        
    c1->EvaluateDeriv (mir, result, deriv);
    Matrix<> tmp (dims[0], dims[1]);

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(dims[1], dims[0], &result(i,0));  // source matrix format
        tmp = Trans(reshape);
        FlatMatrix<> reshape2(dims[0], dims[1], &result(i,0));  // range matrix format
        reshape2 = tmp;
      }
    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(dims[1], dims[0], &deriv(i,0));  // source matrix format
        tmp = Trans(reshape);
        FlatMatrix<> reshape2(dims[0], dims[1], &deriv(i,0));  // range matrix format
        reshape2 = tmp;
      }
  }
  
  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    FlatArray<int> dims = Dimensions();            
    c1->EvaluateDDeriv (mir, result, deriv, dderiv);
    Matrix<> tmp (dims[0], dims[1]);

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(dims[1], dims[0], &result(i,0));  // source matrix format
        tmp = Trans(reshape);
        FlatMatrix<> reshape2(dims[0], dims[1], &result(i,0));  // range matrix format
        reshape2 = tmp;
      }
    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(dims[1], dims[0], &deriv(i,0));  // source matrix format
        tmp = Trans(reshape);
        FlatMatrix<> reshape2(dims[0], dims[1], &deriv(i,0));  // range matrix format
        reshape2 = tmp;
      }
    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(dims[1], dims[0], &dderiv(i,0));  // source matrix format
        tmp = Trans(reshape);
        FlatMatrix<> reshape2(dims[0], dims[1], &dderiv(i,0));  // range matrix format
        reshape2 = tmp;
      }
  }


  virtual void Evaluate (const BaseMappedIntegrationRule & mir,
                         FlatArray<FlatMatrix<>*> input,
                         FlatMatrix<> result) const
  {
    FlatArray<int> dims = Dimensions();
    FlatMatrix<> v1 = *input[0];
    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(dims[1], dims[0], &v1(i,0));  // source matrix format
        FlatMatrix<> reshape2(dims[0], dims[1], &result(i,0));  // range matrix format
        reshape2 = Trans (reshape);
      }
  }  
  
  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const
  {
    FlatArray<int> dims = Dimensions();            
    FlatMatrix<> v1 = *input[0];
    FlatMatrix<> dv1 = *dinput[0];

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(dims[1], dims[0], &v1(i,0));  // source matrix format
        FlatMatrix<> reshape2(dims[0], dims[1], &result(i,0));  // range matrix format
        reshape2 = Trans (reshape);
      }
    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(dims[1], dims[0], &dv1(i,0));  // source matrix format
        FlatMatrix<> reshape2(dims[0], dims[1], &deriv(i,0));  // range matrix format
        reshape2 = Trans (reshape);
      }
  }  


  
  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                               FlatArray<FlatMatrix<>*> input,
                               FlatArray<FlatMatrix<>*> dinput,
                               FlatArray<FlatMatrix<>*> ddinput,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const
  {
    FlatArray<int> dims = Dimensions();        
    FlatMatrix<> v1 = *input[0];
    FlatMatrix<> dv1 = *dinput[0];
    FlatMatrix<> ddv1 = *ddinput[0];

    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(dims[1], dims[0], &v1(i,0));  // source matrix format
        FlatMatrix<> reshape2(dims[0], dims[1], &result(i,0));  // range matrix format
        reshape2 = Trans (reshape);
      }
    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(dims[1], dims[0], &dv1(i,0));  // source matrix format
        FlatMatrix<> reshape2(dims[0], dims[1], &deriv(i,0));  // range matrix format
        reshape2 = Trans (reshape);
      }
    for (int i = 0; i < mir.Size(); i++)
      {
        FlatMatrix<> reshape(dims[1], dims[0], &ddv1(i,0));  // source matrix format
        FlatMatrix<> reshape2(dims[0], dims[1], &dderiv(i,0));  // range matrix format
        reshape2 = Trans (reshape);
      }
    
  }

    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & mir,
                                AFlatMatrix<> result,
                                AFlatMatrix<> deriv) const
    {
      FlatArray<int> dims = Dimensions();
      size_t dim0 = dims[0], dim1 = dims[1];
      STACK_ARRAY(SIMD<double>, hmem, dims[0]*dims[1]*mir.Size());
      AFlatMatrix<double> in0 (dims[0]*dims[1], mir.IR().GetNIP(), &hmem[0]);
      STACK_ARRAY(SIMD<double>, hdmem, dims[0]*dims[1]*mir.Size());
      AFlatMatrix<double> din0 (dims[0]*dims[1], mir.IR().GetNIP(), &hdmem[0]);

      c1->EvaluateDeriv (mir, in0, din0);
      size_t s = mir.Size();

      for (size_t j = 0; j < dim0; j++)
        for (size_t k = 0; k < dim1; k++)
          for (size_t i = 0; i < s; i++)
            result.Get(j*dim1+k, i) = in0.Get(k*dim0+j, i);

      for (int j = 0; j < dim0; j++)
        for (int k = 0; k < dim1; k++)
          for (size_t i = 0; i < s; i++)
            deriv.Get(j*dim1+k, i) = din0.Get(k*dim0+j, i);
    }

    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & mir,
                                FlatArray<AFlatMatrix<>*> input,
                                FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result,
                                AFlatMatrix<> deriv) const
    {
      FlatArray<int> dims = Dimensions();
      size_t dim0 = dims[0], dim1 = dims[1];
      auto in0 = *input[0];
      auto din0 = *dinput[0];
      size_t s = mir.Size();

      for (size_t j = 0; j < dim0; j++)
        for (size_t k = 0; k < dim1; k++)
          for (size_t i = 0; i < s; i++)
            result.Get(j*dim1+k, i) = in0.Get(k*dim0+j, i);

      for (int j = 0; j < dim0; j++)
        for (int k = 0; k < dim1; k++)
          for (size_t i = 0; i < s; i++)
            deriv.Get(j*dim1+k, i) = din0.Get(k*dim0+j, i);
    }

    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & mir,
                                 FlatArray<AFlatMatrix<>*> input,
                                 FlatArray<AFlatMatrix<>*> dinput,
                                 FlatArray<AFlatMatrix<>*> ddinput,
                                 AFlatMatrix<> result,
                                 AFlatMatrix<> deriv,
                                 AFlatMatrix<> dderiv) const
    {
      FlatArray<int> dims = Dimensions();              
      auto in0 = *input[0];
      auto din0 = *dinput[0];
      auto ddin0 = *ddinput[0];
      for (size_t i = 0; i < mir.Size(); i++)
        {
          for (int j = 0; j < dims[0]; j++)
            for (int k = 0; k < dims[1]; k++)
              result.Get(j*dims[1]+k, i) = in0.Get(k*dims[0]+j, i);
        }
      for (size_t i = 0; i < mir.Size(); i++)
        {
          for (int j = 0; j < dims[0]; j++)
            for (int k = 0; k < dims[1]; k++)
              deriv.Get(j*dims[1]+k, i) = din0.Get(k*dims[0]+j, i);
        }
      for (size_t i = 0; i < mir.Size(); i++)
        {
          for (int j = 0; j < dims[0]; j++)
            for (int k = 0; k < dims[1]; k++)
              dderiv.Get(j*dims[1]+k, i) = ddin0.Get(k*dims[0]+j, i);
        }
    }
  
  
  };  




  

  
  // ///////////////////////////// operators  /////////////////////////

  struct GenericPlus {
    template <typename T> T operator() (T x, T y) const { return x+y; }
  };
  struct GenericMinus {
    template <typename T> T operator() (T x, T y) const { return x-y; }
  };
  struct GenericMult {
    template <typename T> T operator() (T x, T y) const { return x*y; }
  };
  struct GenericDiv {
    template <typename T> T operator() (T x, T y) const { return x/y; }
  };
  GenericPlus gen_plus;
  GenericMinus gen_minus;
  GenericMult gen_mult;
  GenericDiv gen_div;
  
  shared_ptr<CoefficientFunction> operator+ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2)
  {
    return BinaryOpCF (c1, c2, 
                       gen_plus, // [](double a, double b) { return a+b; },
                       // [](Complex a, Complex b) { return a+b; },
                       // [](double a, double b, double & dda, double & ddb) { dda = 1; ddb = 1; },
                       // [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                       // { ddada = 0; ddadb = 0; ddbdb = 0; },
                       [](bool a, bool b) { return a||b; }, '+'
                       );
  }
  
  shared_ptr<CoefficientFunction> operator- (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2)
  {
    return BinaryOpCF (c1, c2, 
                       gen_minus, // [](double a, double b) { return a-b; },
                       /*
                       [](Complex a, Complex b) { return a-b; },
                       [](double a, double b, double & dda, double & ddb) { dda = 1; ddb = -1; },
                       [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                       { ddada = 0; ddadb = 0; ddbdb = 0; },
                       */
                       [](bool a, bool b) { return a||b; }, '-'
                       );
  }
  shared_ptr<CoefficientFunction> operator* (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2)
  {
    if (c1->Dimensions().Size() == 2 && c2->Dimensions().Size() == 2)
      return make_shared<MultMatMatCoefficientFunction> (c1, c2);
    if (c1->Dimensions().Size() == 2 && c2->Dimensions().Size() == 1)
      return make_shared<MultMatVecCoefficientFunction> (c1, c2);
    if (c1->Dimension() > 1 && c2->Dimension() > 1)
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
      return make_shared<MultScalVecCoefficientFunction> (c1, c2);
    if (c1->Dimension() > 1 && c2->Dimension() == 1)
      return make_shared<MultScalVecCoefficientFunction> (c2, c1);
    
    return BinaryOpCF (c1, c2, 
                       gen_mult, // [](double a, double b) { return a*b; },
                       /*
                       [](Complex a, Complex b) { return a*b; },
                       [](double a, double b, double & dda, double & ddb) { dda = b; ddb = a; },
                       [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                       { ddada = 0; ddadb = 1; ddbdb = 0; },
                       */
                       [](bool a, bool b) { return a&&b; }, '*'
                       );
  }

  shared_ptr<CoefficientFunction> operator* (double v1, shared_ptr<CoefficientFunction> c2)
  {
    return make_shared<ScaleCoefficientFunction> (v1, c2); 
  }
  
  shared_ptr<CoefficientFunction> operator* (Complex v1, shared_ptr<CoefficientFunction> c2)
  {
    return make_shared<ScaleCoefficientFunctionC> (v1, c2); 
  }

  shared_ptr<CoefficientFunction> InnerProduct (shared_ptr<CoefficientFunction> c1,
                                                shared_ptr<CoefficientFunction> c2)
  {
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

  shared_ptr<CoefficientFunction> TransposeCF (shared_ptr<CoefficientFunction> coef)
  {
    return make_shared<TransposeCoefficientFunction> (coef);
  }

  shared_ptr<CoefficientFunction> NormCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsComplex())
      return make_shared<NormCoefficientFunction<Complex>> (coef);
    else
      return make_shared<NormCoefficientFunction<double>> (coef);
  }
  
  
  shared_ptr<CoefficientFunction> operator/ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2)
  {
    return BinaryOpCF (c1, c2,
                       gen_div, // [](double a, double b) { return a/b; },
                       /*
                       [](Complex a, Complex b) { return a/b; },
                       [](double a, double b, double & dda, double & ddb) { dda = 1.0/b; ddb = -a/(b*b); },
                       [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                       { ddada = 0; ddadb = -1.0/(b*b); ddbdb = 2*a/(b*b*b); },
                       */
                       [](bool a, bool b) { return a; }, '/'
                       );
  }







  
class ComponentCoefficientFunction : public T_CoefficientFunction<ComponentCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  int dim1;
  int comp;
  typedef T_CoefficientFunction<ComponentCoefficientFunction> BASE;
public:
  ComponentCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                int acomp)
    : BASE(1, ac1->IsComplex()), c1(ac1), comp(acomp)
  {
    dim1 = c1->Dimension();
  }
  
  // virtual bool IsComplex() const { return c1->IsComplex(); }
  // virtual int Dimension() const { return 1; }
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    auto dims = c1->Dimensions();
    int i,j;
    GetIndex(dims, comp, i, j);
    code.body += Var(index).Assign( Var(inputs[0], i, j ));
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  { return Array<CoefficientFunction*>({ c1.get() }); }

  using BASE::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    VectorMem<20> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    return v1(comp);
  }

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<> result) const
  {
    VectorMem<20> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    result(0) = v1(comp);
  }  

  virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                         FlatVector<Complex> result) const
  {
    Vector<Complex> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    result(0) = v1(comp);
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         FlatMatrix<> result) const
  {
    // int dim1 = c1->Dimension();
    STACK_ARRAY(double, hmem, ir.Size()*dim1);
    FlatMatrix<> temp(ir.Size(), dim1, hmem);
    // Matrix<> m1(ir.Size(), c1->Dimension());
    
    c1->Evaluate (ir, temp);
    result.Col(0) = temp.Col(comp);
  }  

  virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                         FlatMatrix<Complex> result) const
  {
    // int dim1 = c1->Dimension();
    STACK_ARRAY(double, hmem, 2*ir.Size()*dim1);
    FlatMatrix<Complex> temp(ir.Size(), dim1, (Complex*)hmem);
    c1->Evaluate (ir, temp);
    result.Col(0) = temp.Col(comp);
  }  

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<T> values) const
  {
    STACK_ARRAY(T, hmem, ir.Size()*dim1);
    FlatMatrix<T> temp(dim1, ir.Size(), &hmem[0]);
    
    c1->Evaluate (ir, temp);
    size_t nv = ir.Size();
    __assume(nv > 0);
    for (size_t i = 0; i < nv; i++)
      values(0,i) = temp(comp, i);
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   FlatArray<BareSliceMatrix<T>> input,                       
                   BareSliceMatrix<T> values) const
  {
    auto in0 = input[0];    
    values.Row(0).AddSize(ir.Size()) = in0.Row(comp);
  }
  
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    auto in0 = *input[0];    
    values.Row(0) = in0.Row(comp);
  }

  
  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    STACK_ARRAY(double, hmem, mir.Size()*dim1);
    FlatMatrix<> v1(mir.Size(), dim1, hmem);
    STACK_ARRAY(double, hdmem, mir.Size()*dim1);
    FlatMatrix<> dv1(mir.Size(), dim1, hdmem);

    c1->EvaluateDeriv (mir, v1, dv1);
    result.Col(0) = v1.Col(comp);
    deriv.Col(0) = dv1.Col(comp);
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    STACK_ARRAY(double, hmem, mir.Size()*dim1);
    FlatMatrix<> v1(mir.Size(), dim1, hmem);
    STACK_ARRAY(double, hdmem, mir.Size()*dim1);
    FlatMatrix<> dv1(mir.Size(), dim1, hdmem);
    STACK_ARRAY(double, hddmem, mir.Size()*dim1);
    FlatMatrix<> ddv1(mir.Size(), dim1, hddmem);

    c1->EvaluateDDeriv (mir, v1, dv1, ddv1);
    result.Col(0) = v1.Col(comp);
    deriv.Col(0) = dv1.Col(comp);
    dderiv.Col(0) = ddv1.Col(comp);
  }




  
  virtual void Evaluate (const BaseMappedIntegrationRule & mir,
                         FlatArray<FlatMatrix<>*> input,
                         FlatMatrix<> result) const
  {
    FlatMatrix<> v1 = *input[0];
    result.Col(0) = v1.Col(comp);
  }  
  
  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const
  {
    FlatMatrix<> v1 = *input[0];
    FlatMatrix<> dv1 = *dinput[0];
    
    result.Col(0) = v1.Col(comp);
    deriv.Col(0) = dv1.Col(comp);
   }  

  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                               FlatArray<FlatMatrix<>*> input,
                               FlatArray<FlatMatrix<>*> dinput,
                               FlatArray<FlatMatrix<>*> ddinput,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const
  {
    FlatMatrix<> v1 = *input[0];
    FlatMatrix<> dv1 = *dinput[0];
    FlatMatrix<> ddv1 = *ddinput[0];
    
    result.Col(0) = v1.Col(comp);
    deriv.Col(0) = dv1.Col(comp);
    dderiv.Col(0) = ddv1.Col(comp);
   }  


  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & mir, 
                              AFlatMatrix<double> values, AFlatMatrix<double> deriv) const
  {
    STACK_ARRAY(SIMD<double>, hmem, mir.Size()*dim1);
    AFlatMatrix<> v1(dim1, mir.IR().GetNIP(), hmem);
    STACK_ARRAY(SIMD<double>, hdmem, mir.Size()*dim1);
    AFlatMatrix<> dv1(dim1, mir.IR().GetNIP(), hdmem);
    
    c1->EvaluateDeriv (mir, v1, dv1);
    values.Row(0) = v1.Row(comp);
    deriv.Row(0) = dv1.Row(comp);
  }
  

  
  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & mir, 
                               AFlatMatrix<double> values, AFlatMatrix<double> deriv,
                               AFlatMatrix<double> dderiv) const
  {
    STACK_ARRAY(SIMD<double>, hmem, mir.Size()*dim1);
    AFlatMatrix<> v1(dim1, mir.Size(), hmem);
    STACK_ARRAY(SIMD<double>, hdmem, mir.Size()*dim1);
    AFlatMatrix<> dv1(dim1, mir.Size(), hdmem);
    STACK_ARRAY(SIMD<double>, hddmem, mir.Size()*dim1);
    AFlatMatrix<> ddv1(dim1, mir.Size(), hddmem);
    
    c1->EvaluateDDeriv (mir, v1, dv1, ddv1);
    values.Row(0) = v1.Row(comp);
    deriv.Row(0) = dv1.Row(comp);
    dderiv.Row(0) = ddv1.Row(comp);
  }
  
  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                              FlatArray<AFlatMatrix<>*> input,
                              FlatArray<AFlatMatrix<>*> dinput,
                              AFlatMatrix<> result,
                              AFlatMatrix<> deriv) const
  {
    result.Row(0) = input[0] -> Row(comp);
    deriv.Row(0) = dinput[0] -> Row(comp);
  }
  
  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                               FlatArray<AFlatMatrix<>*> input,
                               FlatArray<AFlatMatrix<>*> dinput,
                               FlatArray<AFlatMatrix<>*> ddinput,
                               AFlatMatrix<> result,
                               AFlatMatrix<> deriv,
                               AFlatMatrix<> dderiv) const
  {
    result.Row(0) = input[0] -> Row(comp);
    deriv.Row(0) = dinput[0] -> Row(comp);
    dderiv.Row(0) = ddinput[0] -> Row(comp);
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
  {
    Vector<bool> v1(c1->Dimension()), d1(c1->Dimension()), dd1(c1->Dimension());
    c1->NonZeroPattern (ud, v1, d1, dd1);
    nonzero(0) = v1(comp);
    nonzero_deriv(0) = d1(comp);
    nonzero_dderiv(0) = dd1(comp);
  }  
};

  shared_ptr<CoefficientFunction>
  MakeComponentCoefficientFunction (shared_ptr<CoefficientFunction> c1, int comp)
  {
    return make_shared<ComponentCoefficientFunction> (c1, comp);
  }




// ************************ DomainWiseCoefficientFunction *************************************

class DomainWiseCoefficientFunction : public T_CoefficientFunction<DomainWiseCoefficientFunction>
{
  Array<shared_ptr<CoefficientFunction>> ci;
  typedef T_CoefficientFunction<DomainWiseCoefficientFunction> BASE;
  using BASE::Evaluate;
public:
  DomainWiseCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
    : BASE(1, false), ci(aci) 
  { 
    for (auto & cf : ci)
      if (cf && cf->IsComplex()) is_complex = true;
    for (auto & cf : ci)
      if (cf) SetDimensions(cf->Dimensions());
  }

  /*
  virtual bool IsComplex() const 
  { 
    for (auto cf : ci)
      if (cf && cf->IsComplex()) return true;
    return false;
  }

  virtual int Dimension() const
  {
    for (auto cf : ci)
      if (cf) return cf->Dimension();
    return 0;
  }
  */
  
  virtual bool DefinedOn (const ElementTransformation & trafo)
  {
    int matindex = trafo.GetElementIndex();
    return (matindex < ci.Size() && ci[matindex]);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    code.body += "// DomainWiseCoefficientFunction:\n";
    string type = "decltype(0.0";
    for(int in : inputs)
        type += "+decltype("+Var(in).S()+")()";
    type += ")";
    TraverseDimensions( Dimensions(), [&](int ind, int i, int j) {
        code.body += Var(index,i,j).Declare(type);
    });
    code.body += "switch(domain_index) {\n";
    for(int domain : Range(inputs))
    {
        code.body += "case " + ToLiteral(domain) + ": \n";
        TraverseDimensions( Dimensions(), [&](int ind, int i, int j) {
            code.body += "  "+Var(index, i, j).Assign(Var(inputs[domain], i, j), false);
        });
        code.body += "  break;\n";
    }
    code.body += "default: \n";
    TraverseDimensions( Dimensions(), [&](int ind, int i, int j) {
        code.body += "  "+Var(index, i, j).Assign(string("0.0"), false);
    });
    code.body += "  break;\n";
    code.body += "}\n";
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)   
  {
    for (auto & cf : ci)
      if (cf)
        cf->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  {
    Array<CoefficientFunction*> cfa;
    for (auto cf : ci)
      cfa.Append (cf.get());
    return Array<CoefficientFunction*>(cfa);
  } 
  
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    result = 0;
    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ip, result);
  }


  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ir, values);
    else
      values = 0.0;
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ir, values);
    else
      values = 0.0;
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<T> values) const
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ir, values);
    else
      values.AddSize(Dimension(), ir.Size()) = T(0.0);
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   FlatArray<BareSliceMatrix<T>> input,                       
                   BareSliceMatrix<T> values) const
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      values.AddSize(Dimension(), ir.Size()) = input[matindex];
    else
      values.AddSize(Dimension(), ir.Size()) = T(0.0);
  }  

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    int matindex = ir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      values = *input[matindex];
    else
      values = 0.0;
  }
  
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    result = 0;
    int matindex = ip.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> Evaluate (ip, result);
  }
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1,Complex> res;
    Evaluate (ip, res);
    return res(0);
  }
    
  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    result = 0;
    deriv = 0;

    int matindex = mir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> EvaluateDeriv (mir, result, deriv);
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    result = 0;
    deriv = 0;
    dderiv = 0;

    int matindex = mir.GetTransformation().GetElementIndex();
    if (matindex < ci.Size() && ci[matindex])
      ci[matindex] -> EvaluateDDeriv (mir, result, deriv, dderiv);
  }
};

  shared_ptr<CoefficientFunction>
  MakeDomainWiseCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
  {
    return make_shared<DomainWiseCoefficientFunction> (move (aci));
  }






// ************************ OtherCoefficientFunction *************************************

class OtherCoefficientFunction : public T_CoefficientFunction<OtherCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  typedef T_CoefficientFunction<OtherCoefficientFunction> BASE;
  using BASE::Evaluate;
public:
  OtherCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : BASE(ac1->Dimension(), ac1->IsComplex()), c1(ac1)
  { ; }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    throw Exception ("OtherCF::GenerateCode not available");
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)   
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  {
    Array<CoefficientFunction*> cfa;
    cfa.Append (c1.get());
    return Array<CoefficientFunction*>(cfa);
  } 
  
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    throw Exception ("OtherCF::Evaluated (mip) not available");    
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    throw Exception ("OtherCF::Evaluated (mip) not available");        
  }


  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
    if (!ir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");
    c1->Evaluate (*ir.GetOtherMIR(), values);
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const
  {
    if (!ir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");    
    c1->Evaluate (*ir.GetOtherMIR(), values);    
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<T> values) const
  {
    if (!ir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");    
    c1->Evaluate (*ir.GetOtherMIR(), values);    
  }

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   FlatArray<BareSliceMatrix<T>> input,                       
                   BareSliceMatrix<T> values) const
  {
    if (!ir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");    
    c1->Evaluate (*ir.GetOtherMIR(), values);    
  }

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    // compile not available
    if (!ir.GetOtherMIR()) throw Exception ("other mir not set, pls report to developers");    
    c1->Evaluate (*ir.GetOtherMIR(), values);        
  }
  
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    throw Exception ("OtherCF::Evaluated (mip) not available");        
  }
  
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const
  {
    throw Exception ("OtherCF::Evaluated (mip) not available");            
  }
    
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
};

shared_ptr<CoefficientFunction>
MakeOtherCoefficientFunction (shared_ptr<CoefficientFunction> me)
{
  me->TraverseTree
    ( [&] (CoefficientFunction & nodecf)
      {
        if (dynamic_cast<const ProxyFunction*> (&nodecf))
          throw Exception ("You cannot create an other - CoefficientFunction from a tree involving a ProxyFunction\n  ---> use the Other()-operator on sub-trees");
      }
      );
  return make_shared<OtherCoefficientFunction> (me);
}







  

  // ///////////////////////////// IfPos   ////////////////////////////////  

  
  class IfPosCoefficientFunction : public CoefficientFunction
  {
    shared_ptr<CoefficientFunction> cf_if;
    shared_ptr<CoefficientFunction> cf_then;
    shared_ptr<CoefficientFunction> cf_else;
  public:
    IfPosCoefficientFunction (shared_ptr<CoefficientFunction> acf_if,
                              shared_ptr<CoefficientFunction> acf_then,
                              shared_ptr<CoefficientFunction> acf_else)
      : CoefficientFunction(acf_then->Dimension(),
                            acf_then->IsComplex() || acf_else->IsComplex()),
                            cf_if(acf_if), cf_then(acf_then), cf_else(acf_else)
    {
      SetDimensions(cf_then->Dimensions());
    }

    virtual ~IfPosCoefficientFunction () { ; }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
    {
      if (cf_if->Evaluate(ip) > 0)
        return cf_then->Evaluate(ip);
      else
        return cf_else->Evaluate(ip);      
    }

    virtual void Evaluate (const BaseMappedIntegrationPoint& ip, FlatVector<double> values) const
    {
      if(cf_if->Evaluate(ip) > 0)
        cf_then->Evaluate(ip,values);
      else
        cf_else->Evaluate(ip,values);
    }
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
    {
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

    virtual void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<Complex> values) const
    {
      if(cf_if->Evaluate(ip)>0)
        cf_then->Evaluate(ip,values);
      else
        cf_else->Evaluate(ip,values);
    }


    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
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
                           AFlatMatrix<double> values) const
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

    // virtual bool IsComplex() const { return cf_then->IsComplex() | cf_else->IsComplex(); }
    // virtual int Dimension() const { return cf_then->Dimension(); }

    void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
    {
      auto var_if = Var(inputs[0]);
      TraverseDimensions( cf_then->Dimensions(), [&](int ind, int i, int j) {
          code.body += Var(index,i,j).Declare("decltype("+Var(inputs[1]).S()+")");
      });
      if(code.is_simd) {
        TraverseDimensions( cf_then->Dimensions(), [&](int ind, int i, int j) {
            code.body += Var(index,i,j).Assign("IfPos("+Var(inputs[0]).S()+','+Var(inputs[1],i,j).S()+','+Var(inputs[2],i,j).S()+")", false);
        });
      } else {
        code.body += "if (" + var_if.S() + ">0.0) {\n";
        TraverseDimensions( cf_then->Dimensions(), [&](int ind, int i, int j) {
            code.body += Var(index,i,j).Assign( Var(inputs[1],i,j), false );
        });
        code.body += "} else {\n";
        TraverseDimensions( cf_then->Dimensions(), [&](int ind, int i, int j) {
            code.body += Var(index,i,j).Assign( Var(inputs[2],i,j), false );
        });
        code.body += "}\n";
      }
    }

    /*
    virtual Array<int> Dimensions() const
    {
      return cf_then->Dimensions();
    }
    */
    
    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                FlatMatrix<> values,
                                FlatMatrix<> deriv) const
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
      /*
      *testout << "IfPos::std" << endl
               << "if = " << endl << Trans(if_values)
               << "then = " << endl << Trans(then_values)
               << "then_deriv" << endl << Trans(then_deriv)
               << "else = " << endl << Trans(else_values)
               << "else_deriv" << endl << Trans(else_deriv)
               << "val = " << endl << Trans(values)
               << "deriv = " << endl << Trans(deriv);
      */
    }

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
      /*
      *testout << "IfPos::simd" << endl
               << "ifval = " << endl << if_values
               << "then-val = " << endl << then_values
               << "then-dval = " << endl << then_deriv
               << "else-val = " << endl << else_values
               << "else-dval = " << endl << else_deriv
               << "val = " << endl << values
               << "deriv = " << endl << deriv;
      */
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
    
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
    {
      cf_if->TraverseTree (func);
      cf_then->TraverseTree (func);
      cf_else->TraverseTree (func);
      func(*this);
    }
    
    virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
    {
      return Array<CoefficientFunction*>( { cf_if.get(), cf_then.get(), cf_else.get() } );
    }
    
    virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                                 FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
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
  };
  
  extern
  shared_ptr<CoefficientFunction> IfPos (shared_ptr<CoefficientFunction> cf_if,
                                         shared_ptr<CoefficientFunction> cf_then,
                                         shared_ptr<CoefficientFunction> cf_else)
  {
    return make_shared<IfPosCoefficientFunction> (cf_if, cf_then, cf_else);
  }


class VectorialCoefficientFunction : public T_CoefficientFunction<VectorialCoefficientFunction>
{
  Array<shared_ptr<CoefficientFunction>> ci;
  Array<size_t> dimi;  // dimensions of components
  Array<std::tuple<CoefficientFunction*, size_t>> both;
  typedef T_CoefficientFunction<VectorialCoefficientFunction> BASE;
public:
  VectorialCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
    : BASE(0, false), ci(aci), dimi(aci.Size()), both(aci.Size()+1)
  {
    int hdim = 0;
    for (int i : Range(ci))
      {
        dimi[i] = ci[i]->Dimension();
        both[i] = make_tuple(aci[i].get(), hdim);
        hdim += dimi[i];
      }
    both[ci.Size()] = make_tuple(nullptr, hdim);
    
    for (auto cf : ci)
      if (cf && cf->IsComplex())
        is_complex = true;

    SetDimension(hdim);
    // dims = Array<int> ( { dimension } ); 
  }


  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
  {
    for (auto cf : ci)
      cf->TraverseTree (func);
    func(*this);
  }

  virtual Array<CoefficientFunction*> InputCoefficientFunctions() const
  {
    Array<CoefficientFunction*> cfa;
    for (auto cf : ci)
      cfa.Append (cf.get());
    return Array<CoefficientFunction*>(cfa);
  } 


  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
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

  virtual bool DefinedOn (const ElementTransformation & trafo)
  {
    for (auto & cf : ci)
      if (!cf->DefinedOn(trafo)) return false;
    return true;
  }
  

  using BASE::Evaluate;  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
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
                        FlatVector<Complex> result) const
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

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
    int base = 0;
    for (auto & cf : ci)
      {
        int dimi = cf->Dimension();
        STACK_ARRAY(double, hmem, ir.Size()*dimi);
        FlatMatrix<> temp(ir.Size(), dimi, hmem);
        cf->Evaluate(ir, temp);
        result.Cols(base,base+dimi) = temp;
        base += dimi;
      }
  }

  /*
  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
  {
    int base = 0;
    for (int i : Range(ci))
      {
        ci[i]->Evaluate(ir, values.Rows(base,base+dimi[i]));
        base += dimi[i];
      }
  }
  */

  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<T> values) const
  {
    FlatArray<std::tuple<CoefficientFunction*, size_t>> hboth = both;
    for (size_t i = 0; i < hboth.Size()-1; i++)
      get<0>(hboth[i])->Evaluate(ir, values.Rows(get<1>(hboth[i]), get<1>(hboth[i+1])));
  }
 
  template <typename T>
  void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                   FlatArray<BareSliceMatrix<T>> input,                       
                   BareSliceMatrix<T> values) const
  {
    size_t base = 0;
    size_t np = ir.Size();
    for (size_t i : Range(ci))
      {
        values.Rows(base,base+dimi[i]).AddSize(dimi[i], np) = input[i];
        base += dimi[i];
      }
  }

  virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                         AFlatMatrix<double> values) const
  {
    size_t base = 0;
    for (size_t i : Range(ci))
      {
        values.Rows(base,base+dimi[i]) = *input[i];
        base += dimi[i];
      }
  }

  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                              AFlatMatrix<> result,
                              AFlatMatrix<> deriv) const
  {
    size_t base = 0;
    for (size_t i : Range(ci))
      {
        ci[i]->EvaluateDeriv (ir, result.Rows(base,base+dimi[i]), deriv.Rows(base,base+dimi[i]));
        base += dimi[i];
      }
  }


  
  virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                              FlatArray<AFlatMatrix<>*> input,
                              FlatArray<AFlatMatrix<>*> dinput,
                              AFlatMatrix<> result,
                              AFlatMatrix<> deriv) const
  {
    size_t base = 0;
    for (size_t i : Range(ci))
      {
        result.Rows(base,base+dimi[i]) = *input[i];
        deriv.Rows(base,base+dimi[i]) = *dinput[i];
        base += dimi[i];
      }
  }

  virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                               FlatArray<AFlatMatrix<>*> input,
                               FlatArray<AFlatMatrix<>*> dinput,
                               FlatArray<AFlatMatrix<>*> ddinput,
                               AFlatMatrix<> result,
                               AFlatMatrix<> deriv,
                               AFlatMatrix<> dderiv) const
  {
    size_t base = 0;
    for (size_t i : Range(ci))
      {
        result.Rows(base,base+dimi[i]) = *input[i];
        deriv.Rows(base,base+dimi[i]) = *dinput[i];
        dderiv.Rows(base,base+dimi[i]) = *ddinput[i];
        base += dimi[i];
      }
  }
  
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<Complex> result) const
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        STACK_ARRAY(double, hmem, 2*ir.Size()*dimi);
        FlatMatrix<Complex> temp(ir.Size(), dimi, (Complex*)hmem);
        cf->Evaluate(ir, temp);
        result.Cols(base,base+dimi) = temp;
        base += dimi;
      }
  }


  
  virtual void EvaluateDeriv(const BaseMappedIntegrationRule & mir,
                             FlatMatrix<> result,
                             FlatMatrix<> deriv) const
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        Matrix<> hval(mir.Size(), dimi);
        Matrix<> hderiv(mir.Size(), dimi);
        cf->EvaluateDeriv(mir, hval, hderiv);
        result.Cols(base, base+dimi) = hval;
        deriv.Cols(base, base+dimi) = hderiv;
        base += dimi;
      }
      // ci[i]->EvaluateDeriv(ip, result.Range(i,i+1), deriv.Range(i,i+1));
  }

  virtual void EvaluateDDeriv(const BaseMappedIntegrationRule & mir,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv,
                              FlatMatrix<> dderiv) const
  {
    int base = 0;
    for (auto cf : ci)
      {
        int dimi = cf->Dimension();
        Matrix<> hval(mir.Size(), dimi);
        Matrix<> hderiv(mir.Size(), dimi);
        Matrix<> hdderiv(mir.Size(), dimi);
        cf->EvaluateDDeriv(mir, hval, hderiv, hdderiv);
        result.Cols(base, base+dimi) = hval;
        deriv.Cols(base, base+dimi) = hderiv;
        dderiv.Cols(base, base+dimi) = hdderiv;
        base += dimi;
      }
  }




  virtual void Evaluate (const BaseMappedIntegrationRule & mir,
                         FlatArray<FlatMatrix<>*> input,
                         FlatMatrix<> result) const
  {
    int base = 0;
    for (int i : Range(ci))
      {
        int d = dimi[i];
        result.Cols(base, base+d) = *input[i];
        base += d;
      }
  }
  
  virtual void EvaluateDeriv (const BaseMappedIntegrationRule & mir,
                              FlatArray<FlatMatrix<>*> input,
                              FlatArray<FlatMatrix<>*> dinput,
                              FlatMatrix<> result,
                              FlatMatrix<> deriv) const
  {
    int base = 0;
    for (int i : Range(ci))
      {
        int d = dimi[i];        
        result.Cols(base,base+d) = *input[i];
        deriv.Cols(base, base+d) = *dinput[i];        
        base += d;
      }

    /*
    for (int i : Range(ci))
      {
        int dimi = ci[i]->Dimension();
        result.Cols(base, base+dimi) = *input[i];
        deriv.Cols(base, base+dimi) = *dinput[i];
        base += dimi;
      }
    */
  }

  
  virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & mir,
                               FlatArray<FlatMatrix<>*> input,
                               FlatArray<FlatMatrix<>*> dinput,
                               FlatArray<FlatMatrix<>*> ddinput,
                               FlatMatrix<> result,
                               FlatMatrix<> deriv,
                               FlatMatrix<> dderiv) const
  {
    int base = 0;
    for (int i : Range(ci))
      {
        int dimi = ci[i]->Dimension();
        result.Cols(base, base+dimi) = *input[i];
        deriv.Cols(base, base+dimi) = *dinput[i];
        dderiv.Cols(base, base+dimi) = *ddinput[i];
        base += dimi;
      }
  }

  
};

  void VectorialCoefficientFunction::GenerateCode(Code &code, FlatArray<int> inputs, int index) const
  {
    int input = 0;
    int input_index = 0;
    FlatArray<int> dims = Dimensions();
    TraverseDimensions( dims, [&](int ind, int i, int j) {
	auto cfi = ci[input];
        int i1, j1;
        GetIndex( cfi->Dimensions(), input_index, i1, j1 );
        code.body += Var(index,i,j).Assign( Var(inputs[input], i1, j1) );
        input_index++;
        if (input_index == cfi->Dimension() )
        {
            input++;
            input_index = 0;
        }
    });
  }


  shared_ptr<CoefficientFunction>
  MakeVectorialCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci)
  {
    return make_shared<VectorialCoefficientFunction> (move (aci));
  }


// ////////////////////////// Coordinate CF ////////////////////////

  class CoordCoefficientFunction
    : public T_CoefficientFunction<CoordCoefficientFunction, CoefficientFunctionNoDerivative>
  {
    int dir;
    typedef T_CoefficientFunction<CoordCoefficientFunction, CoefficientFunctionNoDerivative> BASE;
  public:
    CoordCoefficientFunction (int adir) : BASE(1, false), dir(adir) { ; }
    virtual string GetDescription () const
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
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
    {
      if (!ip.IsComplex())
        return ip.GetPoint()(dir);
      else
        return ip.GetPointComplex()(dir).real();
    }
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
			  FlatMatrix<Complex> result) const
    {
      result.Col(0) = ir.GetPoints().Col(dir);
    }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
        auto v = Var(index);
        // code.body += v.Assign(CodeExpr(string("mir.GetPoints()(i,")+ToLiteral(dir)+")"));
        code.body += v.Assign(CodeExpr(string("points(i,")+ToLiteral(dir)+")"));
    }

    template <typename T>
    void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<T> values) const
    {
      auto points = ir.GetPoints();
      size_t nv = ir.Size();
      __assume (nv > 0);
      for (size_t i = 0; i < nv; i++)
        values(i) = points(i, dir);
    }

    template <typename T>
    void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                     FlatArray<BareSliceMatrix<T>> input,                       
                     BareSliceMatrix<T> values) const
    { T_Evaluate (ir, values); }
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    {
      Evaluate (ir, values);
    }
    
  };


shared_ptr<CoefficientFunction> MakeCoordinateCoefficientFunction (int comp)
{
  return make_shared<CoordCoefficientFunction> (comp);
}


  // ///////////////////////////// Compiled CF /////////////////////////
// int myglobalvar;
// int myglobalvar_eval;
  class CompiledCoefficientFunction : public CoefficientFunction
  {
    typedef void (*lib_function)(const ngfem::BaseMappedIntegrationRule &, ngbla::FlatMatrix<double>);
    typedef void (*lib_function_simd)(const ngfem::SIMD_BaseMappedIntegrationRule &, BareSliceMatrix<SIMD<double>>);
    typedef void (*lib_function_deriv)(const ngfem::BaseMappedIntegrationRule &, ngbla::FlatMatrix<double>, ngbla::FlatMatrix<double>);
    typedef void (*lib_function_simd_deriv)(const ngfem::SIMD_BaseMappedIntegrationRule &, BareSliceMatrix<AutoDiff<1,SIMD<double>>>);
    typedef void (*lib_function_dderiv)(const ngfem::BaseMappedIntegrationRule &, ngbla::FlatMatrix<double>, ngbla::FlatMatrix<double>, ngbla::FlatMatrix<double>);
    typedef void (*lib_function_simd_dderiv)(const ngfem::SIMD_BaseMappedIntegrationRule &, BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>>);

    typedef void (*lib_function_complex)(const ngfem::BaseMappedIntegrationRule &, ngbla::FlatMatrix<Complex>);
    typedef void (*lib_function_simd_complex)(const ngfem::SIMD_BaseMappedIntegrationRule &, BareSliceMatrix<SIMD<Complex>>);

    shared_ptr<CoefficientFunction> cf;
    Array<CoefficientFunction*> steps;
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

  public:
    CompiledCoefficientFunction (shared_ptr<CoefficientFunction> acf, bool realcompile, int maxderiv, bool wait )
      : CoefficientFunction(acf->Dimension(), acf->IsComplex()), cf(acf) // , compiled_function(nullptr), compiled_function_simd(nullptr)
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
               Array<CoefficientFunction*> in = stepcf.InputCoefficientFunctions();
               max_inputsize = max2(in.Size(), max_inputsize);
               for (auto incf : in)
                 inputs.Add (mypos, steps.Pos(incf));
             }
         });
      cout << IM(3) << "inputs = " << endl << inputs << endl;

      if(realcompile)
      {
        std::vector<string> link_flags;
        if(cf->IsComplex())
            maxderiv = 0;
        stringstream s;
        string pointer_code;
        string top_code = ""
             "#include<fem.hpp>\n"
             "using namespace ngfem;\n"
             "extern \"C\" {\n"
             ;

        string parameters[3] = {"results", "deriv", "dderiv"};

        for (int deriv : Range(maxderiv+1))
        for (auto simd : {false,true}) {
            cout << IM(3) << "Compiled CF:" << endl;
            Code code;
            code.is_simd = simd;
            code.deriv = deriv;
            for (auto i : Range(steps)) {
              cout << IM(3) << "step " << i << ": " << typeid(*steps[i]).name() << endl;
              steps[i]->GenerateCode(code, inputs[i],i);
            }

            pointer_code += code.pointer;
            top_code += code.top;

            // set results
            string scal_type = cf->IsComplex() ? "Complex" : "double";
            string res_type = scal_type;
            if(simd) res_type = "SIMD<" + res_type + ">";
            if(deriv==1) res_type = "AutoDiff<1," + res_type + ">";
            if(deriv==2) res_type = "AutoDiffDiff<1," + res_type + ">";
            int ii = 0;
            TraverseDimensions( cf->Dimensions(), [&](int ind, int i, int j) {
                 code.body += Var(steps.Size(),i,j).Declare(res_type);
                 code.body += Var(steps.Size(),i,j).Assign(Var(steps.Size()-1,i,j),false);
                 string sget = "(i," + ToLiteral(ii) + ") =";
                 if(simd) sget = "(" + ToLiteral(ii) + ",i) =";

                 for (auto ideriv : Range(simd ? 1 : deriv+1))
                 {
                   code.body += parameters[ideriv] + sget + Var(steps.Size(),i,j).code;
                   if(deriv>=1 && !simd)
                   {
                     code.body += ".";
                     if(ideriv==2) code.body += "D";
                     if(ideriv>=1) code.body += "DValue(0)";
                     else code.body += "Value()";
                   }
                   // if(simd) code.body +=".Data()";
                   code.body += ";\n";
                 }
                 ii++;
            });

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
              {
                s << "(SIMD_BaseMappedIntegrationRule & mir, BareSliceMatrix<" << res_type << "> results";
              }
            else
              {
                string param_type = simd ? "BareSliceMatrix<SIMD<"+scal_type+">> " : "FlatMatrix<"+scal_type+"> ";
                if (simd && deriv == 0) param_type = "BareSliceMatrix<SIMD<"+scal_type+">> ";
                s << "( " << (simd?"SIMD_":"") << "BaseMappedIntegrationRule &mir";
                for(auto i : Range(deriv+1))
                  s << ", " << param_type << parameters[i];
              }
            s << " ) {" << endl;
            s << code.header << endl;
            s << "auto points = mir.GetPoints();" << endl;
            s << "auto domain_index = mir.GetTransformation().GetElementIndex();" << endl;
            s << "for ( auto i : Range(mir)) {" << endl;
            s << "auto & ip = mir[i];" << endl;
            s << code.body << endl;
            s << "}\n}" << endl << endl;

            for(const auto &lib : code.link_flags)
                if(std::find(std::begin(link_flags), std::end(link_flags), lib) == std::end(link_flags))
                    link_flags.push_back(lib);

        }
        s << "}" << endl;
        string file_code = top_code + s.str();
        std::vector<string> codes;
        codes.push_back(file_code);
        if(pointer_code.size()) {
          pointer_code = "extern \"C\" {\n" + pointer_code;
          pointer_code += "}\n";
          codes.push_back(pointer_code);
        }

        auto compile_func = [this, codes, link_flags, maxderiv] () {
              library = CompileCode( codes, link_flags );
              if(cf->IsComplex())
              {
                  compiled_function_simd_complex = library->GetFunction<lib_function_simd_complex>("CompiledEvaluateSIMD");
                  compiled_function_complex = library->GetFunction<lib_function_complex>("CompiledEvaluate");
              }
              else
              {
                  compiled_function_simd = library->GetFunction<lib_function_simd>("CompiledEvaluateSIMD");
                  compiled_function = library->GetFunction<lib_function>("CompiledEvaluate");
                  if(maxderiv>0)
                  {
                      compiled_function_simd_deriv = library->GetFunction<lib_function_simd_deriv>("CompiledEvaluateDerivSIMD");
                      compiled_function_deriv = library->GetFunction<lib_function_deriv>("CompiledEvaluateDeriv");
                  }
                  if(maxderiv>1)
                  {
                      compiled_function_simd_dderiv = library->GetFunction<lib_function_simd_dderiv>("CompiledEvaluateDDerivSIMD");
                      compiled_function_dderiv = library->GetFunction<lib_function_dderiv>("CompiledEvaluateDDeriv");
                  }
              }
              cout << IM(7) << "Compilation done" << endl;
        };
        if(wait)
            compile_func();
        else
        {
          try {
            std::thread( compile_func ).detach();
          } catch (const std::exception &e) {
              cerr << IM(3) << "Compilation of CoefficientFunction failed: " << e.what() << endl;
          }
        }
      }
    }

    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func)
    {
      cf -> TraverseTree (func);
      func(*cf);
    }

    // virtual bool IsComplex() const { return cf->IsComplex(); }
    // virtual int Dimension() const { return cf->Dimension(); }
    // virtual Array<int> Dimensions() const  { return cf->Dimensions(); } 
    
    
    virtual bool ElementwiseConstant () const { return false; }
    virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                                 FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const
    { cf->NonZeroPattern (ud, nonzero, nonzero_deriv, nonzero_dderiv); }

    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
    {
      return cf -> Evaluate(ip);
    }

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<> result) const
    {
      cf->Evaluate (ip, result);      
    }

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<Complex> result) const
    {
      cf->Evaluate (ip, result);
    }

    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
    {
      if(compiled_function)
      {
        compiled_function(ir,values);
        return;
      }

      // static Timer t1("CompiledCF::Evaluate 1");
      // static Timer t2("CompiledCF::Evaluate 2");
      // static Timer t3("CompiledCF::Evaluate 3");

      // t1.Start();
      // int totdim = 0;
      // for (int d : dim) totdim += d;
      ArrayMem<double, 10000> hmem(ir.Size()*totdim);
      int mem_ptr = 0;
      ArrayMem<FlatMatrix<>,100> temp(steps.Size());
      ArrayMem<FlatMatrix<>*, 100> in(max_inputsize);
      for (int i = 0; i < steps.Size(); i++)
        {
          temp[i].AssignMemory(ir.Size(), dim[i], &hmem[mem_ptr]);
          mem_ptr += ir.Size()*dim[i];
        }
      // t1.Stop();
      // t2.Start();
      for (int i = 0; i < steps.Size(); i++)
        {
          // myglobalvar ++;
          // timers[i]->Start();
          auto inputi = inputs[i];
          for (int nr : Range(inputi))
            in[nr] = &temp[inputi[nr]];
          steps[i] -> Evaluate (ir, in.Range(0, inputi.Size()), temp[i]);



          // timers[i]->Stop();
        }
      values = temp.Last();
      // t2.Stop();
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const
    {
      if(compiled_function_simd_deriv)
      {
        compiled_function_simd_deriv(ir, values);
        return;
      }
      
      typedef AutoDiff<1,SIMD<double>> T;
      STACK_ARRAY(T, hmem, ir.Size()*totdim);      
      int mem_ptr = 0;
      ArrayMem<BareSliceMatrix<T>,100> temp(steps.Size());
      ArrayMem<BareSliceMatrix<T>,100> in(max_inputsize);

      for (int i = 0; i < steps.Size(); i++)
        {
          new (&temp[i]) BareSliceMatrix<T> (ir.Size(), &hmem[mem_ptr], DummySize(dim[i], ir.Size()));
          mem_ptr += ir.Size()*dim[i];
        }

      for (int i = 0; i < steps.Size(); i++)
        {
          auto inputi = inputs[i];
          for (int nr : Range(inputi))
            new (&in[nr]) BareSliceMatrix<T> (temp[inputi[nr]]);
          // in[nr] = &temp[inputi[nr]];

          steps[i] -> Evaluate (ir, in.Range(0, inputi.Size()), temp[i]);
        }

      values.AddSize(Dimension(), ir.Size()) = temp.Last();
    }


    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const
    {
      if(compiled_function_simd_dderiv)
      {
        compiled_function_simd_dderiv(ir, values);
        return;
      }
      
      typedef AutoDiffDiff<1,SIMD<double>> T;
      STACK_ARRAY(T, hmem, ir.Size()*totdim);      
      int mem_ptr = 0;
      ArrayMem<BareSliceMatrix<T>,100> temp(steps.Size());
      ArrayMem<BareSliceMatrix<T>,100> in(max_inputsize);

      for (int i = 0; i < steps.Size(); i++)
        {
          new (&temp[i]) BareSliceMatrix<T> (ir.Size(), &hmem[mem_ptr], DummySize(dim[i], ir.Size()));
          mem_ptr += ir.Size()*dim[i];
        }

      for (int i = 0; i < steps.Size(); i++)
        {
          auto inputi = inputs[i];
          for (int nr : Range(inputi))
            new (&in[nr]) BareSliceMatrix<T> (temp[inputi[nr]]);
          // in[nr] = &temp[inputi[nr]];

          steps[i] -> Evaluate (ir, in.Range(0, inputi.Size()), temp[i]);
        }

      values.AddSize(Dimension(), ir.Size()) = temp.Last();
    }

    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
    {
      if(compiled_function_simd)
      {
        compiled_function_simd(ir, values);
        return;
      }

      STACK_ARRAY(SIMD<double>, hmem, ir.Size()*totdim);      
      int mem_ptr = 0;
      ArrayMem<AFlatMatrix<double>,100> temp(steps.Size());
      ArrayMem<AFlatMatrix<double>*,100> in(max_inputsize);

      for (int i = 0; i < steps.Size(); i++)
        {
          new (&temp[i]) AFlatMatrix<double> (dim[i], ir.IR().GetNIP(), &hmem[mem_ptr]);
          mem_ptr += ir.Size()*dim[i];
        }

      for (int i = 0; i < steps.Size(); i++)
        {
          // timers[i]->Start();          
          auto inputi = inputs[i];
          for (int nr : Range(inputi))
            in[nr] = &temp[inputi[nr]];

          steps[i] -> Evaluate (ir, in.Range(0, inputi.Size()), temp[i]);
          // timers[i]->Stop();                    
        }

      BareSliceMatrix<SIMD<double>> temp_last = temp.Last();
      values.AddSize(Dimension(), ir.Size()) = temp_last;
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const
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

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const
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


    
    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                FlatMatrix<double> values, FlatMatrix<double> deriv) const
    {
      if(compiled_function_deriv)
      {
        compiled_function_deriv(ir, values, deriv);
        return;
      }
      /*
      Array<Matrix<>*> temp;
      Array<Matrix<>*> dtemp;
      for (int i = 0; i < steps.Size(); i++)
        {
          // timers[i]->Start();
          
          temp.Append (new Matrix<>(ir.Size(), dim[i]));
          dtemp.Append (new Matrix<>(ir.Size(), dim[i]));
          
          Array<FlatMatrix<>*> in;
          for (int j : inputs[i])
            in.Append (temp[j]);
          Array<FlatMatrix<>*> din;
          for (int j : inputs[i])
            din.Append (dtemp[j]);
          
          steps[i] -> EvaluateDeriv (ir, in, din, *temp[i], *dtemp[i]);
          // timers[i]->Stop();
        }

      values = *temp.Last();
      deriv = *dtemp.Last();
      
      for (int i = 0; i < steps.Size(); i++)
        delete temp[i];
      for (int i = 0; i < steps.Size(); i++)
        delete dtemp[i];
      */

      // int totdim = 0;
      // for (int d : dim) totdim += d;
      ArrayMem<double, 10000> hmem(ir.Size()*3*totdim);
      int mem_ptr = 0;
      
      ArrayMem<FlatMatrix<>,100> temp(steps.Size());
      ArrayMem<FlatMatrix<>,100> dtemp(steps.Size());
      ArrayMem<FlatMatrix<>*, 20> in, din;

      for (int i = 0; i < steps.Size(); i++)
        {
          // timers[i]->Start();
          temp[i].AssignMemory(ir.Size(), dim[i], &hmem[mem_ptr]);
          mem_ptr += ir.Size()*dim[i];
          dtemp[i].AssignMemory(ir.Size(), dim[i], &hmem[mem_ptr]);          
          mem_ptr += ir.Size()*dim[i];

          in.SetSize(0);
          din.SetSize(0);
          for (int j : inputs[i])
            in.Append (&temp[j]);
          for (int j : inputs[i])
            din.Append (&dtemp[j]);
          steps[i] -> EvaluateDeriv (ir, in, din, temp[i], dtemp[i]);
          // timers[i]->Stop();
        }

      values = temp.Last();
      deriv = dtemp.Last();
    }

    virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatMatrix<double> values, FlatMatrix<double> deriv,
                                 FlatMatrix<double> dderiv) const
    {
      if(compiled_function_dderiv)
      {
        compiled_function_dderiv(ir, values, deriv, dderiv);
        return;
      }
      int totdim = 0;
      for (int d : dim) totdim += d;
      ArrayMem<double, 10000> hmem(ir.Size()*3*totdim);
      int mem_ptr = 0;
      
      Array<FlatMatrix<>> temp(steps.Size());
      Array<FlatMatrix<>> dtemp(steps.Size());
      Array<FlatMatrix<>> ddtemp(steps.Size());
      ArrayMem<FlatMatrix<>*, 20> in, din, ddin;

      for (int i = 0; i < steps.Size(); i++)
        {
          // timers[i]->Start();          
          temp[i].AssignMemory(ir.Size(), dim[i], &hmem[mem_ptr]);
          mem_ptr += ir.Size()*dim[i];
          dtemp[i].AssignMemory(ir.Size(), dim[i], &hmem[mem_ptr]);          
          mem_ptr += ir.Size()*dim[i];
          ddtemp[i].AssignMemory(ir.Size(), dim[i], &hmem[mem_ptr]);          
          mem_ptr += ir.Size()*dim[i];

          in.SetSize(0);
          din.SetSize(0);
          ddin.SetSize(0);
          for (int j : inputs[i])
            in.Append (&temp[j]);
          for (int j : inputs[i])
            din.Append (&dtemp[j]);
          for (int j : inputs[i])
            ddin.Append (&ddtemp[j]);

          steps[i] -> EvaluateDDeriv (ir, in, din, ddin, temp[i], dtemp[i], ddtemp[i]);
          // timers[i]->Stop();                    
        }

      values = temp.Last();
      deriv = dtemp.Last();
      dderiv = ddtemp.Last();
    }

    
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                AFlatMatrix<double> values, AFlatMatrix<double> deriv) const
    {
      /*
      if(compiled_function_simd_deriv)
      {
        compiled_function_simd_deriv(ir, values, deriv);
        return;
      }
      */
      // throw ExceptionNOSIMD ("*************** CompiledCF :: EvaluateDeriv not available without codegeneration");


      STACK_ARRAY(SIMD<double>, hmem, 2*ir.Size()*totdim);      
      int mem_ptr = 0;
      ArrayMem<AFlatMatrix<double>,100> temp(steps.Size());
      ArrayMem<AFlatMatrix<double>,100> dtemp(steps.Size());
      ArrayMem<AFlatMatrix<double>*,100> in(max_inputsize), din(max_inputsize);

      for (int i = 0; i < steps.Size(); i++)
        {
          new (&temp[i]) AFlatMatrix<double> (dim[i], ir.IR().GetNIP(), &hmem[mem_ptr]);
          mem_ptr += ir.Size()*dim[i];
          new (&dtemp[i]) AFlatMatrix<double> (dim[i], ir.IR().GetNIP(), &hmem[mem_ptr]);
          mem_ptr += ir.Size()*dim[i];
        }

      for (int i = 0; i < steps.Size(); i++)
        {
          // timers[i]->Start();          
          auto inputi = inputs[i];
          for (int nr : Range(inputi))
            {
              in[nr] = &temp[inputi[nr]];
              din[nr] = &dtemp[inputi[nr]];
            }

          steps[i] -> EvaluateDeriv (ir, in.Range(0, inputi.Size()), din.Range(0, inputi.Size()),
                                     temp[i], dtemp[i]);

          // timers[i]->Stop();                    
        }
      values = temp.Last();
      deriv = dtemp.Last();
    }

    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                 AFlatMatrix<double> values, AFlatMatrix<double> deriv,
                                 AFlatMatrix<double> dderiv) const
    {
      /*
      if(compiled_function_simd_dderiv)
      {
        compiled_function_simd_dderiv(ir, values, deriv, dderiv);
        return;
      }
      */
      // throw ExceptionNOSIMD ("*************** CompiledCF :: EvaluateDDeriv coming soon");
      STACK_ARRAY(SIMD<double>, hmem, 3*ir.Size()*totdim);      
      int mem_ptr = 0;
      ArrayMem<AFlatMatrix<double>,100> temp(steps.Size());
      ArrayMem<AFlatMatrix<double>,100> dtemp(steps.Size());
      ArrayMem<AFlatMatrix<double>,100> ddtemp(steps.Size());
      ArrayMem<AFlatMatrix<double>*,100> in(max_inputsize), din(max_inputsize), ddin(max_inputsize);

      for (int i = 0; i < steps.Size(); i++)
        {
          new (&temp[i]) AFlatMatrix<double> (dim[i], ir.IR().GetNIP(), &hmem[mem_ptr]);
          mem_ptr += ir.Size()*dim[i];
          new (&dtemp[i]) AFlatMatrix<double> (dim[i], ir.IR().GetNIP(), &hmem[mem_ptr]);
          mem_ptr += ir.Size()*dim[i];
          new (&ddtemp[i]) AFlatMatrix<double> (dim[i], ir.IR().GetNIP(), &hmem[mem_ptr]);
          mem_ptr += ir.Size()*dim[i];
        }

      for (int i = 0; i < steps.Size(); i++)
        {
          // timers[i]->Start();          
          auto inputi = inputs[i];
          for (int nr : Range(inputi))
            {
              in[nr] = &temp[inputi[nr]];
              din[nr] = &dtemp[inputi[nr]];
              ddin[nr] = &ddtemp[inputi[nr]];
            }

          steps[i] -> EvaluateDDeriv (ir, in.Range(0, inputi.Size()), din.Range(0, inputi.Size()),
                                      ddin.Range(0, inputi.Size()),
                                      temp[i], dtemp[i], ddtemp[i]);
          // timers[i]->Stop();                    
        }
      values = temp.Last();
      deriv = dtemp.Last();
      dderiv = ddtemp.Last();
    }



    
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const
    {
      return cf->GenerateCode(code, inputs, index);
    }
  };

  class RealCF : public CoefficientFunction
  {
    shared_ptr<CoefficientFunction> cf;
  public:
    RealCF(shared_ptr<CoefficientFunction> _cf) : cf(_cf), CoefficientFunction(1,false)
    { ; }

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
  };

  class ImagCF : public CoefficientFunction
  {
    shared_ptr<CoefficientFunction> cf;
  public:
    ImagCF(shared_ptr<CoefficientFunction> _cf) : cf(_cf), CoefficientFunction(1,false)
    { ; }

    virtual double Evaluate(const BaseMappedIntegrationPoint& ip) const override
    {
      if(cf->IsComplex())
        {
          Vec<1,Complex> val;
          cf->Evaluate(ip,val);
          return val(0).imag();
        }
      throw Exception("real cf has no imag part!");
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

  shared_ptr<CoefficientFunction> Compile (shared_ptr<CoefficientFunction> c, bool realcompile, int maxderiv, bool wait)
  {
    return make_shared<CompiledCoefficientFunction> (c, realcompile, maxderiv, wait);
  }
  
}

