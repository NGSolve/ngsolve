#ifndef FILE_COEFFICIENT
#define FILE_COEFFICIENT

/*********************************************************************/
/* File:   coefficient.hh                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

#include <bla.hpp>
#include "code_generation.hpp"

#include "intrule.hpp"
#include "elementtransformation.hpp"

namespace ngfem
{
  
  /*
    type to determine (non)zero propagation of arithmetic expressions,
    formerly bool was used instead
   */
  class NonZero
  {
    bool nz;
  public:
    constexpr NonZero () : nz(false) { }
    constexpr NonZero (bool _nz) : nz(_nz) { };
    NonZero & operator= (const NonZero &) = default;
    NonZero & operator= (bool _nz) { nz = _nz; return *this; }
    constexpr operator bool() const { return nz; }

    constexpr NonZero operator+ (const NonZero & nz2) const { return nz || nz2.nz; }
    constexpr NonZero operator- (const NonZero & nz2) const { return nz || nz2.nz; }
    constexpr NonZero operator* (const NonZero & nz2) const { return nz && nz2.nz; }
    NonZero & operator+= (const NonZero & nz2) { nz = nz || nz2.nz; return *this; }
    NonZero & operator*= (const NonZero & nz2) { nz = nz && nz2.nz; return *this; }        
  };

  /** 
      coefficient functions
  */
  class NGS_DLL_HEADER CoefficientFunction : public enable_shared_from_this<CoefficientFunction>
  {
  private:
    size_t dimension = 1;
    Array<int> dims;
  protected:
    bool elementwise_constant = false;
    bool is_complex = false;
    int spacedim = -1;  // needed for grad(x), grad(1), ...
    string description;
    bool is_variable = false;  // variables cannot be optimized away (e.g. for differentiation)
  public:
    static std::true_type shallow_archive;
    typedef std::map<shared_ptr<CoefficientFunction>, shared_ptr<CoefficientFunction>> T_DJC; // DiffJacobi Cache type
    // default constructor for archive
    CoefficientFunction() = default;
    CoefficientFunction (int adimension, bool ais_complex = false)
      : is_complex(ais_complex)
    {
      SetDimension(adimension);
    }

    void SetDimension(int adimension)
    {
      dimension = adimension;
      if (dimension <= 1)
        dims = Array<int> (0);
      else
        dims = Array<int> ( { int(dimension) } );
    }
    
    ///
    virtual ~CoefficientFunction ();

    virtual void DoArchive(Archive& ar) { ar & dimension & dims & is_complex; }
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
    ///
    virtual int NumRegions () { return INT_MAX; }
    virtual bool DefinedOn (const ElementTransformation & trafo) { return true; }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const = 0;
    
    ///
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const;
    void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double,ColMajor> values) const
    { Evaluate (ir, Trans(values)); }
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const;
    void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex,ColMajor> values) const
    { Evaluate (ir, Trans(values)); }
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const;

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const;
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatArray<BareSliceMatrix<double,ColMajor>> input,
                           BareSliceMatrix<double,ColMajor> values) const
    {
      Evaluate (ir, Trans(values));
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<SIMD<double>>> input,
                           BareSliceMatrix<SIMD<double>> values) const
    {
      Evaluate (ir, values);
    }
    

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiff<1,double>> values) const
    {
      throw Exception (string("Evaluate AutoDiff<double> not overloaded, type = ")+typeid(*this).name());
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiff<1,double>,ColMajor>> input,
                           BareSliceMatrix<AutoDiff<1,double>,ColMajor> values) const
    {
      Evaluate (ir, Trans(values));
    }

    void Evaluate (const BaseMappedIntegrationRule & ir, 
                   BareSliceMatrix<AutoDiff<1,double>, ColMajor> values) const
    {
      Evaluate (ir, Trans(values));
    }
    

    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const 
    {
      throw ExceptionNOSIMD (string("cf::Evaluate(AutoDiff<simd>) not overloaded for ")+typeid(*this).name());      
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiff<1,SIMD<double>>>> input,
                           BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const
    {
      Evaluate (ir, values);
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiffDiff<1,double>> values) const
    {
      throw Exception (string("Evaluate AutoDiffDiff<double> not overloaded, type = ")+typeid(*this).name());
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiffDiff<1,double>,ColMajor>> input,
                           BareSliceMatrix<AutoDiffDiff<1,double>,ColMajor> values) const
    {
      Evaluate (ir, Trans(values));
    }

    void Evaluate (const BaseMappedIntegrationRule & ir, 
                   BareSliceMatrix<AutoDiffDiff<1,double>, ColMajor> values) const
    {
      Evaluate (ir, Trans(values));
    }
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const 
    {
      throw ExceptionNOSIMD (string("cf::Evaluate(AutoDiffDiff<simd>) not overloaded for ")+typeid(*this).name());      
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>>> input,
                           BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const
    {
      Evaluate (ir, values);
    }

    ///
    virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
    { 
      return Evaluate (ip);
    }

    template <typename SCAL>
    inline SCAL T_Evaluate (const BaseMappedIntegrationPoint & ip) const
    { 
      return SCAL (Evaluate (ip));    // used by PML : AutoDiff<complex>
    }

    virtual double EvaluateConst () const
    {
      throw Exception (string ("EvaluateConst called for non-const coefficient function ")+
		       typeid(*this).name());
    }

    bool IsComplex() const { return is_complex; }
    size_t Dimension() const { return dimension; }
    FlatArray<int> Dimensions() const { return dims; }
    
    void SetDimensions (FlatArray<int> adims)
    {
      dims = adims;
      dimension = 1;
      for (int d : dims) dimension *= d;
    }

    // creates a wrapper with new shape
    shared_ptr<CoefficientFunction> Reshape (FlatArray<int> adims) const;
    shared_ptr<CoefficientFunction> Reshape (int s) const
    { return Reshape( Array<int>( { s } )); }    
    shared_ptr<CoefficientFunction> Reshape (int h, int w) const
    { return Reshape( Array<int>( { h, w } )); }

    shared_ptr<CoefficientFunction> Transpose () const;    
    shared_ptr<CoefficientFunction> TensorTranspose (int i, int j) const;
    
    int SpaceDim () const { return spacedim; } 
    void SetSpaceDim (int adim);
    
    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<> result) const
    {
      double f = Evaluate (ip);
      result(0) = f; 
    }

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<Complex> result) const
    {
      VectorMem<10,double> dres(result.Size());
      Evaluate(ip, dres);
      result = dres;
    }
    
    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                FlatMatrix<Complex> result,
                                FlatMatrix<Complex> deriv) const
    {
      Evaluate (ir, result);
      deriv = 0;
    }

    bool ElementwiseConstant () const { return elementwise_constant; }
    virtual bool IsZeroCF() const;
    // virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const;

    /*
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<bool> nonzero,
                                 FlatVector<bool> nonzero_deriv,
                                 FlatVector<bool> nonzero_dderiv) const;
    */
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<AutoDiffDiff<1,NonZero>> nonzero) const;

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const
    {
      cout << string("nonzero in-out not overloaded for type")+typeid(*this).name() << endl;
      /*
      Vector<bool> nz(values.Size()), nzd(values.Size()), nzdd(values.Size());
      NonZeroPattern (ud, nz, nzd, nzdd);
      for (size_t i = 0; i < values.Size(); i++)
        {
          values(i).Value() = nz(i);
          values(i).DValue(0) = nzd(i);
          values(i).DDValue(0) = nzdd(i);
        }
      */
      NonZeroPattern (ud, values);
    }
    
    virtual void PrintReport (ostream & ost) const;
    virtual void PrintReportRec (ostream & ost, int level) const;
    virtual string GetDescription () const;
    void SetDescription (string desc) { description = desc; }


    bool IsVariable() const { return is_variable; }
    void SetVariable (bool var = true) { is_variable = var; }
    
    virtual shared_ptr<CoefficientFunction>
      Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const;
    // returns Jacobi-matrix (possible as higher order tensor)
    virtual shared_ptr<CoefficientFunction>
      DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const;


    virtual shared_ptr<CoefficientFunction> Operator (const string & name) const;
    virtual shared_ptr<CoefficientFunction> Operator (shared_ptr<class DifferentialOperator> diffop) const;
    
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func);
    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const
    { return Array<shared_ptr<CoefficientFunction>>(); }
    virtual bool StoreUserData() const { return false; }

  };

  
  inline ostream & operator<< (ostream & ost, const CoefficientFunction & cf)
  {
    cf.PrintReport (ost);
    return ost;
  }


  template <>
  inline double CoefficientFunction :: 
  T_Evaluate<double> (const BaseMappedIntegrationPoint & ip) const
  {
    return Evaluate (ip);
  }

  template <>
  inline Complex CoefficientFunction :: 
  T_Evaluate<Complex> (const BaseMappedIntegrationPoint & ip) const
  {
    return EvaluateComplex (ip);
  }
  
  

  /*
  template <int S, int R>
  inline double Evaluate (const CoefficientFunction & fun,
			  const MappedIntegrationPoint<S,R> & ip) 
  { 
    return fun.Evaluate(ip); 
  }
  */
  
  inline double Evaluate (const CoefficientFunction & fun,
			  const BaseMappedIntegrationPoint & ip) 
  { 
    return fun.Evaluate(ip); 
  }


  class NGS_DLL_HEADER CoefficientFunctionNoDerivative : public CoefficientFunction
  {
  public:
    using CoefficientFunction::CoefficientFunction;
    using CoefficientFunction::Evaluate;

    /*
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                AFlatMatrix<double> values, AFlatMatrix<double> deriv) const override
    {
      Evaluate(ir, values);
      deriv = 0.0;
    }

    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                 AFlatMatrix<double> values, AFlatMatrix<double> deriv,
                                 AFlatMatrix<double> dderiv) const override
    {
      Evaluate (ir, values);
      deriv = 0.0;
      dderiv = 0.0;
    }
    */
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<SIMD<double>>> input,
                           BareSliceMatrix<SIMD<double>> values) const override
    { Evaluate (ir, values); }


    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiff<1,double>> values) const override
    {
      SliceMatrix<double> hvalues(ir.Size(), Dimension(), 2*values.Dist(), &values(0).Value());
      Evaluate (ir, hvalues);
      for (size_t i = 0; i < ir.Size(); i++)
        for (size_t j = Dimension(); j-- > 0; )
          values(i,j) = hvalues(i,j);
    }
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const override
    {
      // BareSliceMatrix<SIMD<double>> hvalues(2*values.Dist(), &values(0).Value(), DummySize(Dimension(), ir.Size()));
      BareSliceMatrix<SIMD<double>> hvalues(Dimension(), ir.Size(), 2*values.Dist(), &values(0).Value());
      Evaluate (ir, hvalues);
      for (size_t i = 0; i < Dimension(); i++)
        for (size_t j = ir.Size(); j-- > 0; )
          values(i,j) = hvalues(i,j);
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiff<1,SIMD<double>>>> input,
                           BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const override
    { Evaluate (ir, values); }


    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiffDiff<1,double>> values) const override
    {
      SliceMatrix<double> hvalues(ir.Size(), Dimension(), 3*values.Dist(), &values(0).Value());
      Evaluate (ir, hvalues);
      for (size_t i = 0; i < ir.Size(); i++)
        for (size_t j = Dimension(); j-- > 0; )
          values(i,j) = hvalues(i,j);
    }
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const override
    {
      // BareSliceMatrix<SIMD<double>> hvalues(3*values.Dist(), &values(0).Value(), DummySize(Dimension(), ir.Size()));
      BareSliceMatrix<SIMD<double>> hvalues(Dimension(), ir.Size(), 3*values.Dist(), &values(0).Value());      
      Evaluate (ir, hvalues);
      for (size_t i = 0; i < Dimension(); i++)
        for (size_t j = ir.Size(); j-- > 0; )
          values(i,j) = hvalues(i,j);
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>>> input,
                           BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const override
    { Evaluate (ir, values); }

    /*
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                FlatArray<AFlatMatrix<>*> input,
                                FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result,
                                AFlatMatrix<> deriv) const
    {
      Evaluate (ir, input, result);
      deriv = 0.0;
    }

    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                 FlatArray<AFlatMatrix<>*> input,
                                 FlatArray<AFlatMatrix<>*> dinput,
                                 FlatArray<AFlatMatrix<>*> ddinput,
                                 AFlatMatrix<> result,
                                 AFlatMatrix<> deriv,
                                 AFlatMatrix<> dderiv) const
    {
      Evaluate (ir, input, result);
      deriv = 0.0;
      dderiv = 0.0;
    }
    */

    /*
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<bool> nonzero,
                                 FlatVector<bool> nonzero_deriv,
                                 FlatVector<bool> nonzero_dderiv) const override
    {
      nonzero = true;
      nonzero_deriv = false;
      nonzero_dderiv = false;
    }
    */

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      values = AutoDiffDiff<1,NonZero> (true);
    }

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      values = AutoDiffDiff<1,NonZero> (true);
    }
    
    virtual shared_ptr<CoefficientFunction>
      Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override;

    virtual shared_ptr<CoefficientFunction>
      DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override;

  };


  

  template <typename TCF, typename BASE = CoefficientFunction>
  class T_CoefficientFunction : public BASE
  {
  public:
    using BASE::IsComplex;
    using BASE::Dimension;
    using BASE::BASE;
      
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <SIMD<double>> */ (ir, values); }
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const override
    {
      if (IsComplex())
        static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <SIMD<Complex>> */ (ir, values);
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

    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      STACK_ARRAY(double, hmem, Dimension());
      FlatMatrix<double,ColMajor> mat(Dimension(), 1, hmem);
      ip.IntegrationRuleFromPoint([this,mat] (const BaseMappedIntegrationRule & ir)
                                  { static_cast<const TCF*>(this)->T_Evaluate (ir, BareSliceMatrix<double,ColMajor>(mat)); });
      return mat(0);
    }


    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<> result) const override
    {
      FlatMatrix<double,ColMajor> mat(Dimension(), 1, &result(0));
      ip.IntegrationRuleFromPoint([this,mat] (const BaseMappedIntegrationRule & ir)
                                  { static_cast<const TCF*>(this)->T_Evaluate (ir, BareSliceMatrix<double,ColMajor>(mat)); });
    }

    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<Complex> result) const override
    {
      FlatMatrix<Complex,ColMajor> mat(Dimension(), 1, &result(0));
      ip.IntegrationRuleFromPoint([this,mat] (const BaseMappedIntegrationRule & ir)
                                  { static_cast<const TCF*>(this)->T_Evaluate (ir, BareSliceMatrix<Complex,ColMajor>(mat)); });
    }
    
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<SIMD<double>>> input,
                           BareSliceMatrix<SIMD<double>> values) const override
    {  static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <SIMD<double>> */ (ir, input, values); }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const override
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate (ir, Trans(values)); }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<double,ColMajor>> input,
                           BareSliceMatrix<double,ColMajor> values) const override
    { static_cast<const TCF*>(this) -> T_Evaluate (ir, input, values); }
    
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const override
    {
      if (!IsComplex())
        {
          /*
          BareSliceMatrix<double> realvalues(2*values.Dist(), (double*)values.Data(),
                                             DummySize(values.Height(), values.Width()));
          */
          BareSliceMatrix<double> realvalues(values.Height(), values.Width(), 2*values.Dist(), (double*)values.Data());
          Evaluate (ir, realvalues);
          for (size_t i = 0; i < ir.Size(); i++)
            for (size_t j = Dimension(); j-- > 0; )
              values(i,j) = realvalues(i,j);
          return;
        } 
      static_cast<const TCF*>(this) -> /* template */ T_Evaluate (ir, Trans(values));
    }
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiff<1,double>> values) const override
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate (ir, Trans(values)); }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiff<1,double>,ColMajor>> input,
                           BareSliceMatrix<AutoDiff<1,double>,ColMajor> values) const override
    { static_cast<const TCF*>(this) -> T_Evaluate (ir, input, values); }
        
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const override
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <AutoDiff<1,SIMD<double>>> */ (ir, values); }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiff<1,SIMD<double>>>> input,
                           BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const override
    {  static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <AutoDiff<1,SIMD<double>>> */ (ir, input, values); }


    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiffDiff<1,double>> values) const override
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate (ir, Trans(values)); }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiffDiff<1,double>,ColMajor>> input,
                           BareSliceMatrix<AutoDiffDiff<1,double>,ColMajor> values) const override
    { static_cast<const TCF*>(this) -> T_Evaluate (ir, input, values); }
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const override
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <AutoDiffDiff<1,SIMD<double>>> */ (ir, values); }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>>> input,
                           BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const override
    {  static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <AutoDiffDiff<1,SIMD<double>>> */ (ir, input, values); }
  };


  /// The coefficient is constant everywhere
  class NGS_DLL_HEADER ConstantCoefficientFunction
    : public T_CoefficientFunction<ConstantCoefficientFunction, CoefficientFunctionNoDerivative>
  {
    ///
    double val;
    typedef T_CoefficientFunction<ConstantCoefficientFunction, CoefficientFunctionNoDerivative> BASE;
    using BASE::T_DJC;
  public:
    ///
    // ConstantCoefficientFunction() = default;
    ConstantCoefficientFunction (double aval);
    ///
    virtual ~ConstantCoefficientFunction ();
    ///

    void DoArchive (Archive & archive) override
    {
      /*
      BASE::DoArchive(archive);
      archive & val;
      */
    }

    auto GetCArgs() const { return tuple { val }; }
    
    using BASE::Evaluate;
    // virtual bool ElementwiseConstant () const override { return true; }
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      return val;
    }

    virtual double EvaluateConst () const override
    {
      return val;
    }
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const override;
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const override;

    template <typename MIR, typename T, ORDERING ORD>
      void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
    {
      size_t np = ir.Size();    
      __assume (np > 0);
      for (size_t i = 0; i < np; i++)
        values(0,i) = val;
    }
      
    template <typename MIR, typename T, ORDERING ORD>
      void T_Evaluate (const MIR & ir,
                       FlatArray<BareSliceMatrix<T,ORD>> input,                       
                       BareSliceMatrix<T,ORD> values) const
    { T_Evaluate (ir, values); }

    virtual void PrintReport (ostream & ost) const override;
    virtual string GetDescription () const override;
    
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override; 
    
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      values = AutoDiffDiff<1,NonZero> (val != 0.0);
    }

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      values = AutoDiffDiff<1,NonZero> (val != 0.0);
    }

    virtual shared_ptr<CoefficientFunction>
      DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override;
  };



  /// The coefficient is constant everywhere
  template<typename SCAL>
  class NGS_DLL_HEADER ParameterCoefficientFunction : public CoefficientFunctionNoDerivative
  {
    ///
    SCAL val;
  public:
    ///
    // ParameterCoefficientFunction() = default;
    ParameterCoefficientFunction (SCAL aval);
    ///
    virtual ~ParameterCoefficientFunction ();
    ///
    void DoArchive (Archive& ar) override;
    auto GetCArgs() const { return tuple { val }; }
    
    using CoefficientFunction::Evaluate;
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override;
    void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const override;
    void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override;
    void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const override;
    void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const override;

    virtual void SetValue (SCAL in) { val = in; }
    virtual SCAL GetValue () { return val; }
    void PrintReport (ostream & ost) const override;
    void GenerateCode (Code &code, FlatArray<int> inputs, int index) const override;
  };

  class NGS_DLL_HEADER PlaceholderCoefficientFunction : public CoefficientFunction
  {
    shared_ptr<CoefficientFunction> cf;
  public:
    PlaceholderCoefficientFunction() = default;
    PlaceholderCoefficientFunction(shared_ptr<CoefficientFunction> _cf)
      : CoefficientFunction(_cf->Dimension(), _cf->IsComplex()), cf(_cf)
      { SetDimensions(cf->Dimensions()); }

    void DoArchive(Archive& ar) override;

    void Set(shared_ptr<CoefficientFunction> _cf);

    double Evaluate(const BaseMappedIntegrationPoint& ip) const override
    { return cf->Evaluate(ip); }
    void Evaluate(const BaseMappedIntegrationRule & ir,
                  BareSliceMatrix<double> values) const override
    { cf->Evaluate(ir, values); }
    void Evaluate(const SIMD_BaseMappedIntegrationRule & ir,
                  BareSliceMatrix<SIMD<double>> values) const override
    { cf->Evaluate(ir, values); }
    void Evaluate(const SIMD_BaseMappedIntegrationRule & ir,
                  BareSliceMatrix<SIMD<Complex>> values) const override
    { cf->Evaluate(ir, values); }
    void Evaluate(const BaseMappedIntegrationRule & ir,
                  BareSliceMatrix<Complex> values) const override
    { cf->Evaluate(ir, values); }
    void Evaluate(const SIMD_BaseMappedIntegrationRule & ir,
                  FlatArray<BareSliceMatrix<SIMD<double>>> input,
                  BareSliceMatrix<SIMD<double>> values) const override
    { cf->Evaluate (ir, values); }
    void Evaluate(const BaseMappedIntegrationRule & ir,
                  BareSliceMatrix<AutoDiff<1,double>> values) const override
    { cf->Evaluate(ir, values); }
    void Evaluate(const SIMD_BaseMappedIntegrationRule & ir,
                  BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const override
    { cf->Evaluate(ir, values); }
    void Evaluate(const SIMD_BaseMappedIntegrationRule & ir,
                  FlatArray<BareSliceMatrix<AutoDiff<1,SIMD<double>>>> input,
                  BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const override
    { cf->Evaluate (ir, values); }
    void Evaluate(const BaseMappedIntegrationRule & ir,
                  BareSliceMatrix<AutoDiffDiff<1,double>> values) const override
    { cf->Evaluate(ir, values); }
    void Evaluate(const SIMD_BaseMappedIntegrationRule & ir,
                  BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const override
    { cf->Evaluate(ir, values); }
    void Evaluate(const SIMD_BaseMappedIntegrationRule & ir,
                  FlatArray<BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>>> input,
                  BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const override
    { cf->Evaluate(ir, input, values); }

    Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const override
    { return cf->EvaluateComplex(ip); }

    double EvaluateConst () const override
    { return cf->EvaluateConst(); }

    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<> result) const override
    { cf->Evaluate(ip, result); }
    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<Complex> result) const override
    { cf->Evaluate(ip, result); }
    void EvaluateDeriv(const BaseMappedIntegrationRule & ir,
                       FlatMatrix<Complex> result,
                       FlatMatrix<Complex> deriv) const override
    { cf->EvaluateDeriv(ir, result, deriv); }

    void NonZeroPattern (const class ProxyUserData & ud,
                         FlatVector<AutoDiffDiff<1,NonZero>> nonzero) const override
    { cf->NonZeroPattern(ud, nonzero); }

    void NonZeroPattern (const class ProxyUserData & ud,
                         FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                         FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    { cf->NonZeroPattern(ud, input, values); }

    void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      cf->TraverseTree(func);
      func(*this);
    }
    Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return { cf }; }
  };
  


#ifdef OLD
  ///
  // template <int DIM>
  class NGS_DLL_HEADER DomainVariableCoefficientFunction : public CoefficientFunction
  {
    Array<shared_ptr<EvalFunction>> fun;
    Array<shared_ptr<CoefficientFunction>> depends_on;
    int numarg;
  public:
    ///
    DomainVariableCoefficientFunction (const EvalFunction & afun);
    DomainVariableCoefficientFunction (const EvalFunction & afun,
				       const Array<shared_ptr<CoefficientFunction>> & adepends_on);
    DomainVariableCoefficientFunction (const Array<shared_ptr<EvalFunction>> & afun);
    DomainVariableCoefficientFunction (const Array<shared_ptr<EvalFunction>> & afun,
				       const Array<shared_ptr<CoefficientFunction>> & adepends_on);

    ///
    virtual ~DomainVariableCoefficientFunction ();
    ///
    virtual int NumRegions () { return (fun.Size() == 1) ? INT_MAX : fun.Size(); }
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;

    virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const;

    EvalFunction & GetEvalFunction(const int index)
    {
      return *(fun[index]);
    }

    virtual bool IsComplex() const;
    /*
    {
      for (int i = 0; i < fun.Size(); i++)
	if (fun[i]->IsResultComplex()) return true;
      return false;
    }
    */
    virtual int Dimension() const; 
    // { return fun[0]->Dimension(); }

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<> result) const;

    virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
			  FlatVector<Complex> result) const;

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
			   BareSliceMatrix<double> values) const;

    virtual void PrintReport (ostream & ost) const;

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
  };
#endif


#ifdef OLD
  ///
  template <int DIM>
  class NGS_DLL_HEADER DomainInternalCoefficientFunction : public CoefficientFunction
  {
    ///
    int matnr;
    ///
    double (*f)(const double*); 
  public:
    ///
    DomainInternalCoefficientFunction (int amatnr, double (*af)(const double*))
      : matnr(amatnr), f(af) { ; }
    ///
    virtual ~DomainInternalCoefficientFunction () { ; }
    ///
    /*
      template <int S, int R>
      double Evaluate (const MappedIntegrationPoint<S,R> & ip)
      {
      int elind = ip.GetTransformation().GetElementIndex();
      if (elind != matnr && matnr > 0) return 0;
  
      return f(&ip.GetPoint()(0));
      }
    */
  
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
    {
      int elind = ip.GetTransformation().GetElementIndex();
      if (elind != matnr && matnr > 0) return 0;

      return f(&static_cast<const DimMappedIntegrationPoint<DIM>&>(ip).GetPoint()(0));
      // return f(&ip.GetPoint().REval(0));
    }
  
  };
#endif


#ifdef OLD

  /**
     coefficient function that is defined in every integration point.
     NOTE: for the constructor, the maximal number of integration
     points per element is required!
  **/
  class IntegrationPointCoefficientFunction : public CoefficientFunction
  {
    int elems, ips_per_elem;
    ///
    Array<double> values;
  public:
    IntegrationPointCoefficientFunction() = default;
    ///
    IntegrationPointCoefficientFunction (int aelems, int size)
      : CoefficientFunction(1, false), elems(aelems), ips_per_elem(size), values(aelems*size) { ; }
    ///
    IntegrationPointCoefficientFunction (int aelems, int size, double val)
      : CoefficientFunction(1, false), elems(aelems), ips_per_elem(size), values(aelems*size)
    {
      values = val;
    } 
    ///
    IntegrationPointCoefficientFunction (int aelems, int size, Array<double> & avalues)
      : CoefficientFunction(1, false), elems(aelems), ips_per_elem(size), values(avalues) 
    { 
      if ( avalues.Size() < aelems * size )
	{
	  cout << "Warning: IntegrationPointCoefficientFunction, constructor: sizes don't match!" << endl;
	  values.SetSize(aelems*size);
	}
    }

    ///
    virtual ~IntegrationPointCoefficientFunction () { ; }

    void DoArchive(Archive& ar) override
    {
      ar & elems & ips_per_elem & values;
    }
    ///
    /*
      template <int S, int R>
      double Evaluate (const MappedIntegrationPoint<S,R> & ip)
      {
      int ipnr = ip.GetIPNr();
      int elnr = ip.GetTransformation().GetElementNr();

      if ( ipnr < 0 || ipnr >= ips_per_elem || elnr < 0 || elnr >= elems ) 
      {
      ostringstream ost;
      ost << "IntegrationPointCoefficientFunction: ip = "
      << ipnr << " / elem = " << elnr << ". Ranges: 0 - " 
      << ips_per_elem << "/ 0 - " << elems << "!" << endl;
      throw Exception (ost.str());
      }
  
      return values[elnr*ips_per_elem+ipnr];
      }
    */
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      int ipnr = ip.GetIPNr();
      int elnr = ip.GetTransformation().GetElementNr();

      if ( ipnr < 0 || ipnr >= ips_per_elem || elnr < 0 || elnr >= elems ) 
	{
	  ostringstream ost;
	  ost << "IntegrationPointCoefficientFunction: ip = "
	      << ipnr << " / elem = " << elnr << ". Ranges: 0 - " 
	      << ips_per_elem << "/ 0 - " << elems << "!" << endl;
	  throw Exception (ost.str());
	}
  
      return values[elnr*ips_per_elem+ipnr];
    }


    // direct access to the values at the integration points
    double & operator() (int elnr, int ipnr)
    {
      if ( ipnr < 0 || ipnr >= ips_per_elem || elnr < 0 || elnr >= elems ) 
	{
	  ostringstream ost;
	  ost << "IntegrationPointCoefficientFunction: ip = "
	      << ipnr << " / elem = " << elnr << ". Ranges: 0 - " 
	      << ips_per_elem-1 << " / 0 - " << elems-1 << "!" << endl;
	  throw Exception (ost.str());
	}
  
      return values[elnr*ips_per_elem+ipnr];
    }

    double operator() (int elnr, int ipnr) const
    {
      if ( ipnr < 0 || ipnr >= ips_per_elem || elnr < 0 || elnr >= elems ) 
	{
	  ostringstream ost;
	  ost << "IntegrationPointCoefficientFunction: ip = "
	      << ipnr << " / elem = " << elnr << ". Ranges: 0 - " 
	      << ips_per_elem-1 << " / 0 - " << elems-1 << "!" << endl;
	  throw Exception (ost.str());
	}
  
      return values[elnr*ips_per_elem+ipnr];
    }



    int GetNumIPs() const { return ips_per_elem; }
    int GetNumElems() const { return elems; }



    void ReSetValues( Array<double> & avalues )
    {
      if ( avalues.Size() < values.Size() )
	{
	  throw Exception("IntegrationPointCoefficientFunction::ReSetValues - sizes don't match!");
	}
      values = avalues;
    }
  
  };
#endif


#ifdef OLD
  /// Coefficient function that depends (piecewise polynomially) on a parameter
  class PolynomialCoefficientFunction : public CoefficientFunction
  {
  private:
    Array < Array< Array<double>* >* > polycoeffs;
    Array < Array<double>* > polybounds;

  private:
    double EvalPoly(const double t, const Array<double> & coeffs) const;
    double EvalPolyDeri(const double t, const Array<double> & coeffs) const;

  public:
    PolynomialCoefficientFunction(const Array < Array<double>* > & polycoeffs_in);
    PolynomialCoefficientFunction(const Array < Array< Array<double>* >* > & polycoeffs_in, const Array < Array<double>* > & polybounds_in);
  
    virtual ~PolynomialCoefficientFunction();
    using CoefficientFunction::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const; 

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip, const double & t) const;

    virtual double EvaluateDeri (const BaseMappedIntegrationPoint & ip, const double & t) const;

    virtual double EvaluateConst () const;
  };
#endif


  //////////////////
#ifdef OLD
  class FileCoefficientFunction : public CoefficientFunction
  {
  private:
    Array < Array < double > * > ValuesAtIps;

    ofstream outfile;

    string valuesfilename;
    string infofilename;
    string ipfilename;

    int maxelnum, maxipnum, totalipnum;

    bool writeips;

  private:
    void EmptyValues(void);

  public:
    FileCoefficientFunction();

    FileCoefficientFunction(const string & filename);

    FileCoefficientFunction(const string & aipfilename,
			    const string & ainfofilename,
			    const string & avaluesfilename,
			    const bool loadvalues = false);

    virtual ~FileCoefficientFunction();

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;

    void LoadValues(const string & filename);
    inline void LoadValues(void){ LoadValues(valuesfilename); }

    void StartWriteIps(const string & filename);
    inline void StartWriteIps(void){ StartWriteIps(ipfilename); }
  
    void StopWriteIps(const string & infofilename);
    inline void StopWriteIps(void){ StopWriteIps(infofilename); }

    void Reset(void);

  };
#endif













  // *************************** CoefficientFunction Algebra ********************************
  template <typename OP>
  class cl_UnaryOpCF : public T_CoefficientFunction<cl_UnaryOpCF<OP>>
{
  shared_ptr<CoefficientFunction> c1;
  OP lam;
  string name;
  typedef  T_CoefficientFunction<cl_UnaryOpCF<OP>> BASE;
  using typename BASE::T_DJC;
public:
  // cl_UnaryOpCF() = default;
  cl_UnaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                OP alam, string aname="undefined")
    : BASE(ac1->Dimension(),
           ac1->IsComplex() && typeid (lam(Complex(0.0))) == typeid(Complex)),
      c1(ac1), lam(alam), name(aname)
  {
    this->SetDimensions (c1->Dimensions());
    this->elementwise_constant = c1->ElementwiseConstant();
    this->SetDescription(string("unary operation '")+name+"'");
  }

  virtual void DoArchive (Archive & archive) override
  {
    BASE::DoArchive(archive);  // for reshape
    // archive.Shallow(c1) & name & lam;
  }
  auto GetCArgs() const { return tuple { c1, lam, name }; }
  /*
  virtual string GetDescription () const override
  {
    return string("unary operation '")+name+"'";
  }
  */
  virtual bool DefinedOn (const ElementTransformation & trafo) override
  { return c1->DefinedOn(trafo); } 

  // virtual bool ElementwiseConstant () const override { return c1->ElementwiseConstant(); }
  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    //code.Declare (code.res_type, index, this->Dimensions());
    code.Declare (index, this->Dimensions(), this->IsComplex());
    if (code_uses_tensors)
      {
        code.body += "for (size_t i = 0; i < "+ToString(this->Dimension())+"; i++)\n";
        code.body += "var_"+ToString(index)+"[i] = "+name+"( var_"+ToString(inputs[0])+"[i]);\n";
      }
    else
      for (int i = 0; i < this->Dimension(); i++)
        code.body += Var(index, i, this->Dimensions())
          .Assign( Var(inputs[0], i, c1->Dimensions()).Func(name), false);
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 }); }

  using BASE::Evaluate;
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    return lam (c1->Evaluate(ip));
  }

  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const override
  {
    return lam (c1->EvaluateComplex(ip));
  }

  virtual double EvaluateConst () const override
  {
    return lam (c1->EvaluateConst());
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const override
  {
    c1->Evaluate (ip, result);
    for (int j = 0; j < result.Size(); j++)
      result(j) = lam(result(j));
  }

  virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> result) const override
  {
    c1->Evaluate (ir, result);
    /*
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lam (result(i));
    */
    size_t np = ir.Size();
    size_t dim = this->Dimension();
    for (size_t i = 0; i < np; i++)
      for (size_t j = 0; j < dim; j++)
        result(i,j) = lam (result(i,j));
  }
  
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const override
  {
    c1->Evaluate (ip, result);
    for (int j = 0; j < result.Size(); j++)
      result(j) = lam(result(j));
  }
  
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        BareSliceMatrix<Complex> result) const override
  {
    c1->Evaluate (ir, result);
    // for (int i = 0; i < result.Height()*result.Width(); i++)
    // result(i) = lam(result(i));
    size_t np = ir.Size();
    size_t dim = this->Dimension();
    for (size_t i = 0; i < np; i++)
      for (size_t j = 0; j < dim; j++)
        result(i,j) = lam (result(i,j));
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    c1->Evaluate (ir, values);
    size_t dim = this->Dimension();
    size_t np = ir.Size();
    for (size_t i = 0; i < dim; i++)
      for (size_t j = 0; j < np; j++)
        values(i,j) = lam (values(i,j));
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    auto in0 = input[0];
    size_t dim = this->Dimension();
    size_t np = ir.Size();
    for (size_t i = 0; i < dim; i++)
      for (size_t j = 0; j < np; j++)
        values(i,j) = lam (in0(i,j));
  }


  virtual shared_ptr<CoefficientFunction>
  Operator (const string & name) const override
  { throw Exception ("unarycf "+name+" does not provide Operator"); }

  
  virtual shared_ptr<CoefficientFunction>
  Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override
  { throw Exception ("unarycf "+name+" does not provide a derivative"); }

  virtual shared_ptr<CoefficientFunction>
  DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  { return BASE::DiffJacobi(var, cache); }

  
  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv,
                               FlatVector<bool> nonzero_dderiv) const override
  {
    size_t dim = this->Dimension();    
    Vector<bool> v1(dim), d1(dim), dd1(dim);
    c1->NonZeroPattern(ud, v1, d1, dd1);
    for (int i = 0; i < nonzero.Size(); i++)
      {
        if (name == "-" || name == " ") // "-" actually not used that way
          {
            nonzero(i) = v1(i);
            nonzero_deriv(i) = d1(i);
            nonzero_dderiv(i) = dd1(i);
          }
        else
          {
            nonzero(i) = v1(i);
            nonzero_deriv(i) = d1(i);
            nonzero_dderiv(i) = d1(i) || dd1(i);
          }
      }
  }
  */

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,NonZero>> values) const override                               
  {
    size_t dim = this->Dimension();    
    Vector<AutoDiffDiff<1,NonZero>> v1(dim);
    c1->NonZeroPattern(ud, v1);
    for (int i = 0; i < values.Size(); i++)
      {
        if (name == "-" || name == " ") // "-" actually not used that way
          {
            values = v1;
          }
        else
          {
            for (size_t i = 0; i < values.Size(); i++)
              {
                values[i].Value() = v1[i].Value();
                values[i].DValue(0) = v1[i].DValue(0);
                values[i].DDValue(0) = v1[i].DValue(0) || v1[i].DDValue(0);
              }
          }
      }
  }

  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                               FlatVector<AutoDiffDiff<1,NonZero>> values) const override
  {
    auto v1 = input[0];
    if (name == "-"  || name == " ") // "-" actually not used that way
      {
        values = v1;
      }
    else
      {
        for (size_t i = 0; i < values.Size(); i++)
          {
            values[i].Value() = v1[i].Value();
            values[i].DValue(0) = v1[i].DValue(0);
            values[i].DDValue(0) = v1[i].DValue(0) || v1[i].DDValue(0);
          }
      }
  }
};



  
  template <typename OP>
  class cl_BinaryOpCF : public T_CoefficientFunction<cl_BinaryOpCF<OP>>
{
  typedef T_CoefficientFunction<cl_BinaryOpCF<OP>> BASE;
  shared_ptr<CoefficientFunction> c1, c2;
  OP lam;
  string opname;
  using BASE::is_complex;
  using BASE::Dimension;
  using BASE::SetDimension;
  using BASE::SetDimensions;
  using BASE::Evaluate;
  using typename BASE::T_DJC;    
public:
  cl_BinaryOpCF() = default;
  cl_BinaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                 shared_ptr<CoefficientFunction> ac2, 
                 OP alam, string aopname)
    : BASE(ac1->Dimension(), ac1->IsComplex() || ac2->IsComplex()),
      c1(ac1), c2(ac2), lam(alam),
      opname(aopname)
  {
    int dim1 = c1->Dimension();
    int dim2 = c2->Dimension();
    // if (!(c1->Dimensions() == c2->Dimensions()))  // too critical ???
    if (dim1 != dim2)
      throw Exception ("Dimensions don't match, op = "+opname + " dims1 = " + ToString(c1->Dimensions()) + ", dims2 = " + ToString(c2->Dimensions()));
    is_complex = c1->IsComplex() || c2->IsComplex();
    this->elementwise_constant = c1->ElementwiseConstant() && c2->ElementwiseConstant();
    SetDimensions (c1->Dimensions());
  }

  virtual void DoArchive (Archive & archive) override
  {
      BASE::DoArchive(archive);
      archive.Shallow(c1).Shallow(c2) & opname;
  }

  virtual string GetDescription () const override
  {
    return string("binary operation '")+opname+"'";
  }
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    // code.Declare (code.res_type, index, this->Dimensions());
    code.Declare (index, this->Dimensions(), this->IsComplex());

    if (code_uses_tensors)
      {
        code.body += "for (int i = 0; i < "+ToString(this->Dimension())+"; i++)\n";
        code.body += "var_" + ToString(index) + "[i] = ";
        if(opname.size()>2) // atan2, pow, etc.
          {
            code.body += opname + '(' + "var_" + ToString(inputs[0]) + "[i],";
            code.body += "var_" + ToString(inputs[1]) + "[i]); \n";            
          }
        else
          {
            code.body += "var_" + ToString(inputs[0]) + "[i]" + opname;
            code.body += "var_" + ToString(inputs[1]) + "[i]; \n";
          }
      }

    else
      {
        for (int i = 0; i < this->Dimension(); i++)
          {
            auto op1 = Var(inputs[0], i, c1->Dimensions()).S();
            auto op2 = Var(inputs[1], i, c2->Dimensions()).S();
            string expr;
            if(opname.size()>2) // atan2, pow, etc.
              expr = opname + '(' + op1 + ',' + op2 + ')';
            else // +,-,*,/, etc.
              expr = op1 + ' ' + opname + ' ' + op2;
            code.body += Var(index,i,this->Dimensions()).Assign( expr, false );
          }
      }
  }

  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    c2->TraverseTree (func);
    func(*this);
  }

  virtual bool DefinedOn (const ElementTransformation & trafo) override
  { return c1->DefinedOn(trafo) && c2->DefinedOn(trafo); }

  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1, c2 }); }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
  {
    return lam (c1->Evaluate(ip), c2->Evaluate(ip));
  }

  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const override
  {
    return lam (c1->EvaluateComplex(ip), c2->EvaluateComplex(ip));
  }

  virtual double EvaluateConst () const override
  {
    return lam (c1->EvaluateConst(), c2->EvaluateConst());
  }


  virtual void Evaluate(const BaseMappedIntegrationPoint & mip,
                        FlatVector<> result) const override
  {
    size_t dim = Dimension();
    STACK_ARRAY(double, hmem, dim);
    FlatVector<> temp(dim, hmem);

    c1->Evaluate (mip, result);
    c2->Evaluate (mip, temp);
    for (int i = 0; i < result.Size(); i++)
      result(i) = lam (result(i), temp(i));
  }


  virtual void Evaluate(const BaseMappedIntegrationPoint & mip,
                        FlatVector<Complex> result) const override
  {
    size_t dim = Dimension();
    if (!is_complex)
      {
        STACK_ARRAY(double, hmem, dim);
        FlatVector<> temp(dim, &hmem[0]);
        Evaluate (mip, temp);
        result = temp;
        return;
      }
    
    STACK_ARRAY(Complex, hmem, dim);
    FlatVector<Complex> temp(dim, hmem);

    c1->Evaluate (mip, result);
    c2->Evaluate (mip, temp);
    for (int i = 0; i < result.Size(); i++)
      result(i) = lam (result(i), temp(i));
  }



  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        BareSliceMatrix<> result) const override
  {
    size_t dim = Dimension();
    size_t np = ir.Size();
    STACK_ARRAY(double, hmem, np*dim);
    FlatMatrix<> temp(np, dim, hmem);

    c1->Evaluate (ir, result);
    c2->Evaluate (ir, temp);
    /*
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lam (result(i), temp(i));
    */
    for (size_t i = 0; i < np; i++)
      for (size_t j = 0; j < dim; j++)
        result(i,j) = lam(result(i,j), temp(i,j));
  }

  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        BareSliceMatrix<Complex> result) const override
  {
    size_t dim = Dimension();    
    if (!is_complex)
      {
        STACK_ARRAY(double, hmem, ir.Size()*dim);
        FlatMatrix<> temp(ir.Size(), dim, &hmem[0]);
        Evaluate (ir, temp);
        result.AddSize(ir.Size(), dim) = temp;
        return;
      }

        
    STACK_ARRAY(double, hmem, 2*ir.Size()*dim);
    FlatMatrix<Complex> temp(ir.Size(), dim, reinterpret_cast<Complex*> (&hmem[0]));

    c1->Evaluate (ir, result);
    c2->Evaluate (ir, temp);
    /*
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lam(result(i), temp(i));
    */
    size_t np = ir.Size();
    for (size_t i = 0; i < np; i++)
      for (size_t j = 0; j < dim; j++)
        result(i,j) = lam (result(i,j), temp(i,j));
    
  }


  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
  {
    size_t np = ir.Size();
    size_t mydim = Dimension();
    STACK_ARRAY(T, hmem, np*mydim);
    FlatMatrix<T,ORD> temp(mydim, np, &hmem[0]);
    c1->Evaluate (ir, values);
    c2->Evaluate (ir, temp);
    for (size_t i = 0; i < mydim; i++)
      for (size_t j = 0; j < np; j++)
        values(i,j) = lam (values(i,j), temp(i,j));
  }

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    size_t np = ir.Size();
    size_t dim = Dimension();
    
    auto in0 = input[0];
    auto in1 = input[1];
    for (size_t i = 0; i < dim; i++)
      for (size_t j = 0; j < np; j++)
        values(i,j) = lam (in0(i,j), in1(i,j));
  }

  using CoefficientFunction::Operator;
  virtual shared_ptr<CoefficientFunction>
  Operator (const string & name) const override
  { throw Exception ("binarycf "+opname+" does not provide Operator"); }
  
  virtual shared_ptr<CoefficientFunction>
  Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override
  { throw Exception ("binarycf "+opname+" does not provide a derivative"); }

  virtual shared_ptr<CoefficientFunction>
  DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  { return BASE::DiffJacobi(var, cache); }

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv,
                               FlatVector<bool> nonzero_dderiv) const override
  {
    size_t dim = Dimension();    
    Vector<bool> v1(dim), v2(dim), d1(dim), d2(dim), dd1(dim), dd2(dim);
    c1->NonZeroPattern(ud, v1, d1, dd1);
    c2->NonZeroPattern(ud, v2, d2, dd2);
    for (int i = 0; i < nonzero.Size(); i++)
      {
        if (opname == "+" || opname == "-")
          {
            nonzero(i) = v1(i) || v2(i);
            nonzero_deriv(i) = d1(i) || d2(i);
            nonzero_dderiv(i) = dd1(i) || dd2(i);
          }
        else if (opname == "*")
          {
            nonzero(i) = v1(i) && v2(i);
            nonzero_deriv(i) = (v1(i) && d2(i)) || (d1(i) && v2(i));
            nonzero_dderiv(i) = (v1(i) && dd2(i)) || (d1(i) && d2(i)) || (dd1(i) && v2(i));
          }
        else
          {
            nonzero(i) = v1(i) || v2(i);
            nonzero_deriv(i) = d1(i) || d2(i);
            nonzero_dderiv(i) = d1(i) || dd1(i) || d2(i) || dd2(i);
          }
      }
  }
  */

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,NonZero>> values) const override
  {
    size_t dim = Dimension();    
    Vector<AutoDiffDiff<1,NonZero>> v1(dim), v2(dim);
    c1->NonZeroPattern(ud, v1);
    c2->NonZeroPattern(ud, v2);
    for (int i = 0; i < values.Size(); i++)
      {
        if (opname == "+" || opname == "-")
          values(i) = v1(i) + v2(i);
        else if (opname == "*")
          values(i) = v1(i) * v2(i);
        else
          {
            values(i).Value() = v1(i).Value() || v2(i).Value();
            values(i).DValue(0) = v1(i).DValue(0) || v2(i).DValue(0);
            values(i).DDValue(0) = v1(i).DValue(0) || v2(i).DValue(0) || v1(i).DDValue(0) || v2(i).DDValue(0);
          }
      }
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                               FlatVector<AutoDiffDiff<1,NonZero>> values) const override
  {
    auto v1 = input[0];
    auto v2 = input[1];
    for (int i = 0; i < values.Size(); i++)
      {
        if (opname == "+" || opname == "-")
          values(i) = v1(i) + v2(i);
        else if (opname == "*")
          values(i) = v1(i) * v2(i);
        else
          {
            values(i).Value() = v1(i).Value() || v2(i).Value();
            values(i).DValue(0) = v1(i).DValue(0) || v2(i).DValue(0);
            values(i).DDValue(0) = v1(i).DValue(0) || v2(i).DValue(0) || v1(i).DDValue(0) || v2(i).DDValue(0);
          }
      }
  }
};


  
  // extern shared_ptr<CoefficientFunction> shape;  // for shape derivative
  // information howto treat GridFunctions (and maybe other stuff) for shape-derivatives
  class DiffShapeCF : public ConstantCoefficientFunction
  {
  public:
    DiffShapeCF() : ConstantCoefficientFunction(1) {
      SetVariable();
    }
    ~DiffShapeCF() override;
    Array<shared_ptr<CoefficientFunction>> Eulerian_gridfunctions;
  };






  template <typename OP>
INLINE shared_ptr<CoefficientFunction> BinaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                                  shared_ptr<CoefficientFunction> c2, 
                                                  OP lam,
                                                  string opname)
{
  return make_shared<cl_BinaryOpCF<OP>>(c1, c2, lam, opname);
}


  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeConstantCoefficientFunction (Complex c);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeComponentCoefficientFunction (shared_ptr<CoefficientFunction> c1, int comp);
  
  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeVectorialCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeSubTensorCoefficientFunction (shared_ptr<CoefficientFunction> c1,
                                    int first, Array<int> num, Array<int> dist);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeExtendDimensionCoefficientFunction (shared_ptr<CoefficientFunction> c1,
                                      Array<int> dims, Array<int> pos, Array<int> stride);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeTensorTransposeCoefficientFunction (shared_ptr<CoefficientFunction> c1, Array<int> ordering);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeTensorTransposeCoefficientFunction (shared_ptr<CoefficientFunction> c1, int i1, int i2);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeTensorTraceCoefficientFunction (shared_ptr<CoefficientFunction> c1, int i1, int i2);

  // cf_ijk v0_i v1_j v2_k
  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeVectorContractionCoefficientFunction (shared_ptr<CoefficientFunction> c1,
                                            Array<shared_ptr<CoefficientFunction>> vectors);

  // contract index'th index with vector
  // cf_i_0, .., i_ind-1, k, i_ind+1, i_dim  v_k
  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeSingleContractionCoefficientFunction (shared_ptr<CoefficientFunction> c1,
                                            shared_ptr<CoefficientFunction> vec,
                                            int index);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeCoordinateCoefficientFunction (int comp);

  // for DG jump terms 
  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeOtherCoefficientFunction (shared_ptr<CoefficientFunction> me);
  NGS_DLL_HEADER bool IsOtherCoefficientFunction (CoefficientFunction & coef);

  
  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeDomainWiseCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci);
  



  /* ******************** matrix operations ********************** */
  


  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator+ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);
  
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator- (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator* (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);
  // coponent-wise multiplication
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> CWMult (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator* (double v1, shared_ptr<CoefficientFunction> c2);
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator* (Complex v1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator- (shared_ptr<CoefficientFunction> c1);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> CreateWrapperCF (shared_ptr<CoefficientFunction> cf);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> ConjCF (shared_ptr<CoefficientFunction> c1);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> InnerProduct (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> CrossProduct (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator/ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator/ (double val, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> IdentityCF (int dim);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> IdentityCF (FlatArray<int> dims);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> ZeroCF (FlatArray<int> dims);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> ConstantCF (double val);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> UnitVectorCF (int dim, int coord);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> LeviCivitaCF(int dimension);

  NGS_DLL_HEADER
  shared_ptr <CoefficientFunction> EinsumCF(const string &index_signature,
                                            const Array <shared_ptr<CoefficientFunction>>& cfs,
                                            const map<string, bool> &options = {});

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> TransposeCF (shared_ptr<CoefficientFunction> coef);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> ReshapeCF (shared_ptr<CoefficientFunction> coef,
                                             FlatArray<int> adims);

  INLINE
  shared_ptr<CoefficientFunction> ReshapeCF (shared_ptr<CoefficientFunction> coef, int s)  
  { return ReshapeCF (std::move(coef), Array<int>{s}); }

  INLINE
  shared_ptr<CoefficientFunction> ReshapeCF (shared_ptr<CoefficientFunction> coef, int h, int w)  
  { return ReshapeCF (std::move(coef), Array<int>{h,w}); }

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> InverseCF (shared_ptr<CoefficientFunction> coef);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> DeterminantCF (shared_ptr<CoefficientFunction> coef);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> CofactorCF (shared_ptr<CoefficientFunction> coef);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> SymmetricCF (shared_ptr<CoefficientFunction> coef);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> SkewCF (shared_ptr<CoefficientFunction> coef);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> TraceCF (shared_ptr<CoefficientFunction> coef);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> NormCF (shared_ptr<CoefficientFunction> coef);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> EigCF (shared_ptr<CoefficientFunction> coef);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> IfPos (shared_ptr<CoefficientFunction> cf_if,
                                         shared_ptr<CoefficientFunction> cf_then,
                                         shared_ptr<CoefficientFunction> cf_else);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> Real(shared_ptr<CoefficientFunction> cf);
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> Imag(shared_ptr<CoefficientFunction> cf);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> Freeze (shared_ptr<CoefficientFunction> cf);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  CreateMinimizationCF(shared_ptr<CoefficientFunction> expression,
                       shared_ptr<CoefficientFunction> startingpoint,
                       std::optional<double> atol, std::optional<double> rtol,
                       std::optional<int> maxiter, std::optional<bool> allow_fail);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  CreateMinimizationCF(shared_ptr<CoefficientFunction> expression,
                       const Array<shared_ptr<CoefficientFunction>> &startingpoints,
                       std::optional<double> tol, std::optional<double> rtol,
                       std::optional<int> maxiter, std::optional<bool> allow_fail);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  CreateNewtonCF (shared_ptr<CoefficientFunction> expression,
                  shared_ptr<CoefficientFunction> startingpoint,
                  std::optional<double> atol,
                  std::optional<double> rtol,
                  std::optional<int> maxiter,
                  std::optional<bool> allow_fail);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  CreateNewtonCF (shared_ptr<CoefficientFunction> expression,
                  const Array<shared_ptr<CoefficientFunction>> &startingpoints,
                  std::optional<double> tol,
                  std::optional<double> rtol,
                  std::optional<int> maxiter,
                  std::optional<bool> allow_fail);


class CompiledCoefficientFunctionInterface : public CoefficientFunction
{
protected:
  Array<CoefficientFunction*> steps;
public:
  using CoefficientFunction::CoefficientFunction;
  virtual Code GenerateProgram (int deriv, bool simd) const = 0;
  const Array<CoefficientFunction*> & Steps() const { return steps; }
};

NGS_DLL_HEADER
shared_ptr<CompiledCoefficientFunctionInterface> Compile (shared_ptr<CoefficientFunction> c, bool realcompile=false, int maxderiv=2, bool wait=false, bool keep_files=false);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> LoggingCF (shared_ptr<CoefficientFunction> func, string logfile="stdout");

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> CacheCF (shared_ptr<CoefficientFunction> func);

  NGS_DLL_HEADER
  void PrecomputeCacheCF (CoefficientFunction & func, SIMD_BaseMappedIntegrationRule & mir,
                          LocalHeap & lh);

  NGS_DLL_HEADER Array<CoefficientFunction*> FindCacheCF (CoefficientFunction & func);
  NGS_DLL_HEADER
  void PrecomputeCacheCF (const Array<CoefficientFunction*> & cachecfs, BaseMappedIntegrationRule & mir,
                          LocalHeap & lh);
  NGS_DLL_HEADER
  void PrecomputeCacheCF (const Array<CoefficientFunction*> & cachecfs, SIMD_BaseMappedIntegrationRule & mir,
                          LocalHeap & lh);

  
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> NormalVectorCF
  (int dim, optional<BitArray> inverted_faces=nullopt);
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> TangentialVectorCF (int dim, bool consistent);
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> JacobianMatrixCF (int dims, int dimr);
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> WeingartenCF (int dim);
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> VertexTangentialVectorsCF (int dim);
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> EdgeFaceTangentialVectorsCF (int dim);
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> EdgeCurvatureCF (int dim);

  template <typename OP /* , typename OPC */>
shared_ptr<CoefficientFunction> UnaryOpCF(shared_ptr<CoefficientFunction> c1,
                                          OP lam, /* OPC lamc, */ string name="undefined")
{
  if (c1->GetDescription() == "ZeroCF" && lam(0.)==0.)
  {
    return ZeroCF(c1->Dimensions());
  }
  return shared_ptr<CoefficientFunction> (new cl_UnaryOpCF<OP /* ,OPC */> (c1, lam/* , lamc */, name));
}

  
}





#endif
