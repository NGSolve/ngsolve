#ifndef FILE_COEFFICIENT
#define FILE_COEFFICIENT

/*********************************************************************/
/* File:   coefficient.hh                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace pybind11 { class module; };

namespace ngfem
{
  /** 
      coefficient functions
  */

  typedef enum {
    CF_Type_undefined,
    CF_Type_constant,
    CF_Type_vectorial,
    CF_Type_coordinate,
    CF_Type_norm,
    CF_Type_trans,
    CF_Type_component,
    CF_Type_real,
    CF_Type_imag,
    CF_Type_ifpos,
    CF_Type_normal_vector,
    CF_Type_tangential_vector,
    CF_Type_mesh_size,
    CF_Type_scale,
    CF_Type_scale_complex,
    CF_Type_add,    
    CF_Type_sub,    
    CF_Type_mult,    
    CF_Type_div,    
    CF_Type_domainconst,    
    CF_Type_domainwise,    
    CF_Type_unary_op,
    CF_Type_binary_op,
    CF_Type_usertype,
    CF_Type_eig,
  } CF_Type;
  
  class NGS_DLL_HEADER CoefficientFunction
  {
  private:
    int dimension;
    Array<int> dims;
  protected:
    bool is_complex;
  public:
    ///
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
        dims = Array<int> ( { dimension } );
    }
    
    ///
    virtual ~CoefficientFunction ();

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
    
    // virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const;    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const;

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const;
    // virtual void EvaluateSoA (const BaseMappedIntegrationRule & ir, AFlatMatrix<Complex> values) const;
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatArray<BareSliceMatrix<double,ColMajor>> input,
                           BareSliceMatrix<double,ColMajor> values) const
    {
      Evaluate (ir, Trans(values));
    }

    /*
    [[deprecated("Use Evaluate (SIMD) instead")]]        
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    {
      throw ExceptionNOSIMD (string("cf::Evaluate(simd, input->output) not overloaded for ")+typeid(*this).name());
    }

    [[deprecated("Use Evaluate (AutoDiff) instead")]]
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                AFlatMatrix<double> values, AFlatMatrix<double> deriv) const 
    {
      throw ExceptionNOSIMD (string("cf::EvaluateDeriv(simd) not overloaded for ")+typeid(*this).name());
    }

    [[deprecated("Use Evaluate (AutoDiffDiff) instead")]]    
    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                 AFlatMatrix<double> values, AFlatMatrix<double> deriv,
                                 AFlatMatrix<double> dderiv) const 
    {
      throw ExceptionNOSIMD (string("cf::EvaluateDDeriv(simd) not overloaded for ")+typeid(*this).name());
    }
    */
    
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

    /*
    [[deprecated("Use Evaluate (AutoDiff) instead")]]    
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                FlatArray<AFlatMatrix<>*> input,
                                FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result,
                                AFlatMatrix<> deriv) const
    {
      throw ExceptionNOSIMD (string("cf::EvaluateDeriv(simd,in-out) not overloaded for ")+typeid(*this).name());
    }

    [[deprecated("Use Evaluate (AutoDiffDiff) instead")]]        
    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                 FlatArray<AFlatMatrix<>*> input,
                                 FlatArray<AFlatMatrix<>*> dinput,
                                 FlatArray<AFlatMatrix<>*> ddinput,
                                 AFlatMatrix<> result,
                                 AFlatMatrix<> deriv,
                                 AFlatMatrix<> dderiv) const
    {
      throw ExceptionNOSIMD (string("cf::EvaluateDDeriv(simd,in-out) not overloaded for ")+typeid(*this).name());
    }
    */
    
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
    int Dimension() const { return dimension; }
    FlatArray<int> Dimensions() const { return dims; }
    
    void SetDimensions (FlatArray<int> adims)
    {
      dims = adims;
      dimension = 1;
      for (int d : dims) dimension *= d;
    }
    
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
      /*
      Complex f = EvaluateComplex (ip);
      result(0) = f; 
      */
    }

    /*
    [[deprecated("Use Evaluate (AutoDiff) instead")]]
    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                FlatMatrix<> result,
                                FlatMatrix<> deriv) const
    {
      Evaluate (ir, result);
      deriv = 0;
    }
    */
    
    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                FlatMatrix<Complex> result,
                                FlatMatrix<Complex> deriv) const
    {
      Evaluate (ir, result);
      deriv = 0;
    }

    /*
    [[deprecated("Use Evaluate (AutoDiffDiff) instead")]]
    virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatMatrix<> result,
                                 FlatMatrix<> deriv,
                                 FlatMatrix<> dderiv) const
    {
      EvaluateDeriv (ir, result, deriv);
      dderiv = 0;
    }
    */
    
    /*
    virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatMatrix<Complex> result,
                                 FlatMatrix<Complex> deriv,
                                 FlatMatrix<Complex> dderiv) const
    {
      EvaluateDeriv (ir, result, deriv);
      dderiv = 0;
    }

    
    [[deprecated("Use Evaluate (AutoDiff) instead")]]
    virtual void EvaluateDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatArray<FlatMatrix<>*> input,
                                 FlatArray<FlatMatrix<>*> dinput,
                                 FlatMatrix<> result,
                                 FlatMatrix<> deriv) const
    {
      EvaluateDeriv (ir, result, deriv);
    }

    [[deprecated("Use Evaluate (AutoDiffDiff) instead")]]    
    virtual void EvaluateDDeriv (const BaseMappedIntegrationRule & ir,
                                 FlatArray<FlatMatrix<>*> input,
                                 FlatArray<FlatMatrix<>*> dinput,
                                 FlatArray<FlatMatrix<>*> ddinput,
                                 FlatMatrix<> result,
                                 FlatMatrix<> deriv,
                                 FlatMatrix<> dderiv) const override
    {
      EvaluateDDeriv (ir, result, deriv, dderiv);
    }
    */

    virtual bool ElementwiseConstant () const { return false; }
    // virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const;
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<bool> nonzero,
                                 FlatVector<bool> nonzero_deriv,
                                 FlatVector<bool> nonzero_dderiv) const;

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                                 FlatVector<AutoDiffDiff<1,bool>> values) const
    {
      cout << string("nonzero in-out not overloaded for type")+typeid(*this).name() << endl;
      Vector<bool> nz(values.Size()), nzd(values.Size()), nzdd(values.Size());
      NonZeroPattern (ud, nz, nzd, nzdd);
      for (size_t i = 0; i < values.Size(); i++)
        {
          values(i).Value() = nz(i);
          values(i).DValue(0) = nzd(i);
          values(i).DDValue(0) = nzdd(i);
        }
    }
    
    virtual void PrintReport (ostream & ost) const;
    virtual void PrintReportRec (ostream & ost, int level) const;
    virtual string GetDescription () const;
    
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func);
    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const
    { return Array<shared_ptr<CoefficientFunction>>(); }
    virtual bool StoreUserData() const { return false; }

    virtual CF_Type GetType() const { return CF_Type_undefined; } 
    virtual void DoArchive (Archive & archive)
    {
      archive & dimension & dims & is_complex;
    } 
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
      BareSliceMatrix<SIMD<double>> hvalues(2*values.Dist(), &values(0).Value(), DummySize(Dimension(), ir.Size()));
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
      FlatMatrix<double> hvalues(ir.Size(), 3*values.Dist(), &values(0).Value());
      Evaluate (ir, hvalues);
      for (size_t i = 0; i < ir.Size(); i++)
        for (size_t j = Dimension(); j-- > 0; )
          values(i,j) = hvalues(i,j);
    }
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const override
    {
      BareSliceMatrix<SIMD<double>> hvalues(3*values.Dist(), &values(0).Value(), DummySize(Dimension(), ir.Size()));
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
    
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<bool> nonzero,
                                 FlatVector<bool> nonzero_deriv,
                                 FlatVector<bool> nonzero_dderiv) const override
    {
      nonzero = true;
      nonzero_deriv = false;
      nonzero_dderiv = false;
    }

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                                 FlatVector<AutoDiffDiff<1,bool>> values) const override
    {
      values = AutoDiffDiff<1,bool> (true);
    }
    
  };


  

  template <typename TCF, typename BASE = CoefficientFunction>
  class T_CoefficientFunction : public BASE
  {
  protected:
    using BASE::IsComplex;
    using BASE::Dimension;
  public:
    using BASE::BASE;
      
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <SIMD<double>> */ (ir, values); }
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const
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

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<SIMD<double>>> input,
                           BareSliceMatrix<SIMD<double>> values) const
    {  static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <SIMD<double>> */ (ir, input, values); }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate (ir, Trans(values)); }
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiff<1,double>> values) const
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate (ir, Trans(values)); }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiff<1,double>,ColMajor>> input,
                           BareSliceMatrix<AutoDiff<1,double>,ColMajor> values) const
    { static_cast<const TCF*>(this) -> T_Evaluate (ir, input, values); }
        
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <AutoDiff<1,SIMD<double>>> */ (ir, values); }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiff<1,SIMD<double>>>> input,
                           BareSliceMatrix<AutoDiff<1,SIMD<double>>> values) const
    {  static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <AutoDiff<1,SIMD<double>>> */ (ir, input, values); }


    virtual void Evaluate (const BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiffDiff<1,double>> values) const
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate (ir, Trans(values)); }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiffDiff<1,double>,ColMajor>> input,
                           BareSliceMatrix<AutoDiffDiff<1,double>,ColMajor> values) const
    { static_cast<const TCF*>(this) -> T_Evaluate (ir, input, values); }
    
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, 
                           BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const
    { static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <AutoDiffDiff<1,SIMD<double>>> */ (ir, values); }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>>> input,
                           BareSliceMatrix<AutoDiffDiff<1,SIMD<double>>> values) const
    {  static_cast<const TCF*>(this) -> /* template */ T_Evaluate /* <AutoDiffDiff<1,SIMD<double>>> */ (ir, input, values); }
  };


  /// The coefficient is constant everywhere
  class NGS_DLL_HEADER ConstantCoefficientFunction
    : public T_CoefficientFunction<ConstantCoefficientFunction, CoefficientFunctionNoDerivative>
  {
    ///
    double val;
    typedef T_CoefficientFunction<ConstantCoefficientFunction, CoefficientFunctionNoDerivative> BASE;
  public:
    ///
    ConstantCoefficientFunction (double aval);
    ///
    virtual ~ConstantCoefficientFunction ();
    ///
    using BASE::Evaluate;
    virtual bool ElementwiseConstant () const override { return true; }
    
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      return val;
    }

    virtual double EvaluateConst () const override
    {
      return val;
    }
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const override;
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const override;

    template <typename MIR, typename T, ORDERING ORD>
      void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const;
    template <typename MIR, typename T, ORDERING ORD>
      void T_Evaluate (const MIR & ir,
                       FlatArray<BareSliceMatrix<T,ORD>> input,                       
                       BareSliceMatrix<T,ORD> values) const
    { T_Evaluate (ir, values); }

    /*
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const override
    { values = val; }

    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                AFlatMatrix<> result, AFlatMatrix<> deriv) const override
    {
      result = val;
      deriv = 0.0;
    }
    
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                FlatArray<AFlatMatrix<>*> input, FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result, AFlatMatrix<> deriv) const override
    {
      result = val;
      deriv = 0.0;
    }

    virtual void EvaluateDDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                 FlatArray<AFlatMatrix<>*> input, FlatArray<AFlatMatrix<>*> dinput,
                                 FlatArray<AFlatMatrix<>*> ddinput,
                                 AFlatMatrix<> result, AFlatMatrix<> deriv,
                                 AFlatMatrix<> dderiv) const override
    {
      result = val;
      deriv = 0.0;
      dderiv = 0.0;
    }
    */
    
    
    virtual void PrintReport (ostream & ost) const override;
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override; 

    /*
    virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero) const
    {
      nonzero(0) = (val != 0.0);
    }
    */
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<bool> nonzero,
                                 FlatVector<bool> nonzero_deriv,
                                 FlatVector<bool> nonzero_dderiv) const override
    {
      nonzero(0) = (val != 0.0);
      nonzero_deriv = 0.0;
      nonzero_dderiv = 0.0;
    }

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                                 FlatVector<AutoDiffDiff<1,bool>> values) const override
    {
      values = AutoDiffDiff<1,bool> (val != 0.0);
    }
    


    
    virtual CF_Type GetType() const override { return CF_Type_constant; } 
    virtual void DoArchive (Archive & archive) override
    {
      CoefficientFunction::DoArchive(archive);
      archive & val;
    } 
  };



  /// The coefficient is constant everywhere
  class NGS_DLL_HEADER ConstantCoefficientFunctionC : public CoefficientFunction
  {
    ///
    Complex val;
  public:
    ConstantCoefficientFunctionC (Complex aval);
    virtual ~ConstantCoefficientFunctionC ();

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;
    virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const;

    virtual void Evaluate (const BaseMappedIntegrationPoint & mip, FlatVector<Complex> values) const;
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values) const;
    
    virtual void PrintReport (ostream & ost) const;
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const;
  };


  /// The coefficient is constant everywhere
  class NGS_DLL_HEADER ParameterCoefficientFunction : public CoefficientFunctionNoDerivative
  {
    ///
    double val;
  public:
    ///
    ParameterCoefficientFunction (double aval);
    ///
    virtual ~ParameterCoefficientFunction ();
    ///
    using CoefficientFunction::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      return val;
    }
    
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const override;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override
    { values.AddSize(Dimension(), ir.Size()) = val; }
    /*
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    { values = val; }
    */
    virtual void SetValue (double in) { val = in; }
    virtual double GetValue () { return val; }
    virtual void PrintReport (ostream & ost) const override;
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override;
  };

  

  /// The coefficient is constant in every sub-domain
  class NGS_DLL_HEADER DomainConstantCoefficientFunction  
    : public T_CoefficientFunction<DomainConstantCoefficientFunction, CoefficientFunctionNoDerivative>
  {
    ///
    Array<double> val;
    typedef T_CoefficientFunction<DomainConstantCoefficientFunction, CoefficientFunctionNoDerivative> BASE;    
  public:
    ///
    DomainConstantCoefficientFunction (const Array<double> & aval);
    ///
    virtual int NumRegions () override { return val.Size(); }
    ///
    virtual ~DomainConstantCoefficientFunction ();
    ///
      using T_CoefficientFunction<DomainConstantCoefficientFunction, CoefficientFunctionNoDerivative>::Evaluate;
      using CoefficientFunction::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override; 
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const override;
    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> values) const override;

    // virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const;
    // virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const;

    template <typename MIR, typename T, ORDERING ORD>
      void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const;
    template <typename MIR, typename T, ORDERING ORD>
      void T_Evaluate (const MIR & ir,
                       FlatArray<BareSliceMatrix<T,ORD>> input,                       
                       BareSliceMatrix<T,ORD> values) const
    { T_Evaluate (ir, values); }
    

    
    virtual double EvaluateConst () const override { return val[0]; }
    double operator[] (int i) const { return val[i]; }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override;
    virtual CF_Type GetType() const override { return CF_Type_domainconst; }
    virtual void DoArchive (Archive & archive) override
    {
        CoefficientFunction::DoArchive(archive);
        archive & val;
    }
    
  protected:
    void CheckRange (int elind) const
    {
      if (elind < 0 || elind >= val.Size())
        {
          ostringstream ost;
          ost << "DomainConstantCoefficientFunction: Element index "
              << elind << " out of range 0 - " << val.Size()-1 << endl;
          throw Exception (ost.str());
        }
    }
  };




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
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
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



  //////////////////

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













  // *************************** CoefficientFunction Algebra ********************************
  template <typename OP>
  class cl_UnaryOpCF : public T_CoefficientFunction<cl_UnaryOpCF<OP>>
{
  shared_ptr<CoefficientFunction> c1;
  OP lam;
  string name;
  typedef  T_CoefficientFunction<cl_UnaryOpCF<OP>> BASE;
public:
  cl_UnaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                OP alam, string aname="undefined")
    : BASE(ac1->Dimension(),
           ac1->IsComplex() && typeid (lam(Complex(0.0))) == typeid(Complex)),
      c1(ac1), lam(alam), name(aname)
  {
    this->SetDimensions (c1->Dimensions());
  }

  virtual string GetDescription () const override
  {
    return string("unary operation '")+name+"'";
  }

  virtual bool ElementwiseConstant () const override { return c1->ElementwiseConstant(); }
  
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    TraverseDimensions( this->Dimensions(), [&](int ind, int i, int j) {
        int i1, j1;
        GetIndex( c1->Dimensions(), ind, i1, j1 );
        code.body += Var(index,i,j).Assign( Var(inputs[0],i1,j1).Func(name) );
        });
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
                        FlatMatrix<Complex> result) const override
  {
    c1->Evaluate (ir, result);
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lam(result(i));
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
        if (name == "-") // actually not used that way
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

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
  {
    auto v1 = input[0];
    if (name == "-") // actually not used that way
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
    


  

  virtual CF_Type GetType() const override { return CF_Type_unary_op; }
  virtual void DoArchive (Archive & archive) override
  {
      archive & name;
      CoefficientFunction::DoArchive(archive);
  }
};

  template <typename OP /* , typename OPC */> 
shared_ptr<CoefficientFunction> UnaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                          OP lam, /* OPC lamc, */ string name="undefined")
{
  return shared_ptr<CoefficientFunction> (new cl_UnaryOpCF<OP /* ,OPC */> (c1, lam/* , lamc */, name));
}



  
  template <typename OP, typename NONZERO> 
  class cl_BinaryOpCF : public T_CoefficientFunction<cl_BinaryOpCF<OP,NONZERO>>
{
  typedef T_CoefficientFunction<cl_BinaryOpCF<OP,NONZERO>> BASE;
  shared_ptr<CoefficientFunction> c1, c2;
  OP lam;
  NONZERO lam_nonzero;
  string opname;
  using BASE::is_complex;
  using BASE::Dimension;
  using BASE::SetDimension;
  using BASE::SetDimensions;
  using BASE::Evaluate;  
public:
  cl_BinaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                 shared_ptr<CoefficientFunction> ac2, 
                 OP alam, NONZERO alam_nonzero, string aopname)
    : BASE(ac1->Dimension(), ac1->IsComplex() || ac2->IsComplex()),
      c1(ac1), c2(ac2), lam(alam),
      lam_nonzero(alam_nonzero),
      opname(aopname)
  {
    int dim1 = c1->Dimension();
    int dim2 = c2->Dimension();
    if (dim1 != dim2) throw Exception ("Dimensions don't match");
    is_complex = c1->IsComplex() || c2->IsComplex();
    SetDimensions (c1->Dimensions());
  }
  virtual string GetDescription () const override
  {
    return string("binary operation '")+opname+"'";
  }
  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
  {
    TraverseDimensions( c1->Dimensions(), [&](int ind, int i, int j) {
        int i2,j2;
        GetIndex(c2->Dimensions(), ind, i2, j2);
        code.body += Var(index,i,j).Assign(   Var(inputs[0],i,j).S()
                                            + opname
                                            + Var(inputs[1],i2,j2).S()
                                          );
    });
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
    
    STACK_ARRAY(double, hmem, 2*dim);
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
                        FlatMatrix<Complex> result) const override
  {
    size_t dim = Dimension();    
    if (!is_complex)
      {
        STACK_ARRAY(double, hmem, ir.Size()*dim);
        FlatMatrix<> temp(ir.Size(), dim, &hmem[0]);
        Evaluate (ir, temp);
        result = temp;
        return;
      }

        
    STACK_ARRAY(double, hmem, 2*ir.Size()*dim);
    FlatMatrix<Complex> temp(ir.Size(), dim, reinterpret_cast<Complex*> (&hmem[0]));

    c1->Evaluate (ir, result);
    c2->Evaluate (ir, temp);
    for (int i = 0; i < result.Height()*result.Width(); i++)
      result(i) = lam(result(i), temp(i));
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

  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                               FlatVector<AutoDiffDiff<1,bool>> values) const override
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

  
  virtual CF_Type GetType() const override {
      if(opname =="+") return CF_Type_add;
      if(opname =="-") return CF_Type_sub;
      if(opname =="*") return CF_Type_mult;
      if(opname =="/") return CF_Type_div;
      return CF_Type_binary_op;
  }
  virtual void DoArchive (Archive & archive) override 
  {
      archive & opname;
      CoefficientFunction::DoArchive(archive);
  }

};

  template <typename OP, typename NONZERO> 
INLINE shared_ptr<CoefficientFunction> BinaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                                  shared_ptr<CoefficientFunction> c2, 
                                                  OP lam,
                                                  NONZERO lam_nonzero,
                                                  string opname)
{
  return shared_ptr<CoefficientFunction> (new cl_BinaryOpCF<OP,NONZERO> 
                                          (c1, c2, lam, lam_nonzero, opname));
}




  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeComponentCoefficientFunction (shared_ptr<CoefficientFunction> c1, int comp);
  
  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeVectorialCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci);

  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeCoordinateCoefficientFunction (int comp);

  // for DG jump terms 
  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeOtherCoefficientFunction (shared_ptr<CoefficientFunction> me);


  
  NGS_DLL_HEADER shared_ptr<CoefficientFunction>
  MakeDomainWiseCoefficientFunction (Array<shared_ptr<CoefficientFunction>> aci);
  



  /* ******************** matrix operations ********************** */
  


  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator+ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);
  
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator- (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator* (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator* (double v1, shared_ptr<CoefficientFunction> c2);
  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator* (Complex v1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> InnerProduct (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> operator/ (shared_ptr<CoefficientFunction> c1, shared_ptr<CoefficientFunction> c2);

  NGS_DLL_HEADER
  shared_ptr<CoefficientFunction> TransposeCF (shared_ptr<CoefficientFunction> coef);

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
  shared_ptr<CoefficientFunction> Compile (shared_ptr<CoefficientFunction> c, bool realcompile=false, int maxderiv=2, bool wait=false);
}





#endif
