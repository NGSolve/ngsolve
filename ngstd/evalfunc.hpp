#ifndef FILE_EVALFUNC
#define FILE_EVALFUNC

/**************************************************************************/
/* File:   evalfunc.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 95                                                    */
/**************************************************************************/


namespace ngstd
{




  class GenericVariable 
  {
    int dim;
    bool iscomplex;
    double * data;
  public:
    GenericVariable (const GenericVariable & v2)
    {
      dim = v2.dim;
      iscomplex = v2.iscomplex;
      int hdim = iscomplex ? 2*dim : dim;
      data = new double[hdim];
      for (int j = 0; j < hdim; j++)
        data[j] = v2.data[j];
    }
    GenericVariable (GenericVariable && v2)
    {
      cout << "move constr" << endl;
      dim = v2.dim;
      iscomplex = v2.iscomplex;
      int hdim = iscomplex ? 2*dim : dim;
      data = new double[hdim];
      for (int j = 0; j < hdim; j++)
        data[j] = v2.data[j];
    }
    GenericVariable (bool acomplex = false, int adim = 1)
      : dim(adim), iscomplex(acomplex)
    {
      int hdim = iscomplex ? 2*dim : dim;
      data = new double[hdim];
    }
    ~GenericVariable ()  
    { 
      delete [] data; 
    }
    GenericVariable & operator= (const GenericVariable & v2)
    {
      dim = v2.dim;
      iscomplex = v2.iscomplex;
      int hdim = iscomplex ? 2*dim : dim;
      delete [] data;
      data = new double[hdim];
      for (int j = 0; j < hdim; j++)
        data[j] = v2.data[j];
      return *this;
    }

    int Dimension() const { return dim; }
    bool IsComplex () const { return iscomplex; }

    double & ValueDouble(int i = 0) { return data[i]; }
    complex<double> & ValueComplex (int i = 0) 
    { return reinterpret_cast<std::complex<double>*>(data)[i]; }
    
    const double & ValueDouble(int i = 0) const { return data[i]; }
    const std::complex<double> & ValueComplex (int i = 0) const 
    { return reinterpret_cast<complex<double>*>(data)[i]; }

    
    template <typename SCAL> SCAL Value (int i) const;
  };

  template<> 
  inline double GenericVariable::Value<double> (int i) const 
  { 
    if (iscomplex) throw Exception ("Value<double> called for complex variable");
    return data[i]; 
  }
  template<> 
  inline std::complex<double> GenericVariable::Value<std::complex<double>> (int i) const 
  { 
    if (iscomplex)
      return complex<double> (data[2*i], data[2*i+1]);
    else
      return complex<double> (data[i]);
  }


  inline ostream & operator<< (ostream & ost, const GenericVariable & var)
  {
    if (var.IsComplex())
      for (int i = 0; i < var.Dimension(); i++)
        ost << var.ValueComplex(i) << ", ";
    else
      for (int i = 0; i < var.Dimension(); i++)
        ost << var.ValueDouble(i) << ", ";
    return ost;
  }









/**
   Numerical expression parser.
   The expression is stored in reverse Polnish notation.
   The evaluation tree can be filled form an external parser,
   see Addxxx methods.
*/
class NGS_DLL_HEADER EvalFunction
{

  ///
  enum EVAL_TOKEN
  {
    ADD = '+', SUB = '-', MULT = '*', DIV = '/', LP ='(', RP = ')',
    COMMA = ',',
    NEG = 100, 
    VEC_ADD, VEC_SUB, VEC_SCAL_MULT, SCAL_VEC_MULT, VEC_VEC_MULT, VEC_SCAL_DIV, VEC_ELEM, VEC_DIM,
    AND, OR, NOT, GREATER, LESS, GREATEREQUAL, LESSEQUAL, EQUAL,
    CONSTANT, IMAG, VARIABLE, FUNCTION, GLOBVAR, GLOBGENVAR, /* COEFF_FUNC,*/ END, STRING,
    SIN, COS, TAN, ATAN, ATAN2, EXP, LOG, ABS, SIGN, SQRT, STEP,
    BESSELJ0, BESSELY0, BESSELJ1, BESSELY1
  };

public:
  /// 
  EvalFunction ();
  /// parse from input stream
  EvalFunction (istream & aist);
  /// parse from string
  EvalFunction (const string & str);
  ///
  EvalFunction (const EvalFunction & eval2);
  /// 
  virtual ~EvalFunction ();

  /// parse from stream
  bool Parse (istream & aist);
  /// define constant 
  void DefineConstant (const string & name, double val);
  /// define constant 
  void DefineGlobalVariable (const string & name, double * var);
  /// define constant 
  void DefineGlobalVariable (const string & name, GenericVariable * var);
  /// define arguments 
  void DefineArgument (const string & name, int num, int vecdim = 1, bool iscomplex = false);

  /// evaluate function
  double Eval (const double * x = NULL) const;
  /// evaluate multi-value function
  void Eval (const double * x, double * y, int ydim) const;

  /// evaluate function
  complex<double> Eval (const complex<double> * x = NULL) const;
  /// evaluate multi-value complex function
  void Eval (const complex<double> * x, complex<double> * y, int ydim) const;
  /// evaluate multi-value complex function with real result
  void Eval (const complex<double> * x, double * y, int ydim) const;

  /*
  /// evaluate multi-value function
  template <typename TIN>
  void Eval (const TIN * x, complex<double> * y, int ydim) const;
  */
  template <typename TIN, typename TCALC>
  void Eval (const TIN * x, TCALC * stack) const;


  /// is expression complex valued ?
  bool IsComplex () const;

  /// is expression complex valued ?
  bool IsResultComplex () const { return res_type.iscomplex; }

  /// is expression a constant ?
  bool IsConstant () const;
  /// evaluate the constant value
  double EvalConstant () const;

  /// vector dimension of result
  int Dimension() const { return res_type.vecdim; }

  /// push constant on stack. 
  void AddConstant (double val)
  { program.Append (step (val)); }

  /// push variable x[varnum-1].
  void AddVariable (int varnum)
  { program.Append (step(varnum)); }

  /// push pointer to global double value.
  void AddGlobVariable (const double * dp)
  { program.Append (step(dp)); }

  /// push pointer to global double value.
  void AddGlobVariable (const GenericVariable * dp)
  { program.Append (step(dp)); }

  /// push operation. 
  void AddOperation (EVAL_TOKEN op)
  { program.Append (step(op)); }

  /// push function call. 
  void AddFunction (double (*fun) (double))
  { program.Append (step(fun)); }

  /// print expression
  void Print (ostream & ost) const;
protected:
   
  /// one step of evaluation
  class step
  {
  public:
    ///
    EVAL_TOKEN op;
    /// the data 
    union UNION_OP
    {
      /// a constant value
      double val;
      /// a pointer to a global variable
      const double *globvar;
      /// a pointer to a global variable
      const GenericVariable *globgenvar;
      /// the input argument number varnum
      int varnum;
      /// a pointer to a unary function
      double (*fun) (double);
    }; 
    ///
    UNION_OP operand;

    /// dimension of vector
    short int vecdim;

    step () { ; }

    step (EVAL_TOKEN hop)
    { 
      op = hop;
      operand.val = 0;
    }

    step (double hval)
    { 
      op = CONSTANT;
      operand.val = hval;
    }

    step (int varnum)
    { 
      op = VARIABLE;
      operand.varnum = varnum;
    }

    step (const double * aglobvar)
    { 
      op = GLOBVAR;
      operand.globvar = aglobvar;
    }

    step (const GenericVariable * aglobvar)
    { 
      op = GLOBGENVAR;
      operand.globgenvar = aglobvar;
    }

    step (double (*fun) (double))
    {
      op = FUNCTION;
      operand.fun = fun;
    }
  };

  /// the evaluation sequence
  Array<step> program;

  class ResultType
  {
  public:
    int vecdim;
    bool isbool;
    bool iscomplex;
    ResultType ()
      : vecdim(1), isbool(false), iscomplex(false)
    { ; }
  };

  ResultType res_type;
  const double eps;

  /// parsing expression (standard parsing grammar)
  ResultType ParseExpression ();
  /// parsing expression (standard parsing grammar)
  ResultType ParseCommaExpression ();
  /// parsing expression (standard parsing grammar)
  ResultType ParseSubExpression ();
  /// parsing expression (standard parsing grammar)
  ResultType ParseTerm ();
  /// parsing expression (standard parsing grammar)
  ResultType ParsePrimary ();

  /// parse from stream
  istream * ist;

  ///
  EVAL_TOKEN token;
  ///
  double num_value;
  ///
  char string_value[1000];
  ///
  int var_num, var_dim;
  ///
  bool var_iscomplex;
  ///
  double * globvar;
  GenericVariable * globgenvar;
  streampos lastpos;

  typedef double(*TFUNP) (double);
  /// registered functions
  static SymbolTable<TFUNP> functions;

  /// registered constants
  SymbolTable<double> constants;

  /// registered variables
  SymbolTable<double*> globvariables;
  /// registered variables
  SymbolTable<GenericVariable*> genericvariables;
  
public:
  /// the arguments passed to the function
  struct argtype
  {
    int argnum;
    int dim;
    bool iscomplex;
  public:
    argtype ()
      : argnum(-1), dim(1), iscomplex(false) { ; }
    argtype (int aanum, int adim = 1, bool acomplex = false)
      : argnum(aanum), dim(adim), iscomplex(acomplex) { ; }
  };
  SymbolTable<argtype> arguments;
  int num_arguments;

  /// returns last token
  EVAL_TOKEN GetToken() const
    { return token; }

  /// returns num_value of last token
  double GetNumValue() const
    { return num_value; }

  /// returns variable number of last token
  int GetVariableNumber() const
    { return var_num; }
  /// returns dimension of variable of last token
  int GetVariableDimension() const
    { return var_dim; }
  bool GetVariableIsComplex() const
    { return var_iscomplex; }

  /// returns identifier of last token
  const char * GetStringValue() const
    { return string_value; }
  
  /// read next token
  void ReadNext(bool optional = true);
  void WriteBack();

  bool ToBool (double x)  const { return x > eps; }
  bool ToBool (complex<double> x) const { return x.real() > eps; }
  double CheckReal (double x)  const { return x; }
  double CheckReal (complex<double> x) const 
  {
    if (x.imag() != 0) cerr << "illegal complex value" << endl; 
    return x.real();
  }

  double Abs (double x) const { return std::fabs(x); }
  double Abs (complex<double> x) const { return abs(x); }
};


}


#endif


