#ifndef FILE_EVALFUNC
#define FILE_EVALFUNC

/**************************************************************************/
/* File:   evalfunc.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 95                                                    */
/**************************************************************************/


namespace ngstd
{


/**
   Numerical expression parser.
   The expression is stored in reverse Polnish notation.
   The evaluatino tree can be filled form an external parser,
   see Addxxx methods.
*/
class EvalFunction
{

  ///
  enum EVAL_TOKEN
  {
    ADD = '+', SUB = '-', MULT = '*', DIV = '/', LP ='(', RP = ')',
    COMMA = ',',
    NEG = 100, 
    AND, OR, NOT, GREATER, LESS, GREATEREQUAL, LESSEQUAL, EQUAL,
    CONSTANT, VARIABLE, FUNCTION, GLOBVAR, COEFF_FUNC, END, STRING,
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
  void Parse (istream & aist);
  /// define constant 
  void DefineConstant (const char * name, double val);
  /// define constant 
  void DefineGlobalVariable (const char * name, double * var);

  /// evaluate function
  double Eval (const double * x = NULL) const;
  /// evaluate multi-value function
  void Eval (const double * x, double * y, int ydim) const;
  /// is expression a constant ?
  bool IsConstant () const;

  /// push constant on stack. Used for external parser.
  void AddConstant (double val);
  /// push variable x[varnum-1]. Used for external parser.
  void AddVariable (int varnum);
  /// push pointer to global double value. Used for external parser.
  void AddGlobVariable (const double * dp);
  /// push operation. Used for external parser.
  void AddOperation (EVAL_TOKEN op);
  /// push function call. Used for external parser.
  void AddFunction (double (*fun) (double));


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
      /// the input argument number varnum
      int varnum;
      /// a pointer to a uniary function
      double (*fun) (double);
    }; 
    ///
    UNION_OP operand;
  };

  /// the evaluation sequence
  Array<step> program;

  const double eps;


  /// parsing expression (standard parsing grammer)
  void ParseExpression ();
  /// parsing expression (standard parsing grammer)
  void ParseSubExpression ();
  /// parsing expression (standard parsing grammer)
  void ParseTerm ();
  /// parsing expression (standard parsing grammer)
  void ParsePrimary ();

  /// parse from stream
  istream * ist;

  ///
  EVAL_TOKEN token;
  ///
  double num_value;
  ///
  char string_value[1000];
  ///
  char var_num;
  ///
  double * globvar;
 
  typedef double(*TFUNP) (double);
  /// registerd functions
  static SymbolTable<TFUNP> functions;

  /// registerd constants
  SymbolTable<double> constants;

  /// registerd variables
  SymbolTable<double*> globvariables;

  /// returns last token
  EVAL_TOKEN GetToken() const
    { return token; }

  /// returns num_value of last token
  double GetNumValue() const
    { return num_value; }

  /// returns variable number of last token
  int GetVariableNumber() const
    { return var_num; }

  /// returns identifier of last token
  const char * GetStringValue() const
    { return string_value; }
  
  /// read next token
  void ReadNext();
};
}


#endif


