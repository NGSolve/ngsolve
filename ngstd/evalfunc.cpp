/**************************************************************************/
/* File:   evalfunc.cpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 95                                                    */
/**************************************************************************/

/* 
   Function parser
*/



#include <ngstd.hpp>
#include "evalfunc.hpp"


namespace ngstd
{
  using std::sqrt;   // for icc 11.0 ???

#include "bessel.hpp"

  using namespace std;
  using namespace ngstd;

  SymbolTable<EvalFunction::TFUNP> EvalFunction::functions;


  EvalFunction :: EvalFunction () : eps(1e-14)
  {
    DefineConstant ("pi", M_PI);
    DefineArgument ("x", 0);
    DefineArgument ("y", 1);
    DefineArgument ("z", 2);
    num_arguments = 3;
  }

  EvalFunction :: EvalFunction (istream & aist) : eps(1e-14)
  {
    DefineConstant ("pi", M_PI);
    DefineArgument ("x", 0);
    DefineArgument ("y", 1);
    DefineArgument ("z", 2);
    num_arguments = 3;
    Parse (aist);
  }


  EvalFunction :: EvalFunction (const string & str) : eps(1e-14)
  {
    DefineConstant ("pi", M_PI);
    DefineArgument ("x", 0);
    DefineArgument ("y", 1);
    DefineArgument ("z", 2);
    num_arguments = 3;
    stringstream strstr(str);
    Parse (strstr);
  }

  EvalFunction :: EvalFunction (const EvalFunction & eval2) : eps(1e-14)
  {
    program = eval2.program;
    res_type = eval2.res_type;
    constants.Update(eval2.constants);
    globvariables.Update(eval2.globvariables);
    arguments.Update(eval2.arguments);
    num_arguments = eval2.num_arguments;
  }

  EvalFunction :: ~EvalFunction ()
  {
    ;
  }

  bool EvalFunction :: Parse (istream & aist)
  {
    ist = &aist;
    ReadNext();
    res_type = ParseExpression ();
    if (GetToken() != END) WriteBack();
    return program.Size() > 0;
  }


  void EvalFunction :: DefineConstant (const string & name, double val)
  {
    constants.Set (name, val);
  }

  void EvalFunction :: DefineGlobalVariable (const string & name, double * var)
  {
    globvariables.Set (name, var);
  }

  void EvalFunction :: DefineGlobalVariable (const string & name, GenericVariable * var)
  {
    genericvariables.Set (name, var);
  }

  void EvalFunction :: DefineArgument (const string & name, int num, int vecdim, bool iscomplex)
  {
    arguments.Set (name, argtype(num, vecdim, iscomplex));
  }





  double EvalFunction :: Eval (const double * x) const
  {
    double y;
    Eval (x, &y, 1);
    return y;
  }

  complex<double> EvalFunction :: Eval (const complex<double> * x) const
  {
    complex<double> y;
    Eval (x, &y, 1);
    return y;
  }

  void EvalFunction :: Eval (const double * x, double * y, int ydim) const
  {
    if (res_type.vecdim != ydim)
      {
	cout << "Eval called with ydim = " << ydim << ", but result.dim = " << res_type.vecdim << endl;
	return;
      }

    ArrayMem<double, 100> stack(program.Size());
    Eval<double,double> (x, &stack[0]);

    for (int i = 0; i < res_type.vecdim; i++)
      y[i] = stack[i];
  }

  void EvalFunction :: Eval (const complex<double> * x, complex<double> * y, int ydim) const
  {
    if (res_type.vecdim != ydim)
      {
	cout << "Eval complex called with ydim = " << ydim << ", but result.dim = " << res_type.vecdim << endl;
	return;
      }

    ArrayMem<complex<double>, 100> stack(program.Size());
    Eval<complex<double>, complex<double> > (x, &stack[0]);

    for (int i = 0; i < res_type.vecdim; i++)
      y[i] = stack[i];
  }

  void EvalFunction :: Eval (const complex<double> * x, double * y, int ydim) const
  {
    if (res_type.vecdim != ydim)
      {
	cout << "Eval complex/double called with ydim = " << ydim << ", but result.dim = " << res_type.vecdim << endl;
	return;
      }

    ArrayMem<complex<double>, 100> stack(program.Size());
    Eval<complex<double>, complex<double> > (x, &stack[0]);

    for (int i = 0; i < res_type.vecdim; i++)
      y[i] = stack[i].real();
  }

  template <typename SCAL> SCAL Imag ();

  template <> double Imag<double> ()
  {
    cerr << "IMAG used for real" << endl;
    return 0;
  }
  template <> complex<double> Imag<complex<double> > ()
  {
    return complex<double> (0,1); 
  }


  template <typename TIN, typename TCALC>
  void EvalFunction :: Eval (const TIN * x, TCALC * stack) const
  {
    int stacksize = -1;
    for (int i = 0; i < program.Size(); i++)
      {
	switch (program[i].op)
	  {
	  case ADD:
	    stack[stacksize-1] += stack[stacksize];
	    stacksize--;
	    break;

	  case SUB:
	    stack[stacksize-1] -= stack[stacksize];
	    stacksize--;
	    break;

	  case MULT:
	    stack[stacksize-1] *= stack[stacksize];
	    stacksize--;
	    break;

	  case DIV:
	    stack[stacksize-1] /= stack[stacksize];
	    stacksize--;
	    break;

	  case VEC_ADD:
	    {
	      int dim = program[i].vecdim;
	      for (int j = 0; j < dim; j++)
		stack[stacksize-2*dim+j+1] += stack[stacksize-dim+j+1];
	      stacksize -= dim;
	      break;
	    }

	  case VEC_SUB:
	    {
	      int dim = program[i].vecdim;
	      for (int j = 0; j < dim; j++)
		stack[stacksize-2*dim+j+1] -= stack[stacksize-dim+j+1];
	      stacksize -= dim;
	      break;
	    }

	  case SCAL_VEC_MULT:
	    {
	      int dim = program[i].vecdim;
	      TCALC scal = stack[stacksize-dim];
	      for (int j = 0; j < dim; j++)
		stack[stacksize-dim+j] = scal * stack[stacksize-dim+j+1];
	      stacksize--;
	      break;
	    }

	  case VEC_VEC_MULT:
	    {
	      int dim = program[i].vecdim;
	      TCALC scal = 0;
	      for (int j = 0; j < dim; j++)
		scal += stack[stacksize-2*dim+j+1] * stack[stacksize-dim+j+1];
	      stacksize-=2*dim-1;
              stack[stacksize] = scal;
	      break;
	    }

	  case VEC_ELEM:
	    {
	      int dim = program[i-1].vecdim;
              int index = int(CheckReal (stack[stacksize]));
              // cout << "vec_elem, dim = " << dim << ", index = " << index << endl;
              stack[stacksize-dim] = stack[stacksize-dim+index-1];
	      stacksize -= dim;
	      break;
	    }

	  case VEC_DIM:
	    {
	      int dim = program[i-1].vecdim;
              // cout << "dim = " << dim << endl;
	      stacksize -= dim-1;
              stack[stacksize]=dim;
	      break;
	    }

	  case NEG:
	    stack[stacksize] = -stack[stacksize];
	    break;

	  case AND:
	    if( ToBool (stack[stacksize-1]) && ToBool (stack[stacksize]) )
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;
	    
	  case OR:
	    if( ToBool(stack[stacksize-1]) || ToBool (stack[stacksize]) )
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;
	    
	  case NOT:
	    if( ToBool (stack[stacksize]) )
	      stack[stacksize] = 0;
	    else
	      stack[stacksize] = 1;
	    break;
	    
	  case GREATER:
	    if( CheckReal (stack[stacksize-1]) > CheckReal (stack[stacksize]))
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case GREATEREQUAL:
	    if( CheckReal (stack[stacksize-1]) >= CheckReal (stack[stacksize]) )
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case EQUAL:
	    if(Abs(stack[stacksize-1] - stack[stacksize]) < eps)
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case LESSEQUAL:
	    if( CheckReal (stack[stacksize-1]) <= CheckReal (stack[stacksize]) )
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case LESS:
	    if( CheckReal (stack[stacksize-1]) < CheckReal (stack[stacksize]) )
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case CONSTANT:
	    stacksize++;
	    stack[stacksize] = program[i].operand.val;
	    break;

	  case VARIABLE:
	    for (int j = 0; j < program[i].vecdim; j++)
	      {
		stacksize++;
		stack[stacksize] = x[program[i].operand.varnum+j];
	      }
	    break;

	  case GLOBVAR:
	    stacksize++;
	    stack[stacksize] = *program[i].operand.globvar;
	    break;

	  case GLOBGENVAR:
            for (int j = 0; j < program[i].operand.globgenvar->Dimension(); j++)
              {
                stacksize++;
                stack[stacksize] = program[i].operand.globgenvar->Value<TCALC>(j);
              }
	    break;

	  case IMAG:
	    // throw Exception ("Real Eval called for Complex EvalFunction\n");

	    stacksize++;
	    stack[stacksize] = Imag<TCALC>(); // complex<double> (0, 1);
	    break;


	  case FUNCTION:
	    stack[stacksize] = (*program[i].operand.fun) ( CheckReal (stack[stacksize]) );
	    break;

	  case SIN:
	    stack[stacksize] = sin (stack[stacksize]);
	    break;
	  case COS:
	    stack[stacksize] = cos (stack[stacksize]);
	    break;
	  case TAN:
	    stack[stacksize] = tan (stack[stacksize]);
	    break;
	  case ATAN:
	    stack[stacksize] = atan (stack[stacksize] );
	    break;
	  case ATAN2:
	    stack[stacksize-1] = atan2( CheckReal (stack[stacksize-1]),
					CheckReal (stack[stacksize]) );
	    stacksize--;
	    break;
	  case EXP:
	    stack[stacksize] = exp (stack[stacksize]);
	    break;
	  case LOG:
	    stack[stacksize] = log (stack[stacksize]);
	    break;
	  case ABS:
	    {
	      int dim = program[i].vecdim;
	      if (dim == 1)
		stack[stacksize] = Abs (stack[stacksize]);
	      else
		{
		  double sum = 0.0; 
		  for (int j = 0; j < dim; j++)
		    sum += sqr (Abs (stack[stacksize-j]));
		  stacksize -= dim-1;
		  stack[stacksize] = sqrt(sum);
		}
	      break;
	    }
	  case SIGN: 
	    if( CheckReal (stack[stacksize]) > 0)
	      stack[stacksize] = 1;
	    else if( CheckReal (stack[stacksize]) < 0)
	      stack[stacksize] = -1;
	    else
	      stack[stacksize] = 0;
	    break;
	  case SQRT:
	    stack[stacksize] = sqrt (stack[stacksize]);
	    break;
	  case STEP:
	    stack[stacksize] = ( CheckReal (stack[stacksize]) >= 0) ? 1 : 0;
	    break;
	    
	  case COMMA:
	    break;

	  case BESSELJ0:
	    stack[stacksize] = bessj0 ( CheckReal (stack[stacksize]) );
	    break;
	  case BESSELJ1:
	    stack[stacksize] = bessj1 ( CheckReal (stack[stacksize]) );
	    break;
	  case BESSELY0:
	    stack[stacksize] = bessy0 ( CheckReal (stack[stacksize]) );
	    break;
	  case BESSELY1:
	    stack[stacksize] = bessy1 ( CheckReal (stack[stacksize]) );
	    break;

	  default:
	    cerr << "undefined operation for EvalFunction" << endl;
	  }
      }
  }




  bool EvalFunction :: IsConstant () const
  {
    if (res_type.iscomplex) return false;
    if (res_type.vecdim > 1) return false;
    
    for (int i = 0; i < program.Size(); i++)
      {
	EVAL_TOKEN op = program[i].op;
	if (op == VARIABLE || op == GLOBVAR)
	  return false;
      }
    return true;
  }

  double EvalFunction :: EvalConstant () const
  {
    return Eval ((double*)NULL);
  }


  bool EvalFunction :: IsComplex () const
  {
    for (int i = 0; i < program.Size(); i++)
      {
	EVAL_TOKEN op = program[i].op;
	if (op == IMAG) return true;
      }
    for (int i = 0; i < arguments.Size(); i++)
      if (arguments[i].argnum != -1 && arguments[i].iscomplex)
	return true;
    return false;
  }

  /*
  int EvalFunction :: Dimension() const
  {
    int dim = 1;
    for (int i = 0; i < program.Size(); i++)
      if (program[i].op == ',') dim++;

    return dim;
  }
  */

  void EvalFunction :: Print (ostream & ost) const
  {
    for (int i = 0; i < program.Size(); i++)
      {
	EVAL_TOKEN op = program[i].op;
	ost << "Step " << i << ": " << (int)op << " = ";
	switch (op)
	  {
	  case CONSTANT: ost << " const, val = " << program[i].operand.val; break;
	  case VARIABLE: ost << " input var " << program[i].operand.varnum; break;
	  case GLOBGENVAR: ost << " global var " << *program[i].operand.globgenvar; break;
	  case SIN: ost << " sin"; break;
	  case COS: ost << " cos"; break;
	  case TAN: ost << " tan"; break;
	  case ATAN: ost << " atan"; break;
	  case ATAN2: ost << " atan2"; break;
	  case EXP: ost << " exp"; break;
	  case LOG: ost << " log"; break;
	  case ABS: ost << " abs"; break;
	  case SIGN: ost << " sign"; break;
	  case SQRT: ost << " sqrt"; break;
	  case STEP: ost << " step"; break;
	  default:
	    ost << (char) op;
	  }
	ost << " vdim = " << program[i].vecdim;
	ost << endl;
      }
  }





  EvalFunction::ResultType EvalFunction :: ParseCommaExpression ()
  {
    ResultType result = ParseExpression ();
    if (GetToken() == COMMA)
      {
	ReadNext();   // ','
	result = ParseCommaExpression ();
	result.vecdim++;
	AddOperation (COMMA);
      }
    return result;

    /*
    while (1)
      {
	switch (GetToken())
	  {
	  case COMMA:
	    {
	      ReadNext();
	      ParseExpression ();
	      AddOperation (COMMA);
	      break;
	    }
	  default:
	    return ResultType();
	  }
      }
    */
  }




  EvalFunction::ResultType EvalFunction :: ParseExpression ()
  {
    ResultType result = ParseSubExpression ();

    while (1)
      {
	switch (GetToken())
	  {
	  case GREATER:
	    {
	      ReadNext();
	      ParseSubExpression ();
	      AddOperation (GREATER);
	      break;
	    }
	  case GREATEREQUAL:
	    {
	      ReadNext();
	      ParseSubExpression ();
	      AddOperation (GREATEREQUAL);
	      break;
	    }
	  case EQUAL:
	    {
	      ReadNext();
	      ParseSubExpression ();
	      AddOperation (EQUAL);
	      break;
	    }
	  case LESSEQUAL:
	    {
	      ReadNext();
	      ParseSubExpression ();
	      AddOperation (LESSEQUAL);
	      break;
	    }
	  case LESS:
	    {
	      ReadNext();
	      ParseSubExpression ();
	      AddOperation (LESS);
	      break;
	    }
	  default:
            {
              return result;
            }
          }
      }
    return result;
  }
  


  EvalFunction::ResultType EvalFunction :: ParseSubExpression ()
  {
    ResultType result = ParseTerm ();

    while (1)
      {
	switch (GetToken())
	  {
	  case ADD:
	    {
	      ReadNext();   // '+'
	      ResultType result2 = ParseTerm ();
	      if (result.vecdim != result2.vecdim) cerr << "vec error" << endl;
	      result.iscomplex |= result2.iscomplex;
	      if (result.vecdim == 1)
		AddOperation (ADD);
	      else
		{
		  AddOperation (VEC_ADD);
		  program.Last().vecdim = result.vecdim;
		}
	      break;
	    }
	  case SUB:
	    {
	      ReadNext();   // '+'
	      ResultType result2 = ParseTerm ();
	      if (result.vecdim != result2.vecdim) cerr << "vec error" << endl;
	      result.iscomplex |= result2.iscomplex;
	      if (result.vecdim == 1)
		AddOperation (SUB);
	      else
		{
		  AddOperation (VEC_SUB);
		  program.Last().vecdim = result.vecdim;
		}
	      break;
	    }
	  case OR:
	    {
	      ReadNext();
	      ParseTerm ();
	      AddOperation (OR);
	      break;
	    }
	  default:
	    return result;
	  }
      }
  }


  EvalFunction::ResultType EvalFunction :: ParseTerm ()
  {
    ResultType result = ParsePrimary();
  
    while (1)
      {
	switch (GetToken())
	  {
	  case MULT:
	    {
	      ReadNext();    // '*'
	      ResultType result2 = ParsePrimary();
	      result.iscomplex |= result2.iscomplex;
	      if (result.vecdim == 1 && result2.vecdim == 1)
		{
		  AddOperation (MULT);
		}
	      else if (result.vecdim == 1 && result2.vecdim > 1)
		{
		  AddOperation (SCAL_VEC_MULT);
		  program.Last().vecdim = result2.vecdim;
		  result.vecdim = result2.vecdim;
		}
	      else if (result.vecdim > 1 && result2.vecdim > 1)
		{
		  AddOperation (VEC_VEC_MULT);
		  program.Last().vecdim = result.vecdim;
                  result.vecdim = 1;
		}
	      break;
	    }
	  case DIV:
	    {
	      ReadNext();
	      ParsePrimary();
	      AddOperation (DIV);
	      break;
	    }
	  case AND:
	    {
	      ReadNext();
	      ParsePrimary();
	      AddOperation (AND);
	      break;
	    }
	  default:
	    return result;
	  }
      }
  }


  
  EvalFunction::ResultType EvalFunction :: ParsePrimary()
  {
    ResultType result;

    switch (GetToken())
      {
      case CONSTANT:
	{
	  AddConstant (GetNumValue());
	  ReadNext();
	  break;
	}
      case SUB:
	{
	  ReadNext();   // '-'
	  AddConstant (-1);
	  result = ParsePrimary();
          if (result.vecdim == 1)
            AddOperation (MULT);
          else
            {
              AddOperation (SCAL_VEC_MULT);
              program.Last().vecdim = result.vecdim;
            }
	  break;
	}
      case LP:
	{
	  ReadNext();  // '('
	  result = ParseCommaExpression();
	  ReadNext();  // ')'
	  break;
	}
      case VARIABLE:
	{
	  AddVariable (GetVariableNumber());
	  program.Last().vecdim = GetVariableDimension();
	  result.vecdim = GetVariableDimension();
	  result.iscomplex = GetVariableIsComplex();
	  ReadNext();
	  break;
	}
      case GLOBVAR:
	{
	  AddGlobVariable (globvar);
	  ReadNext();
	  break;
	}
      case GLOBGENVAR:
	{
	  ReadNext();
	  AddGlobVariable (globgenvar);
          program.Last().vecdim = globgenvar -> Dimension();
          result.iscomplex = globgenvar -> IsComplex();
          result.vecdim = globgenvar -> Dimension();

          if (GetToken() == '(')
            {
              // cout << "array access !!!" << endl;
              /* ResultType result2 = */ 
              ParsePrimary();
              program.Last().vecdim = result.vecdim;
              AddOperation(VEC_ELEM);
              result.vecdim = 1;
            }

	  break;
	}
        
      case VEC_DIM:
        {
          ReadNext();
          ParsePrimary();
          AddOperation(VEC_DIM);
          program.Last().vecdim = result.vecdim;
          break;
        }

      case IMAG:
	{
	  ReadNext();
	  AddOperation (IMAG);
	  result.iscomplex = true;
	  break;
	}
      case FUNCTION:
	{
	  ReadNext();
	  double (*funp)(double) = functions[string_value];
	  result = ParsePrimary();
	  AddFunction (funp);
	  break;
	}
      case SIN:
      case COS:
      case TAN:
      case ATAN:
      case EXP:
      case LOG:
      case SIGN:
      case SQRT:
      case STEP:
      case BESSELJ0:
      case BESSELJ1:
      case BESSELY0:
      case BESSELY1:
	{
	  EVAL_TOKEN op = GetToken();
	  ReadNext();
	  result = ParsePrimary();
	  AddOperation (op);
	  break;
	}
      case ATAN2:  // needs 2 arguments
	{
	  EVAL_TOKEN op = GetToken();
	  ReadNext();
	  ParsePrimary();       // a COMMA expr is an expression 
	  AddOperation (op);
	  break;

	  /*
	  EVAL_TOKEN op = GetToken();
	  ReadNext();
	  ReadNext();      
	  ParseExpression();
	  ReadNext();        //  ','
	  ParseExpression();
	  ReadNext();        //  ')'
	  AddOperation (op);
	  break;
	  */
	}
      case ABS:
	{
	  EVAL_TOKEN op = GetToken();
	  ReadNext();
	  result = ParsePrimary();
	  result.iscomplex = false;
	  AddOperation (op);
	  program.Last().vecdim = GetVariableDimension();	  
	  result.vecdim = 1;
	  break; 
	}
      default:
	cout << "EvalFunction: why did I get here  ???" << endl;
      }
    return result;
  }


  void EvalFunction :: WriteBack ()
  {
    ist -> seekg(lastpos);
  }

  void EvalFunction :: ReadNext (bool optional)
  {
    lastpos = ist -> tellg();
    char ch;

    // skip whitespaces
    do
      {
	if (!ist->good() || ist->eof())
	  {
	    token = END;
	    return;
	  }
	(*ist).get(ch);
	if (!ist->good() || ist->eof())
	  {
	    token = END;
	    return;
	  }
      }
    while (isspace(ch));
    
    switch (ch)
      {
      case '*': case '/': 
      case '+': case '-': 
      case '(': case ')':
      case ',':
	{
	  token = EVAL_TOKEN (ch);

          // check for double --
          if (ch == '-')
            {
              char ch2;
              ist -> get(ch2);
              ist -> putback(ch2);
              if (ch2 == '-')
                {
                  token = END;
                  return;
                }
            }
	  break;
	}
      
      default:
	{
	  if (isdigit (ch) || ch == '.')
	    {
	      (*ist).putback (ch);
	      (*ist) >> num_value;
	      token = CONSTANT;
	    }
	  else
	    {
	      int cnt = 0;
	      while ((*ist) && ((isalnum (ch) || ch == '_' || ch == '>' || ch == '<' || 
				 ch == '=') || ch == '.') )
		{
		  string_value[cnt] = ch;
		  cnt++;
		  (*ist).get(ch);
		}
	      (*ist).putback (ch);
	      string_value[cnt] = 0;

	      if (strcmp (string_value, "and") == 0)
		{
		  token = AND;
		  return;
		}

	      if (strcmp (string_value, "or") == 0)
		{
		  token = OR;
		  return;
		}

	      if (strcmp (string_value, "not") == 0)
		{
		  token = NOT;
		  return;
		}

	      if (strcmp (string_value, ">") == 0)
		{
		  token = GREATER;
		  return;
		}

	      if (strcmp (string_value, ">=") == 0)
		{
		  token = GREATEREQUAL;
		  return;
		}

	      if (strcmp (string_value, "=") == 0)
		{
		  token = EQUAL;
		  return;
		}

	      if (strcmp (string_value, "<=") == 0)
		{
		  token = LESSEQUAL;
		  return;
		}

	      if (strcmp (string_value, "<") == 0)
		{
		  token = LESS;
		  return;
		}

	      if (strcmp (string_value, "sin") == 0)
		{
		  token = SIN;
		  return;
		}

	      if (strcmp (string_value, "cos") == 0)
		{
		  token = COS;
		  return;
		}

	      if (strcmp (string_value, "tan") == 0)
		{
		  token = TAN;
		  return;
		}

	      if (strcmp (string_value, "atan") == 0)
		{
		  token = ATAN;
		  return;
		}

	      if (strcmp (string_value, "atan2") == 0)
		{
		  token = ATAN2;
		  return;
		}

	      if (strcmp (string_value, "exp") == 0)
		{
		  token = EXP;
		  return;
		}

	      if (strcmp (string_value, "log") == 0)
		{
		  token = LOG;
		  return;
		}

	      if (strcmp (string_value, "abs") == 0)
		{
		  token = ABS;
		  return;
		}

	      if (strcmp (string_value, "sign") == 0)
		{
		  token = SIGN;
		  return;
		}

	      if (strcmp (string_value, "sqrt") == 0)
		{
		  token = SQRT;
		  return;
		}

	      if (strcmp (string_value, "step") == 0)
		{
		  token = STEP;
		  return;
		}

	      if (strcmp (string_value, "besselj0") == 0)
		{
		  token = BESSELJ0;
		  return;
		}

	      if (strcmp (string_value, "besselj1") == 0)
		{
		  token = BESSELJ1;
		  return;
		}

	      if (strcmp (string_value, "bessely0") == 0)
		{
		  token = BESSELY0;
		  return;
		}

	      if (strcmp (string_value, "bessely1") == 0)
		{
		  token = BESSELY1;
		  return;
		}

	      if (strcmp (string_value, "dim") == 0)
		{
		  token = VEC_DIM;
		  return;
		}

	      if (strcmp (string_value, "I") == 0)
		{
		  token = IMAG;
		  return;
		}

	      if (functions.Used (string_value))
		{
		  token = FUNCTION;
		  return;
		}

	      if (constants.Used (string_value))
		{
		  token = CONSTANT;
		  num_value = constants[string_value];
		  return;
		}

	      if (globvariables.Used (string_value))
		{
		  token = GLOBVAR;
		  globvar = globvariables[string_value];
		  return;
		}

	      if (genericvariables.Used (string_value))
		{
		  token = GLOBGENVAR;
		  globgenvar = genericvariables[string_value];
		  return;
		}

	      //	    (*ist) >> string_value;
	      //	    cout << "string = " << string_value << endl;
	 
// 		if (keywords.Used (string_value))
// 		{
// 		token = keywords.Get (string_value);
// 		return;
// 		}
	

	      if (arguments.Used (string_value))
		{
		  var_num = arguments[string_value].argnum;
		  var_dim = arguments[string_value].dim;
		  var_iscomplex = arguments[string_value].iscomplex;
		  if (var_num == -1)
		    {
		      var_num = arguments[string_value].argnum = num_arguments;
		      num_arguments += arguments[string_value].dim;
		      /*
		      cout << "argument " << string_value 
			   << " becomes arg " << var_num 
			   << " vecdim = " << arguments[string_value].dim
			   << " complex = " << arguments[string_value].iscomplex
			   << endl;
		      */
		    }
		  token = VARIABLE;
		  return;
		}

	      /*
	      if (strcmp (string_value, "x1") == 0 ||
		  strcmp (string_value, "x") == 0)
		{
		  var_num = 0;
		  token = VARIABLE;
		  return;
		}
		  
	      if (strcmp (string_value, "x2") == 0 ||
		  strcmp (string_value, "y") == 0)
		{
		  var_num = 1;
		  token = VARIABLE;
		  return;
		}
		  
	      if (strcmp (string_value, "x3") == 0 ||
		  strcmp (string_value, "z") == 0)
		{
		  var_num = 2;
		  token = VARIABLE;
		  return;
		}
	      */  

	      /*	    
	      if (cnt == 2 &&
		  string_value[0] == 'x' && 
		  isdigit (string_value[1]))
		{
		  //		  cout << "found coordinate" << endl;
		  var_num = atoi (string_value+1)-1;
		  token = VARIABLE;
		}
	      else
	      */
		token = STRING;
	    }
	}
      }

    if(!optional && token == STRING)
      cerr << "WARNING: Please check function, didn't know what to do with \"" << string_value << "\"" << endl;

    //  cout << "token = " << token << " = " << char(token) << " numval = " << num_value << endl;
  }








}






