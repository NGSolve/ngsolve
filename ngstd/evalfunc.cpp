/**************************************************************************/
/* File:   evalfunc.cpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 95                                                    */
/**************************************************************************/

/* 
   Function parser
*/



#include <ngstd.hpp>


namespace ngstd
{
  using std::sqrt;   // for icc 11.0 ???

#include "bessel.hpp"

  using namespace std;
  using namespace ngstd;

  SymbolTable<EvalFunction::TFUNP> EvalFunction::functions;


  EvalFunction :: EvalFunction () : eps(1e-14)
  {
  }

  EvalFunction :: EvalFunction (istream & aist) : eps(1e-14)
  {
    Parse (aist);
  }


  EvalFunction :: EvalFunction (const string & str) : eps(1e-14)
  {
    stringstream strstr(str);
    Parse (strstr);
  }

  EvalFunction :: EvalFunction (const EvalFunction & eval2) : eps(1e-14)
  {
    program = eval2.program;
    constants = eval2.constants;
    globvariables = eval2.globvariables;
  }

  EvalFunction :: ~EvalFunction ()
  {
    ;
  }

  void EvalFunction :: Parse (istream & aist)
  {
    ist = &aist;
    ReadNext();
    ParseExpression ();
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


  void EvalFunction :: DefineConstant (const char * name, double val)
  {
    constants.Set (name, val);
  }

  void EvalFunction :: DefineGlobalVariable (const char * name, double * var)
  {
    globvariables.Set (name, var);
  }

  void EvalFunction :: Eval (const double * x, double * y, int ydim) const
  {
    ArrayMem<double, 100> stack(program.Size());

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

	  case NEG:
	    stack[stacksize] = -stack[stacksize];
	    break;

	  case AND:
	    if(stack[stacksize-1] > eps && stack[stacksize] > eps)
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;
	    
	  case OR:
	    if(stack[stacksize-1] > eps || stack[stacksize] > eps)
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;
	    
	  case NOT:
	    if(stack[stacksize] > eps)
	      stack[stacksize] = 0;
	    else
	      stack[stacksize] = 1;
	    break;

	  case GREATER:
	    if(stack[stacksize-1] > stack[stacksize])
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case GREATEREQUAL:
	    if(stack[stacksize-1] >= stack[stacksize])
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case EQUAL:
	    if(std::fabs(stack[stacksize-1] - stack[stacksize]) < eps)
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case LESSEQUAL:
	    if(stack[stacksize-1] <= stack[stacksize])
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case LESS:
	    if(stack[stacksize-1] < stack[stacksize])
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
	    stacksize++;
	    stack[stacksize] = x[program[i].operand.varnum];
	    break;

	  case GLOBVAR:
	    stacksize++;
	    stack[stacksize] = *program[i].operand.globvar;
	    break;

	  case IMAG:
	    throw Exception ("Real Eval called for Complex EvalFunction");

	  case FUNCTION:
	    stack[stacksize] = (*program[i].operand.fun) (stack[stacksize]);
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
	    stack[stacksize] = atan (stack[stacksize]);
	    break;
	  case ATAN2:
	    stack[stacksize-1] = atan2(stack[stacksize-1],
				       stack[stacksize]);
	    stacksize--;
	    break;
	  case EXP:
	    stack[stacksize] = exp (stack[stacksize]);
	    break;
	  case LOG:
	    stack[stacksize] = log (stack[stacksize]);
	    break;
	  case ABS:
	    stack[stacksize] = std::fabs (stack[stacksize]);
	    break;
	  case SIGN:
	    if(stack[stacksize] > 0)
	      stack[stacksize] = 1;
	    else if(stack[stacksize] < 0)
	      stack[stacksize] = -1;
	    else
	      stack[stacksize] = 0;
	    break;
	  case SQRT:
	    stack[stacksize] = sqrt (stack[stacksize]);
	    break;
	  case STEP:
	    stack[stacksize] = (stack[stacksize] >= 0) ? 1 : 0;
	    break;
	    
	  case COMMA:
	    break;
	    /*
	  case BESSELJ0:
	    stack[stacksize] = bessj0 (stack[stacksize]);
	    break;
	  case BESSELJ1:
	    stack[stacksize] = bessj1 (stack[stacksize]);
	    break;
	  case BESSELY0:
	    stack[stacksize] = bessy0 (stack[stacksize]);
	    break;
	  case BESSELY1:
	    stack[stacksize] = bessy1 (stack[stacksize]);
	    break;
	    */

	  default:
	    cerr << "undefined operation for EvalFunction" << endl;
	  }
      }

    if (stacksize != ydim-1)
      {
	cout << "final stacksize not matching ydim" << endl;
	return;
      }

    for (int i = 0; i < ydim; i++)
      y[i] = stack[i];
  }



  void EvalFunction :: Eval (const complex<double> * x, complex<double> * y, int ydim) const
  {
    ArrayMem<complex<double>, 100> stack(program.Size());

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

	  case NEG:
	    stack[stacksize] = -stack[stacksize];
	    break;

	    /*
	  case AND:
	    if(stack[stacksize-1] > eps && stack[stacksize] > eps)
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;
	    
	  case OR:
	    if(stack[stacksize-1] > eps || stack[stacksize] > eps)
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;
	    
	  case NOT:
	    if(stack[stacksize] > eps)
	      stack[stacksize] = 0;
	    else
	      stack[stacksize] = 1;
	    break;

	  case GREATER:
	    if(stack[stacksize-1] > stack[stacksize])
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case GREATEREQUAL:
	    if(stack[stacksize-1] >= stack[stacksize])
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case EQUAL:
	    if(std::fabs(stack[stacksize-1] - stack[stacksize]) < eps)
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case LESSEQUAL:
	    if(stack[stacksize-1] <= stack[stacksize])
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;

	  case LESS:
	    if(stack[stacksize-1] < stack[stacksize])
	      stack[stacksize-1] = 1;
	    else
	      stack[stacksize-1] = 0;
	    stacksize--;
	    break;
	    */
	  case CONSTANT:
	    stacksize++;
	    stack[stacksize] = program[i].operand.val;
	    break;

	  case VARIABLE:
	    stacksize++;
	    stack[stacksize] = x[program[i].operand.varnum];
	    break;

	  case GLOBVAR:
	    stacksize++;
	    stack[stacksize] = *program[i].operand.globvar;
	    break;

	  case IMAG:
	    stacksize++;
	    stack[stacksize] = complex<double> (0, 1);
	    break;


	    /*
	  case FUNCTION:
	    stack[stacksize] = (*program[i].operand.fun) (stack[stacksize]);
	    break;
	    */

	  case SIN:
	    stack[stacksize] = sin (stack[stacksize]);
	    break;
	  case COS:
	    stack[stacksize] = cos (stack[stacksize]);
	    break;
	  case TAN:
	    stack[stacksize] = tan (stack[stacksize]);
	    break;
	    /*
	  case ATAN:
	    stack[stacksize] = atan (stack[stacksize]);
	    break;
	  case ATAN2:
	    stack[stacksize-1] = atan2(stack[stacksize-1],
				       stack[stacksize]);
	    stacksize--;
	    break;
	    */
	  case EXP:
	    stack[stacksize] = exp (stack[stacksize]);
	    break;
	  case LOG:
	    stack[stacksize] = log (stack[stacksize]);
	    break;
	    /*
	  case ABS:
	    stack[stacksize] = std::fabs (stack[stacksize]);
	    break;
	    */
	    /*
	  case SIGN:
	    if(stack[stacksize] > 0)
	      stack[stacksize] = 1;
	    else if(stack[stacksize] < 0)
	      stack[stacksize] = -1;
	    else
	      stack[stacksize] = 0;
	    break;
	    */
	  case SQRT:
	    stack[stacksize] = sqrt (stack[stacksize]);
	    break;
	    /*
	  case STEP:
	    stack[stacksize] = (stack[stacksize] >= 0) ? 1 : 0;
	    break;
	    */
	  case COMMA:
	    break;
	    /*
	  case BESSELJ0:
	    stack[stacksize] = bessj0 (stack[stacksize]);
	    break;
	  case BESSELJ1:
	    stack[stacksize] = bessj1 (stack[stacksize]);
	    break;
	  case BESSELY0:
	    stack[stacksize] = bessy0 (stack[stacksize]);
	    break;
	  case BESSELY1:
	    stack[stacksize] = bessy1 (stack[stacksize]);
	    break;
	    */

	  default:
	    cerr << "undefined operation for EvalFunction" << endl;
	  }
      }

    if (stacksize != ydim-1)
      {
	cout << "final stacksize not matching ydim" << endl;
	return;
      }

    for (int i = 0; i < ydim; i++)
      y[i] = stack[i];
  }








  bool EvalFunction :: IsConstant () const
  {
    for (int i = 0; i < program.Size(); i++)
      {
	EVAL_TOKEN op = program[i].op;

	if (op == VARIABLE || op == GLOBVAR || op == COEFF_FUNC  || op == COMMA)
	  return false;
      }
    return true;
  }

  bool EvalFunction :: IsComplex () const
  {
    cout << "evalfunction::iscomplex:" << endl;
    for (int i = 0; i < program.Size(); i++)
      {
	EVAL_TOKEN op = program[i].op;
	if (op == IMAG) return true;
      }
    cout << "return false" << endl;
    return false;
  }

  int EvalFunction :: Dimension() const
  {
    int dim = 1;
    for (int i = 0; i < program.Size(); i++)
      if (program[i].op == ',') dim++;

    return dim;
  }


  void EvalFunction :: Print (ostream & ost) const
  {
    for (int i = 0; i < program.Size(); i++)
      {
	EVAL_TOKEN op = program[i].op;
	ost << "Step " << i << ": " << (int)op << " = " << (char) op 
	    << ", val = " << program[i].operand.val << endl;
      }
  }





  void EvalFunction :: ParseExpression ()
  {
    ParseExpression2 ();

    while (1)
      {
	//      cout << "parseexpr, goken = " << GetToken() << endl;
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
	    return;
	  }
      }
  }




  void EvalFunction :: ParseExpression2 ()
  {
    ParseSubExpression ();

    while (1)
      {
	//      cout << "parseexpr, goken = " << GetToken() << endl;
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
	    return;
	  }
      }
  }
  


  void EvalFunction :: ParseSubExpression ()
  {
    ParseTerm ();

    while (1)
      {
	//      cout << "parseexpr, goken = " << GetToken() << endl;
	switch (GetToken())
	  {
	  case ADD:
	    {
	      ReadNext();
	      ParseTerm ();
	      AddOperation (ADD);
	      break;
	    }
	  case SUB:
	    {
	      ReadNext();
	      ParseTerm ();
	      AddOperation (SUB);
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
	    return;
	  }
      }
  }


  void EvalFunction :: ParseTerm ()
  {
    ParsePrimary();
  
    while (1)
      {
	//      cout << "parseterm, goken = " << GetToken() << endl;
	switch (GetToken())
	  {
	  case MULT:
	    {
	      ReadNext();
	      ParsePrimary();
	      AddOperation (MULT);
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
	    return;
	  }
      }
  }



  void EvalFunction :: ParsePrimary()
  {
    switch (GetToken())
      {
      case CONSTANT:
	{
	  ReadNext();
	  AddConstant (GetNumValue());
	  break;
	}
      case SUB:
	{
	  ReadNext();
	  ParsePrimary();
	  AddConstant (-1);
	  AddOperation (MULT);
	  break;
	}
      case LP:
	{
	  ReadNext();
	  ParseExpression();
	  ReadNext();
	  break;
	}
      case VARIABLE:
	{
	  ReadNext();
	  AddVariable (GetVariableNumber());
	  break;
	}
      case GLOBVAR:
	{
	  ReadNext();
	  AddGlobVariable (globvar);
	  break;
	}
      case IMAG:
	{
	  ReadNext();
	  AddOperation (IMAG);
	  break;
	}
      case FUNCTION:
	{
	  ReadNext();
	  double (*funp)(double) = functions[string_value];
	  ParsePrimary();
	  AddFunction (funp);
	  break;
	}
      case SIN:
      case COS:
      case TAN:
      case ATAN:
      case EXP:
      case LOG:
      case ABS:
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
	  ParsePrimary();
	  AddOperation (op);
	  break;
	}
      case ATAN2:
	{
	  EVAL_TOKEN op = GetToken();
	  ReadNext();
	  ReadNext();
	  ParseExpression();
	  ReadNext();
	  ParseExpression();
	  ReadNext(); 
	  AddOperation (op);
	  break;
	}
      default:
        return;
      }
  }



  void EvalFunction :: ReadNext ()
  {
    char ch;

    // skip whitespaces
    do
      {
	(*ist).get(ch);
	if ((*ist).eof())
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
	      while ((*ist) && (isalnum (ch) || ch == '_' || ch == '>' || ch == '<' || ch == '='))
		{
		  string_value[cnt] = ch;
		  cnt++;
		  (*ist).get(ch);
		}
	      (*ist).putback (ch);
	      string_value[cnt] = 0;

	      //	      cout << "parse string " << string_value << endl;


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
		  //		  cout << "scanner found constant" << endl;
		  token = CONSTANT;
		  num_value = constants[string_value];
		  return;
		}

	      if (globvariables.Used (string_value))
		{
		  token = GLOBVAR;
		  globvar = globvariables[string_value];
		  //		  cout << "scanner found glob. variable: " << *globvar << endl;
		  return;
		}
	      //	    cout << "found string " << string_value << endl;

	      //	    (*ist) >> string_value;
	      //	    cout << "string = " << string_value << endl;
	 
// 		if (keywords.Used (string_value))
// 		{
// 		token = keywords.Get (string_value);
// 		return;
// 		}
	

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

    if(token == STRING)
      cerr << "WARNING: Please check function, didn't know what to do with \"" << string_value << "\"" << endl;

    //  cout << "token = " << token << " = " << char(token) << " numval = " << num_value << endl;
  }








}






