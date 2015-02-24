#include <solve.hpp>
using namespace ngsolve;

class NumProcDemo : public NumProc
{

public:
    
  NumProcDemo (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    ;
  }

  virtual void Do(LocalHeap & lh)
  {
    ;
  }


  virtual string GetClassName () const
  {
    return "NumProcDemo";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc Parabolic:\n"                
	<< endl;
  }
};


static RegisterNumProc<NumProcDemo> npinit1("demo");

