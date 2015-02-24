#include <solve.hpp>
using namespace ngsolve;


class NumProcTEMPLATE : public NumProc
{
protected:
  // BilinearForm * bfa;
  // LinearForm * lff;
  // GridFunction * gfu;
  // double dt;

public:

  NumProcTEMPLATE (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearforma", "a"));
    // lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", "f"));
    // gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));
    // dt = flags.GetNumFlag ("dt", 0.001);
  }


  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "solve TEMPLATE" << endl;
  }


  virtual string GetClassName () const
  {
    return "Numproc TEMPLATE";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc TEMPLATE:\n" 
	<< endl;
  }
};



// register the numproc 'TEMPLATE' 
static RegisterNumProc<NumProcTEMPLATE> npinit1("TEMPLATE");

