/*
  
We create a simple dialogbox to interact with the solver

 */


#include <solve.hpp>
using namespace ngsolve;

#include <tcl.h>


static shared_ptr<PDE> global_pde;


int NGS_TclDemo (ClientData clientData,
                 Tcl_Interp * interp,
                 int argc, const char *argv[])
{
  cout << "argv[0] = " << argv[0] << endl;
  cout << "argv[1] = " << argv[1] << endl;

  Array<shared_ptr<EvalFunction>> ae(1); 
  ae[0] = make_shared<EvalFunction> (argv[1]);

  auto coef = make_shared<DomainVariableCoefficientFunction>(ae);

  LocalHeap lh(1000000, "tcldemo-heap");
  SetValues (coef, *global_pde->GetGridFunction("u"), VOL, NULL, lh);
  return TCL_OK;
}



class NumProcTclDemo : public NumProc
{
protected:
public:
  NumProcTclDemo (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    global_pde = apde;
    apde->Tcl_Eval ("source tcldemo.tcl");
    
    Tcl_CreateCommand (apde->GetTclInterpreter(), 
                       "NGS_TclDemo", NGS_TclDemo,
                       (ClientData)NULL,
                       (Tcl_CmdDeleteProc*) NULL);
  }

  virtual void Do(LocalHeap & lh)
  {
    ;
  }  

};

static RegisterNumProc<NumProcTclDemo> npinit1("tcldemo");
