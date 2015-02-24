/*
  
We create a simple dialogbox to interact with the solver

 */


#include <solve.hpp>
using namespace ngsolve;

#include <tcl.h>


static PDE * global_pde;


int NGS_TclDemo (ClientData clientData,
                 Tcl_Interp * interp,
                 int argc, const char *argv[])
{
  cout << "argv[0] = " << argv[0] << endl;
  cout << "argv[1] = " << argv[1] << endl;

  cout << "parse function ... " << flush;
  EvalFunction eval(argv[1]);
  cout << "success" << endl;
  Array<EvalFunction*> ae(1); 
  ae[0] = &eval;

  DomainVariableCoefficientFunction<2> coef(ae);

  LocalHeap lh(1000000, "tcldemo-heap");
  SetValues (coef, *global_pde->GetGridFunction("u"), VOL, NULL, lh);
  return TCL_OK;
}



class NumProcTclDemo : public NumProc
{
protected:
public:
  NumProcTclDemo (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    global_pde = &pde;
    pde.Tcl_Eval ("source tcldemo.tcl");

    Tcl_CreateCommand (pde.GetTclInterpreter(), 
                       "NGS_TclDemo", NGS_TclDemo,
                       (ClientData)NULL,
                       (Tcl_CmdDeleteProc*) NULL);
  }

  virtual void Do(LocalHeap & lh)
  {
    ;
  }  

};





void TclCallBack ()
{
  ;
}


static RegisterNumProc<NumProcTclDemo> npinit1("tcldemo");
