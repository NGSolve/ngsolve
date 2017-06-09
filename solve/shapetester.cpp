
#include <solve.hpp>

#include <tcl.h>
#if TCL_MAJOR_VERSION==8 && TCL_MINOR_VERSION>=4
#define tcl_const const
#else
#define tcl_const
#endif



namespace ngsolve
{
  int NGS_DrawShape (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[]);

  class NumProcShapeTester;


  
  static NumProcShapeTester * shapetester;
  ///
  class NumProcShapeTester : public NumProc
  {
  protected:
    ///
    shared_ptr<GridFunction> gfu;
    ///
    BilinearForm * bfa;
    ///
    int dof;
  public:
    ///
    NumProcShapeTester (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    {
      gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
      dof = int(flags.GetNumFlag ("dof", 0));
      apde->Tcl_Eval (
        "set w .shapetester;"
        "toplevel $w;" 
        "wm withdraw $w\n"
        "wm geom $w +100+100;"
        "wm deiconify $w;"
        "wm title $w \"Shape Tester\"\n"
        "set dofnr 0\n"
        "ttk::frame $w.frame;"
        "ttk::label $w.frame.l -text \"Dof number\"\n;"
        "ttk::spinbox $w.frame.dofnr -from 0 -to 1e9 -textvariable dofnr -width 5 -command { NGS_DrawShape $dofnr }\n"
        "pack $w.frame.dofnr $w.frame.l -side left -anchor w;"
        "pack $w.frame -padx 10 -pady 10\n"
        "focus .options_dlg\n"
        );
      Tcl_CreateCommand (apde->GetTclInterpreter(), 
            "NGS_DrawShape", NGS_DrawShape,
            (ClientData)NULL,
            (Tcl_CmdDeleteProc*) NULL);
      shapetester = this;
    }

    ///
    virtual ~NumProcShapeTester() { ; }
    
    void SetDof (int adof)
    {
      dof = adof;
      DoIt();
    }
    
    ///
    virtual void Do(LocalHeap & lh)
    {
      DoIt();
    }

    void DoIt ()
    {
      BaseVector & vecu = gfu->GetVector();
      vecu = 0;
      Array<int> dnums(1);
      if (dof >= vecu.Size())
	dof = vecu.Size()-1;
      dnums[0] = dof;
       
      Vector<> elu(1);
      elu(0) = 1;
      vecu.SetIndirect (dnums, elu);
      //gfu->Visualize(gfu->GetName());
      Ng_Redraw ();
    }

    ///
    virtual string GetClassName () const
    {
      return "Shape tester";
    }

    virtual void PrintReport (ostream & ost) const
    {
      ;
    }

    static void PrintDoc (ostream & ost)
    {
      ;
    }
  };



  int NGS_DrawShape (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
    cout << "draw shape nr " << argv[1] << endl;
    shapetester->SetDof (atoi(argv[1]));
    return TCL_OK;
  }


  static RegisterNumProc<NumProcShapeTester> npinitshapetester("shapetester");
}
