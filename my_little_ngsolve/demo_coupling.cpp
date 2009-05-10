/*
  

In this demo, we take the solution of one PDE, and use it as right
hand side of a second one.

For simplicity, we want to solve the first PDE,

-\Delta u = f

the second one is

-\Delta w = u


Similar problems occur when coupling different physical fields.

Please include this file to the src files given in netgen/ngsolve/Makefile
*/



#include <solve.hpp>

using namespace ngsolve;



/*
  Every solver is a class derived from the class NumProc.
  It collects objects (such as bilinear-forms, gridfunctions) and parameters.
*/

class NumProcCouplingDemo : public NumProc
{
protected:
  // grid function provides the solution vector of the first PDE
  S_GridFunction<double> * gfu;

  // linear-form providing the right hand side for the second PDE
  S_LinearForm<double> * lff;

public:
    
  /*
    In the constructor, the solver class gets the flags from the pde - input file.
    the PDE class apde constains all bilinear-forms, etc...
  */
  NumProcCouplingDemo (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the bilinear-forms for the stiffness and for the mass-term
    // like  "-linearform=f". Default arguments are 'f' and 'u'

    // we know that we work with real-valued gridfunctions:
    lff = dynamic_cast<S_LinearForm<double> *> 
      (pde.GetLinearForm (flags.GetStringFlag ("linearform", "f")));
    gfu = dynamic_cast<S_GridFunction<double> *> 
      (pde.GetGridFunction (flags.GetStringFlag ("gridfunction", "u")));
  }

  virtual ~NumProcCouplingDemo() 
  { ; }


  // creates an solver object
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NumProcCouplingDemo (pde, flags);
  }

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Compute coupling terms" << endl;
      
    // reference to the matrices provided by the bi-forms.
    // will be of type SparseSymmetricMatrix<double> for scalar problems

    BaseVector & vecf = lff->GetVector();
    BaseVector & vecu = gfu->GetVector();

    const FESpace & fesu = gfu -> GetFESpace();
    const FESpace & fesf = lff -> GetFESpace();

    Array<int> dnumsu, dnumsf;  // the dof-numbes for u and f
    ElementTransformation eltrans;

    lff -> GetVector() = 0.0;

    int ne = ma.GetNE();
    for (int i = 0; i < ne; i++)   // loop over elements
      {  
        HeapReset hr(lh);    // reset the local heap memory at the end of the loop
		
        ma.GetElementTransformation (i, eltrans, lh);

        const ScalarFiniteElement<2> & felu = dynamic_cast<const ScalarFiniteElement<2>&> (fesu.GetFE (i, lh));
        const ScalarFiniteElement<2> & felf = dynamic_cast<const ScalarFiniteElement<2>&> (fesf.GetFE (i, lh));
		      
        fesu.GetDofNrs (i, dnumsu);
        fesf.GetDofNrs (i, dnumsf);

        FlatVector<> elu (dnumsu.Size(), lh);
        FlatVector<> helf (dnumsf.Size(), lh);
        FlatVector<> elf (dnumsf.Size(), lh);

        gfu -> GetElementVector (dnumsu, elu);

        const IntegrationRule & ir = 
          SelectIntegrationRule (eltrans.GetElementType(), felu.Order()+felf.Order() );

        
        elf = 0.0;
        for (int j = 0; j < ir.GetNIP(); j++)   // integration points
          {
            SpecificIntegrationPoint<2, 2> sip(ir[j], eltrans, lh);  // computes Jacobi matrix etc

            Vec<1> ui;   // value of u in point
            DiffOpId<2>::Apply (felu, sip, elu, ui, lh);   // compute value in point
            
            // could use also other differential operators such as
            // DiffOpGradient<2>, DiffOpCurl<2>, ....


            Vec<1> fi;
            fi = ir[j].Weight() * fabs (sip.GetJacobiDet())  * ui;   // local math

            DiffOpId<1>::ApplyTrans (felf, sip, fi, helf, lh);   // coefficient times test functions
            elf += helf;
          }

        lff -> AddElementVector (dnumsf, elf);
      }

    // *testout << "lff = " << endl << lff->GetVector() << endl;
  }






  virtual string GetClassName () const
  {
    return "Coupling (Demo)";
  }

  virtual void PrintReport (ostream & ost)
  {
    ost << GetClassName() << endl
	<< "Linear-form     = " << lff->GetName() << endl
	<< "Gridfunction    = " << gfu->GetName() << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc Coupling - Demo:\n"                                  
      "------------------------\n"                                      
      "Takes the solution of on PDE as right hand side of a second PDE\n" 
      "Required flags:\n" 
      "-linearform=<lfname>\n"                          \
      "    linear-form providing the right hand side\n" \
      "-gridfunction=<gfname>\n" \
      "    grid-function to store the solution vector\n" 
	<< endl;
  }
};





namespace demo_coupling_cpp
{
  class Init
  { 
  public: 
    Init ();
  };
    
  Init::Init()
  {
    GetNumProcs().AddNumProc ("democoupling", NumProcCouplingDemo::Create, NumProcCouplingDemo::PrintDoc);
  }
    
  Init init;
}
  
