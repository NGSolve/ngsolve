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
  shared_ptr<GridFunction> gfu;

  // linear-form providing the right hand side for the second PDE
  shared_ptr<LinearForm> lff;

public:
    
  /*
    In the constructor, the solver class gets the flags from the pde - input file.
    the PDE class apde constains all bilinear-forms, etc...
  */
  NumProcCouplingDemo (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the linear-form and the gridfunction as
    // -linearform=f -gridfunction=u

    lff = apde->GetLinearForm (flags.GetStringFlag ("linearform", "f"));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));
  }

  virtual ~NumProcCouplingDemo() 
  { ; }
  

  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "Compute coupling terms" << endl;
      
    shared_ptr<FESpace> fesu = gfu -> GetFESpace();
    shared_ptr<FESpace> fesf = lff -> GetFESpace();
    
    lff -> GetVector() = 0.0;

    int ne = ma->GetNE();
    for (int i = 0; i < ne; i++)   // loop over elements
      {  
        HeapReset hr(lh);    // reset the local heap memory at the end of the loop
		
        const ElementTransformation & eltrans = ma->GetTrafo (i, 0, lh);

        const ScalarFiniteElement<2> & felu = 
	  dynamic_cast<const ScalarFiniteElement<2>&> (fesu->GetFE (i, lh));
        const ScalarFiniteElement<2> & felf = 
	  dynamic_cast<const ScalarFiniteElement<2>&> (fesf->GetFE (i, lh));

	int nd_u = felu.GetNDof();
	int nd_f = felf.GetNDof();

	Array<int> dnumsu(nd_u, lh), dnumsf(nd_f, lh);  // the dof-numbes for u and f
		      
        fesu->GetDofNrs (i, dnumsu);
        fesf->GetDofNrs (i, dnumsf);

        FlatVector<> elu (nd_u, lh), shape (nd_f, lh), elf (nd_f, lh);

        gfu -> GetElementVector (dnumsu, elu);

        IntegrationRule ir(eltrans.GetElementType(), felu.Order()+felf.Order() );

        
        elf = 0.0;
        for (int j = 0; j < ir.GetNIP(); j++)   // integration points
          {
            MappedIntegrationPoint<2, 2> mip(ir[j], eltrans);  // computes Jacobi matrix etc

            Vec<1> ui;   // value of u in point
            DiffOpId<2>::Apply (felu, mip, elu, ui, lh);   // compute value in point (= elu * shape)
            
            // could use also other differential operators such as
            // DiffOpGradient<2>, DiffOpCurl<2>, ....
	    
 
            double fac = ir[j].Weight() * mip.GetMeasure();   // integration weights
	    
	    felf.CalcShape (ir[j], shape);
            elf += (fac*ui(0)) * shape; 
          }

        lff -> AddElementVector (dnumsf, elf);
      }
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




// register the numproc 'democoupling' 
static RegisterNumProc<NumProcCouplingDemo> npinit1("democoupling");

