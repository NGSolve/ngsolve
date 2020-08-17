#include "../include/solve.hpp"
#include <parallelngs.hpp>


int dummy_bvp = -1;

namespace ngsolve
{

  /* *************************** Numproc BVP ********************** */

  ///
  class NumProcBVP : public NumProc
  {
  protected:
    shared_ptr<BilinearForm> bfa;
    // std::result_of<PDE.GetBilinearForm(const string&,bool)>::type bfa;
    // decltype ( GetReturnValue (&PDE::GetBilinearForm) ) bfa;
    ///
    shared_ptr<LinearForm> lff;
    ///
    shared_ptr<GridFunction> gfu;
    ///
    shared_ptr<Preconditioner> pre;
    ///
    int maxsteps;
    ///
    double prec;
    ///
    double tau, taui;
    ///
    bool print;
    ///
    enum SOLVER { CG, GMRES, QMR/*, NCG */, SIMPLE, DIRECT, BICGSTAB };
    ///
    enum IP_TYPE { SYMMETRIC, HERMITEAN, CONJ_HERMITEAN };
    ///
    SOLVER solver;

    IP_TYPE ip_type;
    ///
    bool useseedvariant;
  public:
    ///
    NumProcBVP (shared_ptr<PDE> apde, const Flags & flags);
    
    NumProcBVP (shared_ptr<BilinearForm> abfa, 
                shared_ptr<LinearForm> alff,
                shared_ptr<GridFunction> agfu, 
                shared_ptr<Preconditioner> apre, 
                int amaxsteps, 
                double aprec)
      : bfa(abfa), lff(alff), gfu(agfu), pre(apre), 
        maxsteps(amaxsteps), prec(aprec)
    {
      print = false;
      solver = CG;
      ip_type = SYMMETRIC;
    }

    ///
    virtual ~NumProcBVP();

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual string GetClassName () const
    {
      return "Boundary Value Problem";
    }

    virtual void PrintReport (ostream & ost) const
    {
      ost << GetClassName() << endl
	  << "Bilinear-form = " << bfa->GetName() << endl
	  << "Linear-form   = " << lff->GetName() << endl
	  << "Gridfunction  = " << gfu->GetName() << endl
	  << "Preconditioner = " << ((pre) ? pre->ClassName() : "None") << endl
	  << "solver        = "; 
      
      switch(solver)
        {
        case CG:
          ost << "CG" << endl; break;
        case QMR:
          ost << "QMR" << endl; break;
        case GMRES:
          ost << "GMRES" << endl; break;
          // case NCG:
          // ost << "NCG" << endl; break;
        case SIMPLE:
          ost << "Simple" << endl; break;
        case DIRECT:
          ost << "DIRECT" << endl; break;
	case BICGSTAB:
	  ost << "BiCGStab" << endl; break;
        default:
          ost << "Unknown solver-type" << endl;
        }
      ost << "precision     = " << prec << endl
	  << "maxsteps      = " << maxsteps << endl;
    }

    ///
    static void PrintDoc (ostream & ost);
  };



  NumProcBVP :: NumProcBVP (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    lff = apde->GetLinearForm (flags.GetStringFlag ("linearform", ""));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
    if (flags.StringFlagDefined("preconditioner"))
      pre = apde->GetPreconditioner (flags.GetStringFlag ("preconditioner", ""));
    else
      pre = NULL;
    maxsteps = int(flags.GetNumFlag ("maxsteps", 200));
    prec = flags.GetNumFlag ("prec", 1e-12);
    tau = flags.GetNumFlag ("tau", 1);
    taui = flags.GetNumFlag ("taui", 0);
    solver = CG;

    if (flags.GetDefineFlag ("qmr")) solver = QMR;
// old flag-style: no longer supported. only qmr for compatibility 
    if (flags.GetDefineFlag ("gmres")) 
      cout << "*** warning: flag -gmres deprecated: use -solver=gmres instead" << endl;
    if (flags.GetDefineFlag ("ncg"))
      cout << "*** warning: flag -ncg deprecated: use -solver=ncg instead" << endl;
    if (flags.GetDefineFlag ("direct"))
      cout << "*** warning: flag -direct deprecated: use -solver=direct instead" << endl;
    
    // new style: -solver=cg|qmr|gmres|direct|bicgstab
    {
      string solvername = flags.GetStringFlag("solver","cg");
      if (solvername == "cg")     solver = CG;
      if (solvername == "qmr")    solver = QMR;
      if (solvername == "gmres")  solver = GMRES;
      if (solvername == "simple") solver = SIMPLE;
      if (solvername == "direct") solver = DIRECT;
      if (solvername == "bicgstab") solver = BICGSTAB;
    }       
    
    string ipflag = flags.GetStringFlag("innerproduct","symmetric");
    ip_type = SYMMETRIC;
    if (ipflag == "symmetric") ip_type = SYMMETRIC;
    if (ipflag == "hermitean") ip_type = HERMITEAN;
    if (ipflag == "hermitian") ip_type = HERMITEAN;
    if (ipflag == "conj_hermitean") ip_type = CONJ_HERMITEAN;
    if (ipflag == "conj_hermitian") ip_type = CONJ_HERMITEAN;

    print = flags.GetDefineFlag ("print");
    useseedvariant = flags.GetDefineFlag ("seed");

    if (solver != DIRECT)
      apde->AddVariable (string("bvp.")+flags.GetStringFlag ("name",NULL)+".its", 0.0, 6);
  }

  NumProcBVP :: ~NumProcBVP()
  {
    ;
  }

  void NumProcBVP :: PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc BVP:\n" \
      "------------\n" \
      "Solves the linear system resulting from a boundary value problem\n\n" \
      "Required flags:\n" 
      "-bilinearform=<bfname>\n" 
      "    bilinear-form providing the matrix\n" \
      "-linearform=<lfname>\n" \
      "    linear-form providing the right hand side\n" \
      "-gridfunction=<gfname>\n" \
      "    grid-function to store the solution vector\n" 
      "\nOptional flags:\n"\
      "\n-solver=<solvername> (cg|qmr|gmres|direct|bicgstab)\n"\
      "-seed\n"\
      "    use seed variant for multiple rhs\n"\
      "-preconditioner=<prename>\n"
      "-maxsteps=n\n"
      "-prec=eps\n"
      "-print\n"
      "    write matrices and vectors into logfile\n"
	<< endl;
  }


  void NumProcBVP :: Do(LocalHeap & lh)
  {
    static Timer timer("Equation solving");
    timer.Start();

    if (!lff->IsAssembled()) lff->Assemble(lh);

    cout << IM(1) << "solve bvp" << endl;

    shared_ptr<BaseMatrix> mat = bfa->GetMatrixPtr();
    const BaseVector & vecf = lff->GetVector();
    BaseVector & vecu = gfu->GetVector();

    bool eliminate_internal = bfa -> UsesEliminateInternal();

    if (print)
      {
	(*testout) << "MatrixHeight = " << endl << mat->VHeight() << endl;
	(*testout) << "MatrixWidth = " << endl << mat->VWidth() << endl;
	(*testout) << "Matrix = " << endl << mat << endl;
	(*testout) << "RHS-Vector = " << endl << vecf << endl;
      }

    shared_ptr<BaseMatrix> premat;
    if (pre) premat = pre->GetMatrixPtr();
    
    KrylovSpaceSolver * invmat = NULL;
    shared_ptr<BaseMatrix> invmat2;

    if (!bfa->GetFESpace()->IsComplex())
      {
	switch (solver)
	  {
          case CG:
	    cout << IM(1) << "cg solve for real system" << endl;
	    invmat = new CGSolver<double>(mat, premat);
	    break;
          case BICGSTAB:
	    cout << IM(1) << "bicgstab solve for real system" << endl;
	    invmat = new BiCGStabSolver<double>(mat, premat);
	    break;
          case QMR:
            cout << IM(1) << "qmr solve for real system" << endl;
            invmat = new QMRSolver<double>(mat, premat);
            break;
	  case GMRES:
            cout << IM(1) << "gmres solve for real system" << endl;
            invmat = new GMRESSolver<double>(mat, premat);
	    break;
	  case SIMPLE:
            {
              cout << IM(1) << "simple solve for real system" << endl;
              SimpleIterationSolver<double> * hinv = new SimpleIterationSolver<double>(mat, premat);
              hinv -> SetTau (tau);
              invmat = hinv;
              break;
            }
          case DIRECT:
            cout << IM(1) << "direct solve for real system" << endl;
            invmat2 = mat->InverseMatrix(bfa->GetFESpace()->GetFreeDofs(eliminate_internal)); 
            // invmat2 = dynamic_cast<const BaseSparseMatrix&> (mat) . InverseMatrix(bfa->GetFESpace()->GetFreeDofs(eliminate_internal)); 
            break;
	  }
      }
    else if ( ip_type == SYMMETRIC )
      {
 	switch (solver)
	  {
	  case CG:
            cout << IM(1) << "cg solve for complex system" << endl;
            invmat = new CGSolver<Complex>(mat, premat);
	    break;
          case BICGSTAB:
	    cout << IM(1) << "bicgstab solve for complex system" << endl;
	    invmat = new BiCGStabSolver<Complex>(mat, premat);
	    break;
	  case QMR:
            cout << IM(1) << "qmr solve for complex system" << endl;
            invmat = new QMRSolver<Complex>(mat, premat);
	    break;
	  case GMRES:
            cout << IM(1) << "gmres solve for complex system" << endl;
            invmat = new GMRESSolver<Complex>(mat, premat);
	    break;
	  case SIMPLE:
            {
              cout << IM(1) << "simple solve for complex system" << endl;
              SimpleIterationSolver<Complex> * hinv = new SimpleIterationSolver<Complex>(mat, premat);
              hinv -> SetTau (Complex (tau, taui));
              invmat = hinv;
              break;
            }
          case DIRECT:
            cout << IM(1) << "direct solve for complex system" << endl;
            // invmat2 = dynamic_cast<const BaseSparseMatrix&> (mat) 
            //   . InverseMatrix(bfa->GetFESpace()->GetFreeDofs(eliminate_internal)); 
	    invmat2 = mat->InverseMatrix(bfa->GetFESpace()->GetFreeDofs(eliminate_internal)); 
	    break;
          }
      }
    else if ( ip_type == HERMITEAN )
      {
 	switch (solver)
	  {
	  case CG:
            cout << IM(1) << "cg solve for complex system" << endl;
            invmat = new CGSolver<ComplexConjugate>(mat, premat);
	    break;
          case BICGSTAB:
	    cout << IM(1) << "bicgstab solve for complex system" << endl;
	    invmat = new BiCGStabSolver<ComplexConjugate>(mat, premat);
	    break;
	  case QMR:
            cout << IM(1) << "qmr solve for complex system" << endl;
            invmat = new QMRSolver<ComplexConjugate>(mat, premat);
	    break;
	  case GMRES:
            cout << IM(1) << "gmres solve for complex system" << endl;
            invmat = new GMRESSolver<ComplexConjugate>(mat, premat);
	    break;
	  case SIMPLE:
            {
              cout << IM(1) << "simple solve for complex system" << endl;
              SimpleIterationSolver<ComplexConjugate> * hinv = new SimpleIterationSolver<ComplexConjugate>(mat, premat);
              hinv -> SetTau (Complex (tau, taui));
              invmat = hinv;
              break;
            }
          case DIRECT:
            cout << IM(1) << "direct solve for complex system" << endl;
	    invmat2 = mat->InverseMatrix(bfa->GetFESpace()->GetFreeDofs(eliminate_internal)); 
	    // invmat2 = dynamic_cast<const BaseSparseMatrix&> (mat) . InverseMatrix(bfa->GetFESpace()->GetFreeDofs(eliminate_internal)); 
            break;
          }
      }
    else if ( ip_type == CONJ_HERMITEAN )
      {
 	switch (solver)
	  {
	  case CG:
            cout << IM(1) << "cg solve for complex system" << endl;
            invmat = new CGSolver<ComplexConjugate2>(mat, premat);
	    break;
          case BICGSTAB:
	    cout << IM(1) << "bicgstab solve for complex system" << endl;
	    invmat = new BiCGStabSolver<ComplexConjugate2>(mat, premat);
	    break;
	  case QMR:
            cout << IM(1) << "qmr solve for complex system" << endl;
            invmat = new QMRSolver<ComplexConjugate2>(mat, premat);
	    break;
	  case GMRES:
            cout << IM(1) << "gmres solve for complex system" << endl;
            invmat = new GMRESSolver<ComplexConjugate2>(mat, premat);
	    break;
	  case SIMPLE:
            {
              cout << IM(1) << "simple solve for complex system" << endl;
              SimpleIterationSolver<ComplexConjugate2> * hinv = new SimpleIterationSolver<ComplexConjugate2>(mat, premat);
              hinv -> SetTau (Complex (tau, taui));
              invmat = hinv;
              break;
            }
          case DIRECT:
            cout << IM(1) << "direct solve for complex system" << endl;
            invmat2 = dynamic_pointer_cast<BaseSparseMatrix> (mat) -> InverseMatrix(bfa->GetFESpace()->GetFreeDofs(eliminate_internal)); 
            break;
          }
      }

    if (solver != DIRECT)  
      {
        if (ma) ma->PushStatus ("Iterative solver");
        invmat->SetMaxSteps (maxsteps);
        invmat->SetPrecision (prec);
        invmat->SetPrintRates ();
        invmat->SetInitialize (0);
        invmat->SetStatusHandler(ma);
	invmat->UseSeed(useseedvariant);
      }
    else
      if (ma) ma->PushStatus ("Direct solver");

      
    double starttime, endtime;
    starttime = WallTime(); // clock();


    auto hv = vecu.CreateVector();
    if (solver != DIRECT)
      {
	hv = vecf;
        bfa->ModifyRHS (hv);
        
	RunWithTaskManager 
	  ( [invmat,&hv,&vecu] ()
	    {
	      invmat -> Mult(hv, vecu);
	    } );

        // invmat->Mult (hv, vecu);
      }
    else 
      {
	hv = vecf - *mat * vecu;
        bfa->ModifyRHS (hv);
	vecu += *invmat2 * hv;
      }

    if (ma) ma->PopStatus ();
    
    if (print)
      (*testout) << "Solution = " << endl << vecu << endl;

    endtime = WallTime(); // clock();
    
    cout << IM(1) << "Solution time = " << endtime - starttime << " sec wall time" << endl;
    if (solver != DIRECT)
      {
	cout << IM(1) << "Iterations: " << invmat->GetSteps() << endl;

	auto sp = pde.lock();
	if (sp)
	  sp -> AddVariable (string("bvp.")+GetName()+".its", invmat->GetSteps(), 6);
      }

    if (solver != DIRECT) 
      delete invmat;

    timer.Stop();

    bfa -> ComputeInternal (vecu, vecf, lh);

    if (print)
      (*testout) << "Solution = " << endl << vecu << endl;
  }















  /* *************************** Numproc ConstrainedBVP ********************** */


  ///
  class NumProcConstrainedBVP : public NumProc
  {
  protected:
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    shared_ptr<LinearForm> lff;
    ///
    shared_ptr<GridFunction> gfu;
    ///
    shared_ptr<Preconditioner> pre;
    ///
    int maxsteps;
    ///
    double prec;
    ///
    bool print;
    ///
    enum SOLVER { CG, QMR /* , NCG */ };
    ///
    SOLVER solver;
    ///
    Array<shared_ptr<LinearForm>> constraints;
    
  public:
    ///
    NumProcConstrainedBVP (shared_ptr<PDE> apde, const Flags & flags);
    ///
    virtual ~NumProcConstrainedBVP();

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual string GetClassName () const
    {
      return "Boundary Value Problem";
    }

    virtual void PrintReport (ostream & ost) const
    {
      ost << GetClassName() << endl
	  << "Bilinear-form = " << bfa->GetName() << endl
	  << "Linear-form   = " << lff->GetName() << endl
	  << "Gridfunction  = " << gfu->GetName() << endl
	  << "Preconditioner = " << ((pre) ? pre->ClassName() : "None") << endl
	  << "solver        = " << ((solver==CG) ? "CG" : "QMR") << endl
	  << "precision     = " << prec << endl
	  << "maxsteps      = " << maxsteps << endl;
    }

    ///
    static void PrintDoc (ostream & ost);
  };



  class ConstrainedPrecondMatrix : public BaseMatrix
  {
    shared_ptr<BaseMatrix> c1;
    Array<shared_ptr<BaseVector>> constraints;
    Array<shared_ptr<BaseVector>> c1constraints;
    Matrix<double> projection, invprojection;
    int ncnt;
  public:
    ConstrainedPrecondMatrix (shared_ptr<BaseMatrix> ac1)
      : c1(ac1) { ncnt = 0; }
    
    virtual ~ConstrainedPrecondMatrix () { ; }

    int VHeight() const override { return c1->VHeight(); }
    int VWidth() const override { return c1->VWidth(); }
    bool IsComplex() const override { return c1->IsComplex(); }

    AutoVector CreateRowVector() const override { return c1->CreateRowVector(); }
    AutoVector CreateColVector() const override { return c1->CreateColVector(); }
    
    void AddConstraint (shared_ptr<BaseVector> hv)
    {
      constraints.Append (hv);
      c1constraints.Append (hv->CreateVector());
      *c1constraints.Last() = (*c1) * *constraints.Last();

      ncnt = constraints.Size();
      projection.SetSize(ncnt);
      invprojection.SetSize(ncnt);

      for (int i = 0; i < ncnt; i++)
	for (int j = 0; j < ncnt; j++)
	  projection(i,j) = InnerProduct (*constraints[i], *c1constraints[j]);
      for (int i = 0; i < ncnt; i++)
	projection(i,i) += 1;

      CalcInverse (projection, invprojection);
      // cout << "projection = " << endl << projection << endl;
      // cout << "invprojection = " << endl << invprojection << endl;
    }

    void Mult (const BaseVector & x, BaseVector & y) const override
    {
      c1 -> Mult (x, y);
      Vector<double> hv1(ncnt), hv2(ncnt);
      
      for (int i = 0; i < ncnt; i++)
	hv1(i) = InnerProduct (x, *c1constraints[i]);

      hv2 = invprojection * hv1;

      //      cout << "hv1 = " << hv1 << ", hv2 = " << hv2 << endl;

      for (int i = 0; i < ncnt; i++)
	y -= hv2(i) * *c1constraints[i];
    }
  };



  class ConstrainedMatrix : public BaseMatrix
  {
    const BaseMatrix * a1;
    Array<const BaseVector*> constraints;
    int ncnt;
  public:
    ConstrainedMatrix (const BaseMatrix * aa1)
      : a1(aa1) { ncnt = 0; }
    
    virtual ~ConstrainedMatrix () { ; }

    bool IsComplex() const override { return a1 -> IsComplex(); } 
    void AddConstraint (const BaseVector * hv)
    {
      constraints.Append (hv);
      ncnt = constraints.Size();
    }

    AutoVector CreateRowVector() const override { return a1->CreateRowVector(); }
    AutoVector CreateColVector() const override { return a1->CreateColVector(); }

    int VHeight() const override { return a1->VHeight(); }
    int VWidth() const override { return a1->VWidth(); }

    void Mult (const BaseVector & x, BaseVector & y) const override
    {
      a1 -> Mult (x, y);
      Vector<double> hv1(ncnt);
      
      for (int i = 0; i < ncnt; i++)
	hv1(i) = InnerProduct (x, *constraints[i]);
      for (int i = 0; i < ncnt; i++)
	y += hv1(i) * *constraints[i];
    }

    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      a1 -> MultAdd (s, x, y);
      Vector<double> hv1(ncnt);
      
      for (int i = 0; i < ncnt; i++)
	hv1(i) = InnerProduct (x, *constraints[i]);
      for (int i = 0; i < ncnt; i++)
	y += (s*hv1(i)) * *constraints[i];
    }
  };





  NumProcConstrainedBVP :: NumProcConstrainedBVP (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    lff = apde->GetLinearForm (flags.GetStringFlag ("linearform", ""));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
    pre = apde->GetPreconditioner (flags.GetStringFlag ("preconditioner", ""), 1);
    maxsteps = int(flags.GetNumFlag ("maxsteps", 200));
    prec = flags.GetNumFlag ("prec", 1e-12);
    solver = CG;
    if (flags.GetDefineFlag ("qmr")) solver = QMR;
    // if (flags.GetDefineFlag ("ncg")) solver = NCG;
    print = flags.GetDefineFlag ("print");

    const Array<string> & cnts = flags.GetStringListFlag ("constraints");
    for (int i = 0; i < cnts.Size(); i++)
      constraints.Append (apde->GetLinearForm (cnts[i]));

    apde->AddVariable (string("constrbvp.")+flags.GetStringFlag ("name",NULL)+".its", 0.0, 6);
  }

  NumProcConstrainedBVP :: ~NumProcConstrainedBVP()
  {
    ;
  }

  void NumProcConstrainedBVP :: PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc ConstrainedBVP:\n" \
      "------------\n" \
      "Solves the linear system resulting from a boundary value problem\n\n" \
      "Required flags:\n" 
      "-bilinearform=<bfname>\n" 
      "    bilinear-form providing the matrix\n" \
      "-linearform=<lfname>\n" \
      "    linear-form providing the right hand side\n" \
      "-gridfunction=<gfname>\n" \
      "    grid-function to store the solution vector\n" 
      "\nOptional flags:\n"
      "-preconditioner=<prename>\n"
      "-maxsteps=n\n"
      "-prec=eps\n"
      "-qmr\n"
      "-print\n"
      "    write matrices and vectors into logfile\n"
	<< endl;
  }


  void NumProcConstrainedBVP :: Do(LocalHeap & lh)
  {
    cout << "solve constrained bvp" << endl;

    const BaseMatrix & mat = bfa->GetMatrix();
    const BaseVector & vecf = lff->GetVector();
    auto & vecu = gfu->GetVector();

    if (print)
      {
	(*testout) << "MatrixHeight = " << endl << mat.VHeight() << endl;
	(*testout) << "MatrixWidth = " << endl << mat.VWidth() << endl;
	(*testout) << "Matrix = " << endl << mat << endl;
	(*testout) << "RHS-Vector = " << endl << vecf << endl;
      }

    shared_ptr<BaseMatrix> premat;
    if (pre)  
      {
	premat = pre->GetMatrixPtr();

	auto hpre = make_shared<ConstrainedPrecondMatrix> (premat);
	premat = hpre;
	for (int i = 0; i < constraints.Size(); i++)
	  hpre->AddConstraint(constraints[i]->GetVectorPtr());
      }

    auto hmat = make_shared<ConstrainedMatrix> (&mat);
    for (int i = 0; i < constraints.Size(); i++)
      hmat->AddConstraint(&constraints[i]->GetVector());
    

    KrylovSpaceSolver * invmat = NULL;

    if (!bfa->GetFESpace()->IsComplex())
      {
	switch (solver)
	  {
	  case CG:
	    invmat = new CGSolver<double>(hmat, premat);
	    break;
	  case QMR:
	    invmat = new QMRSolver<double>(hmat, premat);
	    break;
	  }
      }
    else
      {
	switch (solver)
	  {
	  case CG:
	    invmat = new CGSolver<Complex>(hmat, premat);
	    break;
	  case QMR:
	    invmat = new QMRSolver<Complex>(hmat, premat);
	    break;
	  }
      }

    if (ma) ma->PushStatus ("Iterative solver");

    invmat->SetMaxSteps (maxsteps);
    invmat->SetPrecision (prec);
    invmat->SetPrintRates ();
    invmat->SetInitialize (0);

    clock_t starttime, endtime;
    starttime = clock();

    invmat->Mult (vecf, vecu);

    if (ma) ma->PopStatus ();
    
    if (print)
      (*testout) << "Solution = " << endl << vecu << endl;

    endtime = clock();
    cout << "Solution time = " << double(endtime - starttime)/CLOCKS_PER_SEC << endl;
    cout << "Iterations: " << invmat->GetSteps() << endl;
    *testout << "Solution time = " << double(endtime - starttime)/CLOCKS_PER_SEC << endl;
    *testout << "Iterations: " << invmat->GetSteps() << endl;


    {
      try
      {
        GetPDE() -> AddVariable (string("constrbvp.")+GetName()+".its", invmat->GetSteps(), 6);
      }
      catch (std::exception e) { ; }
    }
    
    
    delete invmat;


    bfa -> ComputeInternal (vecu, vecf, lh);
  }



  static RegisterNumProc<NumProcBVP> npinitbvp("bvp");
  static RegisterNumProc<NumProcConstrainedBVP> npinitbvp2("constrainedbvp");
}


#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"

using namespace ngsolve;
void ExportBVP(py::module &m)
{
  m.def ("BVP", [](shared_ptr<BilinearForm> bfa,
               shared_ptr<LinearForm> lff,
               shared_ptr<GridFunction> gfu,
               shared_ptr<Preconditioner> pre,
               int maxsteps,
               double prec) -> shared_ptr<NumProc>
            
            {
              return make_shared<NumProcBVP> (bfa, lff, gfu, pre, maxsteps, prec);
            },
            py::arg("bf"), py::arg("lf"), py::arg("gf"), 
            py::arg("pre"), py::arg("maxsteps")=100, py::arg("prec")=1e-8
	   , docu_string(R"raw_string(
Solves the given boundary value problem: bf * gf = lf, non homogeneous boundary conditions
on gf are respected (they must be set in advance). If eliminate_internal is set for the
bf, then static condensation of inner bubbles is used.

Parameters:

bf : ngsolve.comp.BilinearForm
  input bilinear form as the right hand side of the equation

lf : ngsolve.comp.LinearForm
  input linear form as the left hand side of the equation

gf : ngsolve.comp.GridFunction
  input GridFunction where the solution is saved

pre : ngsolve.comp.Preconditioner
  input Preconditioner for the problem

maxsteps : int
  input maximal steps. After the maximal step is reached, the computations stop.

prec : float
  input precision of the residuum. if it is reached the computations stop.

)raw_string"));




}
#endif
