#ifdef PARALLEL


#include "../include/solve.hpp"
#include <parallelngs.hpp>
//#include <mpe.h>

namespace ngsolve
{
  using namespace ngsolve;
  using namespace ngparallel;
  // using namespace ngbla;
  // using namespace ngla;






  /* *************************** Numproc BVP ********************** */


  ///
  class NumProcBVPParallel : public NumProc
  {
  protected:
    ///
    BilinearForm * bfa;
    ///
    LinearForm * lff;
    ///
    GridFunction * gfu;
    ///
    Preconditioner * pre;
    ///
    int maxsteps;
    ///
    double prec;
    ///
    double tau, taui;
    ///
    bool print;
    ///
    enum SOLVER { CG, GMRES, QMR, NCG, SIMPLE, DIRECT, BICGSTAB };
    ///
    enum IP_TYPE { SYMMETRIC, HERMITEAN, CONJ_HERMITEAN };
    ///
    SOLVER solver;
    ///
    IP_TYPE ip_type;

  public:
    ///
    NumProcBVPParallel (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcBVPParallel();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcBVPParallel (pde, flags);
    }


    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual string GetClassName () const
    {
      return "Boundary Value Problem";
    }

    virtual void PrintReport (ostream & ost)
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
        case NCG:
          ost << "NCG" << endl; break;
        case SIMPLE:
          ost << "Simple" << endl; break;
        case DIRECT:
          ost << "DIRECT" << endl; break;
	case BICGSTAB:
	  ost << "BiCGStab" << endl; break;
        default:
          ost << "Unkown solver-type" << endl;
        }
      ost << "precision     = " << prec << endl
	  << "maxsteps      = " << maxsteps << endl;
    }

    ///
    static void PrintDoc (ostream & ost);
  };




  NumProcBVPParallel :: NumProcBVPParallel (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", ""));
    gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
    pre = pde.GetPreconditioner (flags.GetStringFlag ("preconditioner", ""), 1);
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
      const char* sol = flags.GetStringFlag("solver","cg");
      if (strcmp(sol,"cg")==0) solver = CG;
      if (strcmp(sol,"qmr")==0) solver = QMR;
      if (strcmp(sol,"gmres")==0) solver = GMRES;
      if (strcmp(sol,"simple")==0) solver = SIMPLE;
      if (strcmp(sol,"direct")==0) solver = DIRECT;
      if (strcmp(sol,"bicgstab")==0) solver = BICGSTAB;
    }       
    
    const char* ipflag = flags.GetStringFlag("innerproduct","symmetric");
    if (strcmp(ipflag,"symmetric")==0) ip_type = SYMMETRIC;
    if (strcmp(ipflag,"hermitean")==0) ip_type = HERMITEAN;
    if (strcmp(ipflag,"conj_hermitean")==0) ip_type = CONJ_HERMITEAN;

    print = flags.GetDefineFlag ("print");
  }

  NumProcBVPParallel :: ~NumProcBVPParallel()
  {
    ;
  }

  void NumProcBVPParallel :: PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc BVPParallel:\n" \
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
      "\n-solver=<solvername> (cg|qmr|gmres|direct)\n"\
      "-preconditioner=<prename>\n"
      "-maxsteps=n\n"
      "-prec=eps\n"
      "-print\n"
      "    write matrices and vectors into logfile\n"
	<< endl;
  }


  void NumProcBVPParallel :: Do(LocalHeap & lh)
  {
    static int solvetimer = NgProfiler::CreateTimer ("Equation solving");
    NgProfiler::RegionTimer reg (solvetimer);

    if (!lff->IsAssembled()) lff->Assemble(lh);

    cout << "solve bvp parallel" << endl;
    *testout << "solve bvp parallel" << endl;


    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
 
    const BaseMatrix & mat = bfa->GetMatrix();
    const ParallelBaseMatrix * pmat = dynamic_cast<const ParallelBaseMatrix*> (&mat);
    const BaseMatrix * consistentmat;

    if ( pmat ) consistentmat = pmat -> ConsistentMat();

    BaseVector & vecf = lff->GetVector();
    BaseVector & vecu = gfu->GetVector();

    const FESpace & fespace = bfa-> GetFESpace();
    ParallelDofs * paralleldofs = & (fespace.GetParallelDofs() );

//     mat . IsParallelMatrix();
//     (*testout) << vecf << endl;
 
     if (print)
      {
	paralleldofs->Print();
	(*testout) << "MatrixHeight = " << endl << mat.VHeight() << endl;
	(*testout) << "MatrixWidth = " << endl << mat.VWidth() << endl;
	(*testout) << "Matrix = " << endl << mat << endl;
	if ( consistentmat )
	// if ( id > 0 )
	(*testout) << "Consistent mat = " << endl << *consistentmat << endl;
	(*testout) << "RHS-Vector = " << endl << vecf << endl;
      }


 
    const BaseMatrix * premat = NULL;
    if (pre)  premat = &(pre->GetMatrix());
    
    KrylovSpaceSolver * invmat;
    BaseMatrix * invmat2; 

    if (!bfa->GetFESpace().IsComplex())
      {
	switch (solver)
	  {
          case CG:
	    cout << "cg solve for real problem" << endl;
	    invmat = new CGSolver<double>(mat, *premat);
	    break;
          case BICGSTAB:
	    cout << "bicgstab solve for real problem" << endl;
	    invmat = new BiCGStabSolver<double>(mat, *premat);
	    break;
          case QMR:
            cout << "qmr solve for real problem" << endl;
            invmat = new QMRSolver<double>(mat, *premat);
            break;
	  case GMRES:
            cout << "gmres solve for real problem" << endl;
            invmat = new GMRESSolver<double>(mat, *premat);
	    break;
	  case SIMPLE:
            {
              cout << "simple solve for real problem" << endl;
              SimpleIterationSolver<double> * hinv = new SimpleIterationSolver<double>(mat, *premat);
              hinv -> SetTau (tau);
              invmat = hinv;
              break;
            }
          case DIRECT:
            cout << "direct solve for real problem" << endl;
            invmat2 = dynamic_cast<const BaseSparseMatrix&> (mat) . InverseMatrix(); 
            break;
	  }
      }
    else if ( ip_type == SYMMETRIC )
      {
 	switch (solver)
	  {
	  case CG:
            cout << "cg solve for complex problem" << endl;
            invmat = new CGSolver<Complex>(mat, *premat);
	    break;
          case BICGSTAB:
	    cout << "bicgstab solve for real problem" << endl;
	    invmat = new BiCGStabSolver<Complex>(mat, *premat);
	    break;
	  case QMR:
            cout << "qmr solve for complex problem" << endl;
            invmat = new QMRSolver<Complex>(mat, *premat);
	    break;
	  case GMRES:
            cout << "gmres solve for complex problem" << endl;
            invmat = new GMRESSolver<Complex>(mat, *premat);
	    break;
	  case SIMPLE:
            {
              cout << "simple solve for complex problem" << endl;
              SimpleIterationSolver<Complex> * hinv = new SimpleIterationSolver<Complex>(mat, *premat);
              hinv -> SetTau (Complex (tau, taui));
              invmat = hinv;
              break;
            }
          case DIRECT:
            cout << "direct solve for real problem" << endl;
            invmat2 = dynamic_cast<const BaseSparseMatrix&> (mat) . InverseMatrix(); 
            break;
          }
      }
    else if ( ip_type == HERMITEAN )
      {
 	switch (solver)
	  {
	  case CG:
            cout << "cg solve for complex problem" << endl;
            invmat = new CGSolver<ComplexConjugate>(mat, *premat);
	    break;
          case BICGSTAB:
	    cout << "bicgstab solve for real problem" << endl;
	    invmat = new BiCGStabSolver<ComplexConjugate>(mat, *premat);
	    break;
	  case QMR:
            cout << "qmr solve for complex problem" << endl;
            invmat = new QMRSolver<ComplexConjugate>(mat, *premat);
	    break;
	  case GMRES:
            cout << "gmres solve for complex problem" << endl;
            invmat = new GMRESSolver<ComplexConjugate>(mat, *premat);
	    break;
	  case SIMPLE:
            {
              cout << "simple solve for complex problem" << endl;
              SimpleIterationSolver<ComplexConjugate> * hinv = new SimpleIterationSolver<ComplexConjugate>(mat, *premat);
              hinv -> SetTau (Complex (tau, taui));
              invmat = hinv;
              break;
            }
          case DIRECT:
            cout << "direct solve for real problem" << endl;
            invmat2 = dynamic_cast<const BaseSparseMatrix&> (mat) . InverseMatrix(); 
            break;
          }
      }
    else if ( ip_type == CONJ_HERMITEAN )
      {
 	switch (solver)
	  {
	  case CG:
            cout << "cg solve for complex problem" << endl;
            invmat = new CGSolver<ComplexConjugate2>(mat, *premat);
	    break;
          case BICGSTAB:
	    cout << "bicgstab solve for real problem" << endl;
	    invmat = new BiCGStabSolver<ComplexConjugate2>(mat, *premat);
	    break;
	  case QMR:
            cout << "qmr solve for complex problem" << endl;
            invmat = new QMRSolver<ComplexConjugate2>(mat, *premat);
	    break;
	  case GMRES:
            cout << "gmres solve for complex problem" << endl;
            invmat = new GMRESSolver<ComplexConjugate2>(mat, *premat);
	    break;
	  case SIMPLE:
            {
              cout << "simple solve for complex problem" << endl;
              SimpleIterationSolver<ComplexConjugate2> * hinv = new SimpleIterationSolver<ComplexConjugate2>(mat, *premat);
              hinv -> SetTau (Complex (tau, taui));
              invmat = hinv;
              break;
            }
          case DIRECT:
            cout << "direct solve for real problem" << endl;
            invmat2 = dynamic_cast<const BaseSparseMatrix&> (mat) . InverseMatrix(); 
            break;
          }
      }

    if (solver != DIRECT)  
      {
        ma.PushStatus ("Iterative solver");
        invmat->SetMaxSteps (maxsteps);        invmat->SetPrecision (prec);
        invmat->SetPrintRates ();
        invmat->SetInitialize (0);
        invmat->SetStatusHandler(ma);
      }
    else
      ma.PushStatus ("Direct solver");
      
    clock_t starttime, endtime, soltime;
    double startwtime, endwtime;
    starttime = clock();

    MPI_Barrier( MPI_COMM_WORLD );
    if ( id == 0 ) {
        startwtime = MPI_Wtime();
    }
    MPI_Barrier( MPI_COMM_WORLD );

    if (solver != DIRECT)
      invmat->Mult (vecf, vecu);
    else 
      invmat2->Mult (vecf, vecu);
    

    ma.PopStatus ();

    ParallelBaseVector * paru = dynamic_cast<ParallelBaseVector*> (&vecu);
    paru->Distribute();
    paru->AllReduce(&hoprocs);
    
    endtime = clock();
    cout << "Solution time = " << double(endtime - starttime)/CLOCKS_PER_SEC << endl;

    if ( id == 0 ) {
      endwtime = MPI_Wtime();
      printf( "wall clock time = %f\n", double(endwtime-startwtime) );
    }

    if (solver != DIRECT)
      {
        cout << "Iterations: " << invmat->GetSteps() << endl;
        delete invmat;
      }
    else
      delete invmat2;

    bfa -> ComputeInternal (vecu, vecf, lh);
    
    if (print)
      (*testout) << "Solution = " << endl << vecu << endl;


   }



















  


  namespace bvpparallel_cpp
  {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetNumProcs().AddNumProc ("bvpparallel", NumProcBVPParallel::Create, NumProcBVPParallel::PrintDoc);

    }
    
    
    Init init;
    
  }
  


}



#endif
