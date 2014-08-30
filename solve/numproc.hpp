#ifndef FILE_NUMPROC
#define FILE_NUMPROC

/*********************************************************************/
/* File:   numproc.hh                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Jul. 2000                                             */
/*********************************************************************/

namespace ngsolve
{

class PDE;
  

/** 
    numerical procedure.  
    
    A numproc is an object which performs some numerical computations,
    e.g. solving a boundary value problem, or doing some
    post-processing. Typically, a numproc calls a linear algebra
    function or finite element operation for doing so.

    A numproc collects some components in the constructor, and
    performs the action in the virtual function Do.
*/

class NGS_DLL_HEADER NumProc : public NGS_Object
{
protected:
  /// reference to the PDE the numproc belongs to
  PDE & pde;
  
  // int callposition;
public:
  /// 
  // NumProc (PDE & apde); // , const int acallposition = 0);
  /**
     Generate a numproc. 
     Use objects of pde container, parameterized by flags
  */
  NumProc (PDE & apde, const Flags & flags = Flags()); // , const int acallposition = 0);
  ///
  virtual ~NumProc();
  ///
  virtual void Do(LocalHeap & lh) = 0;
  ///
  virtual void PrintReport (ostream & ost);
  ///
  static void PrintDoc (ostream & ost);

  // int GetCallPosition (void) const { return callposition;} 
};




/// Registered numprocs container
class NGS_DLL_HEADER NumProcs
{
public:
  class NumProcInfo
  {
  public:
    string name;
    int dim;

    shared_ptr<NumProc> (*creator)(PDE & pde, const Flags & flags);
    void (*printdoc) (ostream & ost);
    
    NumProcInfo (const string & aname, int adim, 
		 shared_ptr<NumProc> (*acreator)(PDE & pde, const Flags & flags),
		 void (*aprintdoc) (ostream & ost));
  };

  Array<NumProcInfo*> npa;
public:
  NumProcs();
  ~NumProcs();  
  void AddNumProc (const string & aname, 
		   shared_ptr<NumProc> (*acreator)(PDE & pde, const Flags & flags),
		   void (*printdoc) (ostream & ost) = NumProc::PrintDoc);


  void AddNumProc (const string & aname, int dim, 
		   shared_ptr<NumProc> (*acreator)(PDE & pde, const Flags & flags),
		   void (*printdoc) (ostream & ost) = NumProc::PrintDoc);

  const Array<NumProcInfo*> & GetNumProcs() { return npa; }
  const NumProcInfo * GetNumProc(const string & name, int dim);

  void Print (ostream & ost) const;
};



  /// global variable of all registered numprocs 
extern NGS_DLL_HEADER NumProcs & GetNumProcs ();

  /**
     A template-mechanism to register numprocs.

     Use:
     RegisterNumProc<MyNumProc> anyname("mynumproc");
     to register you numproc with a given name.
  */
template <typename NP>
class RegisterNumProc
{
public:
  RegisterNumProc (string label, int dim = -1)
  {
    GetNumProcs().AddNumProc (label, dim, Create, NP::PrintDoc);
  }
  
  static shared_ptr<NumProc> Create (PDE & pde, const Flags & flags)
  {
    return shared_ptr<NumProc> (new NP (pde, flags));
  }
};


}

#endif
