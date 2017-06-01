#ifndef FILE_NUMPROC
#define FILE_NUMPROC

/*********************************************************************/
/* File:   numproc.hh                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Jul. 2000                                             */
/*********************************************************************/

namespace ngcomp
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
  weak_ptr<PDE> pde;

public:
  /**
     Generate a numproc. 
     Use objects of pde container, parameterized by flags
  */
  NumProc (weak_ptr<PDE> apde, const Flags & flags = Flags()); 
  ///
  NumProc (const Flags & flags = Flags()); 
  ///
  virtual ~NumProc();
  ///
  virtual void Do(LocalHeap & lh) = 0;
  ///
  shared_ptr<PDE> GetPDE() const { return shared_ptr<PDE>(pde); }
  ///
  virtual void PrintReport (ostream & ost) const;
  ///
  static void PrintDoc (ostream & ost);
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

    shared_ptr<NumProc> (*creator)(shared_ptr<PDE> pde, const Flags & flags);
    void (*printdoc) (ostream & ost);
    
    NumProcInfo (const string & aname, int adim, 
		 shared_ptr<NumProc> (*acreator)(shared_ptr<PDE> pde, const Flags & flags),
		 void (*aprintdoc) (ostream & ost));
  };

  Array<shared_ptr<NumProcInfo>> npa;
public:
  NumProcs();

  void AddNumProc (const string & aname, 
		   shared_ptr<NumProc> (*acreator)(shared_ptr<PDE> pde, const Flags & flags),
		   void (*printdoc) (ostream & ost) = NumProc::PrintDoc);


  void AddNumProc (const string & aname, int dim, 
		   shared_ptr<NumProc> (*acreator)(shared_ptr<PDE> pde, const Flags & flags),
		   void (*printdoc) (ostream & ost) = NumProc::PrintDoc);

  const Array<shared_ptr<NumProcInfo>> & GetNumProcs() { return npa; }
  shared_ptr<NumProcInfo> GetNumProc(const string & name, int dim);

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
  
  static shared_ptr<NumProc> Create (shared_ptr<PDE> pde, const Flags & flags)
  {
    return make_shared<NP> (pde, flags);
  }
};


}
#endif
