#ifndef FILE_NUMPROC
#define FILE_NUMPROC

/*********************************************************************/
/* File:   numproc.hh                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Jul. 2000                                             */
/*********************************************************************/

namespace ngsolve
{

/** 
    numerical procedures
*/

class PDE;
  
///
class NGS_DLL_HEADER NumProc : public NGS_Object
{
protected:
  ///
  PDE & pde;

  int callposition;
public:
  ///
  NumProc (PDE & apde, const int acallposition = 0);
  ///
  virtual ~NumProc();
  ///
  virtual void Do(LocalHeap & lh) = 0;
  ///
  virtual void PrintReport (ostream & ost);
  ///
  static void PrintDoc (ostream & ost);

  int GetCallPosition (void) const { return callposition;} 
};




/// Registered numprocs
class NGS_DLL_HEADER NumProcs
{
public:
  class NumProcInfo
  {
  public:
    string name;
    int dim;

    NumProc* (*creator)(PDE & pde, const Flags & flags);
    void (*printdoc) (ostream & ost);
    
    NumProcInfo (const string & aname, int adim, 
		 NumProc* (*acreator)(PDE & pde, const Flags & flags),
		 void (*aprintdoc) (ostream & ost));
  };

  Array<NumProcInfo*> npa;
public:
  NumProcs();
  ~NumProcs();  
  void AddNumProc (const string & aname, 
		   NumProc* (*acreator)(PDE & pde, const Flags & flags),
		   void (*printdoc) (ostream & ost) = NumProc::PrintDoc);


  void AddNumProc (const string & aname, int dim, 
		   NumProc* (*acreator)(PDE & pde, const Flags & flags),
		   void (*printdoc) (ostream & ost) = NumProc::PrintDoc);

  const Array<NumProcInfo*> & GetNumProcs() { return npa; }
  const NumProcInfo * GetNumProc(const string & name, int dim);

  void Print (ostream & ost) const;
};

 
extern NGS_DLL_HEADER NumProcs & GetNumProcs ();


template <typename NP>
class RegisterNumProc
{
public:
  RegisterNumProc (string label, int dim = -1)
  {
    GetNumProcs().AddNumProc (label, dim, Create, NP::PrintDoc);
    // cout << "register numproc '" << label << "'" << endl;
  }
  
  static NumProc * Create (PDE & pde, const Flags & flags)
  {
    return new NP (pde, flags);
  }
};


}

#endif
