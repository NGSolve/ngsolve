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
class NumProc : public NGS_Object
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
  /// virtual void Do();
  ///
  virtual void Do(LocalHeap & lh) = 0;
  ///
  virtual void PrintReport (ostream & ost);
  ///
  static void PrintDoc (ostream & ost);

  int GetCallPosition (void) const { return callposition;} 
};




/// Registered numprocs
class NumProcs
{
public:
  class NumProcInfo
  {
  public:
    string name;

    NumProc* (*creator)(PDE & pde, const Flags & flags);
    void (*printdoc) (ostream & ost);
    
    NumProcInfo (const string & aname,
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
  
  const Array<NumProcInfo*> & GetNumProcs() { return npa; }
  const NumProcInfo * GetNumProc(const string & name);

  void Print (ostream & ost) const;
};

 
extern NumProcs & GetNumProcs ();


}

#endif
