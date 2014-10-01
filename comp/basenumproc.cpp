#include <solve.hpp>
//#include <ctime>
// #include <parallelngs.hpp>


namespace ngcomp
{


  NumProc :: NumProc (PDE & apde, const Flags & flags) // , const int acallposition)
    : NGS_Object (apde.GetMeshAccess(int(flags.GetNumFlag("mesh",1))-1), "numproc"), 
      pde(apde) // , callposition(acallposition) 
  {
    if (flags.StringFlagDefined ("name"))
      SetName (flags.GetStringFlag ("name",""));
  }
  
  NumProc :: ~NumProc()
  {
    ;
  }

  void NumProc :: PrintReport (ostream & ost)
  {
    ost << "Base-class NumProc" << endl;
  }

  void NumProc :: PrintDoc (ostream & ost)
  {
    ost << "No documentation available" << endl;
  }




  NumProcs::NumProcInfo::
  NumProcInfo (const string & aname, int adim, 
	       shared_ptr<NumProc> (*acreator)(PDE & pde, const Flags & flags),
	       void (*aprintdoc) (ostream & ost) )
    : name(aname), dim(adim), creator(acreator), printdoc(aprintdoc)
  {
    ;
  }
  
  NumProcs :: NumProcs ()
  {
    ;
  }

  NumProcs :: ~NumProcs()
  {
    for (int i = 0; i < npa.Size(); i++)
      delete npa[i];
  }
  
  void NumProcs :: 
  AddNumProc (const string & aname,
	      shared_ptr<NumProc> (*acreator)(PDE & pde, const Flags & flags),
	      void (*printdoc) (ostream & ost) )
  {
    npa.Append (new NumProcInfo(aname, -1, acreator, printdoc));
  }

  void NumProcs :: 
  AddNumProc (const string & aname, int adim, 
	      shared_ptr<NumProc> (*acreator)(PDE & pde, const Flags & flags),
	      void (*printdoc) (ostream & ost) )
  {
    npa.Append (new NumProcInfo(aname, adim, acreator, printdoc));
  }



  const NumProcs::NumProcInfo * 
  NumProcs::GetNumProc(const string & name, int dim)
  {
    for (int i = 0; i < npa.Size(); i++)
      {
	if (name == npa[i]->name &&
	    ( (dim == npa[i]->dim) || (npa[i]->dim==-1) ))
	  return npa[i];
      }
    return 0;
  }

  void NumProcs :: Print (ostream & ost) const
  {
    ost << endl << "NumProcs:" << endl;
    ost <<         "---------" << endl;
    ost << setw(20) << "Name" << endl;
    for (int i = 0; i < npa.Size(); i++)
      ost << setw(20) << npa[i]->name << endl;
  }


 
  NumProcs & GetNumProcs ()
  {
    static NumProcs nps;
    return nps;
  }



}
