#include <solve.hpp>
//#include <ctime>
// #include <parallelngs.hpp>


namespace ngcomp
{


  NumProc :: NumProc (weak_ptr<PDE> apde, const Flags & flags) 
    : NGS_Object (apde.lock()->GetMeshAccess(int(flags.GetNumFlag("mesh",1))-1),flags,
                  "numproc"), pde(apde)
  {
    if (flags.StringFlagDefined ("name"))
      SetName (flags.GetStringFlag ("name",""));
  }

  NumProc :: NumProc (const Flags & flags) 
    : NGS_Object (nullptr, flags, "numproc")
  {
    if (flags.StringFlagDefined ("name"))
      SetName (flags.GetStringFlag ("name",""));
  }

  
  NumProc :: ~NumProc()
  {
    ;
  }

  void NumProc :: PrintReport (ostream & ost) const
  {
    ost << "Base-class NumProc" << endl;
  }

  void NumProc :: PrintDoc (ostream & ost)
  {
    ost << "No documentation available" << endl;
  }




  NumProcs::NumProcInfo::
  NumProcInfo (const string & aname, int adim, 
	       shared_ptr<NumProc> (*acreator)(shared_ptr<PDE> pde, const Flags & flags),
	       void (*aprintdoc) (ostream & ost) )
    : name(aname), dim(adim), creator(acreator), printdoc(aprintdoc)
  {
    ;
  }
  
  NumProcs :: NumProcs ()
  {
    ;
  }

  void NumProcs :: 
  AddNumProc (const string & aname,
	      shared_ptr<NumProc> (*acreator)(shared_ptr<PDE> pde, const Flags & flags),
	      void (*printdoc) (ostream & ost) )
  {
    npa.Append (make_shared<NumProcInfo> (aname, -1, acreator, printdoc));
  }

  void NumProcs :: 
  AddNumProc (const string & aname, int adim, 
	      shared_ptr<NumProc> (*acreator)(shared_ptr<PDE> pde, const Flags & flags),
	      void (*printdoc) (ostream & ost) )
  {
    npa.Append (make_shared<NumProcInfo> (aname, adim, acreator, printdoc));
  }


  shared_ptr<NumProcs::NumProcInfo>
  NumProcs::GetNumProc(const string & name, int dim)
  {
    for (auto & info : npa)
      if ( (name == info->name) && ( (dim == info->dim) || (info->dim==-1) ))
        return info;

    return nullptr;
  }

  void NumProcs :: Print (ostream & ost) const
  {
    ost << endl << "NumProcs:" << endl;
    ost <<         "---------" << endl;
    ost << setw(20) << "Name" << endl;

    for (auto & info : npa)
      ost << setw(20) << info->name << endl;
  }

 
  NumProcs & GetNumProcs ()
  {
    static NumProcs nps;
    return nps;
  }



}
