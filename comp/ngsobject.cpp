#include "ngsobject.hpp"


namespace ngcomp {

  size_t NGS_Object :: global_timestamp = 1;
  
  NGS_Object :: ~NGS_Object () { ; }

  void NGS_Object :: DefineStringFlag(const char* s) //, const char* val) 
  {
    if (flaglist.StringFlagDefined(s))
      {
        cerr << "WARNING in NGS_Object :: DefineStringFlag: stringflag '" << s << "' already defined" << endl;
        return;
      }    
    flaglist.SetFlag(s, "" /* val */ ); 
  }

  void NGS_Object :: DefineNumFlag(const char* s) // , double val) 
{
  if (flaglist.NumFlagDefined(s))
  {
    cerr << "WARNING in NGS_Object :: DefineNumFlag: numflag '" << s << "' already defined" << endl;
    return;
  }    
  flaglist.SetFlag(s, 0.0 /* val */ ); 
}

void NGS_Object :: DefineDefineFlag(const char* s) 
{
  if (flaglist.GetDefineFlag(s))
  {
    cerr << "WARNING in NGS_Object :: DefineFlag: defineflag '" << s << "' already defined" << endl;
    return;
  }    
  flaglist.SetFlag(s); 
}

void NGS_Object :: DefineStringListFlag(const char* s) 
{
  if (flaglist.StringListFlagDefined(s))
  {
    cerr << "WARNING in NGS_Object :: DefineStringListFlag: stringlistflag '" << s << "' already defined" << endl;
    return;
  }    
  Array<string> as(0);
  flaglist.SetFlag(s,as); 
}

void NGS_Object :: DefineNumListFlag(const char* s) 
{
  if (flaglist.NumListFlagDefined(s))
  {
    cerr << "WARNING in NGS_Object :: DefineNumListFlag: numlistflag '" << s << "' already defined" << endl;
    return;
  }    
  Array<double> as(0);
  flaglist.SetFlag(s,as); 
}

int NGS_Object :: CheckFlags(const Flags& flags)
{
  int ret=0; // ok
  
  // parse string flags
  string s;
  for (int i=0; i<flags.GetNStringFlags(); i++)
  {
    flags.GetStringFlag(i, s);
    if (!flaglist.StringFlagDefined(s))
    {
      cerr << IM(1) << "WARNING in NGS_Object :: CheckFlags(): stringflag '" << s << "' not defined for object " << name  << endl;
      ret++;
    }
  }
  // parse num flags
  for (int i=0; i<flags.GetNNumFlags(); i++)
  {
    flags.GetNumFlag(i, s);
    if (!flaglist.NumFlagDefined(s))
    {
      cerr << IM(1) << "WARNING in NGS_Object :: CheckFlags(): numflag '" << s << "' not defined for object " << name  << endl;
      ret++;
    }
  }
  // parse define flags
  for (int i=0; i<flags.GetNDefineFlags(); i++)
  {
    flags.GetDefineFlag(i, s);
    if (!flaglist.GetDefineFlag(s))
    {
      cerr << IM(1) << "WARNING in NGS_Object :: CheckFlags(): defineflag '" << s << "' not defined for object " << name  << endl;
      ret++;
    }
  }
  // parse stringlist flags
  for (int i=0; i<flags.GetNStringListFlags(); i++)
  {
    flags.GetStringListFlag(i, s);
    if (!flaglist.StringListFlagDefined(s))
    {
      cerr << IM(1) << "WARNING in NGS_Object :: CheckFlags(): stringlistflag '" << s << "' not defined for object " << name  << endl;
      ret++;
    }
  }
  // parse numlist flags
  for (int i=0; i<flags.GetNNumListFlags(); i++)
  {
    flags.GetNumListFlag(i, s);
    if (!flaglist.NumListFlagDefined(s))
    {
      cerr << IM(1) << "WARNING in NGS_Object :: CheckFlags(): numlistflag '" << s << "' not defined for object " << name  << endl;
      ret++;
    }
  }

  return ret; // number of undefined flags
}


string NGS_Object :: GetClassName () const
{
  return typeid(*this).name();
}

void NGS_Object :: PrintReport (ostream & ost) const
{
  ost << typeid(*this).name();
}

  
  Array<MemoryUsage> NGS_Object :: GetMemoryUsage () const
  {
    cout << "MemoryUsage not overloaded for class " << GetClassName() << endl;
    return Array<MemoryUsage>();
  }
  

  string DocInfo :: GetPythonDocString() const
  {
    string docstring = short_docu + "\n\n" + long_docu;

    if (arguments.size())
      {
        docstring += "\nKeyword arguments can be:\n\n";
        for (auto & flagdoc : arguments)
          docstring += get<0> (flagdoc) + ": " + get<1> (flagdoc) + "\n";
      }
    return docstring;
  }

  
} // namespace
