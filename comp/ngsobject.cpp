#include <comp.hpp>

namespace ngcomp {

void NGS_Object :: DefineStringFlag(const char* s, const char* val) 
{
  if (flaglist.StringFlagDefined(s))
  {
      cerr << "WARNING in NGS_Object :: DefineStringFlag: stringflag '" << s << "' already defined" << endl;
      return;
  }    
  flaglist.SetFlag(s,val); 
}

void NGS_Object :: DefineNumFlag(const char* s, double val) 
{
  if (flaglist.NumFlagDefined(s))
  {
    cerr << "WARNING in NGS_Object :: DefineNumFlag: numflag '" << s << "' already defined" << endl;
    return;
  }    
  flaglist.SetFlag(s,val); 
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
  ARRAY<char*> as(0);
  flaglist.SetFlag(s,as); 
}

void NGS_Object :: DefineNumListFlag(const char* s) 
{
  if (flaglist.NumListFlagDefined(s))
  {
    cerr << "WARNING in NGS_Object :: DefineNumListFlag: numlistflag '" << s << "' already defined" << endl;
    return;
  }    
  ARRAY<double> as(0);
  flaglist.SetFlag(s,as); 
}

int NGS_Object :: ParseFlags(const Flags& flags)
{
  int ret=0; // ok
  
  // parse string flags
  const char* s;
  for (int i=0; i<flags.GetNStringFlags(); i++)
  {
    flags.GetStringFlag(i, s);
    if (!flaglist.StringFlagDefined(s))
    {
      cerr << "WARNING in NGS_Object :: ParseFlags(): stringflag '" << s << "' not defined for object " << name  << endl;
      ret++;
    }
  }
  // parse num flags
  for (int i=0; i<flags.GetNNumFlags(); i++)
  {
    flags.GetNumFlag(i, s);
    if (!flaglist.NumFlagDefined(s))
    {
      cerr << "WARNING in NGS_Object :: ParseFlags(): numflag '" << s << "' not defined for object " << name  << endl;
      ret++;
    }
  }
  // parse define flags
  for (int i=0; i<flags.GetNDefineFlags(); i++)
  {
    flags.GetDefineFlag(i, s);
    if (!flaglist.GetDefineFlag(s))
    {
      cerr << "WARNING in NGS_Object :: ParseFlags(): defineflag '" << s << "' not defined for object " << name  << endl;
      ret++;
    }
  }
  // parse stringlist flags
  for (int i=0; i<flags.GetNStringListFlags(); i++)
  {
    flags.GetStringListFlag(i, s);
    if (!flaglist.StringListFlagDefined(s))
    {
      cerr << "WARNING in NGS_Object :: ParseFlags(): stringlistflag '" << s << "' not defined for object " << name  << endl;
      ret++;
    }
  }
  // parse numlist flags
  for (int i=0; i<flags.GetNNumListFlags(); i++)
  {
    flags.GetNumListFlag(i, s);
    if (!flaglist.NumListFlagDefined(s))
    {
      cerr << "WARNING in NGS_Object :: ParseFlags(): numlistflag '" << s << "' not defined for object " << name  << endl;
      ret++;
    }
  }

  return ret; // number of undefined flags
}
 
} // namespace
