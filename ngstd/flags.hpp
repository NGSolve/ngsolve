#ifndef FILE_NGS_FLAGS
#define FILE_NGS_FLAGS


/**************************************************************************/
/* File:   flags.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                    */
/**************************************************************************/

namespace ngstd
{

  /** 
      A storage for command-line flags.
      The flag structure maintains string flags, numerical flags,
      define flags, string list flags, num list flags.
  */



  class NGS_DLL_HEADER Flags 
  {
    /// string flags
    SymbolTable<char *> strflags;
    /// numerical flags
    SymbolTable<double> numflags;
    /// define flags
    SymbolTable<int> defflags;
    /// string list flags
    SymbolTable<Array<char*>*> strlistflags;
    /// numerical list flags
    SymbolTable<Array<double>*> numlistflags;
  public:
    /// no flags
    Flags ();
    /// copy flags 
    Flags (const Flags & flags);
    ///
    Flags (string f1, string f2 = "", string f3 = "", string f4 = "", string f5 = "");
    /// delete mem
    ~Flags ();
  
    /// Deletes all flags
    void DeleteFlags ();
    /// Sets string flag, overwrite if exists
    Flags & SetFlag (const char * name, const char * val);
    /// Sets numerical flag, overwrite if exists
    Flags &  SetFlag (const char * name, double val);
    /// Sets boolean flag
    Flags &  SetFlag (const char * name);
    /// Sets string array flag
    Flags &  SetFlag (const char * name, const Array<char*> & val);
    /// Sets double array flag
    Flags &  SetFlag (const char * name, const Array<double> & val);
  
    /// Save flags to file
    void SaveFlags (const char * filename) const;
    /// write flags to stream
    void PrintFlags (ostream & ost) const;
    /// Load flags from file
    void LoadFlags (const char * filename);
    /**
       Set command line flag.
       Flag must be in form: -name=hello -val=0.5 -defflag 
       -names=[Joe,Jim] -values=[1,3,4]
    */
    void SetCommandLineFlag (const char * st);


    /// Returns string flag, default value if not exists
    const char * GetStringFlag (const char * name, const char * def) const;
    /// Returns string flag, default value if not exists
    string GetStringFlag (const char * name, const string & def) const;
    /// Returns numerical flag, default value if not exists
    double GetNumFlag (const char * name, double def) const;
    /// Returns address of numerical flag, null if not exists
    const double * GetNumFlagPtr (const char * name) const;
    /// Returns address of numerical flag, null if not exists
    double * GetNumFlagPtr (const char * name);
    /// Returns boolean flag
    int GetDefineFlag (const char * name) const;
    int GetDefineFlag (const string & name) const;
    /// Returns string list flag, empty array if not exist
    const Array<char*> & GetStringListFlag (const char * name) const;
    /// Returns num list flag, empty array if not exist
    const Array<double> & GetNumListFlag (const char * name) const;


    /// Test, if string flag is defined
    int StringFlagDefined (const char * name) const;
    /// Test, if num flag is defined
    int NumFlagDefined (const char * name) const;
    /// Test, if string list flag is defined
    int StringListFlagDefined (const char * name) const;
    /// Test, if num list flag is defined
    int NumListFlagDefined (const char * name) const;

    /// number of string flags
    int GetNStringFlags () const { return strflags.Size(); }
    /// number of num flags
    int GetNNumFlags () const { return numflags.Size(); }
    /// number of define flags
    int GetNDefineFlags () const { return defflags.Size(); }
    /// number of string-list flags
    int GetNStringListFlags () const { return strlistflags.Size(); }
    /// number of num-list flags
    int GetNNumListFlags () const { return numlistflags.Size(); }

    ///
    const char * GetStringFlag (int i, const char *& name) const
    { name = strflags.GetName(i); return strflags[i]; }
    double GetNumFlag (int i, const char *& name) const
    { name = numflags.GetName(i); return numflags[i]; }
    void GetDefineFlag (int i, const char *& name) const
    { name = defflags.GetName(i); }
    const Array<double> * GetNumListFlag (int i, const char *& name) const
    { name = numlistflags.GetName(i); return numlistflags[i]; }
    const Array<char*> * GetStringListFlag (int i, const char *& name) const
    { name = strlistflags.GetName(i); return strlistflags[i]; }
  };



  /// Print flags
  inline ostream & operator<< (ostream & s, const Flags & flags)
  {
    flags.PrintFlags (s);
    return s;
  }
}

  
#endif

