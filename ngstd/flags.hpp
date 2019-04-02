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
    SymbolTable<string> strflags;
    /// numerical flags
    SymbolTable<double> numflags;
    /// define flags
    SymbolTable<bool> defflags;
    /// string list flags
    SymbolTable<shared_ptr<Array<string>>> strlistflags;
    /// numerical list flags
    SymbolTable<shared_ptr<Array<double>>> numlistflags;
    /// flags list flags
    SymbolTable<Flags> flaglistflags;
  public:
    /// no flags
    Flags ();
    /// copy flags 
    Flags (const Flags & flags);
    /// steal flags
    Flags (Flags && flags);
    ///
    Flags (std::initializer_list<string> list);
    /// 
    Flags (string f1, string f2 = "", string f3 = "", string f4 = "", string f5 = "");
    /// delete mem
    ~Flags ();

    Flags & operator= (const Flags & f2) { *this = f2; return *this; };
    Flags & operator= (Flags && f2) = default;

    void DoArchive(Archive& ar);

    void Update(const Flags& other);
  
    /// Deletes all flags
    void DeleteFlags ();

    /// Sets string flag, overwrite if exists
    Flags & SetFlag (const char * name, const string & val);
    /// Sets string flag, overwrite if exists
    Flags & SetFlag (const char * name, const char * str)
    { return SetFlag (name, string(str)); }
    /// Sets numerical flag, overwrite if exists
    Flags & SetFlag (const char * name, double val) &;
    /// Sets numerical flag, overwrite if exists    
    Flags & SetFlag (const char * name, int val)
    { return SetFlag (name, double(val)); }
    /// Sets boolean flag
    Flags & SetFlag (const char * name, bool b = true) &;
    /// Sets numerical flag, overwrite if exists
    Flags & SetFlag (const char * name, Flags & val) &;

    /// Sets string flag, overwrite if exists
    Flags & SetFlag (const string & name, const string & val);
    Flags & SetFlag (const string & name, const char * str)
    { return SetFlag (name, string(str)); }
    /// Sets numerical flag, overwrite if exists
    Flags &  SetFlag (const string & name, double val);
    /// Sets numerical flag, overwrite if exists    
    Flags &  SetFlag (const string & name, int val)
    { return SetFlag (name, double(val)); }
    /// Sets boolean flag
    Flags &  SetFlag (const string & name, bool b = true);
    /// Sets numerical flag, overwrite if exists
    Flags &  SetFlag (const string & name, Flags & val);
    /// Sets string array flag
    Flags &  SetFlag (const string & name, const Array<string> & val);
    /// Sets double array flag
    Flags &  SetFlag (const string & name, const Array<double> & val);


    Flags SetFlag (const char * name, bool b = true) &&;
    Flags SetFlag (const char * name, double val) &&;



    /// Save flags to file
    void SaveFlags (const char * filename) const;
    void SaveFlags (ostream & str) const;
    /// write flags to stream
    void PrintFlags (ostream & ost) const;
    /// Load flags from file
    void LoadFlags (const char * filename, SymbolTable<Flags> * sf = nullptr);
    void LoadFlags (istream & str, SymbolTable<Flags> * sf = nullptr);
    /**
       Set command line flag.
       Flag must be in form: -name=hello -val=0.5 -defflag 
       -names=[Joe,Jim] -values=[1,3,4] -solverflags=*abc
    */
    void SetCommandLineFlag (const char * st, SymbolTable<Flags> * sf = nullptr);
    /// Returns string flag, default value if not exists
    string GetStringFlag (const string & name, const char * def) const;
    /// Returns string flag, default value if not exists
    string GetStringFlag (const string & name, string def = "") const;
    /// Returns numerical flag, default value if not exists
    double GetNumFlag (const string & name, double def) const;
    /// Returns address of numerical flag, null if not exists
    const double * GetNumFlagPtr (const string & name) const;
    /// Returns address of numerical flag, null if not exists
    double * GetNumFlagPtr (const string & name);
    /// Returns boolean flag
    // int GetDefineFlag (const char * name) const;
    bool GetDefineFlag (const string & name) const throw();
    xbool GetDefineFlagX (const string & name) const throw();
    /// Returns string list flag, empty array if not exist
    const Array<string> & GetStringListFlag (const string & name) const;
    /// Returns num list flag, empty array if not exist
    const Array<double> & GetNumListFlag (const string & name) const;
    /// Returns flag list flag, empty flag if not exist
    const Flags & GetFlagsFlag (const string & name) const;


    /// Test, if string flag is defined
    bool StringFlagDefined (const string & name) const;
    /// Test, if num flag is defined
    bool NumFlagDefined (const string & name) const;
    /// Test, if num flag is defined
    bool FlagsFlagDefined (const string & name) const;
    /// Test, if string list flag is defined
    bool StringListFlagDefined (const string & name) const;
    /// Test, if num list flag is defined
    bool NumListFlagDefined (const string & name) const;

    /// number of string flags
    int GetNStringFlags () const { return strflags.Size(); }
    /// number of num flags
    int GetNNumFlags () const { return numflags.Size(); }
    /// number of num flags
    int GetNFlagsFlags () const { return flaglistflags.Size(); }
    /// number of define flags
    int GetNDefineFlags () const { return defflags.Size(); }
    /// number of string-list flags
    int GetNStringListFlags () const { return strlistflags.Size(); }
    /// number of num-list flags
    int GetNNumListFlags () const { return numlistflags.Size(); }

    ///
    const string & GetStringFlag (int i, string & name) const
    { name = strflags.GetName(i); return strflags[i]; }
    double GetNumFlag (int i, string & name) const
    { name = numflags.GetName(i); return numflags[i]; }
    bool GetDefineFlag (int i, string & name) const
    { name = defflags.GetName(i); return defflags[i]; }
    const shared_ptr<Array<double>> GetNumListFlag (int i, string & name) const
    { name = numlistflags.GetName(i).c_str(); return numlistflags[i]; }
    const shared_ptr<Array<string>> GetStringListFlag (int i, string & name) const
    { name = strlistflags.GetName(i); return strlistflags[i]; }
    const Flags & GetFlagsFlag (int i, string & name) const
    { name = flaglistflags.GetName(i); return flaglistflags[i]; }
  };

  /// Print flags
  inline ostream & operator<< (ostream & s, const Flags & flags)
  {
    flags.PrintFlags (s);
    return s;
  }
}

  
#endif

