#ifndef FILE_FLAGS
#define FILE_FLAGS


/**************************************************************************/
/* File:   flags.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                   */
/**************************************************************************/

namespace netgen
{

/** 
   Flag - Table.
   A flag table maintains string variables, numerical 
   variables and boolean flags.
*/
class Flags 
{
  ///
  SYMBOLTABLE<char *> strflags;
  ///
  SYMBOLTABLE<double> numflags;
  ///
  SYMBOLTABLE<int> defflags;
  ///
  SYMBOLTABLE<Array<char*>*> strlistflags;
  ///
  SYMBOLTABLE<Array<double>*> numlistflags;
public:
  ///
  DLL_HEADER Flags ();
  ///
  DLL_HEADER ~Flags ();
  
  /// Deletes all flags
  DLL_HEADER void DeleteFlags ();
  /// Sets string flag, overwrite if exists
  DLL_HEADER void SetFlag (const char * name, const char * val);
  /// Sets numerical flag, overwrite if exists
  DLL_HEADER void SetFlag (const char * name, double val);
  /// Sets boolean flag
  DLL_HEADER void SetFlag (const char * name);
  /// Sets string arary falg
  DLL_HEADER void SetFlag (const char * name, const Array<char*> & val);
  /// Sets double array flag
  DLL_HEADER void SetFlag (const char * name, const Array<double> & val);
  
  /// Save flags to file
  DLL_HEADER void SaveFlags (const char * filename) const;
  /// write flags to stream
  DLL_HEADER void PrintFlags (ostream & ost) const;
  /// Load flags from file
  DLL_HEADER void LoadFlags (const char * filename);
  /// set flag of form -name=hello -val=0.5 -defined
  DLL_HEADER void SetCommandLineFlag (const char * st);

  /// Returns string flag, default value if not exists
  DLL_HEADER const char * GetStringFlag (const char * name, const char * def) const;
  /// Returns numerical flag, default value if not exists
  DLL_HEADER double GetNumFlag (const char * name, double def) const;
  /// Returns address of numerical flag, null if not exists
  DLL_HEADER const double * GetNumFlagPtr (const char * name) const;
  /// Returns address of numerical flag, null if not exists
  DLL_HEADER double * GetNumFlagPtr (const char * name);
  /// Returns boolean flag
  DLL_HEADER bool GetDefineFlag (const char * name) const;
  /// Returns string list flag, empty array if not exist
  DLL_HEADER const Array<char*> & GetStringListFlag (const char * name) const;
  /// Returns num list flag, empty array if not exist
  DLL_HEADER const Array<double> & GetNumListFlag (const char * name) const;


  /// Test, if string flag is defined
  DLL_HEADER bool StringFlagDefined (const char * name) const;
  /// Test, if num flag is defined
  DLL_HEADER bool NumFlagDefined (const char * name) const;
  /// Test, if string list flag is defined
  DLL_HEADER bool StringListFlagDefined (const char * name) const;
  /// Test, if num list flag is defined
  DLL_HEADER bool NumListFlagDefined (const char * name) const;
};

}
  
#endif

