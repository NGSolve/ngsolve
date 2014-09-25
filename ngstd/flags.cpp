/**************************************************************************/
/* File:   flags.cpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. Oct. 96                                                    */
/**************************************************************************/

#include <ngstd.hpp>

#ifdef WIN32
#include <float.h>
#endif


namespace ngstd
{
  Flags :: Flags ()
  {
    ;
  }

  Flags :: Flags (const Flags & flags)
  {
    string name;
    for (int i = 0; i < flags.GetNStringFlags(); i++)
      {
	string str = flags.GetStringFlag (i, name);
	SetFlag (name, str);
      }
    for (int i = 0; i < flags.GetNNumFlags(); i++)
      {
	double val = flags.GetNumFlag (i, name);
	SetFlag (name, val);
      }
    for (int i = 0; i < flags.GetNDefineFlags(); i++)
      {
	flags.GetDefineFlag (i, name);
	SetFlag (name);
      }
    for (int i = 0; i < flags.GetNNumListFlags(); i++)
      {
	auto numa = flags.GetNumListFlag (i, name);
	SetFlag (name, *numa);
      }
    for (int i = 0; i < flags.GetNStringListFlags(); i++)
      {
	auto stra = flags.GetStringListFlag (i, name);
	SetFlag (name, *stra);
      }
  }

  Flags :: Flags (std::initializer_list<string> list)
  {
    for (auto i = list.begin(); i < list.end(); i++)
      SetCommandLineFlag ((string("-")+*i).c_str());      
    // data[cnt] = *i;
    // for (int i = 0; i < list.size(); i++)
    //   SetCommandLineFlag ((string("-")+list[i]).c_str());      
  }


  Flags :: Flags (string f1, string f2, string f3, string f4, string f5)
  {
    SetCommandLineFlag ((string("-")+f1).c_str());
    if (f2.length()) SetCommandLineFlag ( (string("-")+f2).c_str() );
    if (f3.length()) SetCommandLineFlag ( (string("-")+f3).c_str() );
    if (f4.length()) SetCommandLineFlag ( (string("-")+f4).c_str() );
    if (f5.length()) SetCommandLineFlag ( (string("-")+f5).c_str() );
  }
  
  Flags :: ~Flags ()
  {
    DeleteFlags ();
  }
  
  void Flags :: DeleteFlags ()
  {
    // for (int i = 0; i < strflags.Size(); i++)
    // delete [] strflags[i];
    strflags.DeleteAll();
    numflags.DeleteAll();
    defflags.DeleteAll();
    /*
    for(int i=0; i<strlistflags.Size(); i++)
      {
        // for (int j = 0; j < strlistflags[i]->Size(); j++)
        // delete [] (*strlistflags[i])[j];
        delete strlistflags[i];
      }
    */
    strlistflags.DeleteAll();
    // for(int i=0; i<numlistflags.Size(); i++)
    // delete numlistflags[i];
    numlistflags.DeleteAll();
  }
  
  Flags & Flags :: SetFlag (const string & name, const string & val)
  {
    // char * hval = new char[strlen (val) + 1];
    // strcpy (hval, val);
    strflags.Set (name, val);
    return *this;
  }
  
  Flags & Flags :: SetFlag (const string & name, double val) &
  {
    numflags.Set (name, val);
    return *this;
  }
  
  Flags & Flags :: SetFlag (const string & name)
  {
    defflags.Set (name, 1);
    return *this;
  }


  Flags & Flags :: SetFlag (const string & name, const Array<string> & val)
  {
    auto strarray = make_shared<Array<string>>(val);
      /*
    for (int i = 0; i < val.Size(); i++)
      {
	strarray->Append (new char[strlen(val[i])+1]);
	strcpy (strarray->Last(), val[i]);
      }
      */
    strlistflags.Set (name, strarray);    
    return *this;
  }

  Flags & Flags :: SetFlag (const string & name, const Array<double> & val)
  {
    // Array<double> * numarray = new Array<double>(val);
    auto numarray = make_shared<Array<double>> (val);

    numlistflags.Set (name, numarray);
    return *this;
  }



  string Flags :: GetStringFlag (const string & name, const char * def) const
  {
    if (strflags.Used (name))
      return strflags[name];
    else
      {
        if (!def) return string("");
        return def;
      }
  }

  string Flags :: GetStringFlag (const string & name, string def) const
  {
    if (strflags.Used (name))
      return strflags[name];
    else
      return def;
  }


  double Flags :: GetNumFlag (const string & name, double def) const
  {
    if (numflags.Used (name))
      return numflags[name];
    else
      return def;
  }
  
  const double * Flags :: GetNumFlagPtr (const string & name) const
  {
    if (numflags.Used (name))
      return & ((SymbolTable<double>&)numflags)[name];
    else
      return NULL;
  }
  
  double * Flags :: GetNumFlagPtr (const string & name) 
  {
    if (numflags.Used (name))
      return & ((SymbolTable<double>&)numflags)[name];
    else
      return NULL;
  }

  /*
  int Flags :: GetDefineFlag (const char * name) const
  {
    return defflags.Used (name);
  }
  */
  int Flags :: GetDefineFlag (const string & name) const
  {
    return defflags.Used (name);
  }


  const Array<string> & 
  Flags :: GetStringListFlag (const string & name) const
  {
    if (strlistflags.Used (name))
      return *strlistflags[name];
    else
      {
	static Array<string> hstra(0);
	return hstra;
      }
  }

  const Array<double> & 
  Flags ::GetNumListFlag (const string & name) const
  {
    if (numlistflags.Used (name))
      return *numlistflags[name];
    else
      {
	static Array<double> hnuma(0);
	return hnuma;
      }
  }


  bool Flags :: StringFlagDefined (const string & name) const
  {
    return strflags.Used (name);
  }

  bool Flags :: NumFlagDefined (const string &name) const
  {
    return numflags.Used (name);
  }
  
  bool Flags :: StringListFlagDefined (const string & name) const
  {
    return strlistflags.Used (name);
  }

  bool Flags :: NumListFlagDefined (const string & name) const
  {
    return numlistflags.Used (name);
  }


  void Flags :: SaveFlags (const char * filename) const 
  {
    ofstream outfile (filename);
  
    for (int i = 0; i < strflags.Size(); i++)
      outfile << strflags.GetName(i) << " = " << strflags[i] << endl;
    for (int i = 0; i < numflags.Size(); i++)
      outfile << numflags.GetName(i) << " = " << numflags[i] << endl;
    for (int i = 0; i < defflags.Size(); i++)
      outfile << defflags.GetName(i) << endl;
  }
 


  void Flags :: PrintFlags (ostream & ost) const 
  {
    for (int i = 0; i < strflags.Size(); i++)
      ost << strflags.GetName(i) << " = " << strflags[i] << endl;
    for (int i = 0; i < numflags.Size(); i++)
      ost << numflags.GetName(i) << " = " << numflags[i] << endl;
    for (int i = 0; i < defflags.Size(); i++)
      ost << defflags.GetName(i) << endl;
    for (int i = 0; i < strlistflags.Size(); i++)
      ost << strlistflags.GetName(i) << " = " << strlistflags[i] << endl;
    for (int i = 0; i < numlistflags.Size(); i++)
      ost << numlistflags.GetName(i) << " = " << numlistflags[i] << endl;
  }


  void Flags :: LoadFlags (const char * filename) 
  {
    char name[100], str[100];
    char ch;
    double val;
    ifstream infile(filename);

    while (infile.good())
      {
	infile >> name;
	if (strlen (name) == 0) break;

	if (name[0] == '/' && name[1] == '/')
	  {
	    ch = 0;
	    while (ch != '\n' && infile.good())
	      {
		ch = infile.get();
	      }
	    continue;
	  }

	ch = 0;
	infile >> ch;
	if (ch != '=')
	  {
	    infile.putback (ch);
	    SetFlag (name);
	  }
	else
	  {
	    infile >> val;
	    if (!infile.good())
	      {
		infile.clear();
		infile >> str;
		SetFlag (name, str);
	      }
	    else
	      {
		SetFlag (name, val);
	      }
	  }
      }
  }

  Archive & operator & (Archive & archive, Flags & flags)
  {
    cout << "Archive Flags" << endl;
    cout << flags << endl;
    archive & flags.strflags;
    archive & flags.numflags;
    archive & flags.defflags;
    archive & flags.numlistflags;
    archive & flags.strlistflags;
    return archive;
  }

  void Flags :: SetCommandLineFlag (const char * st)
  {
    //cout << "SetCommandLineFlag: flag = " << st << endl;
    istringstream inst( (char *)st);

    char name[100];
    double val;


    if (st[0] != '-')
      {
	cerr << "flag must start with '-'" << endl;
	return;
      }

    // flag with double --
    if (st[1] == '-') st++;
  
    const char * pos = strchr (st, '=');
    const char * posbrack = strchr (st, '[');

    if (!pos)
      {
	//      (cout) << "Add def flag: " << st+1 << endl;
	SetFlag (st+1);
      }
    else
      {
	//cout << "pos = " << pos << endl;

	strncpy (name, st+1, (pos-st)-1);
	name[pos-st-1] = 0;

	//cout << "name = " << name << endl;

	pos++;
	char * endptr = NULL;
	val = strtod (pos, &endptr);

        /*
        cout << "val = " << val << endl;
        cout << "isfinite = " << std::isfinite (val) << endl;
        cout << "isinf = " << std::isinf (val) << endl;
        cout << "pos = " << pos << ", endpos = " << endptr << endl;
        */
        if (endptr != pos && !isfinite (val))
          endptr = const_cast<char *>(pos);          

        /*
#ifdef WIN32
	if(endptr != pos && !_finite(val))
	  endptr = const_cast<char *>(pos);
#else
#ifdef MACOS
	if(endptr != pos && (__isnand(val) || __isinfd(val)))
	  endptr = const_cast<char *>(pos);
#else
#ifdef SUN
#else
	if(endptr != pos && (std::isnan(val) || std::isinf(val)))
	  endptr = const_cast<char *>(pos);
#endif
#endif
#endif
        */
	
	//cout << "val = " << val << endl;

	if (!posbrack)
	  {
	    if (endptr == pos)
	      {
		// string-flag
		//(cout) << "Add String Flag: " << name << " = " << pos << endl;
		SetFlag (name, pos);
	      }
	    else
	      {
		// num-flag
		//(cout) << "Add Num Flag: " << name << " = " << val << endl;
		SetFlag (name, val);
	      }
	  }
	else
	  {
	    // list-flag
	    char hc;
	    double val;

	    val = strtod (posbrack+1, &endptr);
	    if (endptr != posbrack+1)
	      {
		Array<double> values;
		
		istringstream ist(posbrack);
		ist >> hc;   // '['
		ist >> val;
		while (ist.good())
		  {
		    values.Append (val);
		    ist >> hc;  // ','
		    ist >> val;
		  }
		SetFlag (name, values);
	      }
	    else
	      {
                // to be cleand up ...
		Array<char *> strs;

		posbrack++;
		char * hstr = new char[strlen(posbrack)+1];
		strcpy (hstr, posbrack);
		
		char * chp = hstr;

		bool start = 1;
		while (*chp && *chp != ']')
		  {
		    if (start)
		      strs.Append (chp);
		    start = 0;
		    if (*chp == ',')
		      {
			*chp = 0;
			start = 1;
		      }
		    chp++;
		  }
		*chp = 0;

                Array<string> strings;
                for (int i = 0; i < strs.Size(); i++)
                  strings.Append (string (strs[i]));
		SetFlag (name, strings);
                delete [] hstr;
	      }
	  }
      }
  }
}
