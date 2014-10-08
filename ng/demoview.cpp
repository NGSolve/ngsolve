/*********************************************************************/
/* File:   demoview.cpp                                              */
/* Author: Robert, Joachim                                           */
/* Date:   6. Mar. 2003                                              */
/*********************************************************************/


#include <mystdlib.h>


//#include <iostream.h>
#include <myadt.hpp>
#include <linalg.hpp>
#include <gprim.hpp>
#include <csg.hpp> 
#include <geometry2d.hpp>
#include <stlgeom.hpp>
#include <meshing.hpp>
#include "inctcl.hpp"
#include <visual.hpp>

namespace netgen {  
  extern VisualScene *vs;
#include "demoview.hpp"


  /*
    static demokwstruct defkw[] =
    {
    { TOK_TIME,     "t" },
    { TOK_CAMPOS,   "camerapos" },
    { TOK_CAMPOINT, "camerapointto" },
    { TOK_CAMUP,    "cameraup" }
    };
  */


  static demoview_kwstruct demoview_defkw[] =
    {
      { DTOK_TIME,     "t" },
      { DTOK_CAMPOS,   "camerapos" },
      { DTOK_CAMPOINT, "camerapointto" },
      { DTOK_CAMUP,    "cameraup" }
    };


  DemoScanner :: DemoScanner (ifstream & ascanin)
  {
    scanin = &ascanin;
    token = DTOK_END;
    num_value = 0;
    linenum = 1;
  }



  void DemoScanner :: ReadNext ()
  {
    char ch;
  

    // whitespaces ueberspringen
    do
      { 
	scanin->get(ch);

	if (ch == '\n') 
	  linenum++;

	// end of file reached
	if (scanin->eof())
	  {
	    token = DTOK_END;
	    return;
	  }

	// skip comment line
	if (ch == '#')
	  {
	    while (ch != '\n')
	      {
		scanin->get(ch);
		if (scanin->eof())
		  {
		    token = DTOK_END;
		    return;
		  }
	      }
	    linenum++;
	  }	
      }
    while (isspace(ch));
  
    switch (ch)
      {
      case '(': case ')': 
      case '[': case ']': 
      case '-': case ':':
      case '=': case ',':
      case ';': case '+':
	{
	  token = DEMOVIEW_TOKEN_TYPE (ch);
	  break;
	}
  
      default:
	{
	  if (isdigit (ch) || ch == '.')
	    {
	      scanin->putback (ch);
	      (*scanin) >> num_value;
	      token = DTOK_NUM;
	      return;
	    }

	  if (isalpha (ch))
	    {
	      string_value = string (1, ch);
	      scanin->get(ch);
	      while (isalnum(ch))
		{
		  string_value += ch;
		  scanin->get(ch);
		}
	      scanin->putback (ch);
	    }

	  int nr = 0;
	  while (demoview_defkw[nr].kw)
	    {
	      if (string_value == demoview_defkw[nr].name)
		{
		  token = demoview_defkw[nr].kw;
		  return;
		}
	      nr++;
	    }

	  token = DTOK_STRING;
	}
      }
  }



  void DemoScanner :: Error (const string & err)
  {
    stringstream errstr;
    errstr << "Parsing error in line " << linenum << ": " << endl << err << endl;
    throw string(errstr.str());
  }



  void ParseChar (DemoScanner & scan, char ch)
  {
    char str[2];
    str[0] = ch;
    str[1] = 0;
    if (scan.GetToken() != DEMOVIEW_TOKEN_TYPE(ch)) 
      scan.Error (string ("token '") + string(str) + string("' expected"));
    scan.ReadNext();
  }
  


  double ParseNumber(DemoScanner & scan)
  {
    if (scan.GetToken() == '-')
      {
	scan.ReadNext();
	return -ParseNumber (scan);
      }
    if (scan.GetToken() != DTOK_NUM) scan.Error ("number expected");
    double val = scan.GetNumValue();
    scan.ReadNext();
    return val;	
  }


  Vec<3> ParseVector (DemoScanner & scan)
  {
    Vec<3> s;

    s(0) = ParseNumber (scan);
    ParseChar (scan, ',');

    s(1) = ParseNumber (scan);
    ParseChar (scan, ',');

    s(2) = ParseNumber (scan);

    return s;
  }


  void ParseConstLineOrSpline (DemoScanner & scan, double * t, Vec<3> * s)
  {
    int np = 1;
  
    scan.ReadNext();
    ParseChar (scan, '(');
  
    t[0] = ParseNumber (scan)*1000;
    ParseChar (scan, ':');

    s[0] = ParseVector (scan);
  
    if (scan.GetToken() != DTOK_RP && 
	scan.GetToken() != DTOK_SEMICOLON)
      scan.Error (") or ; expected");	   
  
    if (scan.GetToken() == DTOK_SEMICOLON)
      {
	np++;
      
	scan.ReadNext();

	t[1] = ParseNumber (scan)*1000;
	ParseChar (scan, ':');
      
	s[1] = ParseVector (scan);
      
	if (scan.GetToken() != DTOK_RP && 
	    scan.GetToken() != DTOK_SEMICOLON)
	  scan.Error (") or ; expected");	   
      
	if (scan.GetToken() == DTOK_SEMICOLON)
	  {
	    np++;
	  
	    scan.ReadNext();
	  
	    t[2] = ParseNumber (scan)*1000;
	    ParseChar (scan, ':');
	  
	    s[2] = ParseVector (scan);
	  
	    ParseChar (scan, ')');
	    ParseChar (scan, ';');
	  }
	else if (scan.GetToken() == DTOK_RP)
	  {
	    scan.ReadNext();
	    ParseChar (scan, ';');
	  }
      }
    else if (scan.GetToken() == DTOK_RP)
      {
	scan.ReadNext();
	ParseChar (scan, ';');
      }
  
    if (np == 1) // constant spline
      {
	t[1] = t[2] = t[0];
	s[1] = s[2] = s[0];
      }
    if (np == 2) // linear spline
      {
	t[2] = t[1]; t[1] = 0.5*(t[0] + t[2]);
	s[2] = s[1]; s[1] = 0.5*(s[0] + s[2]);
      }
  }




  template <class S>
  void InterpolationSpline<S> :: AddSpline(double t1, double t2, double t3, S s1, S s2, S s3)
  {
    int pos, i, j;
    // find pos to insert interpotation point
    for (pos = 0; pos < ip.Size() && ip[pos][0].GetT() < t1; pos++) ;

    ip.SetSize( ip.Size()+1 );
    for (i = ip.Size()-2; i >= pos; i--)
      for (j = 0; j < 3; j++)
	ip[i+1][j] = ip[i][j];

    ip[pos][0].SetTS (t1, s1);
    ip[pos][1].SetTS (t2, s2);
    ip[pos][2].SetTS (t3, s3);
  }



  template <class S>
  S InterpolationSpline<S> :: Evaluate (double t)
  { 
    if (t < ip[0][0].GetT())
      return (ip[0][0].GetS());
    
    if (t > ip[ip.Size()-1][2].GetT())
      {
	finished = 1;
	return (ip[ip.Size()-1][2].GetS());
      }

    int pos;
    for (pos = 0; pos < ip.Size() && t >= ip[pos][0].GetT(); pos++) ;
    pos--;
  
    if (t >= ip[pos][0].GetT() && t <= ip[pos][2].GetT())
      {
	double t0 = ip[pos][0].GetT();
	double t1 = ip[pos][2].GetT();

	double t01 = (t-t0)/(t1-t0);

	double b1, b2, b3, w;

	b1 = (1-t01)*(1-t01);
	b2 = sqrt(2.0) * t01 * (1-t01);
	b3 = t01 * t01;
	w = b1 + b2 + b3;
 
	return ( (1/w) * (b1 * ip[pos][0].GetS() +
			  b2 * ip[pos][1].GetS() +
			  b3 * ip[pos][2].GetS()) );
      }
    else
      return (ip[pos][2].GetS());
  }



  DemoView :: DemoView (const char * filename)
    : campos( Vec<3>(5,0,0) ),
      campoint ( Vec<3>(0,0,0) ),
      camup ( Vec<3>(0,0,1) )
  {
    double time = 0;

    ifstream istr;
    istr.open(filename);

    DemoScanner scan(istr);

    double t[3];
    Vec<3> s[3];

    scan.ReadNext();

    try
      {
	while (1)
	  {
	    if (scan.GetToken() == DTOK_END) break;
	    
	    if (scan.GetToken() == DTOK_CAMPOS)
	      {
		ParseConstLineOrSpline (scan, &t[0], &s[0]);
		campos.AddSpline (time+t[0], time+t[1], time+t[2], s[0], s[1], s[2]);
	      }
	    
	    else if (scan.GetToken() == DTOK_CAMUP)
	      {
		ParseConstLineOrSpline (scan, &t[0], &s[0]);
		camup.AddSpline (time+t[0], time+t[1], time+t[2], s[0], s[1], s[2]);
	      }
	    
	    else if (scan.GetToken() == DTOK_CAMPOINT)
	      {
		ParseConstLineOrSpline (scan, &t[0], &s[0]);
		campoint.AddSpline (time+t[0], time+t[1], time+t[2], s[0], s[1], s[2]);
	      }
	    
	    else if (scan.GetToken() == DTOK_TIME)
	      {
		scan.ReadNext();

		if (scan.GetToken() != DTOK_EQU && 
		    scan.GetToken() != DTOK_PLUS)
		  scan.Error ("= or += expected");	   
      
		if (scan.GetToken() == DTOK_EQU)
		  {
		    scan.ReadNext();
		    time = ParseNumber (scan)*1000;
		    ParseChar (scan, ';');
		  }
		else if (scan.GetToken() == DTOK_PLUS)
		  {
		    scan.ReadNext();
		    ParseChar (scan, '=');
		    time += ParseNumber (scan)*1000;
		    ParseChar (scan, ';');
		  }
	      }
	    
	    else
	      {
		cout << "read unidentified token " << scan.GetToken() 
		     << " string = " << scan.GetStringValue() << endl;
		scan.ReadNext();
	      }
	  }
      }
    catch (string errstr)
      {
	cout << "caught error " << errstr << endl;
      }


  }


  DemoView :: ~DemoView ()
  {
    ;
  }

  int DemoView :: SetTime (double time)
  {
    /*
      cout << "time = " << time << endl;
  
      cout << "campos  : " << campos.Evaluate (time) << endl;
      cout << "campoint: " << campoint.Evaluate (time) << endl;
      cout << "camup   : " << camup.Evaluate (time) << endl;
    */


    vs -> LookAt ( Point<3>(  campos.Evaluate (time)), 
		   Point<3>(campoint.Evaluate (time)),
		   Point<3>(   camup.Evaluate (time)) );

    if (campos.IsFinished() &&
	campoint.IsFinished() &&
	camup.IsFinished())
      {
	return -1;
      }

    return 0;
  }


}
