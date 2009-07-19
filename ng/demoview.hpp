#ifndef FILE_DEMOVIEW
#define FILE_DEMOVIEW

/*********************************************************************/
/* File:   demoview.hpp                                              */
/* Author: Robert, Joachim                                           */
/* Date:   6. Mar. 2003                                              */
/*********************************************************************/

using namespace netgen;


enum DEMOVIEW_TOKEN_TYPE
  { 
    DTOK_MINUS = '-', DTOK_LP = '(', DTOK_RP = ')', DTOK_LSP = '[', DTOK_RSP = ']',
    DTOK_EQU = '=', DTOK_COMMA = ',', DTOK_SEMICOLON = ';', DTOK_COLON = ':', DTOK_PLUS = '+',
    DTOK_NUM = 100, DTOK_STRING, DTOK_TIME, DTOK_CAMPOS, DTOK_CAMPOINT, DTOK_CAMUP,
    DTOK_END 
  };

struct demoview_kwstruct
{
  DEMOVIEW_TOKEN_TYPE kw; 
  const char * name;
};




class DemoScanner
{
  DEMOVIEW_TOKEN_TYPE token;
  double num_value;
  string string_value;
    
  int linenum;
  ifstream * scanin;

public:

  DemoScanner (ifstream & ascanin);

  DEMOVIEW_TOKEN_TYPE GetToken() const
  { return token; }

  double GetNumValue() const
  { return num_value; }

  const string & GetStringValue() const
  { return string_value; }

  void ReadNext();
  void Error (const string & err);
};


void ParseChar (DemoScanner & scan, char ch);

double ParseNumber(DemoScanner & scan);

Vec<3> ParseVector (DemoScanner & scan);



template <class S>
class InterpolationPoint
{
  double t;
  S s;

public:
  InterpolationPoint()
  {};

  ~InterpolationPoint()
  {};

  double GetT() const
  { return t; };

  S GetS() const
  { return s; };

  void SetTS(double at, S as)
  { t = at; s = as; };

  InterpolationPoint & operator= (const InterpolationPoint<S> & ip2)
  {
    SetTS (ip2.t, ip2.s);
    return (*this);
  };

};



template <class S>
class InterpolationSpline
{
protected:
  // Array < InterpolationPoint<S>[3] > ip;

  class intpts
  {
  public:
    InterpolationPoint<S> pts[3];
    InterpolationPoint<S> & operator[](int i) { return pts[i]; }
  };
  Array < intpts > ip;
  
  int finished;

public:
  InterpolationSpline() : finished(0)
  {};

  InterpolationSpline( S s1 ) : finished(0)
  {
    AddSpline (-1e99, -1e99, -1e99, s1, s1, s1);
    // InterpolationSpline();
  }

  ~InterpolationSpline()
  {};

  void AddSpline(double t1, double t2, double t3, S s1, S s2, S s3);

  S Evaluate (double t);

  int IsFinished() const
  {
    return finished;
  }
};





class DemoView
{
  InterpolationSpline< Vec<3> > campos;
  InterpolationSpline< Vec<3> > campoint;
  InterpolationSpline< Vec<3> > camup;

public:
  DemoView (const char * filename);
  ~DemoView ();
  
  int SetTime (double time);
};



#endif
