#ifndef FILE_BASEGEOM
#define FILE_BASEGEOM

/**************************************************************************/
/* File:   basegeom.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   23. Aug. 09                                                    */
/**************************************************************************/


class NetgenGeometry
{
public:
  virtual ~NetgenGeometry () { ; }

  virtual int GenerateMesh (Mesh*& mesh,
			    int perfstepsstart, int perfstepsend, char* optstring);

  virtual const Refinement & GetRefinement () const;
};







#endif
