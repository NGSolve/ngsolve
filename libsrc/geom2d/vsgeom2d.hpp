#ifndef FILE_VSGEOM2D
#define FILE_VSGEOM2D

/**************************************************************************/
/* File:   vsgeom2d.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   05. Jan. 2011                                                  */
/**************************************************************************/

namespace netgen
{

  class VisualSceneGeometry2d : public VisualScene
  {
    const class SplineGeometry2d * geometry2d;
  public:
    VisualSceneGeometry2d ();
    virtual ~VisualSceneGeometry2d ();
    void SetGeometry (const class SplineGeometry2d * ageometry2d) { geometry2d = ageometry2d; }
    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
  };

}



#endif
