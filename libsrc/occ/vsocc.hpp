#ifndef FILE_VSOCC
#define FILE_VSOCC

/**************************************************************************/
/* File:   vsocc.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   05. Jan. 2011                                                  */
/**************************************************************************/

namespace netgen
{

  class VisualSceneOCCGeometry : public VisualScene
  {
    Array<int> trilists;
    Array<int> linelists;
    int selsurf;
    class OCCGeometry * occgeometry;
  public:
    VisualSceneOCCGeometry ();
    virtual ~VisualSceneOCCGeometry ();
    void SetGeometry (class OCCGeometry * ageom) { occgeometry = ageom; }

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
    virtual void MouseDblClick (int px, int py);
  };



}

#endif
