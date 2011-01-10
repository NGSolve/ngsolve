#ifndef FILE_VSCSG
#define FILE_VSCSG

/**************************************************************************/
/* File:   vscsg.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   05. Jan. 2011                                                  */
/**************************************************************************/

namespace netgen
{

  class VisualSceneGeometry : public VisualScene
  {
    class CSGeometry * geometry;
    Array<int> trilists;
    int selsurf;
  public:
    VisualSceneGeometry ();
    virtual ~VisualSceneGeometry ();

    void SetGeometry (class CSGeometry * ageometry) { geometry = ageometry; }
    virtual void SelectSurface (int aselsurf);
    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
  };



}



#endif
