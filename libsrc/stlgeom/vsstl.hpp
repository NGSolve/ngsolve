#ifndef FILE_VSSTL
#define FILE_VSSTL

/**************************************************************************/
/* File:   vsstl.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   05. Jan. 2011                                                  */
/**************************************************************************/

namespace netgen
{

 class VisualSceneSTLGeometry : public VisualScene
  {
    Array<int> trilists;
    class STLGeometry * stlgeometry;

  public:
    DLL_HEADER VisualSceneSTLGeometry ();
	DLL_HEADER virtual ~VisualSceneSTLGeometry ();
	DLL_HEADER void SetGeometry (class STLGeometry * astlgeometry) { stlgeometry = astlgeometry; }

	DLL_HEADER virtual void BuildScene (int zoomall = 0);
	DLL_HEADER virtual void DrawScene ();
  };


  class VisualSceneSTLMeshing : public VisualScene
  {
    Array<int> trilists;
    int selecttrig, nodeofseltrig;
    class STLGeometry * stlgeometry;

  public:
    DLL_HEADER VisualSceneSTLMeshing ();
	DLL_HEADER virtual ~VisualSceneSTLMeshing ();

    void SetGeometry (class STLGeometry * astlgeometry) { stlgeometry = astlgeometry; }

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
    virtual void MouseDblClick (int px, int py);

    int seltria;
  };



}



#endif
