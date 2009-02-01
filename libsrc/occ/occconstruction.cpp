
#ifdef OCCGEOMETRY

#include <mystdlib.h>
#include <occgeom.hpp>  
#include "ShapeAnalysis_ShapeTolerance.hxx"
#include "ShapeAnalysis_ShapeContents.hxx"
#include "ShapeAnalysis_CheckSmallFace.hxx"
#include "ShapeAnalysis_DataMapOfShapeListOfReal.hxx"
#include "BRepAlgoAPI_Fuse.hxx"
#include "BRepCheck_Analyzer.hxx"
#include "BRepLib.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "ShapeFix.hxx"
#include "ShapeFix_FixSmallFace.hxx"
#include "Partition_Spliter.hxx"
//#include "VrmlAPI.hxx"
//#include "StlAPI.hxx"


#include <GC_MakeSegment.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
// #include <BRep_Builder.hxx>
#include <TopoDS_Builder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Section.hxx>
#include <BRepOffsetAPI_Sewing.hxx>
//#include <BRepAlgo_Sewing.hxx>
#include <BRepOffsetAPI_MakeOffsetShape.hxx>
#include <ShapeFix_Shape.hxx>
namespace netgen
{

  void OCCConstructGeometry (OCCGeometry & geom)
  {
#ifdef NOTHING
    cout << "OCC construction" << endl;

    BRep_Builder builder;
    BRepPrimAPI_MakeBox mbox(gp_Pnt(-10e5, -15e5, 0), gp_Pnt(20e5, 15e5, 10e5));


    /*
    TopoDS_Shape air = TopoDS_Solid (mbox);
    air = BRepAlgoAPI_Cut (air, geom.somap(1));
    air = BRepAlgoAPI_Cut (air, geom.somap(2));
    air = BRepAlgoAPI_Cut (air, geom.somap(3));
    air = BRepAlgoAPI_Cut (air, geom.somap(4));
    air = BRepAlgoAPI_Cut (air, geom.somap(5));
    air = BRepAlgoAPI_Cut (air, geom.somap(6));
    air = BRepAlgoAPI_Cut (air, geom.somap(7));
    // air = BRepAlgoAPI_Cut (air, geom.somap(8));
    air = BRepAlgoAPI_Cut (air, geom.somap(9));
    // air = BRepAlgoAPI_Cut (air, geom.somap(10));
    */

    /*
    BRepOffsetAPI_MakeOffsetShape dom8plus (geom.somap(8), 1e4, 1e-6);
    BRepOffsetAPI_MakeOffsetShape dom6plus (geom.somap(6), 1e4, 1e-6);
    dom8plus.Build();
    ShapeFix_Shape fixshape(dom8plus.Shape());
    fixshape.Perform();
    
    ShapeFix_Shape fix_dom2(geom.somap(2));
    fix_dom2.Perform();


    BRepAlgoAPI_Cut dom2m8(fix_dom2.Shape(), fixshape.Shape());
    ShapeFix_Shape fix_dom2m8 (dom2m8);
    fix_dom2m8.Perform();

    builder.Add (geom.shape, 
		 BRepAlgoAPI_Cut 
		 (BRepAlgoAPI_Cut (geom.somap(2), dom6plus),
		  dom8plus));
    // builder.Add (geom.shape, fix_dom2m8.Shape());
    //     builder.Add (geom.shape, fixshape.Shape());
    */

    TopoDS_Shape my_fuse;
    int cnt = 0;
    for (TopExp_Explorer exp_solid(geom.shape, TopAbs_SOLID); exp_solid.More(); exp_solid.Next())
      {
	if (cnt == 0)
	  my_fuse = exp_solid.Current();
	else
	  {
	    cout << "fuse, cnt = " << cnt << endl;
	    if (cnt != 7 && cnt != 9)
	      my_fuse = BRepAlgoAPI_Fuse (my_fuse, exp_solid.Current());
	  }
	cnt++;
      }
    builder.Add (geom.shape, my_fuse);

    /*
    ShapeUpgrade_ShellSewing ss;
    ss.ApplySewing(geom.shape,1e5);
    */

    /*
    BRepAlgo_Sewing sewing(1.e5);
    
    int cnt = 0;
    for (TopExp_Explorer exp_solid(geom.shape, TopAbs_SOLID); exp_solid.More(); exp_solid.Next())
      {
	cout << "swe, cnt = " << cnt << endl;
	if (cnt != 7 && cnt != 9)
	  sewing.Add (exp_solid.Current());
	cnt++;
      }

    sewing.Perform();
    builder.Add (geom.shape, sewing.SewedShape());
    */


    /*
    cout << "build air domain" << endl;
    TopoDS_Shape air = BRepAlgoAPI_Cut (TopoDS_Solid (mbox), my_fuse);

    cnt = 0;
    for (TopExp_Explorer exp_solid(geom.shape, TopAbs_SOLID); exp_solid.More(); exp_solid.Next())
      {
	cout << "section, cnt = " << cnt << endl;
	if (cnt == 7)
	  {
	    builder.Add (geom.shape, 
			 BRepAlgoAPI_Section (air, exp_solid.Current()));
	  }
	cnt++;
      }
    */



    //    builder.Add (geom.shape, air);
    for (int i = 1; i <= 10; i++)
      builder.Remove (geom.shape, geom.somap(i));




    geom.BuildFMap();
    geom.BuildVisualizationMesh();
    geom.changed = 1;
#endif

  }
}


#endif
