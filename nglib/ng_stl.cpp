/*!
   \file ng_stl.cpp
   \author Philippose Rajan
   \date 14 Feb 2009 (Created)

   This sample utility demonstrates the use of the Netgen 
   nglib library for reading, and meshing an STL geometry. 

   The Program takes as input the name of an STL file 
   saved in the STL ASCII Format, and generates a 3D Volume 
   mesh which is saved into the file "test.vol".

   test.vol can be viewed using the usual Netgen Mesher GUI
*/


#include <iostream>
#include <fstream>

using namespace std;

namespace nglib {
#include <nglib.h>
}

int main (int argc, char ** argv)
{
   using namespace nglib;

   cout << "Netgen (nglib) STL Testing" << endl;

   if (argc < 2)
   {
      cerr << "use: ng_stl STL_filename" << endl;
      return 1;
   }

   // Define pointer to a new Netgen Mesh
   Ng_Mesh *mesh;
   
   // Define pointer to STL Geometry
   Ng_STL_Geometry *stl_geom;

   // Result of Netgen Operations
   Ng_Result ng_res;
   
   // Initialise the Netgen Core library
   Ng_Init();

   // Actually create the mesh structure
   mesh = Ng_NewMesh();

   int np, ne;

   // Read in the STL File
   stl_geom = Ng_STL_LoadGeometry(argv[1]);
   if(!stl_geom)
   {
      cout << "Error reading in STL File: " << argv[1] << endl;
	  return 1;
   }
   cout << "Successfully loaded STL File: " << argv[1] << endl;


   // Set the Meshing Parameters to be used
   Ng_Meshing_Parameters mp;
   mp.maxh = 1.0e+6;
   mp.fineness = 0.4;
   mp.second_order = 0;

   cout << "Initialise the STL Geometry structure...." << endl;
   ng_res = Ng_STL_InitSTLGeometry(stl_geom);
   if(ng_res != NG_OK)
   {
      cout << "Error Initialising the STL Geometry....Aborting!!" << endl;
	   return 1;
   }

   cout << "Start Edge Meshing...." << endl;
   ng_res = Ng_STL_MakeEdges(stl_geom, mesh, &mp);
   if(ng_res != NG_OK)
   {
      cout << "Error in Edge Meshing....Aborting!!" << endl;
	   return 1;
   }

   cout << "Start Surface Meshing...." << endl;
   ng_res = Ng_STL_GenerateSurfaceMesh(stl_geom, mesh, &mp);
   if(ng_res != NG_OK)
   {
      cout << "Error in Surface Meshing....Aborting!!" << endl;
	   return 1;
   }
   
   cout << "Start Volume Meshing...." << endl;
   ng_res = Ng_GenerateVolumeMesh (mesh, &mp);
   if(ng_res != NG_OK)
   {
      cout << "Error in Volume Meshing....Aborting!!" << endl;
	  return 1;
   }
   
   cout << "Meshing successfully completed....!!" << endl;

   // volume mesh output
   np = Ng_GetNP(mesh);
   cout << "Points: " << np << endl;

   ne = Ng_GetNE(mesh);
   cout << "Elements: " << ne << endl;

   cout << "Saving Mesh in VOL Format...." << endl;
   Ng_SaveMesh(mesh,"test.vol");


   // refinement without geomety adaption:
   // Ng_Uniform_Refinement (mesh);

   // refinement with geomety adaption:   
   Ng_STL_Uniform_Refinement (stl_geom, mesh);

   cout << "elements after refinement: " << Ng_GetNE(mesh) << endl;
   cout << "points   after refinement: " << Ng_GetNP(mesh) << endl;

   Ng_SaveMesh(mesh,"test_ref.vol");

   return 0;
}
