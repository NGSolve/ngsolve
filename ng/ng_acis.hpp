#include <acisgeom.hpp>
#ifdef ACIS_R17
extern void unlock_spatial_products_661();
#endif


namespace netgen {


ACISGeometry * acisgeometry = NULL;

static VisualSceneACISGeometry vsacisgeom;


extern int ACISGenerateMesh (ACISGeometry & geometry, Mesh*& mesh,
                             int perfstepsstart, int perfstepsend, char* optstring);










  int Ng_ACISCommand (ClientData /* clientData */,
		      Tcl_Interp * interp,
		      int argc, tcl_const char *argv[]) 
  {
    if (argc >= 2)
      {
	if (strcmp (argv[1], "isACISavailable") == 0)
	  {
	    Tcl_SetResult (interp, (char*)"yes", TCL_STATIC);
	    return TCL_OK;
	  }
      }
    

    if (!acisgeometry)
      {
	Tcl_SetResult (interp, (char*)"This operation needs an ACIS geometry", TCL_STATIC);
	return TCL_ERROR;
      }




    if (argc >= 2 && strcmp (argv[1], "getentities") == 0)
      {
        stringstream str;
        const ENTITY_LIST & entlist = acisgeometry -> entlist;


        ENTITY_LIST cellList;
        api_ct_get_all_cells( entlist, cellList );
        for(int j=0; j<cellList.count(); j++)
          {
            str << "cell" << j << " {Cell " << j << "} \n";

            ATTRIB * ap = cellList[j]->attrib();
            int k = 0; 

            while (ap)
              {
                ATTRIB_GEN_INTEGER * aip = dynamic_cast<ATTRIB_GEN_INTEGER *>(ap);
                if (aip)
                  str << "cell" << j << "/attrib" << k 
                      << " { name = " << aip->name() << " val = " << aip->value() << "}  \n";
                else
                  str << "cell" << j << "/attrib" << k 
                      << " { typename = " << ap->type_name() << "}  \n";
                      
                ap = ap->next();
                k++;
              }
          }


        for (int i = 0; i < entlist.count(); i++)
          {
            str << "entity" << i << " {Entity " << i << "} \n";
            ENTITY_LIST faceList;
            ENTITY_LIST edgeList;
            
            api_get_faces_from_all_entities(entlist[i], faceList);

            for(int j=0; j<faceList.count(); j++)
              {
                FACE * face = (FACE*) faceList[j];

                str << "entity" << i << "/face" << j << " {Face " << j << "} \n";
              
                ATTRIB * ap = faceList[j]->attrib();
                int k = 0; 

                while (ap)
                  {
                    ATTRIB_GEN_INTEGER * aip = dynamic_cast<ATTRIB_GEN_INTEGER *>(ap);
                    if (aip)
                      str << "entity" << i << "/face" << j << "/attrib" << k 
                          << " { name = " << aip->name() << " val = " << aip->value() << "}  \n";
                    else
                      str << "entity" << i << "/face" << j << "/attrib" << k 
                          << " { typename = " << ap->type_name() << "}  \n";
                      
                    ap = ap->next();
                    k++;
                  }

                SPAbox * box = face->bound();
                str << "entity" << i << "/face" << j << "/bbox" 
                    << " { BBox (" 
                    << box->low().x() << ", " << box->low().y() << ", " << box->low().z() << ";" 
                    << box->high().x() << ", " << box->high().y() << ", " << box->high().z() << ")" 
                    << " }\n";
              }

        
            api_get_edges_from_all_entities(entlist[i], edgeList);
            
            for(int j=0; j<edgeList.count(); j++)
              {
                str << "entity" << i << "/edge" << j << " {Edge " << j << "} \n";
              }
          }
        
        Tcl_SetResult (interp, (char*)str.str().c_str(), TCL_VOLATILE);
        return TCL_OK;
      }
    
    if (argc >= 2 && strcmp (argv[1], "selectentity") == 0)
      {
        int ientry = -1, iface = -1;

        string select = argv[2];
        cout << "ACIS selectentity: " << select << endl;

        int pos_ent1 = select.find("entity", 0);
        int pos_ent2 = select.find("/", pos_ent1);
        if (pos_ent1 != -1 && pos_ent2 == -1) pos_ent2 = select.length();
        if (pos_ent1 != -1) ientry = atoi (select.substr(pos_ent1+6, pos_ent2).c_str());

        int pos_face1 = select.find("face", 0);
        int pos_face2 = select.find("/", pos_face1);
        if (pos_face1 != -1 && pos_face2 == -1) pos_face2 = select.length();
        if (pos_face1 != -1) iface = atoi (select.substr(pos_face1+4, pos_face2).c_str());

        cout << "entry = " << ientry << ", face = " << iface << endl;

        const ENTITY_LIST & entlist = acisgeometry -> entlist;

        if (ientry != -1 && iface == -1)
          vsacisgeom.SelectEntity (entlist[ientry]);

        if (iface != -1)
          {
            ENTITY_LIST faceList;
            api_get_faces_from_all_entities(entlist[ientry], faceList);
            vsacisgeom.SelectEntity (faceList[iface]);
          }
      }


    if (argc >= 2 && strcmp (argv[1], "createct") == 0)
      {
	acisgeometry -> CreateCT();
      }

    if (argc >= 2 && strcmp (argv[1], "combineall") == 0)
      {
	cout << "combineall " << endl;
	acisgeometry -> Combine();
      }



    if (argc >= 4)
      {
	if (strcmp (argv[1], "subtract") == 0)
	  {
	    cout << "subtract " << argv[2] << " minus " << argv[3] << endl;
	    acisgeometry -> Subtract (atoi (argv[2])-1, atoi (argv[3])-1);
	  }

      }

    return TCL_OK;
  }


}
