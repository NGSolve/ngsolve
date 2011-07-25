#ifdef OCCGEOMETRY

#include <mystdlib.h>
#include <occgeom.hpp>
#include <meshing.hpp>


namespace netgen
{

#include "occmeshsurf.hpp"

#define TCL_OK 0
#define TCL_ERROR 1

#define DIVIDEEDGESECTIONS 1000
#define IGNORECURVELENGTH 1e-4
#define VSMALL 1e-10


   bool merge_solids = 1;


  // can you please explain what you intend to compute here (JS) !!!
   double Line :: Dist (Line l)
   {
      Vec<3> n = p1-p0;
      Vec<3> q = l.p1-l.p0;
      double nq = n*q;

      Point<3> p = p0 + 0.5*n;
      double lambda = (p-l.p0)*n / (nq + VSMALL);

      if (lambda >= 0 && lambda <= 1)
      {
         double d = (p-l.p0-lambda*q).Length();
         //        if (d < 1e-3) d = 1e99;
         return d;
      }
      else
         return 1e99;
   }



   double Line :: Length ()
   {
      return (p1-p0).Length();
   }



   inline Point<3> occ2ng (const gp_Pnt & p)
   {
      return  Point<3> (p.X(), p.Y(), p.Z());
   }



   double ComputeH (double kappa)
   {
      double hret;
      kappa *= mparam.curvaturesafety;

      if (mparam.maxh * kappa < 1)
         hret = mparam.maxh;
      else
         hret = 1 / (kappa + VSMALL);

      if (mparam.maxh < hret)
         hret = mparam.maxh;

      return (hret);
   }




   void RestrictHTriangle (gp_Pnt2d & par0, gp_Pnt2d & par1, gp_Pnt2d & par2,
                           BRepLProp_SLProps * prop, Mesh & mesh, int depth, double h = 0)
   {
      int ls = -1;

      gp_Pnt pnt0,pnt1,pnt2;

      prop->SetParameters (par0.X(), par0.Y());
      pnt0 = prop->Value();

      prop->SetParameters (par1.X(), par1.Y());
      pnt1 = prop->Value();

      prop->SetParameters (par2.X(), par2.Y());
      pnt2 = prop->Value();

      double aux;
      double maxside = pnt0.Distance(pnt1);
      ls = 2;
      aux = pnt1.Distance(pnt2);
      if(aux > maxside)
      {
         maxside = aux;
         ls = 0;
      }
      aux = pnt2.Distance(pnt0);
      if(aux > maxside)
      {
         maxside = aux;
         ls = 1;
      }



      gp_Pnt2d parmid;

      parmid.SetX( (par0.X()+par1.X()+par2.X()) / 3 );
      parmid.SetY( (par0.Y()+par1.Y()+par2.Y()) / 3 );

      if (depth%3 == 0)
      {
         double curvature = 0;

         prop->SetParameters (parmid.X(), parmid.Y());
         if (!prop->IsCurvatureDefined())
         {
            (*testout) << "curvature not defined!" << endl;
            return;
         }
         curvature = max(fabs(prop->MinCurvature()),
            fabs(prop->MaxCurvature()));

         prop->SetParameters (par0.X(), par0.Y());
         if (!prop->IsCurvatureDefined())
         {
            (*testout) << "curvature not defined!" << endl;
            return;
         }
         curvature = max(curvature,max(fabs(prop->MinCurvature()),
            fabs(prop->MaxCurvature())));

         prop->SetParameters (par1.X(), par1.Y());
         if (!prop->IsCurvatureDefined())
         {
            (*testout) << "curvature not defined!" << endl;
            return;
         }
         curvature = max(curvature,max(fabs(prop->MinCurvature()),
            fabs(prop->MaxCurvature())));

         prop->SetParameters (par2.X(), par2.Y());
         if (!prop->IsCurvatureDefined())
         {
            (*testout) << "curvature not defined!" << endl;
            return;
         }
         curvature = max(curvature,max(fabs(prop->MinCurvature()),
            fabs(prop->MaxCurvature())));

         //(*testout) << "curvature " << curvature << endl;

         if (curvature < 1e-3)
         {
            //(*testout) << "curvature too small (" << curvature << ")!" << endl;
            return;
            // return war bis 10.2.05 auskommentiert
         }



         h = ComputeH (curvature+1e-10);

         if(h < 1e-4*maxside)
            return;


         if (h > 30) return;
      }

      if (h < maxside && depth < 10)
      {
         //cout << "\r h " << h << flush;
         gp_Pnt2d pm;

         //cout << "h " << h << " maxside " << maxside << " depth " << depth << endl;
         //cout << "par0 " << par0.X() << " " << par0.Y()
         //<< " par1 " << par1.X() << " " << par1.Y()
         //   << " par2 " << par2.X() << " " << par2.Y()<< endl;

         if(ls == 0)
         {
            pm.SetX(0.5*(par1.X()+par2.X())); pm.SetY(0.5*(par1.Y()+par2.Y()));
            RestrictHTriangle(pm, par2, par0, prop, mesh, depth+1, h);
            RestrictHTriangle(pm, par0, par1, prop, mesh, depth+1, h);
         }
         else if(ls == 1)
         {
            pm.SetX(0.5*(par0.X()+par2.X())); pm.SetY(0.5*(par0.Y()+par2.Y()));
            RestrictHTriangle(pm, par1, par2, prop, mesh, depth+1, h);
            RestrictHTriangle(pm, par0, par1, prop, mesh, depth+1, h);
         }
         else if(ls == 2)
         {
            pm.SetX(0.5*(par0.X()+par1.X())); pm.SetY(0.5*(par0.Y()+par1.Y()));
            RestrictHTriangle(pm, par1, par2, prop, mesh, depth+1, h);
            RestrictHTriangle(pm, par2, par0, prop, mesh, depth+1, h);
         }

      }
      else
      {
         gp_Pnt pnt;
         Point3d p3d;

         prop->SetParameters (parmid.X(), parmid.Y());
         pnt = prop->Value();
         p3d = Point3d(pnt.X(), pnt.Y(), pnt.Z());
         mesh.RestrictLocalH (p3d, h);

         p3d = Point3d(pnt0.X(), pnt0.Y(), pnt0.Z());
         mesh.RestrictLocalH (p3d, h);

         p3d = Point3d(pnt1.X(), pnt1.Y(), pnt1.Z());
         mesh.RestrictLocalH (p3d, h);

         p3d = Point3d(pnt2.X(), pnt2.Y(), pnt2.Z());
         mesh.RestrictLocalH (p3d, h);

         //(*testout) << "p = " << p3d << ", h = " << h << ", maxside = " << maxside << endl;

      }
   }



   void DivideEdge (TopoDS_Edge & edge, Array<MeshPoint> & ps,
                    Array<double> & params, Mesh & mesh)
   {
      double s0, s1;
      double maxh = mparam.maxh;
      int nsubedges = 1;
      gp_Pnt pnt, oldpnt;
      double svalue[DIVIDEEDGESECTIONS];

      GProp_GProps system;
      BRepGProp::LinearProperties(edge, system);
      double L = system.Mass();

      Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);

      double hvalue[DIVIDEEDGESECTIONS+1];
      hvalue[0] = 0;
      pnt = c->Value(s0);

      double olddist = 0;
      double dist = 0;

      int tmpVal = (int)(DIVIDEEDGESECTIONS);

      for (int i = 1; i <= tmpVal; i++)
      {
         oldpnt = pnt;
         pnt = c->Value(s0+(i/double(DIVIDEEDGESECTIONS))*(s1-s0));
         hvalue[i] = hvalue[i-1] +
            1.0/mesh.GetH(Point3d(pnt.X(), pnt.Y(), pnt.Z()))*
            pnt.Distance(oldpnt);

         //(*testout) << "mesh.GetH(Point3d(pnt.X(), pnt.Y(), pnt.Z())) " << mesh.GetH(Point3d(pnt.X(), pnt.Y(), pnt.Z()))
         //	   <<  " pnt.Distance(oldpnt) " << pnt.Distance(oldpnt) << endl;

         olddist = dist;
         dist = pnt.Distance(oldpnt);
      }

      //  nsubedges = int(ceil(hvalue[DIVIDEEDGESECTIONS]));
      nsubedges = max (1, int(floor(hvalue[DIVIDEEDGESECTIONS]+0.5)));

      ps.SetSize(nsubedges-1);
      params.SetSize(nsubedges+1);

      int i = 1;
      int i1 = 0;
      do
      {
         if (hvalue[i1]/hvalue[DIVIDEEDGESECTIONS]*nsubedges >= i)
         {
            params[i] = s0+(i1/double(DIVIDEEDGESECTIONS))*(s1-s0);
            pnt = c->Value(params[i]);
            ps[i-1] = MeshPoint (Point3d(pnt.X(), pnt.Y(), pnt.Z()));
            i++;
         }
         i1++;
         if (i1 > DIVIDEEDGESECTIONS)
         {
            nsubedges = i;
            ps.SetSize(nsubedges-1);
            params.SetSize(nsubedges+1);
            cout << "divide edge: local h too small" << endl;
         }
      } while (i < nsubedges);

      params[0] = s0;
      params[nsubedges] = s1;

      if (params[nsubedges] <= params[nsubedges-1])
      {
         cout << "CORRECTED" << endl;
         ps.SetSize (nsubedges-2);
         params.SetSize (nsubedges);
         params[nsubedges] = s1;
      }
   }




   void OCCFindEdges (OCCGeometry & geom, Mesh & mesh)
   {
      const char * savetask = multithread.task;
      multithread.task = "Edge meshing";

      (*testout) << "edge meshing" << endl;

      int nvertices = geom.vmap.Extent();
      int nedges = geom.emap.Extent();

      (*testout) << "nvertices = " << nvertices << endl;
      (*testout) << "nedges = " << nedges << endl;

      double eps = 1e-6 * geom.GetBoundingBox().Diam();

      for (int i = 1; i <= nvertices; i++)
      {
         gp_Pnt pnt = BRep_Tool::Pnt (TopoDS::Vertex(geom.vmap(i)));
         MeshPoint mp( Point<3>(pnt.X(), pnt.Y(), pnt.Z()) );

         bool exists = 0;
         if (merge_solids)
            for (PointIndex pi = 1; pi <= mesh.GetNP(); pi++)
               if ( Dist2 (mesh[pi], Point<3>(mp)) < eps*eps)
               {
                  exists = 1;
                  break;
               }

               if (!exists)
                  mesh.AddPoint (mp);
      }

      (*testout) << "different vertices = " << mesh.GetNP() << endl;


      int first_ep = mesh.GetNP()+1;

      Array<int> face2solid[2];
      for (int i = 0; i<2; i++)
      {
         face2solid[i].SetSize (geom.fmap.Extent());
         face2solid[i] = 0;
      }

      int solidnr = 0;
      for (TopExp_Explorer exp0(geom.shape, TopAbs_SOLID); exp0.More(); exp0.Next())
      {
         solidnr++;
         for (TopExp_Explorer exp1(exp0.Current(), TopAbs_FACE); exp1.More(); exp1.Next())
         {
            TopoDS_Face face = TopoDS::Face(exp1.Current());
            int facenr = geom.fmap.FindIndex(face);

            if (face2solid[0][facenr-1] == 0)
               face2solid[0][facenr-1] = solidnr;
            else
               face2solid[1][facenr-1] = solidnr;
         }
      }


      int total = 0;
      for (int i3 = 1; i3 <= geom.fmap.Extent(); i3++)
         for (TopExp_Explorer exp2(geom.fmap(i3), TopAbs_WIRE); exp2.More(); exp2.Next())
            for (TopExp_Explorer exp3(exp2.Current(), TopAbs_EDGE); exp3.More(); exp3.Next())
               total++;


      int facenr = 0;
      int edgenr = 0;


      (*testout) << "faces = " << geom.fmap.Extent() << endl;
      int curr = 0;

      for (int i3 = 1; i3 <= geom.fmap.Extent(); i3++)
      {
         TopoDS_Face face = TopoDS::Face(geom.fmap(i3));
         facenr = geom.fmap.FindIndex (face);       // sollte doch immer == i3 sein ??? JS

         int solidnr0 = face2solid[0][i3-1];
         int solidnr1 = face2solid[1][i3-1];

         /* auskommentiert am 3.3.05 von robert
         for (exp2.Init (geom.somap(solidnr0), TopAbs_FACE); exp2.More(); exp2.Next())
         {
         TopoDS_Face face2 = TopoDS::Face(exp2.Current());
         if (geom.fmap.FindIndex(face2) == facenr)
         {
         //		      if (face.Orientation() != face2.Orientation()) swap (solidnr0, solidnr1);
         }
         }
         */

         mesh.AddFaceDescriptor (FaceDescriptor(facenr, solidnr0, solidnr1, 0));

         // Philippose - 06/07/2009
         // Add the face colour to the mesh data
         Quantity_Color face_colour;

         if(!(geom.face_colours.IsNull())
            && (geom.face_colours->GetColor(face,XCAFDoc_ColorSurf,face_colour)))
         {
            mesh.GetFaceDescriptor(facenr).SetSurfColour(Vec3d(face_colour.Red(),face_colour.Green(),face_colour.Blue()));
         }
         else
         {
            mesh.GetFaceDescriptor(facenr).SetSurfColour(Vec3d(0.0,1.0,0.0));
         }
         // ACHTUNG! STIMMT NICHT ALLGEMEIN (RG)


         Handle(Geom_Surface) occface = BRep_Tool::Surface(face);

         for (TopExp_Explorer exp2 (face, TopAbs_WIRE); exp2.More(); exp2.Next())
         {
            TopoDS_Shape wire = exp2.Current();

            for (TopExp_Explorer exp3 (wire, TopAbs_EDGE); exp3.More(); exp3.Next())
            {
               curr++;
               (*testout) << "edge nr " << curr << endl;

               multithread.percent = 100 * curr / double (total);
               if (multithread.terminate) return;

               TopoDS_Edge edge = TopoDS::Edge (exp3.Current());
               if (BRep_Tool::Degenerated(edge))
               {
                  //(*testout) << "ignoring degenerated edge" << endl;
                  continue;
               }

               if (geom.vmap.FindIndex(TopExp::FirstVertex (edge)) ==
                  geom.vmap.FindIndex(TopExp::LastVertex (edge)))
               {
                  GProp_GProps system;
                  BRepGProp::LinearProperties(edge, system);

                  if (system.Mass() < eps)
                  {
                     cout << "ignoring edge " << geom.emap.FindIndex (edge)
                        << ". closed edge with length < " << eps << endl;
                     continue;
                  }
               }


               Handle(Geom2d_Curve) cof;
               double s0, s1;
               cof = BRep_Tool::CurveOnSurface (edge, face, s0, s1);

               int geomedgenr = geom.emap.FindIndex(edge);

               Array <MeshPoint> mp;
               Array <double> params;

               DivideEdge (edge, mp, params, mesh);
 
               Array <int> pnums;
               pnums.SetSize (mp.Size()+2);

               if (!merge_solids)
               {
                  pnums[0] = geom.vmap.FindIndex (TopExp::FirstVertex (edge));
                  pnums[pnums.Size()-1] = geom.vmap.FindIndex (TopExp::LastVertex (edge));
               }
               else
               {
                  Point<3> fp = occ2ng (BRep_Tool::Pnt (TopExp::FirstVertex (edge)));
                  Point<3> lp = occ2ng (BRep_Tool::Pnt (TopExp::LastVertex (edge)));

                  pnums[0] = -1;
                  pnums.Last() = -1;
                  for (PointIndex pi = 1; pi < first_ep; pi++)
                  {
                     if (Dist2 (mesh[pi], fp) < eps*eps) pnums[0] = pi;
                     if (Dist2 (mesh[pi], lp) < eps*eps) pnums.Last() = pi;
                  }
               }


               for (int i = 1; i <= mp.Size(); i++)
               {
                  bool exists = 0;
                  int j;
                  for (j = first_ep; j <= mesh.GetNP(); j++)
                     if ((mesh.Point(j)-Point<3>(mp[i-1])).Length() < eps)
                     {
                        exists = 1;
                        break;
                     }

                     if (exists)
                        pnums[i] = j;
                     else
                     {
                        mesh.AddPoint (mp[i-1]);
                        (*testout) << "add meshpoint " << mp[i-1] << endl;
                        pnums[i] = mesh.GetNP();
                     }
               }
               (*testout) << "NP = " << mesh.GetNP() << endl;

               //(*testout) << pnums[pnums.Size()-1] << endl;

               for (int i = 1; i <= mp.Size()+1; i++)
               {
                  edgenr++;
                  Segment seg;

                  seg[0] = pnums[i-1];
                  seg[1] = pnums[i];
                  seg.edgenr = edgenr;
                  seg.si = facenr;
                  seg.epgeominfo[0].dist = params[i-1];
                  seg.epgeominfo[1].dist = params[i];
                  seg.epgeominfo[0].edgenr = geomedgenr;
                  seg.epgeominfo[1].edgenr = geomedgenr;

                  gp_Pnt2d p2d;
                  p2d = cof->Value(params[i-1]);
                  //			if (i == 1) p2d = cof->Value(s0);
                  seg.epgeominfo[0].u = p2d.X();
                  seg.epgeominfo[0].v = p2d.Y();
                  p2d = cof->Value(params[i]);
                  //			if (i == mp.Size()+1) p2d = cof -> Value(s1);
                  seg.epgeominfo[1].u = p2d.X();
                  seg.epgeominfo[1].v = p2d.Y();

                  /*
                  if (occface->IsUPeriodic())
                  {
                  cout << "U Periodic" << endl;
                  if (fabs(seg.epgeominfo[1].u-seg.epgeominfo[0].u) >
                  fabs(seg.epgeominfo[1].u-
                  (seg.epgeominfo[0].u-occface->UPeriod())))
                  seg.epgeominfo[0].u = p2d.X()+occface->UPeriod();

                  if (fabs(seg.epgeominfo[1].u-seg.epgeominfo[0].u) >
                  fabs(seg.epgeominfo[1].u-
                  (seg.epgeominfo[0].u+occface->UPeriod())))
                  seg.epgeominfo[0].u = p2d.X()-occface->UPeriod();
                  }

                  if (occface->IsVPeriodic())
                  {
                  cout << "V Periodic" << endl;
                  if (fabs(seg.epgeominfo[1].v-seg.epgeominfo[0].v) >
                  fabs(seg.epgeominfo[1].v-
                  (seg.epgeominfo[0].v-occface->VPeriod())))
                  seg.epgeominfo[0].v = p2d.Y()+occface->VPeriod();

                  if (fabs(seg.epgeominfo[1].v-seg.epgeominfo[0].v) >
                  fabs(seg.epgeominfo[1].v-
                  (seg.epgeominfo[0].v+occface->VPeriod())))
                  seg.epgeominfo[0].v = p2d.Y()-occface->VPeriod();
                  }
                  */

                  if (edge.Orientation() == TopAbs_REVERSED)
                  {
                     swap (seg[0], seg[1]);
                     swap (seg.epgeominfo[0].dist, seg.epgeominfo[1].dist);
                     swap (seg.epgeominfo[0].u, seg.epgeominfo[1].u);
                     swap (seg.epgeominfo[0].v, seg.epgeominfo[1].v);
                  }

                  mesh.AddSegment (seg);

                  //edgesegments[geomedgenr-1]->Append(mesh.GetNSeg());

               }
            }
         }
      }

      //	for(i=1; i<=mesh.GetNSeg(); i++)
      //		(*testout) << "edge " << mesh.LineSegment(i).edgenr << " face " << mesh.LineSegment(i).si
      //				<< " p1 " << mesh.LineSegment(i)[0] << " p2 " << mesh.LineSegment(i)[1] << endl;
      //	exit(10);

      mesh.CalcSurfacesOfNode();
      multithread.task = savetask;
   }




   void OCCMeshSurface (OCCGeometry & geom, Mesh & mesh, int perfstepsend)
   {
      int i, j, k;
      int changed;

      const char * savetask = multithread.task;
      multithread.task = "Surface meshing";

      geom.facemeshstatus = 0;

      int noldp = mesh.GetNP();

      double starttime = GetTime();

      Array<int> glob2loc(noldp);

      //int projecttype = PARAMETERSPACE;

      int projecttype = PARAMETERSPACE;

      int notrys = 1;

      int surfmesherror = 0;

      for (k = 1; k <= mesh.GetNFD(); k++)
      {
         if(1==0 && !geom.fvispar[k-1].IsDrawable())
         {
            (*testout) << "ignoring face " << k << endl;
            cout << "ignoring face " << k << endl;
            continue;
         }

         (*testout) << "mesh face " << k << endl;
         multithread.percent = 100 * k / (mesh.GetNFD() + VSMALL);
         geom.facemeshstatus[k-1] = -1;


         /*
         if (k != 42)
         {
         cout << "skipped" << endl;
         continue;
         }
         */


         FaceDescriptor & fd = mesh.GetFaceDescriptor(k);

         int oldnf = mesh.GetNSE();

         Box<3> bb = geom.GetBoundingBox();

         //      int projecttype = PLANESPACE;

         Meshing2OCCSurfaces meshing(TopoDS::Face(geom.fmap(k)), bb, projecttype);

         if (meshing.GetProjectionType() == PLANESPACE)
            PrintMessage (2, "Face ", k, " / ", mesh.GetNFD(), " (plane space projection)");
         else
            PrintMessage (2, "Face ", k, " / ", mesh.GetNFD(), " (parameter space projection)");

         if (surfmesherror)
            cout << "Surface meshing error occured before (in " << surfmesherror << " faces)" << endl;

         //      Meshing2OCCSurfaces meshing(f2, bb);
         meshing.SetStartTime (starttime);

         //(*testout) << "Face " << k << endl << endl;


         if (meshing.GetProjectionType() == PLANESPACE)
         {
            int cntp = 0;
            glob2loc = 0;
            for (i = 1; i <= mesh.GetNSeg(); i++)
            {
               Segment & seg = mesh.LineSegment(i);
               if (seg.si == k)
               {
                  for (j = 1; j <= 2; j++)
                  {
                     int pi = (j == 1) ? seg[0] : seg[1];
                     if (!glob2loc.Get(pi))
                     {
                        meshing.AddPoint (mesh.Point(pi), pi);
                        cntp++;
                        glob2loc.Elem(pi) = cntp;
                     }
                  }
               }
            }

            for (i = 1; i <= mesh.GetNSeg(); i++)
            {
               Segment & seg = mesh.LineSegment(i);
               if (seg.si == k)
               {
                  PointGeomInfo gi0, gi1;
                  gi0.trignum = gi1.trignum = k;
                  gi0.u = seg.epgeominfo[0].u;
                  gi0.v = seg.epgeominfo[0].v;
                  gi1.u = seg.epgeominfo[1].u;
                  gi1.v = seg.epgeominfo[1].v;

                  meshing.AddBoundaryElement (glob2loc.Get(seg[0]), glob2loc.Get(seg[1]), gi0, gi1);
                  //(*testout) << gi0.u << " " << gi0.v << endl;
                  //(*testout) << gi1.u << " " << gi1.v << endl;
               }
            }
         }
         else
         {
            int cntp = 0;

            for (i = 1; i <= mesh.GetNSeg(); i++)
               if (mesh.LineSegment(i).si == k)
                  cntp+=2;


            Array< PointGeomInfo > gis;

            gis.SetAllocSize (cntp);
            gis.SetSize (0);

            for (i = 1; i <= mesh.GetNSeg(); i++)
            {
               Segment & seg = mesh.LineSegment(i);
               if (seg.si == k)
               {
                  PointGeomInfo gi0, gi1;
                  gi0.trignum = gi1.trignum = k;
                  gi0.u = seg.epgeominfo[0].u;
                  gi0.v = seg.epgeominfo[0].v;
                  gi1.u = seg.epgeominfo[1].u;
                  gi1.v = seg.epgeominfo[1].v;

                  int locpnum[2] = {0, 0};

                  for (j = 0; j < 2; j++)
                  {
                     PointGeomInfo gi = (j == 0) ? gi0 : gi1;

                     int l;
                     for (l = 0; l < gis.Size() && locpnum[j] == 0; l++)
                     {
                        double dist = sqr (gis[l].u-gi.u)+sqr(gis[l].v-gi.v);

                        if (dist < 1e-10)
                           locpnum[j] = l+1;
                     }

                     if (locpnum[j] == 0)
                     {
                        int pi = (j == 0) ? seg[0] : seg[1];
                        meshing.AddPoint (mesh.Point(pi), pi);

                        gis.SetSize (gis.Size()+1);
                        gis[l] = gi;
                        locpnum[j] = l+1;
                     }
                  }

                  meshing.AddBoundaryElement (locpnum[0], locpnum[1], gi0, gi1);
                  //(*testout) << gi0.u << " " << gi0.v << endl;
                  //(*testout) << gi1.u << " " << gi1.v << endl;

               }
            }
         }





         // Philippose - 15/01/2009
         double maxh = geom.face_maxh[k-1];
         //double maxh = mparam.maxh;
         mparam.checkoverlap = 0;
         //      int noldpoints = mesh->GetNP();
         int noldsurfel = mesh.GetNSE();

         GProp_GProps sprops;
         BRepGProp::SurfaceProperties(TopoDS::Face(geom.fmap(k)),sprops);
         meshing.SetMaxArea(2.*sprops.Mass());

         MESHING2_RESULT res;

         try {
	   res = meshing.GenerateMesh (mesh, mparam, maxh, k);
         }

         catch (SingularMatrixException)
         {
            (*myerr) << "Singular Matrix" << endl;
            res = MESHING2_GIVEUP;
         }

         catch (UVBoundsException)
         {
            (*myerr) << "UV bounds exceeded" << endl;
            res = MESHING2_GIVEUP;
         }

         projecttype = PARAMETERSPACE;

         if (res != MESHING2_OK)
         {
            if (notrys == 1)
            {
               for (int i = noldsurfel+1; i <= mesh.GetNSE(); i++)
                  mesh.DeleteSurfaceElement (i);

               mesh.Compress();

               cout << "retry Surface " << k << endl;

               k--;
               projecttype*=-1;
               notrys++;
               continue;
            }
            else
            {
               geom.facemeshstatus[k-1] = -1;
               PrintError ("Problem in Surface mesh generation");
               surfmesherror++;
               //	      throw NgException ("Problem in Surface mesh generation");
            }
         }
         else
         {
            geom.facemeshstatus[k-1] = 1;
         }

         notrys = 1;

         for (i = oldnf+1; i <= mesh.GetNSE(); i++)
            mesh.SurfaceElement(i).SetIndex (k);

      }

//      ofstream problemfile("occmesh.rep");

//      problemfile << "SURFACEMESHING" << endl << endl;

      if (surfmesherror)
      {
         cout << "WARNING! NOT ALL FACES HAVE BEEN MESHED" << endl;
         cout << "SURFACE MESHING ERROR OCCURED IN " << surfmesherror << " FACES:" << endl;
         for (int i = 1; i <= geom.fmap.Extent(); i++)
            if (geom.facemeshstatus[i-1] == -1)
            {
               cout << "Face " << i << endl;
//               problemfile << "problem with face " << i << endl;
//               problemfile << "vertices: " << endl;
               TopExp_Explorer exp0,exp1,exp2;
               for ( exp0.Init(TopoDS::Face (geom.fmap(i)), TopAbs_WIRE); exp0.More(); exp0.Next() )
               {
                  TopoDS_Wire wire = TopoDS::Wire(exp0.Current());
                  for ( exp1.Init(wire,TopAbs_EDGE); exp1.More(); exp1.Next() )
                  {
                     TopoDS_Edge edge = TopoDS::Edge(exp1.Current());
                     for ( exp2.Init(edge,TopAbs_VERTEX); exp2.More(); exp2.Next() )
                     {
                        TopoDS_Vertex vertex = TopoDS::Vertex(exp2.Current());
                        gp_Pnt point = BRep_Tool::Pnt(vertex);
//                        problemfile << point.X() << " " << point.Y() << " " << point.Z() << endl;
                     }
                  }
               }
//               problemfile << endl;

            }
            cout << endl << endl;
            cout << "for more information open IGES/STEP Topology Explorer" << endl;
//            problemfile.close();
            throw NgException ("Problem in Surface mesh generation");
      }
      else
      {
//         problemfile << "OK" << endl << endl;
//         problemfile.close();
      }




      if (multithread.terminate || perfstepsend < MESHCONST_OPTSURFACE)
         return;

      multithread.task = "Optimizing surface";

      static int timer_opt2d = NgProfiler::CreateTimer ("Optimization 2D");
      NgProfiler::StartTimer (timer_opt2d);

      for (k = 1; k <= mesh.GetNFD(); k++)
      {
         //      if (k != 42) continue;
         //      if (k != 36) continue;

         //      (*testout) << "optimize face " << k << endl;
         multithread.percent = 100 * k / (mesh.GetNFD() + VSMALL);

         FaceDescriptor & fd = mesh.GetFaceDescriptor(k);

         PrintMessage (1, "Optimize Surface ", k);
         for (i = 1; i <= mparam.optsteps2d; i++)
         {
            //	  (*testout) << "optstep " << i << endl;
            if (multithread.terminate) return;

            {
               MeshOptimize2dOCCSurfaces meshopt(geom);
               meshopt.SetFaceIndex (k);
               meshopt.SetImproveEdges (0);
               meshopt.SetMetricWeight (mparam.elsizeweight);
               //meshopt.SetMetricWeight (0.2);
               meshopt.SetWriteStatus (0);

               //	    (*testout) << "EdgeSwapping (mesh, (i > mparam.optsteps2d/2))" << endl;
               meshopt.EdgeSwapping (mesh, (i > mparam.optsteps2d/2));
            }

            if (multithread.terminate) return;
            {
               MeshOptimize2dOCCSurfaces meshopt(geom);
               meshopt.SetFaceIndex (k);
               meshopt.SetImproveEdges (0);
               //meshopt.SetMetricWeight (0.2);
               meshopt.SetMetricWeight (mparam.elsizeweight);
               meshopt.SetWriteStatus (0);

               //	    (*testout) << "ImproveMesh (mesh)" << endl;
               meshopt.ImproveMesh (mesh, mparam);
            }

            {
               MeshOptimize2dOCCSurfaces meshopt(geom);
               meshopt.SetFaceIndex (k);
               meshopt.SetImproveEdges (0);
               //meshopt.SetMetricWeight (0.2);
               meshopt.SetMetricWeight (mparam.elsizeweight);
               meshopt.SetWriteStatus (0);

               //	    (*testout) << "CombineImprove (mesh)" << endl;
               meshopt.CombineImprove (mesh);
            }

            if (multithread.terminate) return;
            {
               MeshOptimize2dOCCSurfaces meshopt(geom);
               meshopt.SetFaceIndex (k);
               meshopt.SetImproveEdges (0);
               //meshopt.SetMetricWeight (0.2);
               meshopt.SetMetricWeight (mparam.elsizeweight);
               meshopt.SetWriteStatus (0);

               //	    (*testout) << "ImproveMesh (mesh)" << endl;
               meshopt.ImproveMesh (mesh, mparam);
            }
         }

      }


      mesh.CalcSurfacesOfNode();
      mesh.Compress();

      NgProfiler::StopTimer (timer_opt2d);

      multithread.task = savetask;
   }



   void OCCSetLocalMeshSize(OCCGeometry & geom, Mesh & mesh)
   {
      mesh.SetGlobalH (mparam.maxh);
      mesh.SetMinimalH (mparam.minh);

      Array<double> maxhdom;
      maxhdom.SetSize (geom.NrSolids());
      maxhdom = mparam.maxh;

      mesh.SetMaxHDomain (maxhdom);

      Box<3> bb = geom.GetBoundingBox();
      bb.Increase (bb.Diam()/10);

      mesh.SetLocalH (bb.PMin(), bb.PMax(), 0.5);

      if (mparam.uselocalh)
      {
         const char * savetask = multithread.task;
         multithread.percent = 0;

         mesh.SetLocalH (bb.PMin(), bb.PMax(), mparam.grading);

         int nedges = geom.emap.Extent();

         double maxedgelen = 0;
         double minedgelen = 1e99;

         multithread.task = "Setting local mesh size (elements per edge)";

         // setting elements per edge

         for (int i = 1; i <= nedges && !multithread.terminate; i++)
         {
            TopoDS_Edge e = TopoDS::Edge (geom.emap(i));
            multithread.percent = 100 * (i-1)/double(nedges);
            if (BRep_Tool::Degenerated(e)) continue;

            GProp_GProps system;
            BRepGProp::LinearProperties(e, system);
            double len = system.Mass();

            if (len < IGNORECURVELENGTH)
            {
               (*testout) << "ignored" << endl;
               continue;
            }

            double localh = len/mparam.segmentsperedge;
            double s0, s1;

            // Philippose - 23/01/2009
            // Find all the parent faces of a given edge
            // and limit the mesh size of the edge based on the
            // mesh size limit of the face
            TopTools_IndexedDataMapOfShapeListOfShape edge_face_map;
            edge_face_map.Clear();

            TopExp::MapShapesAndAncestors(geom.shape, TopAbs_EDGE, TopAbs_FACE, edge_face_map);
            const TopTools_ListOfShape& parent_faces = edge_face_map.FindFromKey(e);

            TopTools_ListIteratorOfListOfShape parent_face_list;

            for(parent_face_list.Initialize(parent_faces); parent_face_list.More(); parent_face_list.Next())
            {
               TopoDS_Face parent_face = TopoDS::Face(parent_face_list.Value());

               int face_index = geom.fmap.FindIndex(parent_face);

               if(face_index >= 1) localh = min(localh,geom.face_maxh[face_index - 1]);
            }

            Handle(Geom_Curve) c = BRep_Tool::Curve(e, s0, s1);

            maxedgelen = max (maxedgelen, len);
            minedgelen = min (minedgelen, len);

            // Philippose - 23/01/2009
            // Modified the calculation of maxj, because the
            // method used so far always results in maxj = 2,
            // which causes the localh to be set only at the
            // starting, mid and end of the edge.
            // Old Algorithm:
            // int maxj = 2 * (int) ceil (localh/len);
            int maxj = max((int) ceil(len/localh), 2);

            for (int j = 0; j <= maxj; j++)
            {
               gp_Pnt pnt = c->Value (s0+double(j)/maxj*(s1-s0));
               mesh.RestrictLocalH (Point3d(pnt.X(), pnt.Y(), pnt.Z()), localh);
            }
         }

         multithread.task = "Setting local mesh size (edge curvature)";

         // setting edge curvature

         int nsections = 20;

         for (int i = 1; i <= nedges && !multithread.terminate; i++)
         {
            double maxcur = 0;
            multithread.percent = 100 * (i-1)/double(nedges);
            TopoDS_Edge edge = TopoDS::Edge (geom.emap(i));
            if (BRep_Tool::Degenerated(edge)) continue;
            double s0, s1;
            Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);
            BRepAdaptor_Curve brepc(edge);
            BRepLProp_CLProps prop(brepc, 2, 1e-5);

            for (int j = 1; j <= nsections; j++)
            {
               double s = s0 + j/(double) nsections * (s1-s0);
               prop.SetParameter (s);
               double curvature = prop.Curvature();
               if(curvature> maxcur) maxcur = curvature;

               if (curvature >= 1e99)
                  continue;

               gp_Pnt pnt = c->Value (s);

               mesh.RestrictLocalH (Point3d(pnt.X(), pnt.Y(), pnt.Z()), ComputeH (fabs(curvature)));
            }
            // (*testout) << "edge " << i << " max. curvature: " << maxcur << endl;
         }

         multithread.task = "Setting local mesh size (face curvature)";

         // setting face curvature

         int nfaces = geom.fmap.Extent();

         for (int i = 1; i <= nfaces && !multithread.terminate; i++)
         {
            multithread.percent = 100 * (i-1)/double(nfaces);
            TopoDS_Face face = TopoDS::Face(geom.fmap(i));
            TopLoc_Location loc;
            Handle(Geom_Surface) surf = BRep_Tool::Surface (face);
            Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation (face, loc);

            if (triangulation.IsNull()) continue;

            BRepAdaptor_Surface sf(face, Standard_True);
            BRepLProp_SLProps prop(sf, 2, 1e-5);

            int ntriangles = triangulation -> NbTriangles();
            for (int j = 1; j <= ntriangles; j++)
            {
               gp_Pnt p[3];
               gp_Pnt2d par[3];

               for (int k = 1; k <=3; k++)
               {
                  int n = triangulation->Triangles()(j)(k);
                  p[k-1] = triangulation->Nodes()(n).Transformed(loc);
                  par[k-1] = triangulation->UVNodes()(n);
               }

               //double maxside = 0;
               //maxside = max (maxside, p[0].Distance(p[1]));
               //maxside = max (maxside, p[0].Distance(p[2]));
               //maxside = max (maxside, p[1].Distance(p[2]));
               //cout << "\rFace " << i << " pos11 ntriangles " << ntriangles << " maxside " << maxside << flush;

               RestrictHTriangle (par[0], par[1], par[2], &prop, mesh, 0);
               //cout << "\rFace " << i << " pos12 ntriangles " << ntriangles << flush;
            }
         }

         // setting close edges

         if (occparam.resthcloseedgeenable)
         {
            multithread.task = "Setting local mesh size (close edges)";

            int sections = 100;

            Array<Line> lines(sections*nedges);

            Box3dTree* searchtree =
               new Box3dTree (bb.PMin(), bb.PMax());

            int nlines = 0;
            for (int i = 1; i <= nedges && !multithread.terminate; i++)
            {
               TopoDS_Edge edge = TopoDS::Edge (geom.emap(i));
               if (BRep_Tool::Degenerated(edge)) continue;

               double s0, s1;
               Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);
               BRepAdaptor_Curve brepc(edge);
               BRepLProp_CLProps prop(brepc, 1, 1e-5);
               prop.SetParameter (s0);

               gp_Vec d0 = prop.D1().Normalized();
               double s_start = s0;
               int count = 0;
               for (int j = 1; j <= sections; j++)
               {
                  double s = s0 + (s1-s0)*(double)j/(double)sections;
                  prop.SetParameter (s);
                  gp_Vec d1 = prop.D1().Normalized();
                  double cosalpha = fabs(d0*d1);
                  if ((j == sections) || (cosalpha < cos(10.0/180.0*M_PI)))
                  {
                     count++;
                     gp_Pnt p0 = c->Value (s_start);
                     gp_Pnt p1 = c->Value (s);
                     lines[nlines].p0 = Point<3> (p0.X(), p0.Y(), p0.Z());
                     lines[nlines].p1 = Point<3> (p1.X(), p1.Y(), p1.Z());

                     Box3d box;
                     box.SetPoint (Point3d(lines[nlines].p0));
                     box.AddPoint (Point3d(lines[nlines].p1));

                     searchtree->Insert (box.PMin(), box.PMax(), nlines+1);
                     nlines++;

                     s_start = s;
                     d0 = d1;
                  }
               }
            }

            Array<int> linenums;

            for (int i = 0; i < nlines; i++)
            {
               multithread.percent = (100*i)/double(nlines);
               Line & line = lines[i];

               Box3d box;
               box.SetPoint (Point3d(line.p0));
               box.AddPoint (Point3d(line.p1));
               double maxhline = max (mesh.GetH(box.PMin()),
                  mesh.GetH(box.PMax()));
               box.Increase(maxhline);

               double mindist = 1e99;
               linenums.SetSize(0);
               searchtree->GetIntersecting(box.PMin(),box.PMax(),linenums);

               for (int j = 0; j < linenums.Size(); j++)
               {
                  int num = linenums[j]-1;
                  if (i == num) continue;
                  if ((line.p0-lines[num].p0).Length2() < 1e-15) continue;
                  if ((line.p0-lines[num].p1).Length2() < 1e-15) continue;
                  if ((line.p1-lines[num].p0).Length2() < 1e-15) continue;
                  if ((line.p1-lines[num].p1).Length2() < 1e-15) continue;
                  mindist = min (mindist, line.Dist(lines[num]));
               }

               mindist /= (occparam.resthcloseedgefac + VSMALL);

               if (mindist < 1e-3)
               {
                  (*testout) << "extremely small local h: " << mindist
                     << " --> setting to 1e-3" << endl;
                  (*testout) << "somewhere near " << line.p0 << " - " << line.p1 << endl;
                  mindist = 1e-3;
               }

               mesh.RestrictLocalHLine(line.p0, line.p1, mindist);
            }
         }

         multithread.task = savetask;

      }

      // Philippose - 09/03/2009
      // Added the capability to load the mesh size from a 
      // file also for OpenCascade Geometry
      // Note: 
      // ** If the "uselocalh" option is ticked in 
      // the "mesh options...insider" menu, the mesh 
      // size will be further modified by the topology 
      // analysis routines.
      // ** To use the mesh size file as the sole source 
      // for defining the mesh size, uncheck the "uselocalh"
      // option.
      mesh.LoadLocalMeshSize (mparam.meshsizefilename);
   }



  int OCCGenerateMesh (OCCGeometry & geom, Mesh *& mesh, MeshingParameters & mparam,
		       int perfstepsstart, int perfstepsend)
   {
      multithread.percent = 0;

      if (perfstepsstart <= MESHCONST_ANALYSE)
      {
         delete mesh;
         mesh = new Mesh();
         mesh->geomtype = Mesh::GEOM_OCC;

         OCCSetLocalMeshSize(geom,*mesh);
      }

      if (multithread.terminate || perfstepsend <= MESHCONST_ANALYSE)
         return TCL_OK;

      if (perfstepsstart <= MESHCONST_MESHEDGES)
      {
         OCCFindEdges (geom, *mesh);

         /*
         cout << "Removing redundant points" << endl;

         int i, j;
         int np = mesh->GetNP();
         Array<int> equalto;

         equalto.SetSize (np);
         equalto = 0;

         for (i = 1; i <= np; i++)
         {
         for (j = i+1; j <= np; j++)
         {
         if (!equalto[j-1] && (Dist2 (mesh->Point(i), mesh->Point(j)) < 1e-12))
         equalto[j-1] = i;
         }
         }

         for (i = 1; i <= np; i++)
         if (equalto[i-1])
         {
         cout << "Point " << i << " is equal to Point " << equalto[i-1] << endl;
         for (j = 1; j <= mesh->GetNSeg(); j++)
         {
         Segment & seg = mesh->LineSegment(j);
         if (seg[0] == i) seg[0] = equalto[i-1];
         if (seg[1] == i) seg[1] = equalto[i-1];
         }
         }

         cout << "Removing degenerated segments" << endl;
         for (j = 1; j <= mesh->GetNSeg(); j++)
         {
         Segment & seg = mesh->LineSegment(j);
         if (seg[0] == seg[1])
         {
         mesh->DeleteSegment(j);
         cout << "Deleting Segment " << j << endl;
         }
         }

         mesh->Compress();
         */

         /*
         for (int i = 1; i <= geom.fmap.Extent(); i++)
         {
         Handle(Geom_Surface) hf1 =
         BRep_Tool::Surface(TopoDS::Face(geom.fmap(i)));
         for (int j = i+1; j <= geom.fmap.Extent(); j++)
         {
         Handle(Geom_Surface) hf2 =
         BRep_Tool::Surface(TopoDS::Face(geom.fmap(j)));
         if (hf1 == hf2) cout << "face " << i << " and face " << j << " lie on same surface" << endl;
         }
         }
         */

#ifdef LOG_STREAM
         (*logout) << "Edges meshed" << endl
            << "time = " << GetTime() << " sec" << endl
            << "points: " << mesh->GetNP() << endl;
#endif
      }

      if (multithread.terminate || perfstepsend <= MESHCONST_MESHEDGES)
         return TCL_OK;

      if (perfstepsstart <= MESHCONST_MESHSURFACE)
      {
         OCCMeshSurface (geom, *mesh, perfstepsend);
         if (multithread.terminate) return TCL_OK;

#ifdef LOG_STREAM
         (*logout) << "Surfaces meshed" << endl
            << "time = " << GetTime() << " sec" << endl
            << "points: " << mesh->GetNP() << endl;
#endif

#ifdef STAT_STREAM
         (*statout) << mesh->GetNSeg() << " & "
            << mesh->GetNSE() << " & - &"
            << GetTime() << " & " << endl;
#endif

         //      MeshQuality2d (*mesh);
         mesh->CalcSurfacesOfNode();
      }

      if (multithread.terminate || perfstepsend <= MESHCONST_OPTSURFACE)
         return TCL_OK;

      if (perfstepsstart <= MESHCONST_MESHVOLUME)
      {
         multithread.task = "Volume meshing";

         MESHING3_RESULT res = MeshVolume (mparam, *mesh);

/*
         ofstream problemfile("occmesh.rep",ios_base::app);

         problemfile << "VOLUMEMESHING" << endl << endl;
         if(res != MESHING3_OK)
            problemfile << "ERROR" << endl << endl;
         else
            problemfile << "OK" << endl
            << mesh->GetNE() << " elements" << endl << endl;

         problemfile.close();
*/

         if (res != MESHING3_OK) return TCL_ERROR;

         if (multithread.terminate) return TCL_OK;

         RemoveIllegalElements (*mesh);
         if (multithread.terminate) return TCL_OK;

         MeshQuality3d (*mesh);

#ifdef STAT_STREAM
         (*statout) << GetTime() << " & ";
#endif

#ifdef LOG_STREAM
         (*logout) << "Volume meshed" << endl
            << "time = " << GetTime() << " sec" << endl
            << "points: " << mesh->GetNP() << endl;
#endif
      }

      if (multithread.terminate || perfstepsend <= MESHCONST_MESHVOLUME)
         return TCL_OK;

      if (perfstepsstart <= MESHCONST_OPTVOLUME)
      {
         multithread.task = "Volume optimization";

         OptimizeVolume (mparam, *mesh);
         if (multithread.terminate) return TCL_OK;

#ifdef STAT_STREAM
         (*statout) << GetTime() << " & "
            << mesh->GetNE() << " & "
            << mesh->GetNP() << " " << '\\' << '\\' << " \\" << "hline" << endl;
#endif

#ifdef LOG_STREAM
         (*logout) << "Volume optimized" << endl
            << "time = " << GetTime() << " sec" << endl
            << "points: " << mesh->GetNP() << endl;
#endif

         // cout << "Optimization complete" << endl;

      }

      (*testout) << "NP: " << mesh->GetNP() << endl;
      for (int i = 1; i <= mesh->GetNP(); i++)
         (*testout) << mesh->Point(i) << endl;

      (*testout) << endl << "NSegments: " << mesh->GetNSeg() << endl;
      for (int i = 1; i <= mesh->GetNSeg(); i++)
         (*testout) << mesh->LineSegment(i) << endl;

      return TCL_OK;
   }
}

#endif
