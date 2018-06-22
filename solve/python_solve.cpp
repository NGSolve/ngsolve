#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <solve.hpp>
using namespace ngsolve;


typedef GridFunction GF;

extern void ExportBVP(py::module &m);
extern void ExportDrawFlux(py::module &m);

void NGS_DLL_HEADER ExportNgsolve(py::module &m ) {

    m.def ("Tcl_Eval", &Ng_TclCmd);

    m.def ("_Redraw",
           ([](bool blocking, double fr)
             {
               static auto last_time = std::chrono::system_clock::now()-std::chrono::seconds(10);
               auto now = std::chrono::system_clock::now();
               double elapsed = std::chrono::duration<double>(now-last_time).count();
               if (elapsed * fr > 1)
                 {
                   Ng_Redraw(blocking);
                   last_time = std::chrono::system_clock::now();
                 }
             }),
           py::arg("blocking")=false, py::arg("fr") = 25, docu_string(R"raw_string(
Redraw all

Parameters:

blocking : bool
  input blocking

fr : double
  input framerate

)raw_string")
             );


    m.def ("Draw",[](shared_ptr<MeshAccess> mesh) 
                                     {
                                       mesh->SelectMesh();
                                       Ng_TclCmd ("set ::selectvisual mesh;\n");
                                     }
             );


    m.def("SetVisualization",
            [](py::object deformation, py::object min, py::object max,
                /* py::object clippnt, */ py::object clipnormal, py::object clipping)
             {
               bool need_redraw = false;
               if (py::extract<bool>(deformation).check())
                 {
                   bool def = py::extract<bool>(deformation)();
                   Ng_TclCmd ("set ::visoptions.deformation "+ToString(def)+";\n");
                   Ng_TclCmd ("Ng_Vis_Set parameters;\n");
                   need_redraw = true;
                 }
               if (py::extract<double>(min).check())
                 {
                   Ng_TclCmd ("set ::visoptions.autoscale 0\n");
                   Ng_TclCmd ("set ::visoptions.mminval "+ToString(py::extract<double>(min)())+";\n");
                   Ng_TclCmd ("Ng_Vis_Set parameters;\n");                                      
                   need_redraw = true;
                 }
               if (py::extract<double>(max).check())
                 {
                   Ng_TclCmd ("set ::visoptions.autoscale 0\n");
                   Ng_TclCmd ("set ::visoptions.mmaxval "+ToString(py::extract<double>(max)())+";\n");
                   Ng_TclCmd ("Ng_Vis_Set parameters;\n");                   
                   need_redraw = true;
                 }
               if (py::extract<py::tuple>(clipnormal).check())
                 {
                   py::tuple norm = py::extract<py::tuple>(clipnormal)();
                   if (py::len(norm)==3)
                     {
                       // cout << "setting clipping normal" << endl;
                       // tclstring << "set ::viewoptions.clipping.enable 1" << endl
                       Ng_TclCmd ("set ::viewoptions.clipping.nx "+ToString(py::extract<double>(norm[0])())+";\n");
                       Ng_TclCmd ("set ::viewoptions.clipping.ny "+ToString(py::extract<double>(norm[1])())+";\n");
                       Ng_TclCmd ("set ::viewoptions.clipping.nz "+ToString(py::extract<double>(norm[2])())+";\n");
                       // << "set ::viewoptions.clipping.dist " << clipdist << endl;
                       need_redraw = true;
                     }
                 }
               if (py::extract<bool>(clipping).check())
                 {
                   bool clip = py::extract<bool>(clipping)();
                   Ng_TclCmd ("set ::viewoptions.clipping.enable "+ToString(int(clip))+";\n");
                   Ng_TclCmd ("Ng_SetVisParameters");
                   
                   need_redraw = true;
                 }
               if (need_redraw)
                 Ng_Redraw(true);
             },
             py::arg("deformation")=DummyArgument(),
             py::arg("min")=DummyArgument(),
             py::arg("max")=DummyArgument(),
             // py::arg("clippnt")=DummyArgument(),
             py::arg("clipnormal")=DummyArgument(),
          py::arg("clipping")=DummyArgument(), docu_string(R"raw_string(
Set visualization options

Parameters:

deformation : object
            input deformation

min : object
  input min

max : object
  input max

clipnormal : object
  input clipnormal

clipping : object
  input clipping

)raw_string")
            )
      ;
    
    m.def ("Draw",
           [](shared_ptr<GridFunction> gf, int sd, bool autoscale, double min, double max)
              {
                // cout << "in draw" << endl;
                // cout << "gf in draw = " << *gf << endl;
                // cout << "dims of gf =  " << gf->Dimensions() << endl;
                gf->GetMeshAccess()->SelectMesh();
                // Visualize (gf, gf->GetName());

                // netgen::SolutionData * vis = new VisualizeCoefficientFunction (gf->GetMeshAccess(), gf);
                auto gfcf = make_shared<GridFunctionCoefficientFunction> (gf);
                netgen::SolutionData * vis = new VisualizeCoefficientFunction (gf->GetMeshAccess(), gfcf);                

                Ng_SolutionData soldata;
                Ng_InitSolutionData (&soldata);
  
                soldata.name = gf->GetName().c_str();
                soldata.data = 0;
                soldata.components = gf -> Dimension();
                if (gf->IsComplex()) soldata.components *= 2;
                soldata.iscomplex = gf -> IsComplex();
                soldata.draw_surface = true; // draw_surf;
                soldata.draw_volume  = true; // draw_vol; 
                /* 
                if (flags.GetDefineFlag("volume"))
                  soldata.draw_surface = false;
                if (flags.GetDefineFlag("boundary"))
                  soldata.draw_volume = false;
                */
                soldata.dist = 1;
                soldata.soltype = NG_SOLUTION_VIRTUAL_FUNCTION;
                soldata.solclass = vis;
                Ng_SetSolutionData (&soldata);
                
                if (gf->Dimension() == 1)
                  Ng_TclCmd (string("set ::visoptions.scalfunction ")+gf->GetName()+":1;\n");
                else
                  Ng_TclCmd (string("set ::visoptions.vecfunction ")+gf->GetName()+";\n");
                Ng_TclCmd (string("set ::visoptions.subdivisions ")+ToString(sd)+";\n");
		Ng_TclCmd ("set ::visoptions.autoscale "+ToString(autoscale)+";\n");
		if(!autoscale){
		  Ng_TclCmd ("set ::visoptions.mminval "+ToString(min)+";\n");
		  Ng_TclCmd ("set ::visoptions.mmaxval "+ToString(max)+";\n");
		}
		Ng_TclCmd ("Ng_Vis_Set parameters;\n");
                Ng_TclCmd ("set ::selectvisual solution;\n");
              },
             py::arg("gf"),py::arg("sd")=2,py::arg("autoscale")=true,
           py::arg("min")=0.0,py::arg("max")=1.0, docu_string(R"raw_string(
Parameters:

gf : ngsolve.comp.GridFunction
  input GridFunction to draw

sd : int
  input subdivisions

autoscale : bool
  input autscale

min : float
  input minimum value. Need autscale = false

max : float
  input maximum value. Need autscale = false

)raw_string")
             );

    
    m.def ("Draw", [](shared_ptr<CoefficientFunction> cf, shared_ptr<MeshAccess> ma, string name,
                 int sd, bool autoscale, double min, double max,
                 bool draw_vol, bool draw_surf) 
              {
                ma->SelectMesh();
                netgen::SolutionData * vis;
                if(dynamic_cast<ProlongateCoefficientFunction *>(cf.get()))
                {
                  shared_ptr<CoefficientFunction> wrapper(new ProlongateCoefficientFunctionVisualization(*static_cast<ProlongateCoefficientFunction *>(cf.get())));
                  vis = new VisualizeCoefficientFunction (ma, wrapper);
                }
                else
                  vis = new VisualizeCoefficientFunction (ma, cf);
                Ng_SolutionData soldata;
                Ng_InitSolutionData (&soldata);
  
                soldata.name = (char*)name.c_str();
                soldata.data = 0;
                soldata.components = cf -> Dimension();
                if (cf->IsComplex()) soldata.components *= 2;
                soldata.iscomplex = cf -> IsComplex();
                soldata.draw_surface = draw_surf;
                soldata.draw_volume  = draw_vol; 
                /* 
                if (flags.GetDefineFlag("volume"))
                  soldata.draw_surface = false;
                if (flags.GetDefineFlag("boundary"))
                  soldata.draw_volume = false;
                */
                soldata.dist = 1;
                soldata.soltype = NG_SOLUTION_VIRTUAL_FUNCTION;
                soldata.solclass = vis;
                Ng_SetSolutionData (&soldata);
                

                if (cf->Dimension() == 1)
                  Ng_TclCmd (string("set ::visoptions.scalfunction ")+name+":1;\n");
                else
                  if (cf->Dimension() == 3 || cf->Dimension() == ma->GetDimension())
                    Ng_TclCmd (string("set ::visoptions.vecfunction ")+name+";\n");

                Ng_TclCmd (string("set ::visoptions.subdivisions ")+ToString(sd)+";\n");
		Ng_TclCmd ("set ::visoptions.autoscale "+ToString(autoscale)+";\n");
		if(!autoscale){
		  Ng_TclCmd ("set ::visoptions.mminval "+ToString(min)+";\n");
		  Ng_TclCmd ("set ::visoptions.mmaxval "+ToString(max)+";\n");
		}
                Ng_TclCmd ("Ng_Vis_Set parameters;\n");
                Ng_TclCmd ("set ::selectvisual solution;\n");

              },
              py::arg("cf"),py::arg("mesh"),py::arg("name"),
              py::arg("sd")=2,py::arg("autoscale")=true,
	      py::arg("min")=0.0,py::arg("max")=1.0,
              py::arg("draw_vol")=true,py::arg("draw_surf")=true, docu_string(R"raw_string(
Parameters:

cf : ngsolve.comp.CoefficientFunction
  input CoefficientFunction to draw

mesh : ngsolve.comp.Mesh
  input mesh

name : string
  input name

sd : int
  input subdivisions

autoscale : bool
  input autscale

min : float
  input minimum value. Need autscale = false

max : float
  input maximum value. Need autscale = false

draw_vol : bool
  input draw volume

draw_surf : bool
  input draw surface

)raw_string")
             );





    ExportBVP(m);
    ExportDrawFlux(m);
}


PYBIND11_MODULE(libngsolve, m) {
  m.attr("__name__") = "ngsolve";
  ExportNgsolve(m);
}


#endif
