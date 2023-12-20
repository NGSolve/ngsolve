/*********************************************************************/
/* File:   vtkoutput.cpp                                             */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   1. June 2014                                              */
/*********************************************************************/

#include <comp.hpp>
#include "vtkoutput.hpp"

namespace ngcomp
{

  ValueField::ValueField(int adim, string aname) : Array<double>(), dim(adim), name(aname) { ; }

  template <int D>
  VTKOutput<D>::VTKOutput(const Array<shared_ptr<CoefficientFunction>> &a_coefs,
                          const Flags &flags,
                          shared_ptr<MeshAccess> ama)
      : VTKOutput(ama, a_coefs,
                  flags.GetStringListFlag("fieldnames"),
                  flags.GetStringFlag("filename", "output"),
                  (int)flags.GetNumFlag("subdivision", 0),
                  (int)flags.GetNumFlag("only_element", -1),
                  flags.GetStringFlag("floatsize", "double"),
                  flags.GetDefineFlag("legacy"),
                  (int)flags.GetNumFlag("order", 1))
  {
    ;
  }

  template <int D>
  VTKOutput<D>::VTKOutput(shared_ptr<MeshAccess> ama,
                          const Array<shared_ptr<CoefficientFunction>> &a_coefs,
                          const Array<string> &a_field_names,
                          string a_filename, int a_subdivision, int a_only_element, 
                          string a_floatsize, bool a_legacy, int a_order)
      : ma(ama), coefs(a_coefs), fieldnames(a_field_names),
        filename(a_filename), subdivision(a_subdivision), order(a_order), only_element(a_only_element), floatsize(a_floatsize), legacy(a_legacy)
  {
    r = 1 << (subdivision + order -1);
    h = 1.0/r;
    if ((floatsize != "double") && (floatsize != "float") && (floatsize != "single"))
      cout << IM(1) << "VTKOutput: floatsize is not int {\"double\",\"single\",\"float\"}. Using \"float|single\".";
    value_field.SetSize(a_coefs.Size());
    for (int i = 0; i < a_coefs.Size(); i++)
      if (fieldnames.Size() > i)
        value_field[i] = make_shared<ValueField>(coefs[i]->Dimension(), fieldnames[i]);
      else
        value_field[i] = make_shared<ValueField>(coefs[i]->Dimension(), "dummy" + to_string(i));
  }

  /// Empty all field
  template <int D>
  void VTKOutput<D>::ResetArrays()
  {
    points.SetSize(0);
    cells.SetSize(0);
    for (auto field : value_field)
      field->SetSize(0);
  }

  int VTKCell :: GetVtkType(ELEMENT_TYPE et, int order)
  {
    if(order==1)
    {
      switch(et) {
        case ET_POINT: return VTK_VERTEX;
        case ET_SEGM: return VTK_LINE;
        case ET_TRIG: return VTK_TRIANGLE;
        case ET_QUAD: return VTK_QUAD;
        case ET_TET: return VTK_TETRA;
        case ET_HEX: return VTK_HEXAHEDRON;
        case ET_PRISM: return VTK_WEDGE;
        case ET_PYRAMID: return VTK_PYRAMID;
        case ET_HEXAMID: throw Exception("hexamid not handled in vtk output");
      }
    }
    else if(order==2)
    {
      switch(et) {
        case ET_POINT: throw Exception("Have no second order points");
        case ET_SEGM: return VTK_QUADRATIC_EDGE;
        case ET_TRIG: return VTK_QUADRATIC_TRIANGLE;
        case ET_QUAD: return VTK_BIQUADRATIC_QUAD; // 9 nodes
        case ET_TET: return VTK_QUADRATIC_TETRA;
        // case ET_HEX: return VTK_TRIQUADRATIC_HEXAHEDRON; // 27 nodes (not supported by paraview tessellation)
        case ET_HEX: return VTK_QUADRATIC_HEXAHEDRON; // 20 nodes
        case ET_PRISM: return VTK_BIQUADRATIC_QUADRATIC_WEDGE; // 18 nodes
        // case ET_PRISM: return VTK_QUADRATIC_WEDGE; // 18 nodes (not supported by paraview tessellation)
        case ET_PYRAMID: return VTK_PYRAMID; //VTK_TRIQUADRATIC_PYRAMID; // 19 nodes
        case ET_HEXAMID: throw Exception("hexamid not handled in vtk output");          
      }
    }
    else
      throw Exception("Invalid element order: " + ToString(order));
  }
 
  VTKCell :: VTKCell(ELEMENT_TYPE et, int order, const std::map<tuple<int,int,int>, int> &m,
      int i, int j, int k, Vec<3,int> vi, Vec<3,int> vj, Vec<3,int> vk)
  {
    type = GetVtkType(et, order);
    auto f = [&m, i,j,k, vi, vj, vk](int i0,int j0, int k0) {
      auto ijk = Vec<3,int>{i,j,k};
      ijk += i0*vi+j0*vj+k0*vk;
      return m.at({ijk[0], ijk[1], ijk[2]});
    };

    switch(type) {
      case VTK_TRIANGLE: pi = {f(0,0,0), f(1,0,0), f(0,1,0)}; break;
      case VTK_QUAD: pi = {f(0,0,0), f(0,1,0), f(1,1,0), f(1,0,0)}; break;
      case VTK_TETRA: pi = {f(0,0,0), f(0,0,1), f(0,1,0), f(1,0,0)}; break;
      case VTK_WEDGE: pi = {f(0,0,0), f(1,0,0), f(0,1,0), f(0,0,1), f(1,0,1), f(0,1,1)}; break;
      case VTK_HEXAHEDRON: pi = {f(0,0,0), f(0,0,1), f(0,1,1), f(0,1,0), f(1,0,0), f(1,0,1), f(1,1,1), f(1,1,0)}; break;
      case VTK_PYRAMID: pi = {f(0,0,0), f(0,1,0), f(1,1,0), f(1,0,0), f(0,0,1)}; break;

      case VTK_QUADRATIC_TRIANGLE: pi = {f(0,0,0), f(2,0,0), f(0,2,0),
                                         f(1,0,0), f(1,1,0), f(0,1,0)}; break;
      case VTK_BIQUADRATIC_QUAD: pi = {f(0,0,0), f(0,2,0), f(2,2,0), f(2,0,0),
                                       f(0,1,0), f(1,2,0), f(2,1,0), f(1,0,0), f(1,1,0)}; break;

      case VTK_QUADRATIC_TETRA: pi = { f(0,0,0), f(0,0,2), f(0,2,0), f(2,0,0),
                                       f(0,0,1), f(0,1,1), f(0,1,0), f(1,0,0),
                                       f(1,0,1), f(1,1,0) }; break;
      case VTK_TRIQUADRATIC_HEXAHEDRON: pi = { f(0,0,0), f(0,0,2), f(0,2,2), f(0,2,0), f(2,0,0), f(2,0,2), f(2,2,2), f(2,2,0),
                                               f(0,0,1), f(0,1,2), f(0,2,1), f(0,1,0), f(2,0,1), f(2,1,2), f(2,2,1), f(2,1,0),
                                               f(1,0,0), f(1,0,2), f(1,2,2), f(1,2,0),
                                               f(1,1,0), f(1,1,2), f(1,0,1), f(1,2,1), f(0,1,1), f(2,1,1), f(1,1,1)}; break;
      case VTK_QUADRATIC_HEXAHEDRON: pi = { f(0,0,0), f(0,0,2), f(0,2,2), f(0,2,0), f(2,0,0), f(2,0,2), f(2,2,2), f(2,2,0),
                                               f(0,0,1), f(0,1,2), f(0,2,1), f(0,1,0), f(2,0,1), f(2,1,2), f(2,2,1), f(2,1,0),
                                               f(1,0,0), f(1,0,2), f(1,2,2), f(1,2,0)}; break;

      case VTK_BIQUADRATIC_QUADRATIC_WEDGE: pi = {f(0,0,0), f(2,0,0), f(0,2,0), f(0,0,2), f(2,0,2), f(0,2,2),
                                                  f(1,0,0), f(1,1,0), f(0,1,0), f(1,0,2), f(1,1,2), f(0,1,2),
                                                  f(0,0,1), f(2,0,1), f(0,2,1), f(1,0,1), f(1,1,1), f(0,1,1)}; break;

      case VTK_QUADRATIC_WEDGE: pi = { f(0,0,0), f(2,0,0), f(0,2,0), f(0,0,2), f(2,0,2), f(0,2,2),
                                       f(1,0,0), f(1,1,0), f(0,1,0), f(1,0,2), f(1,1,2), f(0,1,2),
                                       f(0,0,1), f(2,0,1), f(0,2,1)}; break;

      default: throw Exception("VTK type not implemented: "+ToString(type));
    }
  }

  /// Fill principil lattices (points and connections on subdivided reference simplex) in 2D
  template <int D>
  void VTKOutput<D>::FillReferenceTrig(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems)
  {
    std::map<tuple<int,int,int>, int> m;

    for (int i = 0; i <= r; ++i)
      for (int j = 0; i + j <= r; ++j)
      {
        m[{i,j,0}] = ref_coords.Size();
        ref_coords.Append(IntegrationPoint(j * h, i * h));
      }

    for (int i = 0; i < r; i+=order)
      for (int j = 0; i + j < r; j+=order)
      {
        if(i+j+order<r)
          ref_elems.Append({ET_QUAD, order, m, i, j, 0});
        else
          ref_elems.Append({ET_TRIG, order, m, i, j, 0});
      }
  }

  /// Fill principil lattices (points and connections on subdivided reference simplex) in 2D
  template <int D>
  void VTKOutput<D>::FillReferenceQuad(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems)
  {
    std::map<tuple<int,int,int>, int> m;

    for (int i : Range(r+1))
      for (int j : Range(r+1))
      {
        m[{i,j,0}] = ref_coords.Size();
        ref_coords.Append(IntegrationPoint(j * h, i * h));
      }

    for (int i=0; i<r; i+=order)
      for (int j=0; j<r; j+=order)
        ref_elems.Append({ET_QUAD, order, m, i, j, 0});
  }

  /// Fill principil lattices (points and connections on subdivided reference simplex) in 3D
  template <int D>
  void VTKOutput<D>::FillReferenceTet(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems)
  {
    std::map<tuple<int,int,int>, int> m;

    for (int i = 0; i <= r; ++i)
      for (int j = 0; i + j <= r; ++j)
        for (int k = 0; i + j + k <= r; ++k)
        {
          m[{i,j,k}] = ref_coords.Size();
          ref_coords.Append(IntegrationPoint(i * h, j * h, k * h));
        }

    for (int i = 0; i < r; i+=order)
      for (int j = 0; i + j < r; j+=order)
        for (int k = 0; i + j + k < r; k+=order)
        {
          if(i+j+k+2*order<r)
            ref_elems.Append({ET_HEX, order, m, i, j, k});
          else if(i+j+k+order<r)
          {
            ref_elems.Append({ET_PRISM, order, m, i, j, k});
            ref_elems.Append({ET_TET, order, m, i+order, j, k, {-1,1,0}});
            ref_elems.Append({ET_TET, order, m, i, j+order, k, {1,0,0}, {1,-1,1}});
          }
          else
            ref_elems.Append({ET_TET, order, m, i, j, k});
        }
  }

  /// Fill principil lattices (points and connections on subdivided reference hexahedron) in 3D
  template <int D>
  void VTKOutput<D>::FillReferenceHex(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems)
  {
    std::map<tuple<int,int,int>, int> m;

    for (int i = 0; i <= r; ++i)
      for (int j = 0; j <= r; ++j)
        for (int k = 0; k <= r; ++k)
        {
          m[{i,j,k}] = ref_coords.Size();
          ref_coords.Append(IntegrationPoint(k * h, j * h, i * h));
        }

    for (int i = 0; i < r; i+=order)
      for (int j = 0; j < r; j+=order)
        for (int k = 0; k < r; k+=order)
          ref_elems.Append({ET_HEX, order, m, i, j, k});
  }

  template <int D>
  void VTKOutput<D>::FillReferencePrism(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems)
  {
    std::map<tuple<int,int,int>, int> m;

    for (int k = 0; k <= r; k++)
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i + j <= r; ++j)
        {
          m[{i,j,k}] = ref_coords.Size();
          ref_coords.Append(IntegrationPoint(j * h, i * h, k * h));
        }

    for (int k = 0; k < r; k+=order)
      for (int i = 0; i < r; i+=order)
        for (int j = 0; i + j < r; j+=order)
        {
          if (i + j + order < r)
            ref_elems.Append({ET_HEX, order, m, i, j, k});
          else
            ref_elems.Append({ET_PRISM, order, m, i, j, k});
        }
  }

  template <int D>
  void VTKOutput<D>::FillReferencePyramid(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems)
  {
    std::map<tuple<int,int,int>, int> m;

    for (int k = 0; k + order <= r; k++)
      for (int i = 0; i <= r; ++i)
        for (int j = 0; j <= r; ++j)
        {
          m[{i,j,k}] = ref_coords.Size();
          double h1 = (1.0-h*k)*h;
          ref_coords.Append(IntegrationPoint(j * h1, i * h1, k * h));
        }

    m[{0,0,r}] = ref_coords.Size();
    ref_coords.Append(IntegrationPoint(0, 0, 1.0));
    auto s = r/order;
    if(order==2)
    {
      double h1 = 2*(1.0-h*(r-1))*h;
      double hr = (r-1)*h;
      m[{0,0,r-1}] = ref_coords.Size();
      ref_coords.Append(IntegrationPoint(0, 0, hr));
      m[{s,0,r-1}] = ref_coords.Size();
      ref_coords.Append(IntegrationPoint(0, s*h1, hr));
      m[{0,s,r-1}] = ref_coords.Size();
      ref_coords.Append(IntegrationPoint(s*h1, 0, hr));
      m[{s,s,r-1}] = ref_coords.Size();
      ref_coords.Append(IntegrationPoint(s*h1, s*h1, hr));
    }

    for (int k = 0; k + order < r; k+=order)
      for (int i = 0; i < r; i+=order)
        for (int j = 0; j < r; j+=order)
          ref_elems.Append({ET_HEX, order, m, i, j, k});

    // pyramids not well supported yet, use two tets instead
    ref_elems.Append({ET_TET, order, m, 0, 0, r-order, {s,s,0}, {0,s,0}});
    ref_elems.Append({ET_TET, order, m, 0, 0, r-order, {s,0,0}, {s,s,0}});
  }

  /* ###########################
     ###########################
     # Legacy Files as Fallback#
     ###########################*/
  template <int D>
  void VTKOutput<D>::PrintPointsLegacy()
  {
    *fileout << "POINTS " << points.Size() << " float" << endl;
    for (auto p : points)
    {
      *fileout << p;
      if (D == 2)
        *fileout << "\t 0.0";
      *fileout << endl;
    }
  }

  /// output of cells in form vertices
  template <int D>
  void VTKOutput<D>::PrintCellsLegacy()
  {
    // count number of data for cells, one + number of vertices
    int ndata = 0;
    for (auto c : cells)
    {
      ndata++;
      ndata += c.pi.Size();
    }
    *fileout << "CELLS " << cells.Size() << " " << ndata << endl;
    for (auto c : cells)
    {
      int nv = c.pi.Size();
      *fileout << nv << "\t";
      for (int i = 0; i < nv; i++)
        *fileout << c.pi[i] << "\t";
      *fileout << endl;
    }
  }

  /// output of cell types (here only simplices)
  template <int D>
  void VTKOutput<D>::PrintCellTypesLegacy(VorB vb, const BitArray *drawelems)
  {
    *fileout << "CELL_TYPES " << cells.Size() << endl;
    for(auto & c : cells)
      *fileout << c.type << " ";

    *fileout << "CELL_DATA " << cells.Size() << endl;
    *fileout << "POINT_DATA " << points.Size() << endl;
  }

  /// output of field data (coefficient values)
  template <int D>
  void VTKOutput<D>::PrintFieldDataLegacy()
  {
    for (auto field : value_field)
    {
      *fileout << "SCALARS " << field->Name()
               << " float " << field->Dimension() << endl
               << "LOOKUP_TABLE default" << endl;

      for (auto v : *field)
        *fileout << v << " ";
      *fileout << endl;
    }
  }
  /// output of data points, XML file format
  template <int D>
  void VTKOutput<D>::PrintPoints(int *offset, stringstream *appenddata)
  {
    *fileout << "<Points>" << endl;
    if (floatsize == "double")
    {
      *fileout << "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"" << 3 << "\" format=\"appended\" offset=\"0\">" << endl;
    }
    else
    {
      *fileout << "<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"" << 3 << "\" format=\"appended\" offset=\"0\">" << endl;
    }
    const double valdouble = 0;
    double wvaldouble = 0;
    const float val = 0;
    float wval = 0;
    stringstream data;

    uint32_t count = 0;
    if (floatsize == "double")
    {
      for (auto p : points)
      {
        for (int k = 0; k < D; k++)
        {
          wvaldouble = p[k];
          data.write((char *)&wvaldouble, sizeof(wvaldouble));
          count += sizeof(wvaldouble);
        }
        if (D == 2)
        {
          data.write((char *)&valdouble, sizeof(valdouble));
          count += sizeof(valdouble);
        }
      }
    }
    else
    {
      for (auto p : points)
      {
        for (int k = 0; k < D; k++)
        {
          wval = p[k];
          data.write((char *)&wval, sizeof(wval));
          count += sizeof(wval);
        }
        if (D == 2)
        {
          data.write((char *)&val, sizeof(val));
          count += sizeof(val);
        }
      }
    }
    appenddata->write((char *)&count, sizeof(uint32_t));
    *appenddata << data.str();
    *offset = count + 4;
    *fileout << endl
             << "</DataArray>" << endl;
    *fileout << "</Points>" << endl;
  }
  /// output of cells in form vertices
  template <int D>
  void VTKOutput<D>::PrintCells(int *offset, stringstream *appenddata)
  {
    stringstream connectivity;
    stringstream offsets;
    uint32_t sizecon = 0;
    uint32_t sizeoff = 0;
    // count number of data for cells, one + number of vertices
    // int32_t ndata = 0;
    int32_t offs = 0;
    for (auto c : cells)
    {

      int nv = c.pi.Size();
      offs += nv;
      offsets.write((char *)&offs, sizeof(int32_t));
      sizeoff += sizeof(int32_t);
      for (int i = 0; i < nv; i++)
      {
        connectivity.write((char *)&c.pi[i], sizeof(int));
        sizecon += sizeof(int);
      }
    }
    *fileout << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << *offset << "\">" << endl;
    *fileout
        << "</DataArray>" << endl;
    *fileout << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << *offset + 4 + sizecon << "\">" << endl;

    *fileout
        << "</DataArray>" << endl;
    *offset += 8 + sizecon + sizeoff;

    appenddata->write((char *)&sizecon, sizeof(uint32_t));
    *appenddata << connectivity.str();
    appenddata->write((char *)&sizeoff, sizeof(uint32_t));
    *appenddata << offsets.str();
  }

  /// output of cell types (here only simplices)
  template <int D>
  void VTKOutput<D>::PrintCellTypes(VorB vb, int *offset, stringstream *appenddata, const BitArray *drawelems)
  {
    *fileout << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << *offset << "\">" << endl;
    stringstream data;
    uint32_t sizetypes = 0;
    uint8_t eltype;

    for(auto & cell : cells)
    {
      sizetypes += sizeof(uint8_t);
      eltype = cell.type;
      data.write((char*)&eltype, sizeof(uint8_t));
    }

    appenddata->write((char *)&sizetypes, sizeof(uint32_t));
    *appenddata << data.str();
    *offset += sizetypes + 4;

    *fileout << endl
             << "</DataArray>" << endl;
  }

  /// output of field data (coefficient values)
  template <int D>
  void VTKOutput<D>::PrintFieldData(int *offset, stringstream *appenddata)
  {
    string header = "";
    string content = "";
    stringstream data;
    int32_t fieldsize = 0;
    header += "<PointData>"; // \"" << field->Name() << "\">" << endl;
    *fileout << header << endl;
    for (auto field : value_field)
    {
      if (floatsize == "double")
      {
        *fileout << "<DataArray type=\"Float64\" Name=\"" << field->Name() << "\" NumberOfComponents=\"" << field->Dimension() << "\" format=\"appended\" offset=\"" << *offset << "\">" << endl;
        //      *fileout << "<DataArray type=\"Float64\" Name=\"" << field->Name() << "\" NumberOfComponents=\"" << field->Dimension() << "\" format=\"appended\" offset=\"0\">" << endl;
      }
      else
      {
        *fileout << "<DataArray type=\"Float32\" Name=\"" << field->Name() << "\" NumberOfComponents=\"" << field->Dimension() << "\" format=\"appended\" offset=\"" << *offset << "\">" << endl;
      }
      double temp = 0;
      float temp2 = 0;
      for (auto v : *field)
      {
        if (floatsize == "double")
        {
          temp = v;
          //fileout->write((char *)&temp, sizeof(double));
          data.write((char *)&temp, sizeof(temp));
          fieldsize += sizeof(temp);
        }
        else
        {
          temp2 = v;
          //fileout->write((char *)&temp, sizeof(double));
          data.write((char *)&temp2, sizeof(temp2));
          fieldsize += sizeof(temp2);
        }
      }
      *offset += 4 + fieldsize;
      appenddata->write((char *)&fieldsize, sizeof(int32_t));
      *appenddata << data.str();

      data.str("");
      data.clear();
      fieldsize = 0;
      *fileout << endl;
      *fileout << "</DataArray>" << endl;
    }
    *fileout << "</PointData>" << endl;
  }
  template <int D>
  void VTKOutput<D>::PrintAppended(stringstream *appended)
  {
    *fileout << "<AppendedData encoding=\"raw\">" << endl
             << "_";
    *fileout << appended->str();
    /*   for (auto field : value_field)
    {
      for (auto v : *field)
      {
        fileout->write((char *)&v, sizeof(double));
      }
    }*/
    *fileout << endl
             << "</AppendedData>" << endl;
  }
  template <int D>
  void VTKOutput<D>::PvdFile(string fname, int index)
  {
    ostringstream filenamefinal;
    stringstream contents;
    std::string fnamepart;
    fnamepart = fname.substr(fname.find_last_of("/\\")+1);
    filenamefinal << fname << ".pvd";

    contents
        << "<?xml version=\"1.0\"?>" << endl;
    contents << "<VTKFile type =\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
    contents << "<Collection>" << endl;
    auto comm = ma->GetCommunicator();
    for (int l = comm.Size() == 1 ? 0 : 1; l < comm.Size(); l++)
    { 
      contents << "<DataSet timestep=\"" << times[0] << "\"";
      if (comm.Size() > 1)
        contents << " part=\"" << l << "\"";
      contents << " file=\"" << fnamepart;
      if (comm.Size() > 1)
        contents << "_proc" << l; 
      contents << ".vtu\"/>" << endl;
      for (int k = 1; k < index; k++)
      { 
        contents << "<DataSet timestep=\"" << times[k] << "\"";
        if (comm.Size() > 1)
          contents << " part=\"" << l << "\"";
        contents << " file=\"" << fnamepart;
        if (comm.Size() > 1)
          contents << "_proc" << l; 
        contents << "_step" << setw(5) << setfill('0') << k << ".vtu\"/>" << endl;
      } 
    } 
    contents << "</Collection>" << endl;
    contents << "</VTKFile>";

    ofstream fileout;
    fileout.open(filenamefinal.str(), ofstream::trunc);
    fileout << contents.str();
    fileout.close();
  }
  template <int D>
  void VTKOutput<D>::Do(LocalHeap &lh, double time, VorB vb, const BitArray *drawelems)
  {
    ostringstream filenamefinal;
    stringstream appended;
    vector<int> datalength;
    int offs = 0;

    filenamefinal << filename;

    auto comm = ma->GetCommunicator();
    if (comm.Size() > 1)
      filenamefinal << "_proc" << comm.Rank();
    if (output_cnt > 0)
      filenamefinal << "_step" << setw(5) << setfill('0') << output_cnt;
    lastoutputname = filenamefinal.str();
    if (!legacy)
      filenamefinal << ".vtu";
    else
      filenamefinal << ".vtk";


    if ((comm.Size() == 1) || (comm.Rank() > 0) )
      cout << IM(4) << " Writing VTK-Output (" << lastoutputname << ")";
    if (output_cnt > 0)
    {
      //cout << IM(4) << " ( " << output_cnt << " )";
      if (time == -1)
      {
        times.push_back(output_cnt);
      }
      else
      {
        times.push_back(time);
      }
    }
    else
    {
      if (time != -1)
        times[0] = time;
    }
    output_cnt++;

    if (!legacy)
    { 
      if ((comm.Size()==1 && output_cnt > 1) || (comm.Size() > 1 && comm.Rank()==0))
        PvdFile(filename, output_cnt);
    } 

    if ((comm.Size() == 1) || (comm.Rank() > 0) )
      cout << IM(4) << ":" << flush;
    else
      return;

    if(legacy)
      fileout = make_shared<ofstream>(filenamefinal.str());
    else
      fileout = make_shared<ofstream>(filenamefinal.str(), ofstream::binary);


    ResetArrays();

    std::map<ELEMENT_TYPE, Array<IntegrationPoint>> ref_vertices;
    std::map<ELEMENT_TYPE, Array<VTKCell>> ref_elems;
    // Array<IntegrationPoint> ref_vertices_tet(0), ref_vertices_prism(0), ref_vertices_pyramid(0), ref_vertices_trig(0), ref_vertices_quad(0), ref_vertices_hex(0);
    // Array<VTKCell> ref_tets(0), ref_prisms(0), ref_trigs(0), ref_quads(0), ref_hexes(0), ref_pyramids(0);
    // FlatArray<IntegrationPoint> ref_vertices;
    // FlatArray<VTKCell> ref_elems;
    /*
    if (D==3)
      FillReferenceData3D(ref_vertices,ref_tets);
    else
      FillReferenceData2D(ref_vertices,ref_tets);
    */
    FillReferenceTet(ref_vertices[ET_TET], ref_elems[ET_TET]);
    FillReferencePrism(ref_vertices[ET_PRISM], ref_elems[ET_PRISM]);
    FillReferencePyramid(ref_vertices[ET_PYRAMID], ref_elems[ET_PYRAMID]);
    FillReferenceQuad(ref_vertices[ET_QUAD], ref_elems[ET_QUAD]);
    FillReferenceTrig(ref_vertices[ET_TRIG], ref_elems[ET_TRIG]);
    FillReferenceHex(ref_vertices[ET_HEX], ref_elems[ET_HEX]);

    // header:
    if (!legacy)
    {
      *fileout << "<?xml version=\"1.0\"?>" << endl;

      *fileout << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
      *fileout << "<UnstructuredGrid>" << endl;
    }
    else
    {
      *fileout << "# vtk DataFile Version 3.0" << endl;
      *fileout << "vtk output" << endl;
      *fileout << "ASCII" << endl;
      *fileout << "DATASET UNSTRUCTURED_GRID" << endl;
    }
    int ne = ma->GetNE(vb);

    IntRange range = only_element >= 0 ? IntRange(only_element, only_element + 1) : IntRange(ne);

    for (int elnr : range)
    {
      if (drawelems && !(drawelems->Test(elnr)))
        continue;

      HeapReset hr(lh);

      ElementId ei(vb, elnr);
      ElementTransformation &eltrans = ma->GetTrafo(ei, lh);
      ELEMENT_TYPE eltype = ma->GetElType(ei);

      int offset = points.Size();
      for (auto ip : ref_vertices[eltype])
      {
        if (vb == VOL)
        {
          MappedIntegrationPoint<D, D> mip(ip, eltrans);
          points.Append(mip.GetPoint());
        }
        else
        {
          MappedIntegrationPoint<D - 1, D> mip(ip, eltrans);
          points.Append(mip.GetPoint());
        }
      }

      for (int i = 0; i < coefs.Size(); i++)
      {
        for (auto ip : ref_vertices[eltype])
        {
          if (vb == VOL)
          {
            MappedIntegrationPoint<D, D> mip(ip, eltrans);
            const int dim = coefs[i]->Dimension();
            FlatVector<> tmp(dim, lh);
            coefs[i]->Evaluate(mip, tmp);
            for (int d = 0; d < dim; ++d)
              value_field[i]->Append(tmp(d));
          }
          else
          {
            MappedIntegrationPoint<D - 1, D> mip(ip, eltrans);
            const int dim = coefs[i]->Dimension();
            FlatVector<> tmp(dim, lh);
            coefs[i]->Evaluate(mip, tmp);
            for (int d = 0; d < dim; ++d)
              value_field[i]->Append(tmp(d));
          }
        }
      }

      for (auto new_elem : ref_elems[eltype])
      {
        for (auto & pi : new_elem.pi)
          pi += offset;
        cells.Append(new_elem);
      }
    }
    if (!legacy)
    {
      *fileout << "<Piece NumberOfPoints=\"" << points.Size() << "\" NumberOfCells=\"" << cells.Size() << "\">" << endl;
      PrintPoints(&offs, &appended);
      *fileout << "<Cells>" << endl;
      PrintCells(&offs, &appended);
      PrintCellTypes(vb, &offs, &appended, drawelems);
      *fileout << "</Cells>" << endl;
      PrintFieldData(&offs, &appended);

      // Footer:
      *fileout << "</Piece>" << endl;
      *fileout << "</UnstructuredGrid>" << endl;
      PrintAppended(&appended);
      *fileout << "</VTKFile>" << endl;
    }
    else
    {
      PrintPointsLegacy();
      PrintCellsLegacy();
      PrintCellTypesLegacy(vb, drawelems);
      PrintFieldDataLegacy();
    }
    cout << IM(4) << " Done." << endl;
    //cout << IM(4) << " [ VTKOutput Counter: " << output_cnt << ""]" << endl;
  }

  template class VTKOutput<2>;
  template class VTKOutput<3>;
}
