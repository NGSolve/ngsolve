/*********************************************************************/
/* File:   vtkoutput.cpp                                             */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   1. June 2014                                              */
/*********************************************************************/

#include <comp.hpp>

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
                  flags.GetDefineFlag("legacy"))
  {
    ;
  }

  template <int D>
  VTKOutput<D>::VTKOutput(shared_ptr<MeshAccess> ama,
                          const Array<shared_ptr<CoefficientFunction>> &a_coefs,
                          const Array<string> &a_field_names,
                          string a_filename, int a_subdivision, int a_only_element, 
                          string a_floatsize, bool a_legacy)
      : ma(ama), coefs(a_coefs), fieldnames(a_field_names),
        filename(a_filename), subdivision(a_subdivision), only_element(a_only_element), floatsize(a_floatsize), legacy(a_legacy)
  {
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

  /// Fill principil lattices (points and connections on subdivided reference simplex) in 2D
  template <int D>
  void VTKOutput<D>::FillReferenceTrig(Array<IntegrationPoint> &ref_coords, Array<INT<ELEMENT_MAXPOINTS + 1>> &ref_elems)
  {
    if (subdivision == 0)
    {
      ref_coords.Append(IntegrationPoint(0.0, 0.0, 0.0));
      ref_coords.Append(IntegrationPoint(1.0, 0.0, 0.0));
      ref_coords.Append(IntegrationPoint(0.0, 1.0, 0.0));
      ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(3, 0, 1, 2));
    }
    else
    {
      const int r = 1 << subdivision;
      const int s = r + 1;

      const double h = 1.0 / r;

      int pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i + j <= r; ++j)
        {
          ref_coords.Append(IntegrationPoint(j * h, i * h));
        }

      pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i + j <= r; ++j, pidx++)
        {
          // int pidx_curr = pidx;
          if (i + j == r)
            continue;
          int pidx_incr_i = pidx + 1;
          int pidx_incr_j = pidx + s - i;

          ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(3, pidx, pidx_incr_i, pidx_incr_j));

          int pidx_incr_ij = pidx_incr_j + 1;

          if (i + j + 1 < r)
            ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(3, pidx_incr_i, pidx_incr_ij, pidx_incr_j));
        }
    }
  }
  /// Fill principil lattices (points and connections on subdivided reference simplex) in 2D
  template <int D>
  void VTKOutput<D>::FillReferenceQuad(Array<IntegrationPoint> &ref_coords, Array<INT<ELEMENT_MAXPOINTS + 1>> &ref_elems)
  {
    if (subdivision == 0)
    {
      ref_coords.Append(IntegrationPoint(0.0, 0.0, 0.0));
      ref_coords.Append(IntegrationPoint(1.0, 0.0, 0.0));
      ref_coords.Append(IntegrationPoint(1.0, 1.0, 0.0));
      ref_coords.Append(IntegrationPoint(0.0, 1.0, 0.0));
      INT<ELEMENT_MAXPOINTS + 1> quad;
      quad[0] = 4;
      quad[1] = 0;
      quad[2] = 1;
      quad[3] = 2;
      quad[4] = 3;
      ref_elems.Append(quad);
    }
    else
    {
      const int r = 1 << subdivision;
      // const int s = r + 1;

      const double h = 1.0 / r;

      int pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; j <= r; ++j)
        {
          ref_coords.Append(IntegrationPoint(j * h, i * h));
        }

      for (int i = 0; i < r; ++i)
      {
        int incr_i = r + 1;
        pidx = i * incr_i;
        for (int j = 0; j < r; ++j, pidx++)
        {
          ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(4, pidx, pidx + 1, pidx + incr_i + 1, pidx + incr_i));
        }
      }
    }
  }

  /// Fill principil lattices (points and connections on subdivided reference simplex) in 3D
  template <int D>
  void VTKOutput<D>::FillReferenceTet(Array<IntegrationPoint> &ref_coords, Array<INT<ELEMENT_MAXPOINTS + 1>> &ref_elems)
  {
    if (subdivision == 0)
    {
      ref_coords.Append(IntegrationPoint(0.0, 0.0, 0.0));
      ref_coords.Append(IntegrationPoint(1.0, 0.0, 0.0));
      ref_coords.Append(IntegrationPoint(0.0, 1.0, 0.0));
      ref_coords.Append(IntegrationPoint(0.0, 0.0, 1.0));
      ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(4, 0, 1, 2, 3));
    }
    else
    {
      const int r = 1 << subdivision;
      const int s = r + 1;

      const double h = 1.0 / r;

      int pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i + j <= r; ++j)
          for (int k = 0; i + j + k <= r; ++k)
          {
            ref_coords.Append(IntegrationPoint(i * h, j * h, k * h));
          }

      for (int i = 0; i <= r; ++i)
        for (int j = 0; i + j <= r; ++j)
          for (int k = 0; i + j + k <= r; ++k, pidx++)
          {
            if (i + j + k == r)
              continue;
            // int pidx_curr = pidx;
            int pidx_incr_k = pidx + 1;
            int pidx_incr_j = pidx + s - i - j;
            int pidx_incr_i = pidx + (s - i) * (s + 1 - i) / 2 - j;

            int pidx_incr_kj = pidx_incr_j + 1;

            int pidx_incr_ij = pidx + (s - i) * (s + 1 - i) / 2 - j + s - (i + 1) - j;
            int pidx_incr_ki = pidx + (s - i) * (s + 1 - i) / 2 - j + 1;
            int pidx_incr_kij = pidx + (s - i) * (s + 1 - i) / 2 - j + s - (i + 1) - j + 1;

            ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(4, pidx, pidx_incr_k, pidx_incr_j, pidx_incr_i));
            if (i + j + k + 1 == r)
              continue;

            ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(4, pidx_incr_k, pidx_incr_kj, pidx_incr_j, pidx_incr_i));
            ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(4, pidx_incr_k, pidx_incr_kj, pidx_incr_ki, pidx_incr_i));

            ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(4, pidx_incr_j, pidx_incr_i, pidx_incr_kj, pidx_incr_ij));
            ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(4, pidx_incr_i, pidx_incr_kj, pidx_incr_ij, pidx_incr_ki));

            if (i + j + k + 2 != r)
              ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(4, pidx_incr_kj, pidx_incr_ij, pidx_incr_ki, pidx_incr_kij));
          }
    }
  }

  /// Fill principil lattices (points and connections on subdivided reference hexahedron) in 3D
  template <int D>
  void VTKOutput<D>::FillReferenceHex(Array<IntegrationPoint> &ref_coords, Array<INT<ELEMENT_MAXPOINTS + 1>> &ref_elems)
  {
    if (subdivision == 0)
    {
      ref_coords.Append(IntegrationPoint(0.0, 0.0, 0.0));
      ref_coords.Append(IntegrationPoint(1.0, 0.0, 0.0));
      ref_coords.Append(IntegrationPoint(1.0, 1.0, 0.0));
      ref_coords.Append(IntegrationPoint(0.0, 1.0, 0.0));
      ref_coords.Append(IntegrationPoint(0.0, 0.0, 1.0));
      ref_coords.Append(IntegrationPoint(1.0, 0.0, 1.0));
      ref_coords.Append(IntegrationPoint(1.0, 1.0, 1.0));
      ref_coords.Append(IntegrationPoint(0.0, 1.0, 1.0));
      INT<ELEMENT_MAXPOINTS + 1> hex;
      hex[0] = 8;
      hex[1] = 0;
      hex[2] = 1;
      hex[3] = 2;
      hex[4] = 3;
      hex[5] = 4;
      hex[6] = 5;
      hex[7] = 6;
      hex[8] = 7;
      ref_elems.Append(hex);
    }
    else
    {
      const int r = 1 << subdivision;
      // const int s = r + 1;

      const double h = 1.0 / r;

      int pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; j <= r; ++j)
          for (int k = 0; k <= r; ++k)
          {
            ref_coords.Append(IntegrationPoint(k * h, j * h, i * h));
          }

      for (int i = 0; i < r; ++i)
      {
        int incr_i = (r + 1) * (r + 1);
        for (int j = 0; j < r; ++j)
        {
          int incr_j = r + 1;
          pidx = i * incr_i + j * incr_j;
          for (int k = 0; k < r; ++k, pidx++)
          {
            ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(8, pidx, pidx + 1, pidx + incr_j + 1, pidx + incr_j,
                                                        pidx + incr_i, pidx + incr_i + 1, pidx + incr_i + incr_j + 1, pidx + incr_j + incr_i));
          }
        }
      }
    }
  }

  template <int D>
  void VTKOutput<D>::FillReferencePrism(Array<IntegrationPoint> &ref_coords, Array<INT<ELEMENT_MAXPOINTS + 1>> &ref_elems)
  {
    if (subdivision == 0)
    {
      ref_coords.Append(IntegrationPoint(0.0, 0.0, 0.0));
      ref_coords.Append(IntegrationPoint(1.0, 0.0, 0.0));
      ref_coords.Append(IntegrationPoint(0.0, 1.0, 0.0));
      ref_coords.Append(IntegrationPoint(0.0, 0.0, 1.0));
      ref_coords.Append(IntegrationPoint(1.0, 0.0, 1.0));
      ref_coords.Append(IntegrationPoint(0.0, 1.0, 1.0));
      INT<ELEMENT_MAXPOINTS + 1> elem;
      elem[0] = 6;
      for (int i = 0; i < ElementTopology::GetNVertices(ET_PRISM); i++)
        elem[i + 1] = i;
      ref_elems.Append(elem);
    }
    else
    {
      const int r = 1 << subdivision;
      const int s = r + 1;

      const double h = 1.0 / r;

      int pidx = 0;
      for (int k = 0; k <= r; k++)
        for (int i = 0; i <= r; ++i)
          for (int j = 0; i + j <= r; ++j)
          {
            ref_coords.Append(IntegrationPoint(j * h, i * h, k * h));
          }

      pidx = 0;
      for (int k = 0; k < r; k++)
      {
        int incr_k = (r + 2) * (r + 1) / 2;
        pidx = k * incr_k;
        for (int i = 0; i <= r; ++i)
          for (int j = 0; i + j <= r; ++j, pidx++)
          {
            // int pidx_curr = pidx;
            if (i + j == r)
              continue;
            int pidx_incr_i = pidx + 1;
            int pidx_incr_j = pidx + s - i;

            ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(6, pidx, pidx_incr_i, pidx_incr_j, pidx + incr_k, pidx_incr_i + incr_k, pidx_incr_j + incr_k, 0, 0));

            int pidx_incr_ij = pidx_incr_j + 1;

            if (i + j + 1 < r)
              ref_elems.Append(INT<ELEMENT_MAXPOINTS + 1>(6, pidx_incr_i, pidx_incr_ij, pidx_incr_j, pidx_incr_i + incr_k, pidx_incr_ij + incr_k, pidx_incr_j + incr_k, 0, 0));
          }
      }
    }
  }
  /* ###########################
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
      ndata += c[0];
    }
    *fileout << "CELLS " << cells.Size() << " " << ndata << endl;
    for (auto c : cells)
    {
      int nv = c[0];
      *fileout << nv << "\t";
      for (int i = 0; i < nv; i++)
        *fileout << c[i + 1] << "\t";
      *fileout << endl;
    }
  }

  /// output of cell types (here only simplices)
  template <int D>
  void VTKOutput<D>::PrintCellTypesLegacy(VorB vb, const BitArray *drawelems)
  {
    *fileout << "CELL_TYPES " << cells.Size() << endl;
    int factor = (1 << subdivision) * (1 << subdivision);
    if (D == 3 && vb == VOL)
      factor *= (1 << subdivision);
    for (auto e : ma->Elements(vb))
    {
      if (drawelems && !(drawelems->Test(e.Nr())))
        continue;

      switch (ma->GetElType(e))
      {
      case ET_TET:
        for (int i = 0; i < factor; i++)
          *fileout << "10 " << endl; //(void)c;
        break;
      case ET_QUAD:
        for (int i = 0; i < factor; i++)
          *fileout << "9 " << endl;
        break;
      case ET_TRIG:
        for (int i = 0; i < factor; i++)
          *fileout << "5 " << endl;
        break;
      case ET_PRISM:
        for (int i = 0; i < factor; i++)
          *fileout << "13 " << endl;
        break;
      case ET_HEX:
        for (int i = 0; i < factor; i++)
          *fileout << "12 " << endl;
        break;
      default:
        cout << "VTKOutput Element Type " << ma->GetElType(e) << " not supported!" << endl;
      }
    }
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
    int32_t ndata = 0;
    int32_t offs = 0;
    for (auto c : cells)
    {

      int nv = c[0];
      offs += nv;
      offsets.write((char *)&offs, sizeof(int32_t));
      sizeoff += sizeof(int32_t);
      for (int i = 0; i < nv; i++)
      {
        connectivity.write((char *)&c[i + 1], sizeof(int));
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
    int factor = (1 << subdivision) * (1 << subdivision);
    stringstream data;
    uint32_t sizetypes = 0;
    uint8_t eltype;

    if (D == 3 && vb == VOL)
      factor *= (1 << subdivision);
    for (auto e : ma->Elements(vb))
    {
      if (drawelems && !(drawelems->Test(e.Nr())))
        continue;

      switch (ma->GetElType(e))
      {
      case ET_TET:
        for (int i = 0; i < factor; i++)
        {
          sizetypes += sizeof(uint8_t);
          eltype = 10;
          data.write((char *)&eltype, sizeof(uint8_t)); //(void)c;
        }
        break;
      case ET_QUAD:
        for (int i = 0; i < factor; i++)
        {
          sizetypes += sizeof(uint8_t);
          eltype = 9;
          data.write((char *)&eltype, sizeof(uint8_t));
        }
        break;
      case ET_TRIG:
        for (int i = 0; i < factor; i++)
        {
          sizetypes += sizeof(uint8_t);
          eltype = 5;
          data.write((char *)&eltype, sizeof(uint8_t));
        }
        break;
      case ET_PRISM:
        for (int i = 0; i < factor; i++)
        {
          sizetypes += sizeof(uint8_t);
          eltype = 13;
          data.write((char *)&eltype, sizeof(uint8_t));
        }
        break;
      case ET_HEX:
        for (int i = 0; i < factor; i++)
        {
          sizetypes += sizeof(uint8_t);
          eltype = 12;
          data.write((char *)&eltype, sizeof(uint8_t));
        }
        break;
      default:
        cout << "VTKOutput Element Type " << ma->GetElType(e) << " not supported!" << endl;
      }
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
      contents << " file=\"" << fname;
      if (comm.Size() > 1)
        contents << "_proc" << l; 
      contents << ".vtu\"/>" << endl;
      for (int k = 1; k < index; k++)
      { 
        contents << "<DataSet timestep=\"" << times[k] << "\"";
        if (comm.Size() > 1)
          contents << " part=\"" << l << "\"";
        contents << " file=\"" << fname;
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
      if ((comm.Size()==1 and output_cnt > 1) || (comm.Size() > 1 && comm.Rank()==0))
        PvdFile(filename, output_cnt);
    } 

    if ((comm.Size() == 1) || (comm.Rank() > 0) )
      cout << IM(4) << ":" << flush;
    else
      return;

    fileout = make_shared<ofstream>(filenamefinal.str());


    ResetArrays();

    Array<IntegrationPoint> ref_vertices_tet(0), ref_vertices_prism(0), ref_vertices_trig(0), ref_vertices_quad(0), ref_vertices_hex(0);
    Array<INT<ELEMENT_MAXPOINTS + 1>> ref_tets(0), ref_prisms(0), ref_trigs(0), ref_quads(0), ref_hexes(0);
    FlatArray<IntegrationPoint> ref_vertices;
    FlatArray<INT<ELEMENT_MAXPOINTS + 1>> ref_elems;
    /*
    if (D==3)
      FillReferenceData3D(ref_vertices,ref_tets);
    else
      FillReferenceData2D(ref_vertices,ref_tets);
    */
    FillReferenceTet(ref_vertices_tet, ref_tets);
    FillReferencePrism(ref_vertices_prism, ref_prisms);
    FillReferenceQuad(ref_vertices_quad, ref_quads);
    FillReferenceTrig(ref_vertices_trig, ref_trigs);
    FillReferenceHex(ref_vertices_hex, ref_hexes);

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

      switch (eltype)
      {
      case ET_TRIG:
        ref_vertices.Assign(ref_vertices_trig);
        ref_elems.Assign(ref_trigs);
        break;
      case ET_QUAD:
        ref_vertices.Assign(ref_vertices_quad);
        ref_elems.Assign(ref_quads);
        break;
      case ET_TET:
        ref_vertices.Assign(ref_vertices_tet);
        ref_elems.Assign(ref_tets);
        break;
      case ET_HEX:
        ref_vertices.Assign(ref_vertices_hex);
        ref_elems.Assign(ref_hexes);
        break;
      case ET_PRISM:
        ref_vertices.Assign(ref_vertices_prism);
        ref_elems.Assign(ref_prisms);
        break;
      default:
        throw Exception("VTK output for element-type" + ToString(eltype) + "not supported");
      }

      int offset = points.Size();
      for (auto ip : ref_vertices)
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
        for (auto ip : ref_vertices)
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

      for (auto elem : ref_elems)
      {
        INT<ELEMENT_MAXPOINTS + 1> new_elem = elem;
        for (int i = 1; i <= new_elem[0]; ++i)
          new_elem[i] += offset;
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

  NumProcVTKOutput::NumProcVTKOutput(shared_ptr<PDE> apde, const Flags &flags)
      : NumProc(apde)
  {
    const Array<string> &coefs_strings = flags.GetStringListFlag("coefficients");

    Array<shared_ptr<CoefficientFunction>> coefs;
    for (int i = 0; i < coefs_strings.Size(); ++i)
      coefs.Append(apde->GetCoefficientFunction(coefs_strings[i]));

    if (apde->GetMeshAccess()->GetDimension() == 2)
      vtkout = make_shared<VTKOutput<2>>(coefs, flags, apde->GetMeshAccess());
    else
      vtkout = make_shared<VTKOutput<3>>(coefs, flags, apde->GetMeshAccess());
  }

  void NumProcVTKOutput::Do(LocalHeap &lh)
  {
    vtkout->Do(lh);
  }

  static RegisterNumProc<NumProcVTKOutput> npvtkout("vtkoutput");

  template class VTKOutput<2>;
  template class VTKOutput<3>;
}
