/*********************************************************************/
/* File:   vtkoutput.cpp                                             */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   1. June 2014                                              */
/*********************************************************************/

#include <comp.hpp>

namespace ngcomp
{ 

  ValueField::ValueField(int adim, string aname) : Array<double>(),  dim(adim), name(aname){;}

  template <int D> 
  VTKOutput<D>::VTKOutput (const Array<shared_ptr<CoefficientFunction>> & a_coefs,
                           const Flags & flags,
                           shared_ptr<MeshAccess> ama )
    : VTKOutput(ama, a_coefs, 
                flags.GetStringListFlag ("fieldnames" ),
                flags.GetStringFlag ("filename","output"),
                (int) flags.GetNumFlag ( "subdivision", 0),
                (int) flags.GetNumFlag ( "only_element", -1))
  {;}


  template <int D> 
  VTKOutput<D>::VTKOutput (shared_ptr<MeshAccess> ama,
                           const Array<shared_ptr<CoefficientFunction>> & a_coefs,
                           const Array<string> & a_field_names,
                           string a_filename, int a_subdivision, int a_only_element)
    : ma(ama), coefs(a_coefs), fieldnames(a_field_names),
      filename(a_filename), subdivision(a_subdivision), only_element(a_only_element)
  {
    value_field.SetSize(a_coefs.Size());
    for (int i = 0; i < a_coefs.Size(); i++)
      if (fieldnames.Size() > i)
        value_field[i] = make_shared<ValueField>(coefs[i]->Dimension(),fieldnames[i]);
      else 
        value_field[i] = make_shared<ValueField>(coefs[i]->Dimension(),"dummy" + to_string(i));
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
  template<>
  void VTKOutput<2>::FillReferenceData(Array<IntegrationPoint> & ref_coords, Array<INT<2+1>> & ref_trigs)
  {
    enum { D = 2 };
    if (subdivision == 0)
    {
      ref_coords.Append(IntegrationPoint(0.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(1.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(0.0,1.0,0.0));
      ref_trigs.Append(INT<D+1>(0,1,2));
    }
    else
    {
      const int r = 1<<subdivision;
      const int s = r + 1;

      const double h = 1.0/r;

      int pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i+j <= r; ++j)
          {
            ref_coords.Append(IntegrationPoint(j*h,i*h));
          }

      pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i+j <= r; ++j, pidx++)
          {
            // int pidx_curr = pidx;
            if (i+j == r) continue;
            int pidx_incr_i = pidx+1;
            int pidx_incr_j = pidx+s-i;

            ref_trigs.Append(INT<3>(pidx,pidx_incr_i,pidx_incr_j));

            int pidx_incr_ij = pidx_incr_j + 1;

            if(i+j+1<r) 
              ref_trigs.Append(INT<3>(pidx_incr_i,pidx_incr_ij,pidx_incr_j));
          }              
    }
  }

  /// Fill principil lattices (points and connections on subdivided reference simplex) in 3D
  template <int D> 
  void VTKOutput<D>::FillReferenceData(Array<IntegrationPoint> & ref_coords, Array<INT<D+1>> & ref_tets)
  {
    if (subdivision == 0)
    {
      ref_coords.Append(IntegrationPoint(0.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(1.0,0.0,0.0));
      ref_coords.Append(IntegrationPoint(0.0,1.0,0.0));
      ref_coords.Append(IntegrationPoint(0.0,0.0,1.0));
      ref_tets.Append(INT<D+1>(0,1,2,3));
    }
    else
    {
      const int r = 1<<subdivision;
      const int s = r + 1;

      const double h = 1.0/r;

      int pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i+j <= r; ++j)
          for (int k = 0; i+j+k <= r; ++k)
          {
            ref_coords.Append(IntegrationPoint(i*h,j*h,k*h));
          }

      pidx = 0;
      for (int i = 0; i <= r; ++i)
        for (int j = 0; i+j <= r; ++j)
          for (int k = 0; i+j+k <= r; ++k, pidx++)
          {
            if (i+j+k == r) continue;
            // int pidx_curr = pidx;
            int pidx_incr_k = pidx+1;
            int pidx_incr_j = pidx+s-i-j;
            int pidx_incr_i = pidx+(s-i)*(s+1-i)/2-j;

            int pidx_incr_kj = pidx_incr_j + 1;

            int pidx_incr_ij = pidx+(s-i)*(s+1-i)/2-j + s-(i+1)-j;
            int pidx_incr_ki = pidx+(s-i)*(s+1-i)/2-j + 1;
            int pidx_incr_kij = pidx+(s-i)*(s+1-i)/2-j + s-(i+1)-j + 1;

            ref_tets.Append(INT<4>(pidx,pidx_incr_k,pidx_incr_j,pidx_incr_i));
            if (i+j+k+1 == r)
              continue;

            ref_tets.Append(INT<4>(pidx_incr_k,pidx_incr_kj,pidx_incr_j,pidx_incr_i));
            ref_tets.Append(INT<4>(pidx_incr_k,pidx_incr_kj,pidx_incr_ki,pidx_incr_i));

            ref_tets.Append(INT<4>(pidx_incr_j,pidx_incr_i,pidx_incr_kj,pidx_incr_ij));
            ref_tets.Append(INT<4>(pidx_incr_i,pidx_incr_kj,pidx_incr_ij,pidx_incr_ki));
           
            if (i+j+k+2 != r)
              ref_tets.Append(INT<4>(pidx_incr_kj,pidx_incr_ij,pidx_incr_ki,pidx_incr_kij));
          }              
    }
  }

  /// output of data points
  template <int D> 
  void VTKOutput<D>::PrintPoints()
  {
    *fileout << "POINTS " << points.Size() << " float" << endl;
    for (auto p : points)
    {
      *fileout << p;
      if (D==2)
        *fileout << "\t 0.0";
      *fileout << endl;
    }
  }

  /// output of cells in form vertices
  template <int D> 
  void VTKOutput<D>::PrintCells()
  {
    *fileout << "CELLS " << cells.Size() << " " << (D+2) * cells.Size() << endl;
    for (auto c : cells)
      *fileout << D+1 <<" " << c << endl;
  }

  /// output of cell types (here only simplices)
  template <int D> 
  void VTKOutput<D>::PrintCellTypes()
  {
    *fileout << "CELL_TYPES " << cells.Size() << endl;
    if (D==3)
      for (auto c : cells)
        { *fileout << "10 " << endl; (void)c; } // no warning 
    else
      for (auto c : cells)
        { *fileout << "5 " << endl; (void)c; } // no warning 
    *fileout << "CELL_DATA " << cells.Size() << endl;
    *fileout << "POINT_DATA " << points.Size() << endl;
  }

  /// output of field data (coefficient values)
  template <int D> 
  void VTKOutput<D>::PrintFieldData()
  {
    *fileout << "FIELD FieldData " << value_field.Size() << endl;

    for (auto field : value_field)
    {
      *fileout << field->Name() << " "
               << field->Dimension() << " "
               << field->Size()/field->Dimension() << " float" << endl;
      for (auto v : *field)
        *fileout << v << " ";
      *fileout << endl;
    }
    
  }
    

  template <int D> 
  void VTKOutput<D>::Do (LocalHeap & lh, const BitArray * drawelems)
  {
    ostringstream filenamefinal;
    filenamefinal << filename;
    if (output_cnt > 0)
      filenamefinal << "_" << output_cnt;
    filenamefinal << ".vtk";
    fileout = make_shared<ofstream>(filenamefinal.str());
    cout << " Writing VTK-Output";
    if (output_cnt > 0)
      cout << " ( " << output_cnt << " )";
    cout << ":" << flush;
    
    output_cnt++;

    ResetArrays();

    Array<IntegrationPoint> ref_vertices(0);
    Array<INT<D+1>> ref_tets(0);
    /*
    if (D==3)
      FillReferenceData3D(ref_vertices,ref_tets);
    else
      FillReferenceData2D(ref_vertices,ref_tets);
    */
    FillReferenceData(ref_vertices,ref_tets);
      
    // header:
    *fileout << "# vtk DataFile Version 3.0" << endl;
    *fileout << "vtk output" << endl;
    *fileout << "ASCII" << endl;
    *fileout << "DATASET UNSTRUCTURED_GRID" << endl;

    int ne = ma->GetNE();

    IntRange range = only_element >= 0 ? IntRange(only_element,only_element+1) : IntRange(ne);
    
    for ( int elnr : range)
    {
      if (drawelems && !(drawelems->Test(elnr)))
          continue;
      
      HeapReset hr(lh);

      ElementId ei(VOL, elnr);
      ElementTransformation & eltrans = ma->GetTrafo (ei, lh);

      int offset = points.Size();
      for ( auto ip : ref_vertices)
      {
        MappedIntegrationPoint<D,D> mip(ip, eltrans);
        points.Append(mip.GetPoint());
      }
      
      for (int i = 0; i < coefs.Size(); i++)
      {
        for ( auto ip : ref_vertices)
        {
          MappedIntegrationPoint<D,D> mip(ip, eltrans);
          const int dim = coefs[i]->Dimension();
          FlatVector<> tmp(dim,lh);
          coefs[i]->Evaluate(mip,tmp);
          for (int d = 0; d < dim; ++d)
            value_field[i]->Append(tmp(d));
        }
      }
      
      for ( auto tet : ref_tets)
      {
        INT<D+1> new_tet = tet;
        for (int i = 0; i < D+1; ++i)
          new_tet[i] += offset;
        cells.Append(new_tet);
      }

    }

    PrintPoints();
    PrintCells();
    PrintCellTypes();
    PrintFieldData();
      
    cout << " Done." << endl;
  }    

  NumProcVTKOutput::NumProcVTKOutput (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    const Array<string> & coefs_strings = flags.GetStringListFlag ("coefficients");
    
    Array<shared_ptr<CoefficientFunction>> coefs;
    for (int i = 0; i < coefs_strings.Size(); ++i)
      coefs.Append(apde->GetCoefficientFunction (coefs_strings[i]));

    if (apde->GetMeshAccess()->GetDimension() == 2)
      vtkout = make_shared<VTKOutput<2>>(coefs, flags, apde->GetMeshAccess());
    else 
      vtkout = make_shared<VTKOutput<3>>(coefs, flags, apde->GetMeshAccess());
  }


  void NumProcVTKOutput::Do (LocalHeap & lh)
  {
    vtkout->Do(lh);
  }

  
  static RegisterNumProc<NumProcVTKOutput> npvtkout("vtkoutput");

  template class VTKOutput<2>;
  template class VTKOutput<3>;
}

