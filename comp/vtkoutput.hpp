#pragma once

/*********************************************************************/
/* File:   vtkoutput.hpp                                             */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   1. June 2014                                              */
/*********************************************************************/

namespace ngcomp
{

  class ValueField : public Array<double>
  {
    int dim = 1;
    string name = "none";

  public:
    ValueField() { ; };
    ValueField(int adim, string aname);
    void SetDimension(int adim) { dim = adim; }
    int Dimension() { return dim; }
    void SetName(string aname) { name = aname; }
    string Name() { return name; }
  };

  class BaseVTKOutput
  {
  public:
    virtual ~BaseVTKOutput() { ; }
    virtual void Do(LocalHeap &lh, double time = -1, VorB vb = VOL, const BitArray *drawelems = 0) = 0;
    string lastoutputname = "";
  };

  struct VTKCell {
    int type = 0;
    ArrayMem<int, ELEMENT_MAXPOINTS> pi;

    VTKCell() = default;

    VTKCell & operator=(const VTKCell&other) {
      type = other.type;
      pi.SetSize(other.pi.Size());
      for(auto i : Range(other.pi))
        pi[i] = other.pi[i];
      return *this;
    }

    VTKCell(ELEMENT_TYPE et, int order, const std::map<tuple<int,int,int>, int> &m,
        int i, int j, int k, Vec<3,int> vi={1,0,0}, Vec<3,int> vj={0,1,0}, Vec<3,int> vk={0,0,1});

    static int GetVtkType(ELEMENT_TYPE et, int order=1);

    enum VTK_CELL_TYPE
    {
      // Linear cells
      VTK_EMPTY_CELL = 0,
      VTK_VERTEX = 1,
      VTK_POLY_VERTEX = 2,
      VTK_LINE = 3,
      VTK_POLY_LINE = 4,
      VTK_TRIANGLE = 5,
      VTK_TRIANGLE_STRIP = 6,
      VTK_POLYGON = 7,
      VTK_PIXEL = 8,
      VTK_QUAD = 9,
      VTK_TETRA = 10,
      VTK_VOXEL = 11,
      VTK_HEXAHEDRON = 12,
      VTK_WEDGE = 13,
      VTK_PYRAMID = 14,
      VTK_PENTAGONAL_PRISM = 15,
      VTK_HEXAGONAL_PRISM = 16,

      // Quadratic, isoparametric cells
      VTK_QUADRATIC_EDGE = 21,
      VTK_QUADRATIC_TRIANGLE = 22,
      VTK_QUADRATIC_QUAD = 23,
      VTK_QUADRATIC_POLYGON = 36,
      VTK_QUADRATIC_TETRA = 24,
      VTK_QUADRATIC_HEXAHEDRON = 25,
      VTK_QUADRATIC_WEDGE = 26,
      VTK_QUADRATIC_PYRAMID = 27,
      VTK_BIQUADRATIC_QUAD = 28,
      VTK_TRIQUADRATIC_HEXAHEDRON = 29,
      VTK_TRIQUADRATIC_PYRAMID = 37,
      VTK_QUADRATIC_LINEAR_QUAD = 30,
      VTK_QUADRATIC_LINEAR_WEDGE = 31,
      VTK_BIQUADRATIC_QUADRATIC_WEDGE = 32,
      VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33,
      VTK_BIQUADRATIC_TRIANGLE = 34,

      // Cubic, isoparametric cell
      VTK_CUBIC_LINE = 35,

      // Special class of cells formed by convex group of points
      VTK_CONVEX_POINT_SET = 41,

      // Polyhedron cell (consisting of polygonal faces)
      VTK_POLYHEDRON = 42,

      // Higher order cells in parametric form
      VTK_PARAMETRIC_CURVE = 51,
      VTK_PARAMETRIC_SURFACE = 52,
      VTK_PARAMETRIC_TRI_SURFACE = 53,
      VTK_PARAMETRIC_QUAD_SURFACE = 54,
      VTK_PARAMETRIC_TETRA_REGION = 55,
      VTK_PARAMETRIC_HEX_REGION = 56,

      // Higher order cells
      VTK_HIGHER_ORDER_EDGE = 60,
      VTK_HIGHER_ORDER_TRIANGLE = 61,
      VTK_HIGHER_ORDER_QUAD = 62,
      VTK_HIGHER_ORDER_POLYGON = 63,
      VTK_HIGHER_ORDER_TETRAHEDRON = 64,
      VTK_HIGHER_ORDER_WEDGE = 65,
      VTK_HIGHER_ORDER_PYRAMID = 66,
      VTK_HIGHER_ORDER_HEXAHEDRON = 67,

      // Arbitrary order Lagrange elements (formulated separated from generic higher order cells)
      VTK_LAGRANGE_CURVE = 68,
      VTK_LAGRANGE_TRIANGLE = 69,
      VTK_LAGRANGE_QUADRILATERAL = 70,
      VTK_LAGRANGE_TETRAHEDRON = 71,
      VTK_LAGRANGE_HEXAHEDRON = 72,
      VTK_LAGRANGE_WEDGE = 73,
      VTK_LAGRANGE_PYRAMID = 74,

      // Arbitrary order Bezier elements (formulated separated from generic higher order cells)
      VTK_BEZIER_CURVE = 75,
      VTK_BEZIER_TRIANGLE = 76,
      VTK_BEZIER_QUADRILATERAL = 77,
      VTK_BEZIER_TETRAHEDRON = 78,
      VTK_BEZIER_HEXAHEDRON = 79,
      VTK_BEZIER_WEDGE = 80,
      VTK_BEZIER_PYRAMID = 81,

      VTK_NUMBER_OF_CELL_TYPES
    };
  };

  template <int D>
  class VTKOutput : public BaseVTKOutput
  {
  protected:
    shared_ptr<MeshAccess> ma = nullptr;
    Array<shared_ptr<CoefficientFunction>> coefs;
    Array<string> fieldnames;

    string filename;
    int subdivision, r;
    double h;
    int order = 1;
    int only_element = -1;
    string floatsize = "double";
    bool legacy = false;
    Array<shared_ptr<ValueField>>
        value_field;
    Array<Vec<D>> points;
    Array<VTKCell> cells;

    int output_cnt = 0;
    std::vector<double> times = {0};
    shared_ptr<ofstream> fileout;

  public:
    VTKOutput(const Array<shared_ptr<CoefficientFunction>> &,
              const Flags &, shared_ptr<MeshAccess>);

    VTKOutput(shared_ptr<MeshAccess>, const Array<shared_ptr<CoefficientFunction>> &,
              const Array<string> &, string, int, int, string, bool, int);
    virtual ~VTKOutput() { ; }

    void ResetArrays();

    void FillReferenceTrig(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems);
    void FillReferenceQuad(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems);
    void FillReferenceTet(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems);
    void FillReferenceHex(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems);
    void FillReferencePrism(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems);
    void FillReferencePyramid(Array<IntegrationPoint> &ref_coords, Array<VTKCell> &ref_elems);
    // void FillReferenceData3D(Array<IntegrationPoint> & ref_coords, Array<INT<D+1>> & ref_tets);
    // XML Methods
    void PrintPoints(int *offset, stringstream *appenddata);
    void PrintCells(int *offset, stringstream *appenddata);
    void PrintCellTypes(VorB vb, int *offset, stringstream *appenddata, const BitArray *drawelems = nullptr);
    void PrintFieldData(int *offset, stringstream *appenddata);

    void PrintAppended(stringstream *appenddata);
    void PvdFile(string filename, int index);
    // Legacy Methods
    void PrintPointsLegacy();
    void PrintCellsLegacy();
    void PrintCellTypesLegacy(VorB vb, const BitArray *drawelems = nullptr);
    void PrintFieldDataLegacy();
    virtual void Do(LocalHeap &lh, double time = -1, VorB vb = VOL, const BitArray *drawelems = 0);
  };

}
