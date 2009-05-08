#include <fem.hpp>

namespace ngfem {
  
  /****************************************************************************
   * FacetVolumeElement 
   ****************************************************************************/

  template<int D>
  FacetVolumeFiniteElement<D>::FacetVolumeFiniteElement (ELEMENT_TYPE aeltype) 
    : FiniteElement(D, aeltype,-1,-1)
  {
    for (int i=0; i<8; i++)
      vnums[i] = i;
    for (int i = 0; i < 6; i++)
      facet_order[i] = -1;
    for (int i=0; i < 7; i++)
      first_facet_dof[i] = 0;
  }
 
  template<int D>
  void FacetVolumeFiniteElement<D>::SetVertexNumbers (FlatArray<int> & avnums)
  {
    for (int i = 0; i < avnums.Size(); i++)
      vnums[i] = avnums[i];
  }
 
  template<int D>
  void FacetVolumeFiniteElement<D>::SetOrder(int ao)
  {
    order = ao;
    for (int i = 0; i < 6; i++)
      facet_order[i] = ao;
  }

  template<int D>
  void FacetVolumeFiniteElement<D>::SetOrder(FlatArray<int> & ao)
  {
    for (int i=0; i<ao.Size(); i++)
      facet_order[i] = ao[i];
  
    order = facet_order[0];        // integration order (JS, Spet 07)
    for (int i = 1; i < ao.Size(); i++)
      order = max(order, ao[i]);
  }

  template<int D>
  void FacetVolumeFiniteElement<D>::GetFacetDofNrs(int fnr, Array<int>& dnums) const
  {
    dnums.SetSize(0);
    for (int i=first_facet_dof[fnr]; i<first_facet_dof[fnr+1]; i++)
      dnums.Append(i);
  }

  template<int D>
  void FacetVolumeFiniteElement<D> :: 
  CalcFacetShape(int fnr, const IntegrationPoint & ip, FlatVector<> shape) const
  {
    shape = 0.0;
    facets[fnr]->CalcShape(ip, shape.Range(first_facet_dof[fnr], first_facet_dof[fnr+1]));
  }


 

  /****************************************************************************
   * FacetVolumeTrig
   ****************************************************************************/

  FacetVolumeTrig :: FacetVolumeTrig()
    : FacetVolumeFiniteElement<2> (ET_TRIG) 
  { 
    for (int i = 0; i < 3; i++)
      facets[i] = &facets2[i];
  };


  void FacetVolumeTrig::ComputeNDof() 
  {
    ndof = 0;
    for (int i=0; i < 3; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += facet_order[i]+1;
      }
    first_facet_dof[3] = ndof;

    const EDGE * edges = ElementTopology::GetEdges (eltype);
    ArrayMem<int,2> fvnums(2);

    for (int fnr = 0; fnr < 3; fnr++)
      {
        fvnums[0] = vnums[edges[fnr][0]]; 
        fvnums[1] = vnums[edges[fnr][1]]; 
   
        facets2[fnr].SetVertexNumbers(fvnums);
        facets2[fnr].SetOrder(facet_order[fnr]);
        facets2[fnr].ComputeNDof();
      }
  }


  /****************************************************************************
   * FacetVolumeQuad
   ****************************************************************************/

  FacetVolumeQuad :: FacetVolumeQuad()
    : FacetVolumeFiniteElement<2> (ET_QUAD) 
  { 
    for (int i = 0; i < 4; i++)
      facets[i] = &facets2[i];
  };


  // -------------------------------------------------------------------------------
  /*
    void FacetVolumeQuad::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
    {
    // topology: points = 0:(0,0), 1:(1,0), 2:(1,1), 3:(0,1)
    //           edges = 0:(0,1), 1:(2,3), 2:(3,0), 4:(1,2)
    IntegrationPoint ip1d;  
    double x=ip(0), y=ip(1);
    int fnr=-1;

    shape = 0.0;
    if (fabs(y)<1e-12)
    {
    fnr=0;
    ip1d(0) = x;
    //     CalcFacetShape(0, ip1d, shape.Range(first_facet_dof[0], first_facet_dof[1]));
    CalcFacetShape(0, ip1d, shape);
    }
    else if (fabs(y)>1-(1e-12))
    {
    fnr=1;
    ip1d(0) = (1-x);
    //     CalcFacetShape(1, ip1d, shape.Range(first_facet_dof[1], first_facet_dof[2]));
    CalcFacetShape(1, ip1d, shape);
    }
    else if (fabs(x)<1e-12)
    {
    fnr=2;
    ip1d(0) = 1-y;
    //     CalcFacetShape(2, ip1d, shape.Range(first_facet_dof[2], first_facet_dof[3]));
    CalcFacetShape(2, ip1d, shape);
    }
    else if (fabs(x)>1-(1e-12))
    {
    fnr=3;
    ip1d(0) = y;
    //     CalcFacetShape(3, ip1d, shape.Range(first_facet_dof[3], first_facet_dof[4]));
    CalcFacetShape(3, ip1d, shape);
    }

    }
  */
 
  /*
  void FacetVolumeQuad::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
  {
    SetFacet(afnr);
    shape=0.0;
    facet.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
  }
   
  void FacetVolumeQuad::SetFacet(int afnr) const
  {
    if (fnr == afnr) return;
  
    FacetVolumeQuad * quad=const_cast<FacetVolumeQuad*>(this);
    quad->fnr = afnr;
    const EDGE * edges = ElementTopology::GetEdges (eltype);
    Array<int> fvnums(2);
    fvnums[0] = vnums[edges[fnr][0]]; 
    fvnums[1] = vnums[edges[fnr][1]]; 
   
    quad->facet.SetVertexNumbers(fvnums);
    quad->facet.SetOrder(facet_order[fnr]);
  }  
  */

  void FacetVolumeQuad::ComputeNDof() 
  {
    ndof = 0;
    for (int i=0; i < 4; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += facet_order[i]+1;
      }
    first_facet_dof[4] = ndof;


    const EDGE * edges = ElementTopology::GetEdges (ET_QUAD);
    Array<int> fvnums(2);

    for (int fnr = 0; fnr < 4; fnr++)
      {
        fvnums[0] = vnums[edges[fnr][0]]; 
        fvnums[1] = vnums[edges[fnr][1]]; 
   
        facets2[fnr].SetVertexNumbers(fvnums);
        facets2[fnr].SetOrder(facet_order[fnr]);
        facets2[fnr].ComputeNDof();
      }
  }




  /****************************************************************************
   * FacetVolumeTet
   ****************************************************************************/

  FacetVolumeTet :: FacetVolumeTet()
    : FacetVolumeFiniteElement<3> (ET_TET) 
  {
    for (int i = 0; i < 4; i++)
      facets[i] = &facets2[i];
  };


  void FacetVolumeTet::ComputeNDof() 
  {
    ndof = 0;
    for (int i = 0; i < 4; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += ( (facet_order[i]+1) * (facet_order[i]+2) ) / 2;
      }
    first_facet_dof[4] = ndof;


    Array<int> fvnums(3);  
    const FACE * faces = ElementTopology::GetFaces (eltype);

    for (int fnr = 0; fnr < 4; fnr++)
      {
        fvnums[0] = vnums[faces[fnr][0]]; 
        fvnums[1] = vnums[faces[fnr][1]]; 
        fvnums[2] = vnums[faces[fnr][2]]; 
   
        facets2[fnr].SetVertexNumbers(fvnums);
        facets2[fnr].SetOrder(facet_order[fnr]);
        facets2[fnr].ComputeNDof();
      }
  }




  /****************************************************************************
   * FacetVolumeHex
   ****************************************************************************/
  // -------------------------------------------------------------------------------

  FacetVolumeHex :: FacetVolumeHex()
    : FacetVolumeFiniteElement<3> (ET_HEX) 
  {
    for (int i = 0; i < 6; i++)
      facets[i] = &facetsq[i];
  };


  /*
    void FacetVolumeHex::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
    {
    cout << "error: CalcShape not implemented yet for FacetVolumeHex" << endl;
    exit(0);
    }
  */

  /*
  void FacetVolumeHex::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
  {
    shape=0.0;
    SetFacet(afnr);
    facet.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
  }
   
  void FacetVolumeHex::SetFacet(int afnr) const
  {
    if (fnr == afnr) return;
  
    FacetVolumeHex * hex=const_cast<FacetVolumeHex*>(this);
    hex->fnr = afnr;
    Array<int> fvnums(4);  
    const FACE * faces = ElementTopology::GetFaces (eltype);

    fvnums[0] = vnums[faces[fnr][0]]; 
    fvnums[1] = vnums[faces[fnr][1]]; 
    fvnums[2] = vnums[faces[fnr][2]]; 
    fvnums[3] = vnums[faces[fnr][3]]; 
   
    hex->facet.SetVertexNumbers(fvnums);
    hex->facet.SetOrder(facet_order[fnr]); 
  }  
  */

  void FacetVolumeHex::ComputeNDof()  
  {
    ndof = 0;
    for (int i=0; i < 6; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += (facet_order[i]+1) * (facet_order[i]+1);
      }
    first_facet_dof[6] = ndof;


    Array<int> fvnums(4);  
    const FACE * faces = ElementTopology::GetFaces (ET_HEX);

    for (int fnr = 0; fnr < 6; fnr++)
      {
        fvnums[0] = vnums[faces[fnr][0]]; 
        fvnums[1] = vnums[faces[fnr][1]]; 
        fvnums[2] = vnums[faces[fnr][2]]; 
        fvnums[3] = vnums[faces[fnr][3]]; 
   
        facetsq[fnr].SetVertexNumbers(fvnums);
        facetsq[fnr].SetOrder(facet_order[fnr]);
        facetsq[fnr].ComputeNDof();
      }

  }

  // const FiniteElement & FacetVolumeHex::GetFacetFE(int fnr, LocalHeap& lh) const
  // {
  //   FacetFacetQuad* fe = new(lh.Alloc(sizeof(FacetFacetQuad))) FacetFacetQuad();
  //    
  //   const FACE * faces = ElementTopology::GetFaces (eltype);
  //   Array<int> fvnums(4);
  //   for (int i=0; i<4; i++)
  //     fvnums[i] = vnums[faces[fnr][i]]; 
  //   
  //   fe->SetVertexNumbers(fvnums);
  //   fe->SetOrder(facet_order[fnr]);
  //    
  //   return *fe;
  // }



  /****************************************************************************
   * FacetVolumePrism
   ****************************************************************************/
  // -------------------------------------------------------------------------------

  FacetVolumePrism :: FacetVolumePrism()
    : FacetVolumeFiniteElement<3> (ET_PRISM) 
  {
    for (int i = 0; i < 2; i++)
      facets[i] = &facetst[i];
    for (int i = 0; i < 3; i++)
      facets[i+2] = &facetsq[i];
  };


  /*
    void FacetVolumePrism::CalcShape (const IntegrationPoint & ip3d, FlatVector<> shape) const
    {  
    //   topology: points = 0:()
    IntegrationPoint ip2d;  
    double x=ip3d(0), y=ip3d(1), z=ip3d(2);

    shape = 0.0;
    if (fabs(y)<1e-12) // front
    {
    ip2d(0) = x; ip2d(1) = z;
    CalcFacetShape(2, ip2d, shape);
    }
    else if (fabs(x)<1e-12) // left
    {
    ip2d(0) = 1-y; ip2d(1)=z;
    CalcFacetShape(4, ip2d, shape);
    }
    else if (fabs(1-x-y)<1e-12) // right
    {
    ip2d(0) = (1-x+y)*0.5; ip2d(1) = z;
    CalcFacetShape(3, ip2d, shape);
    }
    else if (fabs(z)<1e-12) // botom
    {
    ip2d(0) = 1-x-y; ip2d(1) = y;
    CalcFacetShape(0, ip2d, shape);
    }
    else if (fabs(z)>1-1e-12) // top
    {
    ip2d(0) = x; ip2d(1) = y;
    CalcFacetShape(1, ip2d, shape);
    }
    } 
  */

  /*
  void FacetVolumePrism::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
  {
    shape=0.0;
    SetFacet(afnr);
    if (afnr < 2) //
      trig.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
    else // quad
      quad.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
  }
   
  void FacetVolumePrism::SetFacet(int afnr) const
  {
    if (qnr == afnr || tnr == afnr) return;
  
    FacetVolumePrism * prism=const_cast<FacetVolumePrism*>(this);
    if (afnr < 2) // triangles
      {
        prism->tnr = afnr;
        Array<int> fvnums(3);  
        const FACE * faces = ElementTopology::GetFaces (eltype);

        fvnums[0] = vnums[faces[afnr][0]]; 
        fvnums[1] = vnums[faces[afnr][1]]; 
        fvnums[2] = vnums[faces[afnr][2]]; 
   
        prism->trig.SetVertexNumbers(fvnums);
        prism->trig.SetOrder(facet_order[tnr]);
      }
    else // quad face
      {
        prism->qnr = afnr;
        Array<int> fvnums(4);  
        const FACE * faces = ElementTopology::GetFaces (eltype);

        fvnums[0] = vnums[faces[afnr][0]]; 
        fvnums[1] = vnums[faces[afnr][1]]; 
        fvnums[2] = vnums[faces[afnr][2]]; 
        fvnums[3] = vnums[faces[afnr][3]]; 
   
        prism->quad.SetVertexNumbers(fvnums);
        prism->quad.SetOrder(facet_order[qnr]);
      }
  }  
  */


  void FacetVolumePrism::ComputeNDof() 
  {
    ndof = 0;
    // triangles
    for (int i=0; i<2; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += ( (facet_order[i]+1) * (facet_order[i]+2) ) / 2;
      }
    //quads
    for (int i=2; i<5; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += (facet_order[i]+1) * (facet_order[i]+1);
      }
  
    first_facet_dof[5] = ndof;

    Array<int> fvnums(4);  
    const FACE * faces = ElementTopology::GetFaces (ET_PRISM);

    for (int fnr = 0; fnr < 5; fnr++)
      {
        fvnums[0] = vnums[faces[fnr][0]]; 
        fvnums[1] = vnums[faces[fnr][1]]; 
        fvnums[2] = vnums[faces[fnr][2]]; 
        fvnums[3] = vnums[faces[fnr][3]]; 
   
        facets[fnr]->SetVertexNumbers(fvnums);
        facets[fnr]->SetOrder(facet_order[fnr]);
        facets[fnr]->ComputeNDof();
      }


  }

  // const FiniteElement & FacetVolumePrism::GetFacetFE(int fnr, LocalHeap& lh) const
  // {
  //    
  //   const FACE * faces = ElementTopology::GetFaces (eltype);
  //   if (fnr < 2) // triangle
  //   {
  //     FacetFacetTrig* fe = new(lh.Alloc(sizeof(FacetFacetTrig))) FacetFacetTrig();
  //     Array<int> fvnums(3);
  //     for (int i=0; i<3; i++)
  //       fvnums[i] = vnums[faces[fnr][i]]; 
  //     fe->SetVertexNumbers(fvnums);
  //     fe->SetOrder(facet_order[fnr]);
  //     return *fe;
  //   }
  //   else // quad
  //   {
  //     FacetFacetQuad* fe = new(lh.Alloc(sizeof(FacetFacetQuad))) FacetFacetQuad();
  //     Array<int> fvnums(4);
  //     for (int i=0; i<4; i++)
  //       fvnums[i] = vnums[faces[fnr][i]]; 
  //     fe->SetVertexNumbers(fvnums);
  //     fe->SetOrder(facet_order[fnr]);
  //     return *fe;
  //   }
  // }


  /****************************************************************************
   * FacetVolumePyramid
   ****************************************************************************/
  // -------------------------------------------------------------------------------
  /*
    void FacetVolumePyramid::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
    {
    cout << "error: CalcShape not implemented yet for FacetVolumePyramid" << endl;
    exit(0);
    }
  */
 
  void FacetVolumePyramid::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
  {
    shape=0.0;
    SetFacet(afnr);
    if (afnr < 4) // trig
      trig.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
    else
      quad.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
  }
   
  void FacetVolumePyramid::SetFacet(int afnr) const
  {
    if (qnr == afnr || tnr == afnr) return;
  
    FacetVolumePyramid * pyramid=const_cast<FacetVolumePyramid*>(this);
    if (afnr < 4) // triangles
      {
        pyramid->tnr = afnr;
        Array<int> fvnums(3);  const FACE * faces = ElementTopology::GetFaces (eltype);

        fvnums[0] = vnums[faces[tnr][0]]; 
        fvnums[1] = vnums[faces[tnr][1]]; 
        fvnums[2] = vnums[faces[tnr][2]]; 
   
        pyramid->trig.SetVertexNumbers(fvnums);
        pyramid->trig.SetOrder(facet_order[tnr]);
      }
    else // quad face
      {
        pyramid->qnr = afnr;
        Array<int> fvnums(4);  const FACE * faces = ElementTopology::GetFaces (eltype);

        fvnums[0] = vnums[faces[qnr][0]]; 
        fvnums[1] = vnums[faces[qnr][1]]; 
        fvnums[2] = vnums[faces[qnr][2]]; 
        fvnums[3] = vnums[faces[qnr][3]]; 
   
        pyramid->quad.SetVertexNumbers(fvnums);
        pyramid->quad.SetOrder(facet_order[qnr]);
      }
  }  


  void FacetVolumePyramid::ComputeNDof() 
  {
    ndof = 0;
    // triangles
    for (int i=0; i<4; i++)
      {
        first_facet_dof[i] = ndof;
        ndof += ( (facet_order[i]+1) * (facet_order[i]+2) ) / 2;
      }
    //quad - basis
    first_facet_dof[4] = ndof;
    ndof += (facet_order[4]+1) * (facet_order[4]+1);
  
    // final
    first_facet_dof[5] = ndof;
  }

  // const FiniteElement & FacetVolumePyramid::GetFacetFE(int fnr, LocalHeap& lh) const
  // {
  //   const FACE * faces = ElementTopology::GetFaces (eltype);
  //   if (fnr < 4) // triangle
  //   {
  //     FacetFacetTrig* fe = new(lh.Alloc(sizeof(FacetFacetTrig))) FacetFacetTrig();
  //     Array<int> fvnums(3);
  //     for (int i=0; i<3; i++)
  //       fvnums[i] = vnums[faces[fnr][i]]; 
  //     fe->SetVertexNumbers(fvnums);
  //     fe->SetOrder(facet_order[fnr]);
  //     return *fe;
  //   }
  //   else // facet 5: quad
  //   {
  //     FacetFacetQuad* fe = new(lh.Alloc(sizeof(FacetFacetQuad))) FacetFacetQuad();
  //     Array<int> fvnums(4);
  //     for (int i=0; i<4; i++)
  //       fvnums[i] = vnums[faces[fnr][i]]; 
  //     fe->SetVertexNumbers(fvnums);
  //     fe->SetOrder(facet_order[fnr]);
  //     return *fe;
  //   }
  // }


  
  //template class FacetVolumeFiniteElement<1>;

  template class FacetVolumeFiniteElement<2>;
  template class FacetVolumeFiniteElement<3>;

} // namespace


