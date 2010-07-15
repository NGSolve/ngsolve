#include <fem.hpp>
#include "facethofe.hpp"

namespace ngfem {
  
  /****************************************************************************
   * FacetVolumeElement 
   ****************************************************************************/

  template<int D>
  FacetVolumeFiniteElement<D>::FacetVolumeFiniteElement (ELEMENT_TYPE aeltype) 
    : ScalarFiniteElement<D> (aeltype,-1,-1)
  {
    // highest_order_dc = false;
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


  template<int D>
  void FacetVolumeFiniteElement<D> :: 
  EvaluateFacet (int fnr, const IntegrationRule & ir, FlatVector<> coefs, FlatVector<> values) const
  {
    facets[fnr]->Evaluate (ir, coefs.Range(first_facet_dof[fnr], first_facet_dof[fnr+1]), values);
  }

  template<int D>
  void FacetVolumeFiniteElement<D> :: 
  EvaluateFacetTrans (int fnr, const IntegrationRule & ir, FlatVector<> values, FlatVector<> coefs) const
  {
    coefs = 0.0;
    facets[fnr]->EvaluateTrans (ir, values, coefs.Range(first_facet_dof[fnr], first_facet_dof[fnr+1]));
  }
  
  /*
  template <int D>
  void FacetVolumeFiniteElement<D>::
  GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize(0);
    if (highest_order_dc)
      {
	if (D==2)
	  {
	    for (int i=0; i<ElementTopology::GetNFacets(eltype); i++)
	      idofs.Append(first_facet_dof[i+1]-1);   
	  }
	else
	  {
	    for (int i=0; i<ElementTopology::GetNFacets(eltype); i++)
	      {
		int pos = first_facet_dof[i]-2;
		int fac = 4 - ElementTopology::GetNVertices(ElementTopology::GetFaceType(eltype,i));
		for (int k = 0; k <= facet_order[i]; k++){
		  pos += (facet_order[i]+1-fac*k);
		  idofs.Append(pos);
		}
	      }
	  }//end if Dimension
      }
  }
  */
 

  /****************************************************************************
   * FacetVolumeTrig
   ****************************************************************************/

  FacetFE<ET_TRIG> :: FacetFE()
  { 
    for (int i = 0; i < 3; i++)
      facets[i] = &facets2[i];
  };


  void FacetFE<ET_TRIG>::ComputeNDof() 
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

  FacetFE<ET_QUAD> :: FacetFE()
  { 
    for (int i = 0; i < 4; i++)
      facets[i] = &facets2[i];
  };


  // -------------------------------------------------------------------------------
  /*
    void FacetFE<ET_QUAD>::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
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
    void FacetFE<ET_QUAD>::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
    {
    SetFacet(afnr);
    shape=0.0;
    facet.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
    }
   
    void FacetFE<ET_QUAD>::SetFacet(int afnr) const
    {
    if (fnr == afnr) return;
  
    FacetFE<ET_QUAD> * quad=const_cast<FacetFE<ET_QUAD>*>(this);
    quad->fnr = afnr;
    const EDGE * edges = ElementTopology::GetEdges (eltype);
    Array<int> fvnums(2);
    fvnums[0] = vnums[edges[fnr][0]]; 
    fvnums[1] = vnums[edges[fnr][1]]; 
   
    quad->facet.SetVertexNumbers(fvnums);
    quad->facet.SetOrder(facet_order[fnr]);
    }  
  */

  void FacetFE<ET_QUAD>::ComputeNDof() 
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

  FacetFE<ET_TET> :: FacetFE()
  {
    ;
    for (int i = 0; i < 4; i++)
      facets[i] = &facets2[i];
  };


  void FacetFE<ET_TET>::ComputeNDof() 
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

  FacetFE<ET_HEX> :: FacetFE()
  {
    for (int i = 0; i < 6; i++)
      facets[i] = &facetsq[i];
  };


  /*
    void FacetFE<ET_HEX>::CalcShape (const IntegrationPoint & ip, FlatVector<> shape) const
    {
    cout << "error: CalcShape not implemented yet for FacetFE<ET_HEX>" << endl;
    exit(0);
    }
  */

  /*
    void FacetFE<ET_HEX>::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
    {
    shape=0.0;
    SetFacet(afnr);
    facet.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
    }
   
    void FacetFE<ET_HEX>::SetFacet(int afnr) const
    {
    if (fnr == afnr) return;
  
    FacetFE<ET_HEX> * hex=const_cast<FacetFE<ET_HEX>*>(this);
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

  void FacetFE<ET_HEX>::ComputeNDof()  
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

  // const FiniteElement & FacetFE<ET_HEX>::GetFacetFE(int fnr, LocalHeap& lh) const
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

  FacetFE<ET_PRISM> :: FacetFE()
  {
    ;
    /*
    for (int i = 0; i < 2; i++)
      facets[i] = &facetst[i];
    for (int i = 0; i < 3; i++)
      facets[i+2] = &facetsq[i];
    */
  };



  void FacetFE<ET_PRISM>::ComputeNDof() 
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

    /*
    Array<int> fvnums(4);  
    const FACE * faces = ElementTopology::GetFaces (ET_PRISM);

    for (int fnr = 0; fnr < 2; fnr++)
      {
        fvnums[0] = vnums[faces[fnr][0]]; 
        fvnums[1] = vnums[faces[fnr][1]]; 
        fvnums[2] = vnums[faces[fnr][2]]; 
        fvnums[3] = vnums[faces[fnr][3]]; 
   
        facetst[fnr].SetVertexNumbers(fvnums);
        facetst[fnr].SetOrder(facet_order[fnr]);
        facetst[fnr].ComputeNDof();
      }
    for (int fnr = 0; fnr < 3; fnr++)
      {
        fvnums[0] = vnums[faces[fnr+2][0]]; 
        fvnums[1] = vnums[faces[fnr+2][1]]; 
        fvnums[2] = vnums[faces[fnr+2][2]]; 
        fvnums[3] = vnums[faces[fnr+2][3]]; 
   
        facetsq[fnr].SetVertexNumbers(fvnums);
        facetsq[fnr].SetOrder(facet_order[fnr+2]);
        facetsq[fnr].ComputeNDof();
      }
    */




  }

  // const FiniteElement & FacetFE<ET_PRISM>::GetFacetFE(int fnr, LocalHeap& lh) const
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
 
  void FacetFE<ET_PYRAMID>::CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const
  {
    shape=0.0;
    SetFacet(afnr);
    if (afnr < 4) // trig
      trig.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
    else
      quad.CalcShape(ip, shape.Range(first_facet_dof[afnr], first_facet_dof[afnr+1]));
  }
   
  void FacetFE<ET_PYRAMID>::SetFacet(int afnr) const
  {
    if (qnr == afnr || tnr == afnr) return;
  
    FacetFE<ET_PYRAMID> * pyramid=const_cast<FacetFE<ET_PYRAMID>*>(this);
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


  void FacetFE<ET_PYRAMID>::ComputeNDof() 
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




















  template <int DIM>
  EdgeVolumeFiniteElement<DIM> :: EdgeVolumeFiniteElement (ELEMENT_TYPE aeltype, int aorder)
  {
    eltype = aeltype;
    ndof = (aorder+1) * ElementTopology::GetNEdges (aeltype);
    order = aorder;

    edgenr = 0;
  }

  template <int DIM>
  void EdgeVolumeFiniteElement<DIM> :: CalcEdgeShape(int enr, const IntegrationPoint & ip, FlatVector<> shape) const
  {
    switch (eltype)
      {
      case ET_TRIG:
	{
	  double lami[3] = { ip(0), ip(1), 1-ip(0)-ip(1) };
	
	  const EDGE & edge = ElementTopology::GetEdges(eltype)[enr];
	  int vs = edge[0], ve = edge[1];
	  if (vnums[vs] > vnums[ve]) swap (vs, ve); 
	
	  double s = lami[ve] - lami[vs];
	
	  shape = 0;
	  int base = (order+1) * enr;
	  LegendrePolynomial (order, s, shape.Range(base, base+order+1));

	  break;
	}
      case ET_TET:
	{
	  double lami[4] = { ip(0), ip(1), ip(2), 1-ip(0)-ip(1)-ip(2) };
	
	  const EDGE & edge = ElementTopology::GetEdges(eltype)[enr];
	  int vs = edge[0], ve = edge[1];
	  if (vnums[vs] > vnums[ve]) swap (vs, ve);
	
	  double s = lami[ve] - lami[vs];
	
	  shape = 0;
	  int base = (order+1) * enr;
	  LegendrePolynomial (order, s, shape.Range(base, base+order+1));
	  break;
	}
      default:
	throw Exception ("EdgeVolumeFiniteElement not implemnted only for tets");
      }
  }





  
  //template class FacetVolumeFiniteElement<1>;

  template class EdgeVolumeFiniteElement<2>;
  template class EdgeVolumeFiniteElement<3>;

  template class FacetVolumeFiniteElement<2>;
  template class FacetVolumeFiniteElement<3>;

} // namespace


