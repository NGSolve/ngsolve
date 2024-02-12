#ifndef FILE_L2HOFEFO
#define FILE_L2HOFEFO

/*********************************************************************/
/* File:   l2hofefo.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   Apr. 2009                                                 */
/*********************************************************************/

#include "thdivfe.hpp"


namespace ngfem
{

  /**
     High order finite elements for L2 of fixed order
  */

  /*
  template<int DIM>
  class HDivHighOrderFiniteElementFO : virtual public HDivFiniteElement<DIM>
  {
  protected:
    int vnums[8];
    bool ho_div_free;
    bool only_ho_div;
  public:
    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }
    void SetHODivFree (bool aho_div_free) { ho_div_free = aho_div_free; only_ho_div = only_ho_div && !ho_div_free;};  
    void SetOnlyHODiv (bool aonly_ho_div) { only_ho_div = aonly_ho_div; ho_div_free = ho_div_free && !only_ho_div;}; 
    virtual void ComputeNDof () = 0;
  };
  */
  
  template <ELEMENT_TYPE ET, int ORDER> class HDivHighOrderFEFO;


  /*
  template <ELEMENT_TYPE ET, int ORDER>
  class T_HDivHighOrderFiniteElementFO : 
    public HDivHighOrderFiniteElementFO<ET_trait<ET>::DIM>,
    public T_HDivFiniteElement< HDivHighOrderFEFO<ET,ORDER>, ET >,
    public ET_trait<ET>
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    using HDivFiniteElement<DIM>::ndof;
    using HDivFiniteElement<DIM>::order;

    using HDivHighOrderFiniteElementFO<DIM>::vnums;
    using HDivHighOrderFiniteElementFO<DIM>::ho_div_free;
    using HDivHighOrderFiniteElementFO<DIM>::only_ho_div;

    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::N_EDGE;
    using ET_trait<ET>::N_FACE;
    using ET_trait<ET>::FaceType;
    using ET_trait<ET>::GetEdgeSort;
    using ET_trait<ET>::GetFaceSort;



  public:

    T_HDivHighOrderFiniteElementFO () 
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
	vnums[i] = i;
      ho_div_free = false;
      only_ho_div = false;
      order = ORDER;
    }
  };
  */



  /**
     High order triangular finite element
  */
  template <int ORDER>
  class HDivHighOrderFEFO<ET_TRIG, ORDER> : 
    public T_HDivFiniteElement<HDivHighOrderFEFO<ET_TRIG,ORDER>, ET_TRIG>,
    public ET_trait<ET_TRIG>
  {
    using T_HDivFiniteElement<HDivHighOrderFEFO<ET_TRIG,ORDER>, ET_TRIG> :: ndof;
    using T_HDivFiniteElement<HDivHighOrderFEFO<ET_TRIG,ORDER>, ET_TRIG> :: order;
    /*
    using T_HDivFiniteElement<ET_TRIG, ORDER>::ndof;
    using T_HDivFiniteElement<ET_TRIG, ORDER>::vnums; 
    using T_HDivFiniteElement<ET_TRIG, ORDER>::ho_div_free; 
    using T_HDivFiniteElement<ET_TRIG, ORDER>::only_ho_div; 
    */

    int vnums[3];
    bool ho_div_free;
    bool only_ho_div;

  public:

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
    typedef TrigShapesInnerLegendre T_TRIGSHAPES;

  public:
    HDivHighOrderFEFO () 
      : ho_div_free(false), only_ho_div(false)
    { 
      order = ORDER;
    }


    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }
    void SetHODivFree (bool aho_div_free) { ho_div_free = aho_div_free; only_ho_div = only_ho_div && !ho_div_free;};  
    void SetOnlyHODiv (bool aonly_ho_div) { only_ho_div = aonly_ho_div; ho_div_free = ho_div_free && !only_ho_div;}; 



    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }

    virtual void ComputeNDof()
    {
      if (only_ho_div)
        ndof = ORDER*(ORDER+1)/2 - 1;
      else
      {
        ndof = 3 * (ORDER+1);
        ndof += (ho_div_free) ? ORDER*(ORDER-1)/2 : ORDER*ORDER-1;
      }
    }

    virtual void GetInternalDofs (Array<int> & idofs) const
    {
      idofs.SetSize(0);
      idofs += IntRange (3*(ORDER+1), ndof);
    }

    virtual void GetFacetDofs(int i, Array<int> & dnums) const
    {
      dnums.SetSize (0);
      dnums += i;
      dnums += 3+IntRange (i*ORDER,(i+1)*ORDER);
    }


    template<typename Tx, typename TFA>   
    void T_CalcShape (TIP<2,Tx> ip, TFA & shape) const
    {
      if (only_ho_div && (ORDER <= 1)) return;
      Tx x = ip.x, y = ip.y;
      Tx lami[3] = { x, y, 1-x-y };

      Tx adpol1[ORDER], adpol2[ORDER];

      int ii = 3; 

      if (!only_ho_div){
        // const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
        for (int i = 0; i < 3; i++)
        {
          IVec<2> e = this->GetEdgeSort (i, vnums);
          
          //Nedelec low order edge shape function 
          shape[i] = uDv_minus_vDu (lami[e[0]], lami[e[1]]);
          
          //HO-Edge shapes (Gradient Fields)   
          if(ORDER > 0) //  && usegrad_edge[i]) 
            { 
              Tx xi = lami[e[1]] - lami[e[0]]; 
              // LegendrePolynomial::
              IntLegNoBubble::
                EvalScaledMult (ORDER-1, xi, lami[e[0]]+lami[e[1]], 
                                lami[e[0]]*lami[e[1]], 
                                SBLambda([&](int i, Tx v)
                                         {
                                           shape[ii++] = Du(v);
                                         }));
              /*
              AutoDiff<2> eta = 1 - lami[e[1]] - lami[e[0]]; 
              // T_ORTHOPOL::CalcTrigExt(ORDER+1, xi, eta, adpol1); 
              T_ORTHOPOL::CalcScaled<ORDER+1> (xi, 1-eta, adpol1); 
              
              for(int j = 0; j < ORDER; j++) 
              shape[ii++] = Du<2> (adpol1[j]);
              */
            }
        }   
      }
      else
        ii = 0;
      //Inner shapes (Face) 
      if(ORDER > 1) 
	{
	  IVec<4> fav = this->GetFaceSort (0, vnums);

	  Tx xi  = lami[fav[2]]-lami[fav[1]];
	  Tx eta = lami[fav[0]]; 
	  
	  TrigShapesInnerLegendre::CalcSplitted<ORDER+1> (xi, eta, adpol1,adpol2);
	  
        if (!only_ho_div)
        {
          // rotated gradients:
          for (int j = 0; j < ORDER-1; j++)
            for (int k = 0; k < ORDER-1-j; k++, ii++)
              shape[ii] = Du (adpol1[j] * adpol2[k]);
	  }
	  
	  if (!ho_div_free)
	    {
	      // other combination
	      for (int j = 0; j < ORDER-1; j++)
		for (int k = 0; k < ORDER-1-j; k++, ii++)
		  shape[ii] = uDv_minus_vDu (adpol2[k], adpol1[j]);
	      
	      // rec_pol * Nedelec0 
	      for (int j = 0; j < ORDER-1; j++, ii++)
		shape[ii] = wuDv_minus_wvDu (lami[fav[1]], lami[fav[2]], adpol2[j]);
	    }
	}
    }
    
  };

}


#endif
