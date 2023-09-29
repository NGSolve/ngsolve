/*********************************************************************/
/* File:   intrule.cpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Mar. 2000                                             */
/*********************************************************************/

/* 
   Finite Element Integration Rules
*/


#include <fem.hpp>
   
namespace ngfem
{
  static mutex intruletpfacet_mutex;
  static mutex genintrule_mutex;
  static mutex simd_genintrule_mutex[40];

  ostream & operator<< (ostream & ost, const IntegrationPoint & ip)
  {
    ost << " locnr = " << ip.nr << ": (" 
	<< ip.pi[0] << ", " << ip.pi[1] << ", " << ip.pi[2] 
	<< "), weight = " << ip.weight;
    return ost;
  }
  

  BaseMappedIntegrationPoint :: ~BaseMappedIntegrationPoint ()
  {
    if (owns_trafo) delete eltrans;
  }


  FlatVector<> BaseMappedIntegrationPoint :: GetPoint() const
  {
    switch (eltrans->SpaceDim())
      {
      case 1: return static_cast<const DimMappedIntegrationPoint<1>&> (*this).GetPoint();
      case 2: return static_cast<const DimMappedIntegrationPoint<2>&> (*this).GetPoint();
      case 3: return static_cast<const DimMappedIntegrationPoint<3>&> (*this).GetPoint();
      }

    throw Exception("BaseMappedIntegrationPoint::GetPoint, illegal dimension");
  }

  


  
  FlatMatrix<> BaseMappedIntegrationPoint :: GetJacobian() const
  {
    switch(eltrans->VB())
      {
      case VOL:
	switch (eltrans->SpaceDim())
	  {
	  case 1: return static_cast<const MappedIntegrationPoint<1,1>&> (*this).GetJacobian();
	  case 2: return static_cast<const MappedIntegrationPoint<2,2>&> (*this).GetJacobian();
	  case 3: return static_cast<const MappedIntegrationPoint<3,3>&> (*this).GetJacobian();
	  }
        break;
      case BND:
	switch (eltrans->SpaceDim())
	  {
	  case 1: return static_cast<const MappedIntegrationPoint<0,1>&> (*this).GetJacobian();
	  case 2: return static_cast<const MappedIntegrationPoint<1,2>&> (*this).GetJacobian();
	  case 3: return static_cast<const MappedIntegrationPoint<2,3>&> (*this).GetJacobian();
	  }
        break;
      case BBND:
	switch (eltrans->SpaceDim())
	  {
	  case 2: return static_cast<const MappedIntegrationPoint<0,2>&> (*this).GetJacobian();
	  case 3: return static_cast<const MappedIntegrationPoint<1,3>&> (*this).GetJacobian();
	  }
        break;
      case BBBND:
	switch (eltrans->SpaceDim())
	  {
	  case 3: return static_cast<const MappedIntegrationPoint<0,3>&> (*this).GetJacobian();
	  }
        break;
      }
    throw Exception("BaseMappedIntegrationPoint::GetJacobian, illegal dimension");
  }

  /*
  int BaseMappedIntegrationPoint :: DimElement() const
  {
    switch(eltrans->VB())
      {
      case VOL:
        switch (eltrans->SpaceDim())
	  {
          case 1: return 1;
          case 2: return 2;
          case 3: return 3;
          }
        break;
      case BND:
        switch (eltrans->SpaceDim())
	  {
          case 1: return 0;
          case 2: return 1;
          case 3: return 2;
          }
        break;
      case BBND:
        switch (eltrans->SpaceDim())
	  {
          case 2: return 0;
          case 3: return 1;
          }
        break;
      case BBBND:
        switch (eltrans->SpaceDim())
	  {
          case 3: return 0;
          }
        break;
      }
    throw Exception("BaseMappedIntegrationPoint::DimElement, illegal dimension");
  }
  
  int BaseMappedIntegrationPoint :: DimSpace() const
  {
    return eltrans->SpaceDim();
  }
*/
  
  FlatVector<Complex> BaseMappedIntegrationPoint :: GetPointComplex() const
  {
    if (is_complex)
    {
      switch (eltrans->SpaceDim())
        {
        case 1: return static_cast<const DimMappedIntegrationPoint<1,Complex>&> (*this).GetPoint();
        case 2: return static_cast<const DimMappedIntegrationPoint<2,Complex>&> (*this).GetPoint();
        case 3: return static_cast<const DimMappedIntegrationPoint<3,Complex>&> (*this).GetPoint();
        }
      throw Exception("BaseMappedIntegrationPoint::GetPointComplex, illegal dimension");
    }
    throw Exception("BaseMappedIntegrationPoint::GetPointComplex, is not complex");
  }

  


  
  FlatMatrix<Complex> BaseMappedIntegrationPoint :: GetJacobianComplex() const
  {
    if (is_complex)
      {
        switch(eltrans->VB())
          {
          case VOL:
            switch (eltrans->SpaceDim())
              {
              case 1: return static_cast<const MappedIntegrationPoint<1,1,Complex>&> (*this).GetJacobian();
              case 2: return static_cast<const MappedIntegrationPoint<2,2,Complex>&> (*this).GetJacobian();
              case 3: return static_cast<const MappedIntegrationPoint<3,3,Complex>&> (*this).GetJacobian();
              }
            break;
          case BND:
            switch (eltrans->SpaceDim())
              {
              case 1: return static_cast<const MappedIntegrationPoint<0,1,Complex>&> (*this).GetJacobian();
              case 2: return static_cast<const MappedIntegrationPoint<1,2,Complex>&> (*this).GetJacobian();
              case 3: return static_cast<const MappedIntegrationPoint<2,3,Complex>&> (*this).GetJacobian();
              }
            break;
          case BBND:
            switch (eltrans->SpaceDim())
              {
              case 2: return static_cast<const MappedIntegrationPoint<0,2,Complex>&> (*this).GetJacobian();
              case 3: return static_cast<const MappedIntegrationPoint<1,3,Complex>&> (*this).GetJacobian();
              }
            break;
          case BBBND:
            switch (eltrans->SpaceDim())
              {
              case 3: return static_cast<const MappedIntegrationPoint<0,3,Complex>&> (*this).GetJacobian();
              }
            break;
          }
        throw Exception("BaseMappedIntegrationPoint::GetJacobian, illegal dimension");
      }
    throw Exception("BaseMappedIntegrationPoint::GetPointComplex, is not complex");
  }

  /*
  int BaseMappedIntegrationPoint :: Dim() const
  {
    return eltrans->SpaceDim();
  }
  */
  
  /*
  VorB BaseMappedIntegrationPoint :: ElementVB() const
  {
    // cout << "BaseMappedIP, eltrans->VB = " << eltrans->VB() << ", ip.vb = " << ip.VB() << endl;
    return ip.VB();
    // return eltrans->VB();
  }
  */
  

  template <int S, int R, typename SCAL>
  void MappedIntegrationPoint<S,R,SCAL> :: 
  CalcHesse (Mat<2> & ddx1, Mat<2> & ddx2) const
  {
    double eps = 1e-6;

    Mat<2> jacr, jacl;
    for (int dir = 0; dir < 2; dir++)
      {
	IntegrationPoint ipr = this->ip;
	IntegrationPoint ipl = this->ip;
	ipr(dir) += eps;
	ipl(dir) -= eps;
	this->eltrans->CalcJacobian (ipr, jacr);    
	this->eltrans->CalcJacobian (ipl, jacl);    

	for (int j = 0; j < 2; j++)
	  {
	    ddx1(dir,j) = (jacr(0,j) - jacl(0,j) ) / (2*eps);
	    ddx2(dir,j) = (jacr(1,j) - jacl(1,j) ) / (2*eps);
	  }
      }
  }

  template <int S, int R, typename SCAL>
  void MappedIntegrationPoint<S,R,SCAL> ::
  CalcHesse (Mat<1> & ddx1, Mat<1> & ddx2) const
  {
    //if (S == R)
    //  {
    //    this -> CalcHesse(ddx1, ddx2);
    //    return;
    //  }

    // for 1D -> 2D

    double eps = 1e-6;

    Mat<2,1> jacr, jacl;
    for (int dir = 0; dir < 1; dir++)
      {
	IntegrationPoint ipr = this->ip;
	IntegrationPoint ipl = this->ip;
	ipr(dir) += eps;
	ipl(dir) -= eps;
	this->eltrans->CalcJacobian (ipr, jacr);
	this->eltrans->CalcJacobian (ipl, jacl);

	for (int j = 0; j < 1; j++)
	  {
	    ddx1(dir,j) = (jacr(0,j) - jacl(0,j) ) / (2*eps);
	    ddx2(dir,j) = (jacr(1,j) - jacl(1,j) ) / (2*eps);
	  }
      }
  }


  template <int S, int R, typename SCAL>
  void MappedIntegrationPoint<S,R,SCAL> ::
  CalcHesse (Mat<1> & ddx1, Mat<1> & ddx2, Mat<1> & ddx3) const
  {
    //if (S == R)
    //  {
    //    this -> CalcHesse(ddx1, ddx2);
    //    return;
    //  }

    // for 1D -> 3D

    double eps = 1e-6;

    Mat<3,1> jacr, jacl;
    for (int dir = 0; dir < 1; dir++)
      {
	IntegrationPoint ipr = this->ip;
	IntegrationPoint ipl = this->ip;
	ipr(dir) += eps;
	ipl(dir) -= eps;
	this->eltrans->CalcJacobian (ipr, jacr);
	this->eltrans->CalcJacobian (ipl, jacl);

	for (int j = 0; j < 1; j++)
	  {
	    ddx1(dir,j) = (jacr(0,j) - jacl(0,j) ) / (2*eps);
	    ddx2(dir,j) = (jacr(1,j) - jacl(1,j) ) / (2*eps);
            ddx3(dir,j) = (jacr(2,j) - jacl(2,j) ) / (2*eps);
	  }
      }
  }




  template <int S, int R, typename SCAL>
  void MappedIntegrationPoint<S,R,SCAL> :: 
  CalcHesse (Mat<3> & ddx1, Mat<3> & ddx2, Mat<3> & ddx3) const
  {
    double eps = 1e-6;

    Mat<3> jacr, jacl;
    for (int dir = 0; dir < 3; dir++)
      {
	IntegrationPoint ipr = this->IP();
	IntegrationPoint ipl = this->IP();
	ipr(dir) += eps;
	ipl(dir) -= eps;
	this->eltrans->CalcJacobian (ipr, jacr);    
	this->eltrans->CalcJacobian (ipl, jacl);    
	
	for (int j = 0; j < 3; j++)
	  {
	    ddx1(dir,j) = (jacr(0,j) - jacl(0,j) ) / (2*eps);
	    ddx2(dir,j) = (jacr(1,j) - jacl(1,j) ) / (2*eps);
	    ddx3(dir,j) = (jacr(2,j) - jacl(2,j) ) / (2*eps);
	  }
      }
  }

  template <int S, int R, typename SCAL>
  void MappedIntegrationPoint<S,R,SCAL> :: 
  CalcHesse (Mat<2> & ddx1, Mat<2> & ddx2, Mat<2> & ddx3) const
  {
    if (S == R)
      {
        this -> CalcHesse(ddx1, ddx2);
        return;
      }

    // for 2D -> 3D
        
    double eps = 1e-6;

    Mat<3,2> jacr, jacl;
    for (int dir = 0; dir < 2; dir++)
      {
	IntegrationPoint ipr = this->ip;
	IntegrationPoint ipl = this->ip;
	ipr(dir) += eps;
	ipl(dir) -= eps;
	this->eltrans->CalcJacobian (ipr, jacr);    
	this->eltrans->CalcJacobian (ipl, jacl);    

	for (int j = 0; j < 2; j++)
	  {
	    ddx1(dir,j) = (jacr(0,j) - jacl(0,j) ) / (2*eps);
	    ddx2(dir,j) = (jacr(1,j) - jacl(1,j) ) / (2*eps);
            ddx3(dir,j) = (jacr(2,j) - jacl(2,j) ) / (2*eps);
	  }
      }

  }

  template <int DIMS, int DIMR, typename SCAL>
    void MappedIntegrationPoint<DIMS,DIMR,SCAL>::
    CalcHesse (Vec<DIMR,Mat<DIMS,DIMS>> & ddx) const
  {
    double eps = 1e-6;
    Mat<DIMR,DIMS> jacr, jacl;
    for (int dir = 0; dir < DIMS; dir++)
      {
        IntegrationPoint ipr = this->IP();
        IntegrationPoint ipl = this->IP();
        ipr(dir) += eps;
        ipl(dir) -= eps;
        this->eltrans->CalcJacobian (ipr, jacr);    
        this->eltrans->CalcJacobian (ipl, jacl);    
        
        for (int j = 0; j < DIMS; j++)
          for (int k = 0; k < DIMR; k++)
            ddx(k)(dir,j) = (jacr(k,j) - jacl(k,j) ) / (2*eps);
      }
  }



  

  template class MappedIntegrationPoint<0,0>;
  template class MappedIntegrationPoint<0,1>;
  template class MappedIntegrationPoint<0,2>;
  template class MappedIntegrationPoint<0,3>;
  template class MappedIntegrationPoint<1,1>;
  template class MappedIntegrationPoint<1,2>;
  template class MappedIntegrationPoint<2,2>;
  template class MappedIntegrationPoint<1,3>;
  template class MappedIntegrationPoint<2,3>;
  template class MappedIntegrationPoint<3,3>;
  



 
  /*
  IntegrationRule :: IntegrationRule (ELEMENT_TYPE eltype, int order)
  { 
    const IntegrationRule & ir = SelectIntegrationRule (eltype, order);
    size = ir.Size();
    data = &ir[0];
    mem_to_delete = NULL;
  }
  */

  IntegrationRule :: IntegrationRule (size_t asize, double (*pts)[3], double * weights)
  {
    for (size_t i = 0; i < asize; i++)
      {
	IntegrationPoint ip(pts[i], weights[i]);
	ip.SetNr(i);
	AddIntegrationPoint (ip);
      }
    // mem_to_delete = NULL;
  }

  IntegrationRule IntegrationRule :: Copy() const
  {
    IntegrationRule ir2;
    ir2.SetSize (Size());
    for (size_t i = 0; i < Size(); i++)
      ir2[i] = (*this)[i];
    ir2.dimension = dimension;
    return ir2;
  }

  ostream & operator<< (ostream & ost, const IntegrationRule & ir)
  {
    /*
    for (size_t i = 0; i < ir.GetNIP(); i++)
      ost << ir[i] << endl;
    */
    for (auto & ip : ir)
      ost << ip << endl;      
    return ost;
  }

  ostream & operator<< (ostream & ost, const BaseMappedIntegrationRule & mir)
  {
    /*
    for (size_t i = 0; i < mir.Size(); i++)
      ost << mir[i] << endl;
    */
    for (auto & mip : mir)
      ost << mip << endl;
    return ost;
  }


  int BaseMappedIntegrationRule :: DimElement() const
  {
    return eltrans.ElementDim();
  }

  int BaseMappedIntegrationRule :: DimSpace() const
  {
    return eltrans.SpaceDim();
  }

  

  template <int DIM_ELEMENT, int DIM_SPACE, typename SCAL>
  MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE,SCAL> :: 
  MappedIntegrationRule (const IntegrationRule & ir, 
			 const ElementTransformation & aeltrans, 
			 Allocator & lh)
    : BaseMappedIntegrationRule (ir, aeltrans), mips(ir.GetNIP(), lh)
  {
    baseip = (char*)(void*)(BaseMappedIntegrationPoint*)(&mips[0]);
    incr = sizeof (MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE, SCAL>);

    for (int i = 0; i < ir.GetNIP(); i++)
      new (&mips[i]) MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE, SCAL> (ir[i], eltrans, -1);

    eltrans.CalcMultiPointJacobian (ir, *this);

    if (ir.Size())
      if (ir[0].VB() != VOL)
        ComputeNormalsAndMeasure (eltrans.GetElementType(), ir[0].FacetNr());
  }

  template <int DIM_ELEMENT, int DIM_SPACE, typename SCAL>
  void MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE,SCAL> :: 
  ComputeNormalsAndMeasure (ELEMENT_TYPE et, int facetnr)
  {
    auto hmips = mips;
    if (hmips.Size() == 0) return;

    if ((Dim(et) >= 2) && (Dim(et)-int(hmips[0].IP().VB()) == 0))
      { // pnts in dim >= 2
        for (size_t i = 0; i < hmips.Size(); i++)
          hmips[i].SetMeasure(1);
        return;
      }
    if constexpr(DIM_ELEMENT == 3)
      if (hmips[0].IP().VB() == BBND && Dim(et) == 3)
      {
        // throw Exception ("ComputeNormalsAndMeasure not yet available for volume-edges");
	FlatVector<Vec<3>> points(99,(Vec<3>*)ElementTopology::GetVertices (et));
        auto edges = ElementTopology::GetEdges (et);
        Vec<3> tau_ref = points(edges[facetnr][1]) - points(edges[facetnr][0]);
        for (int i = 0; i < hmips.Size(); i++)
          {
            auto & mip = hmips[i];        
            Vec<3> tau = mip.GetJacobian() * tau_ref;
            mip.SetMeasure(L2Norm(tau));
            tau /= L2Norm(tau);
            mip.SetTV(tau);
          }
        return;
      }
      
    
    Vec<DIM_ELEMENT> normal_ref = ElementTopology::GetNormals<DIM_ELEMENT>(et)[facetnr];
    for (int i = 0; i < hmips.Size(); i++)
      {
        auto & mip = hmips[i];
        // Mat<DIM_ELEMENT,DIM_SPACE,SCAL> inv_jac = mip.GetJacobianInverse();
        // double det = fabs (mip.GetJacobiDet()); // GetMeasure();
        // Vec<DIM_SPACE> normal = det * Trans (inv_jac) * normal_ref;
        Vec<DIM_SPACE> normal = mip.GetJacobianCofactor() * normal_ref;
        double sign = mip.GetJacobiDet() > 0 ? 1 : -1;
        double len = L2Norm (normal);       // that's the surface measure
        normal *= sign/len;                  // normal vector on physical element
        mip.SetMeasure (len);

        if constexpr (DIM_ELEMENT == DIM_SPACE)
          {
            mip.SetNV(normal);
            if (DIM_SPACE == 2)
              mip.SetTV(Vec<2>(-normal(1),normal(0)));
          }
        else
          {
            if constexpr (DIM_SPACE == 3 && DIM_ELEMENT == 1)
              {
                mip.SetTV(Vec<3>(normal));
              }
            else if constexpr ( DIM_SPACE == 3 && DIM_ELEMENT == 2)
              {
                Vec<3> tang = Cross(Vec<3> (normal), Vec<3> (mip.GetNV()));
                mip.SetTV(tang);
              }
            else if constexpr ( DIM_SPACE == 2 && DIM_ELEMENT == 1)
              {
                mip.SetTV(Vec<2>(normal));
              }
            else if constexpr ( DIM_ELEMENT == 0 )
               ;                 
            else
              {
                throw Exception("In ComputeNormalsAndMeasure: not implemented yet!");
              }
            
          }
      }
  }


  template <int DIM_ELEMENT, int DIM_SPACE, typename SCAL>
  void MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE,SCAL> :: 
  ComputeNormalsAndMeasure (ELEMENT_TYPE et)
  {
    auto hmips = mips;
    if (hmips.Size() == 0) return;

    if ((Dim(et) >= 2) && (Dim(et)-int(hmips[0].IP().VB()) == 0))
      { // pnts in dim >= 2
        for (size_t i = 0; i < hmips.Size(); i++)
          hmips[i].SetMeasure(1);
        return;
      }
    if constexpr(DIM_ELEMENT == 3)
      if (hmips[0].IP().VB() == BBND && Dim(et) == 3)
      {
        // throw Exception ("ComputeNormalsAndMeasure not yet available for volume-edges");
	FlatVector<Vec<3>> points(99,(Vec<3>*)ElementTopology::GetVertices (et));
        auto edges = ElementTopology::GetEdges (et);
        for (int i = 0; i < hmips.Size(); i++)
          {
            auto & mip = hmips[i];        
            int facetnr = mip.IP().FacetNr();
            Vec<3> tau_ref = points(edges[facetnr][1]) - points(edges[facetnr][0]);
            Vec<3> tau = mip.GetJacobian() * tau_ref;
            mip.SetMeasure(L2Norm(tau));
            tau /= L2Norm(tau);
            mip.SetTV(tau);
          }
        return;
      }
      
    
    for (int i = 0; i < hmips.Size(); i++)
      {
        auto & mip = hmips[i];
        int facetnr = mip.IP().FacetNr();
        Vec<DIM_ELEMENT> normal_ref = ElementTopology::GetNormals<DIM_ELEMENT>(et)[facetnr];
        
        Mat<DIM_ELEMENT,DIM_SPACE,SCAL> inv_jac = mip.GetJacobianInverse();
        double det = fabs (mip.GetJacobiDet()); // GetMeasure();
        Vec<DIM_SPACE> normal = det * Trans (inv_jac) * normal_ref;
        double len = L2Norm (normal);       // that's the surface measure
        normal *= 1.0/len;                  // normal vector on physical element
        mip.SetMeasure (len);

        if constexpr (DIM_ELEMENT == DIM_SPACE)
          {
            mip.SetNV(normal);
            if (DIM_SPACE == 2)
              mip.SetTV(Vec<2>(-normal(1),normal(0)));
          }
        else
          {
            if constexpr (DIM_SPACE == 3 && DIM_ELEMENT == 1)
              {
                mip.SetTV(Vec<3>(normal));
              }
            else if constexpr ( DIM_SPACE == 3 && DIM_ELEMENT == 2)
              {
                Vec<3> tang = Cross(Vec<3> (normal), Vec<3> (mip.GetNV()));
                mip.SetTV(tang);
              }
            else if constexpr ( DIM_SPACE == 2 && DIM_ELEMENT == 1)
              {
                mip.SetTV(Vec<2>(normal));
              }
            else if constexpr ( DIM_ELEMENT == 0 )
               ;                 
            else
              {
                throw Exception("In ComputeNormalsAndMeasure: not implemented yet!");
              }
            
          }
      }
  }


  
  template <int DIM_ELEMENT, int DIM_SPACE>
  MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE,Complex> :: 
  MappedIntegrationRule (const IntegrationRule & ir, 
			 const ElementTransformation & aeltrans, 
			 Allocator & lh)
    : BaseMappedIntegrationRule (ir, aeltrans), mips(ir.GetNIP(), lh)
  {
    baseip = (char*)(void*)(BaseMappedIntegrationPoint*)(&mips[0]);
    incr = sizeof (MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE, SCAL>);

    for (int i = 0; i < ir.GetNIP(); i++)
      new (&mips[i]) MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE, Complex> (ir[i], eltrans, -1);
    eltrans.CalcMultiPointJacobian (ir, *this);

    if (ir.Size())
      if (ir[0].VB() != VOL)
        ComputeNormalsAndMeasure (eltrans.GetElementType(), ir[0].FacetNr());
  }

  template <int DIM_ELEMENT, int DIM_SPACE>
  void MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE,Complex> :: 
  ComputeNormalsAndMeasure (ELEMENT_TYPE et, int facetnr)
  {
    throw Exception ("ComputeNormalsAndMeasure not available for complex IR");
  }
  
  template class MappedIntegrationRule<0,0>;
  template class MappedIntegrationRule<0,1>;
  template class MappedIntegrationRule<0,2>;
  template class MappedIntegrationRule<0,3>;
  template class MappedIntegrationRule<1,1>;
  template class MappedIntegrationRule<2,2>;
  template class MappedIntegrationRule<3,3>;
  template class MappedIntegrationRule<1,2>;
  template class MappedIntegrationRule<1,3>;
  template class MappedIntegrationRule<2,3>;

  template class MappedIntegrationRule<0,0, Complex>;
  template class MappedIntegrationRule<1,1, Complex>;
  template class MappedIntegrationRule<2,2, Complex>;
  template class MappedIntegrationRule<3,3, Complex>;


  template class MappedIntegrationRule<0,1, Complex>;
  template class MappedIntegrationRule<0,2, Complex>;
  template class MappedIntegrationRule<0,3, Complex>;
  template class MappedIntegrationRule<1,2, Complex>;
  template class MappedIntegrationRule<1,3, Complex>;
  template class MappedIntegrationRule<2,3, Complex>;


  
  // computes Gaussean integration formula on (0,1) with n points
  // in: Numerical algs in C (or, was it the Fortran book ?)
  void ComputeGaussRule (int n, 
			 Array<double> & xi, Array<double> & wi)
  {
    //  cout << "compute gauss rule, n = " << n << endl;
    xi.SetSize (n);
    wi.SetSize (n);
    
    int m = (n+1)/2;
    double p1, p2, p3;
    double pp, z, z1;
    for (int i = 1; i <= m; i++)
      {
	z = cos ( M_PI * (i - 0.25) / (n + 0.5));
	
	while(1)
	  {
	    p1 = 1;
	    p2 = 0;
	    for (int j = 1; j <= n; j++)
	      {
		p3 = p2;
		p2 = p1;
		p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
	      }
	    // p1 is legendre polynomial
	    
	    pp = n * (z*p1-p2) / (z*z - 1);
	    z1 = z;
	    z = z1-p1/pp;
	    
	    if (fabs (z - z1) < 1e-14) break;
	  }
	
	xi[i-1] = 0.5 * (1 - z);
	xi[n-i] = 0.5 * (1 + z);
	wi[i-1] = wi[n-i] = 1.0 / ( (1  - z * z) * pp * pp);
      }
    //  (*testout) << "Gauss points with n = " << n << ":" << endl;
    //  for (i = 1; i <= n; i++)
    //    (*testout) << xi.Elem(i) << ",  w= " << wi.Elem(i) << endl;
  }





  double  gammln(double xx)
  {
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
			  -1.231739516,0.120858003e-2,-0.536382e-5};
    
    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (int j=0;j<=5;j++) {
      x += 1.0;
      ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
  }



  const double  EPS = 3.0e-14;
  const int  MAXIT = 10;

  void ComputeGaussJacobiRule (int n, 
			       Array<double> & x, 
			       Array<double> & w,
			       double alf, double bet)
  {
    x.SetSize (n);
    w.SetSize (n);

    int i,its,j;
    double alfbet,an,bn,r1,r2,r3;
    double a,b,c,p1,p2,p3,pp,temp,z=0,z1;
    
    for (i=1;i<=n;i++) {
      if (i == 1) {
	an=alf/n;
	bn=bet/n;
	r1=(1.0+alf)*(2.78/(4.0+n*n)+0.768*an/n);
	r2=1.0+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
	z=1.0-r1/r2;
      } else if (i == 2) {
	r1=(4.1+alf)/((1.0+alf)*(1.0+0.156*alf));
	r2=1.0+0.06*(n-8.0)*(1.0+0.12*alf)/n;
	r3=1.0+0.012*bet*(1.0+0.25*fabs(alf))/n;
	z -= (1.0-z)*r1*r2*r3;
      } else if (i == 3) {
	r1=(1.67+0.28*alf)/(1.0+0.37*alf);
	r2=1.0+0.22*(n-8.0)/n;
	r3=1.0+8.0*bet/((6.28+bet)*n*n);
	z -= (x[0]-z)*r1*r2*r3;
      } else if (i == n-1) {
	r1=(1.0+0.235*bet)/(0.766+0.119*bet);
	r2=1.0/(1.0+0.639*(n-4.0)/(1.0+0.71*(n-4.0)));
	r3=1.0/(1.0+20.0*alf/((7.5+alf)*n*n));
	z += (z-x[n-4])*r1*r2*r3;
      } else if (i == n) {
	r1=(1.0+0.37*bet)/(1.67+0.28*bet);
	r2=1.0/(1.0+0.22*(n-8.0)/n);
	r3=1.0/(1.0+8.0*alf/((6.28+alf)*n*n));
	z += (z-x[n-3])*r1*r2*r3;
      } else {
	z=3.0*x[i-2]-3.0*x[i-3]+x[i-4];
      }
      alfbet=alf+bet;
      for (its=1;its<=MAXIT;its++) {
	temp=2.0+alfbet;
	p1=(alf-bet+temp*z)/2.0;
	p2=1.0;
	for (j=2;j<=n;j++) {
	  p3=p2;
	  p2=p1;
	  temp=2*j+alfbet;
	  a=2*j*(j+alfbet)*(temp-2.0);
	  b=(temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z);
	  c=2.0*(j-1+alf)*(j-1+bet)*temp;
	  p1=(b*p2-c*p3)/a;
	}
	pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z));
	z1=z;
	z=z1-p1/pp;
	if (fabs(z-z1) <= EPS) break;
      }
      if (its > MAXIT) 
	cout << "too many iterations in gaujac";
      x[i-1]=z;


      /*
      // is not double precision accurate !!!
      w[i-1]=exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0)-
      gammln(n+alfbet+1.0))*temp*pow(2.0,alfbet)/(pp*p2);    
      */

      /*
	w[i-1]=
	(alf+n-1)! * (beta+n-1) ! /  (  (n)!  (n+alpha+beta)! )
	*temp*pow(2.0,alfbet)/(pp*p2);    
	*/

      if (bet == 0.0)
	{
	  w[i-1] = 1.0 / ( (n+alf) * n )
	    *temp*pow(2.0,alfbet)/(pp*p2);    
	}
      else
	w[i-1]=
	  exp(  gammln(alf+n) + gammln(bet+n) - gammln(n+1.0) - gammln(n+alfbet+1.0))
	  *temp*pow(2.0,alfbet)/(pp*p2);    
    }
    
    for (int i = 0; i < n; i++)
      {
	w[i] *= 0.5 * pow(1-x[i],-alf) * pow(1+x[i], -bet);
	x[i] = 0.5 * (x[i]+1);
      }
  }
    


  void ComputeGaussLobattoRule (int n, Array<double>& xi, Array<double>& wi)
  {
    Array<double> axi;
    Array<double> awi;
    
    ComputeGaussJacobiRule(n-2,axi,awi,1, 1);
    xi.SetSize(0);
    wi.SetSize(0);
    double sum = 1;
    for (int i = 0; i < n-2; i++)
      sum -= awi[i];

    xi.Append(0);
    wi.Append(sum/2);

    xi.Append (axi);
    wi.Append (awi);

    xi.Append(1);
    wi.Append(sum/2);
  }

  
  void ComputeGaussRadauRule (int n, Array<double>& xi, Array<double>& wi)
  {
    Array<double> axi;
    Array<double> awi;
    
    ComputeGaussJacobiRule(n-1,axi,awi,0, 1);
    xi.SetSize(0);
    wi.SetSize(0);
    double sum = 1;
    for (int i = 0; i < n-1; i++)
      sum -= awi[i];
    
    xi.Append(0);
    wi.Append(sum);
    
    xi.Append (axi);
    wi.Append (awi);
  }



#ifdef LAPACK
  void ComputeHermiteRule (int n, 
			   Array<double> & x,
			   Array<double> & w)
  {
    Matrix<> m(n,n), evecs(n,n);
    m = 0.0;
    for (int i = 0; i < n-1; i++)
      m(i,i+1) = m(i+1, i) = sqrt( (i+1.0)/2);
    
    Vector<> lam(n);
    LapackEigenValuesSymmetric (m, lam, evecs);

    x.SetSize (n);
    w.SetSize (n);

    for (int i = 0; i < n; i++)
      {
	x[i] = lam(i);
	w[i] = sqrt(M_PI) * sqr(evecs(i,0));
      }
  }
#else

  void ComputeHermiteRule (int n, 
                           Array<double> & x,
                           Array<double> & w)
  {
    const double EPS = 3e-14;
    const double PIM4 = 1.0 / pow (M_PI, 0.25);

    x.SetSize (n);
    w.SetSize (n);

    int its, m;
    double p1, p2, p3, pp, z = 0, z1;

    m = (n+1)/2;
    for (int i = 1; i <= m; i++)
      {
        if (i == 1) 
          {
            z = sqrt ( (double) (2*n+1)) - 1.85575*pow( (double)(2*n+1), -0.16667);
          }
        else if (i == 2)
          {
            z -= 1.14*pow( (double)n, 0.426)/z;
          } 
        else if (i == 3)
          {
            z = 1.86*z-0.86*x[0];
          }
        else if (i == 4)
          {
            z = 1.91*z-0.91*x[1];
          } 
        else
          z = 2.0*z-x[i-3];


        for (its = 1; its <= 1000; its++)
          {
            p1 = PIM4;
            p2 = 0;
            for (int j = 1; j <= n; j++)
              {
                p3 = p2;
                p2 = p1;
                p1 = z * sqrt(2.0/j) * p2 - sqrt (((double) (j-1))/j)*p3;
              }
            
            pp = sqrt ( (double)2*n) * p2;
            z1 = z;
            z = z1-p1/pp;
            if (fabs (z-z1) <= EPS) break;

            if (its > 20) cout << "too many steps" << endl;
          }

        x[i-1] = z;
        x[n-i] = -z;
        w[i-1] = 2.0 / (pp*pp);
        w[n-i] = w[i-1];
      }
  }
#endif


  DGIntegrationRule :: DGIntegrationRule(ELEMENT_TYPE eltype, int order)
  {

    switch (eltype)
      {
      case ET_TRIG:
	{
          for (int j = 0; j < 3; j++)
            facetrules.Append (const_cast<IntegrationRule*> (&SelectIntegrationRule (ET_SEGM, order)));


	  int len = floor((order+2.0)/2.0);
	  int lengr = floor((order+5.0)/2.0);
	  Array<double> axi(lengr);
	  Array<double> wi(lengr);
	  Array<double> axib(len);
	  Array<double> wib(len);
		
	  ComputeGaussLobattoRule(lengr, axi , wi);
	  ComputeGaussRule(len,axib,wib);
    
	  for (int i = 0; i < wi.Size(); i++)
	    if (fabs (axi[i]) < 1e-10)
	      boundary_volume_factor = 3/wi[i];
		
	  int ii = 0;
	  for (int j = 0; j < len; j++)
	    for (int k=0; k < lengr;k++)
	      {
                if (fabs (axi[k]-1) < 1e-10) continue;
		IntegrationPoint ip(axib[j]*(1-axi[k]),axi[k],0,wib[j]*wi[k]*(1-axi[k])/3);
		ip.SetNr(ii);
		// ip.FacetNr() = (k == 0) ? 0 : -1;
		AddIntegrationPoint(ip);
		(*this)[ii].SetFacetNr((k == 0) ? 0 : -1);
		ii++;
	      }
    
	  for (int j=0;j<len;j++)
	    for (int k=0;k<lengr;k++)
	      {
                if (fabs (axi[k]-1) < 1e-10) continue;
		IntegrationPoint ip(axi[k], axib[j]*(1-axi[k]), 0, wib[j]*wi[k]*(1-axi[k])/3);
		ip.SetNr(ii);
		ip.SetFacetNr((k == 0) ? 1 : -1);			
		AddIntegrationPoint(ip);
		(*this)[ii].SetFacetNr((k == 0) ? 1 : -1);
		ii++;
	      }

	  for (int j=0;j<len;j++)
	    for (int k=0;k<lengr;k++)
	      {
                if (fabs (axi[k]-1) < 1e-10) continue;
		IntegrationPoint ip( axib[j]*(1-axi[k]),(1-axib[j])*(1-axi[k]), 0, wib[j]*wi[k]*(1-axi[k])/3);
		ip.SetNr(ii);
		ip.SetFacetNr((k == 0) ? 2 : -1);
		AddIntegrationPoint(ip);
		(*this)[ii].SetFacetNr((k == 0) ? 2 : -1);
		ii++;
	      }
	  break;
	}
      case ET_QUAD:
	{
          for (int j = 0; j < 4; j++)
            facetrules.Append (const_cast<IntegrationRule*> (&SelectIntegrationRule (ET_SEGM, order)));


	  int lengl = order/2+2;
	  int len   = order/2+1;
	  Array<double> xi(lengl);
	  Array<double> wi(lengl);
	  Array<double> xig(len);
	  Array<double> wig(len);

	  ComputeGaussLobattoRule(lengl, xi , wi);
	  ComputeGaussRule(len, xig , wig);

	  boundary_volume_factor = 2/wi[0];

	  for (int i = 0, ii = 0; i < len; i++)
	    for (int j = 0; j < lengl; j++) 
	      {
		IntegrationPoint ip1(xig[i],xi[j], 0, 0.5*wig[i]*wi[j]);
		ip1.SetNr(ii); ii++;
		AddIntegrationPoint (ip1);
		IntegrationPoint ip2(xi[j],xig[i], 0, 0.5*wi[j]*wig[i]);
		ip2.SetNr(ii); ii++;
		AddIntegrationPoint (ip2);
	      }
          break;
	}
    

      case ET_TET:
	{

          /*
            o = 0 -> p2 -> 3 p
            o = 1 -> p3 -> 3 p
            o = 2 -> p4 -> 4 p
           */
          Facet2ElementTrafo trafo(ET_TET);
          const POINT3D * pnts = ElementTopology::GetVertices(ET_TET);

          /*
	  int lengr = floor((order+6.0)/2.0);
	  Array<double> axi(lengr);
	  Array<double> wi(lengr);
	  ComputeGaussLobattoRule(lengr, axi , wi);
          */

          int lengr = floor((order+3.0)/2.0);
	  Array<double> axi, wi;
          ComputeGaussJacobiRule(lengr,axi,wi,1,0); 
          double sum = 1.0/2;
          for (int i = 0; i < wi.Size(); i++) sum -= axi[i]*wi[i];
          wi.Append(sum);
          axi.Append(1.0);

          /*
          int lengr = floor((order+1.0)/2.0);
	  Array<double> axi, wi;
          ComputeGaussJacobiRule(lengr,axi,wi,1,2); 
          double sum = 1.0/3;
          for (int i = 0; i < wi.Size(); i++) sum -= axi[i]*axi[i]*wi[i];
          wi.Append(sum);
          axi.Append(1.0);
          */

          int ii = 0;
          for (int j = 0; j < 4; j++)
            {
              facetrules.Append (const_cast<IntegrationRule*> (&SelectIntegrationRule (ET_TRIG, order)));
              IntegrationRule & ir = *facetrules.Last();
              Vec<3> v (pnts[j][0], pnts[j][1], pnts[j][2]);

              for (int k = 0; k < ir.Size(); k++)
                {
                  IntegrationPoint ip3d = trafo(j, ir[k]);
                  for (int l = 0; l < axi.Size(); l++)
                    {
                      Vec<3> pf = ip3d.Point();
                      Vec<3> p = axi[l] * pf + (1-axi[l]) * v;
                      double weight = wi[l]*sqr(axi[l])*ir[k].Weight()/4;
                      if (weight > 1e-10)
                        {
                          IntegrationPoint ip (p, weight);
                          ip.SetNr(ii);
                          ii++;
                          ip.SetFacetNr(-1);
                          if (axi[l] > 1-1e-10) ip.SetFacetNr(j);
                          AddIntegrationPoint (ip);
                        }
                      
                      if (axi[l] > 1-1e-10)
                        boundary_volume_factor = 4/wi[l];
                    }
                }
            }              
              
          break;
        }

      default:
        throw Exception (string("Sorry, no DG intrules for element type ")
                         +ElementTopology::GetElementName(eltype));
      }
  }


  ostream & operator<< (ostream & ost, const DGIntegrationRule & ir)
  {
    ost << "DG-IntegrationRule" << endl;
    ost << "vol-ir: " << endl << (const IntegrationRule&)(ir);
    for (int i = 0; i < ir.GetNFacets(); i++)
      ost << "facet " << i << ": " << endl << ir.GetFacetIntegrationRule(i) << endl;
    ost << "bound-vol-factor = " << ir.BoundaryVolumeFactor() << endl;
    return ost;
  }



#ifdef OLD

  TensorProductIntegrationRule ::   
  TensorProductIntegrationRule (const FiniteElement & el, const IntegrationRule & aseg_rule)
    : type(el.ElementType()), seg_rule(aseg_rule) 
  {
    useho = 0;
    nx = seg_rule.GetNIP();


    if (el.ElementType() == ET_TET)
      {
        const H1HighOrderFiniteElement<3> * h1ho = 
          dynamic_cast<const H1HighOrderFiniteElement<3>*> (&el);

	n = nx*nx*nx;

	if (h1ho)
	  {
	    useho = 1;
	    
	    for (int i = 0; i < 4; i++)
	      vnums[i] = h1ho->GetVNums()[i];
	    
	    for (int i = 0; i < 4; i++)
	      sort[i] = i;
	    
	    for (int i = 0; i < 4; i++)
	      for (int j = 0; j < 3; j++)
		if (vnums[sort[j]] < vnums[sort[j+1]])
		  swap (sort[j], sort[j+1]);
	    
	    for (int i = 0; i < 4; i++)
	      isort[sort[i]] = i;
	    
	    const double refpts[][3] =
	      { { 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 1 },
		{ 0, 0, 0 } };
	    
	    for (int i = 0; i < 4; i++)
	      for (int j = 0; j < 3; j++)
		vertices[i](j) = refpts[sort[i]][j];
	  }
      }

    if (el.ElementType() == ET_TRIG)
      {
        const H1HighOrderFiniteElement<2> * h1ho = 
          dynamic_cast<const H1HighOrderFiniteElement<2>*> (&el);

	n = nx * nx;

	if (h1ho)
	  {
	    useho = 1;
	    
	    for (int i = 0; i < 3; i++)
	      vnums[i] = h1ho->GetVNums()[i];
	    
	    for (int i = 0; i < 3; i++)
	      sort[i] = i;
	
	    for (int i = 0; i < 3; i++)
	      for (int j = 0; j < 2; j++)
		if (vnums[sort[j]] < vnums[sort[j+1]])
		  swap (sort[j], sort[j+1]);
	    
	    for (int i = 0; i < 3; i++)
	      isort[sort[i]] = i;

	    const double refpts[][3] =
	      { { 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 0 } };

	    for (int i = 0; i < 3; i++)
	      for (int j = 0; j < 3; j++)
		vertices[i](j) = refpts[sort[i]][j];
	  }
      }

  }








  TensorProductIntegrationRule ::   
  TensorProductIntegrationRule (ELEMENT_TYPE atype, const int * avnums, const IntegrationRule & aseg_rule)
    : type(atype), seg_rule(aseg_rule) 
  {
    useho = 1;
    nx = seg_rule.GetNIP();

    if (type == ET_TET)
      {
	n = nx * nx * nx;
	useho = 1;

	for (int i = 0; i < 4; i++)
	  vnums[i] = avnums[i];

	for (int i = 0; i < 4; i++)
	  sort[i] = i;
	
	for (int i = 0; i < 4; i++)
	  for (int j = 0; j < 3; j++)
	    if (vnums[sort[j]] < vnums[sort[j+1]])
	      swap (sort[j], sort[j+1]);

	for (int i = 0; i < 4; i++)
	  isort[sort[i]] = i;

	const double refpts[][3] =
	  { { 1, 0, 0 },
	    { 0, 1, 0 },
	    { 0, 0, 1 },
	    { 0, 0, 0 } };

	for (int i = 0; i < 4; i++)
	  for (int j = 0; j < 3; j++)
	    vertices[i](j) = refpts[sort[i]][j];
      }

    if (type == ET_TRIG)
      {
	n = nx * nx;
	useho = 1;

	for (int i = 0; i < 3; i++)
	  vnums[i] = avnums[i];

	for (int i = 0; i < 3; i++)
	  sort[i] = i;
	
	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 2; j++)
	    if (vnums[sort[j]] < vnums[sort[j+1]])
	      swap (sort[j], sort[j+1]);

	for (int i = 0; i < 3; i++)
	  isort[sort[i]] = i;

	const double refpts[][3] =
	  { { 1, 0, 0 },
	    { 0, 1, 0 },
	    { 0, 0, 0 } };

	for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 3; j++)
	    vertices[i](j) = refpts[sort[i]][j];
      }


  }


  IntegrationPoint 
  TensorProductIntegrationRule :: GetIP(int i) const
  {
    int ix = i % nx;
    i /= nx;
    int iy = i % nx;
    i /= nx;
    int iz = i;

    return GetIP (ix, iy, iz);
  }


  IntegrationPoint 
  TensorProductIntegrationRule :: GetIP(int ix, int iy, int iz) const
  {
    IntegrationPoint ip;
    switch (type)
      {
      case ET_TET:
	{
	  double xi = seg_rule[ix](0);
	  double eta = seg_rule[iy](0);
	  double zeta = seg_rule[iz](0);

	  if (useho)
	    {
	      double lami[4];
	      
	      lami[3] = zeta;
	      lami[2] = eta * (1-zeta);
	      lami[1] = xi * (1-eta)*(1-zeta);
	      lami[0] = 1-lami[1]-lami[2]-lami[3];

	      for (int j = 0; j < 3; j++)
		ip(j) = lami[isort[j]];
	    }
	  else
	    {
	      ip(0) = xi * (1-eta)*(1-zeta);
	      ip(1) = eta * (1-zeta);
	      ip(2) = zeta;
	    }

	  ip.SetWeight(seg_rule[ix].Weight() * (1-eta)*(1-zeta) *
		       seg_rule[iy].Weight() * (1-zeta) *
		       seg_rule[iz].Weight());
	  
	  ip.SetNr (ix + nx * iy + nx*nx * iz);
	  break;
	}


      case ET_TRIG:
	{
	  double xi = seg_rule[ix](0);
	  double eta = seg_rule[iy](0);

	  if (useho)
	    {
	      double lami[3];
	      
	      lami[2] = eta;
	      lami[1] = xi * (1-eta);
	      lami[0] = 1-lami[1]-lami[2];

	      for (int j = 0; j < 2; j++)
		ip(j) = lami[isort[j]];
	    }
	  else
	    {
	      ip(0) = xi * (1-eta);
	      ip(1) = eta;
	    }

	  ip.SetWeight(seg_rule[ix].Weight() * (1-eta) *
		       seg_rule[iy].Weight());

	  break;
	}
      };

    // ip.SetPrecomputedGeometry(true);
    return ip;
  }



  void
  TensorProductIntegrationRule :: GetIP_DTet_DHex(int i, Mat<3> & dtet_dhex) const
  {
    int ix = i % nx;
    i /= nx;
    int iy = i % nx;
    i /= nx;
    int iz = i;

    switch (type)
      {
      case ET_TET:
	{
	  double xi = seg_rule[ix](0);
	  double eta = seg_rule[iy](0);
	  double zeta = seg_rule[iz](0);

	  if (useho)
	    {
	      //	      double lami[4];
	      double dlami[4][4];
	      
	      /*
		lami[3] = zeta;
		lami[2] = eta * (1-zeta);
		lami[1] = xi * (1-eta)*(1-zeta);
		lami[0] = 1-lami[1]-lami[2]-lami[3];
	      */

	      dlami[3][0] = 0;
	      dlami[3][1] = 0;
	      dlami[3][2] = 1;
	      
	      dlami[2][0] = 0;
	      dlami[2][1] = 1-zeta;
	      dlami[2][2] = -eta;
	      
	      dlami[1][0] = (1-eta)*(1-zeta);
	      dlami[1][1] = -xi * (1-zeta);
	      dlami[1][2] = -xi * (1-eta);
	      
	      dlami[0][0] = -(1-eta)*(1-zeta);
	      dlami[0][1] = -(1-xi) * (1-zeta);
	      dlami[0][2] = -(1-xi) * (1-eta);
	      
	      
	      for (int j = 0; j < 3; j++)
		for (int k = 0; k < 3; k++)
		  dtet_dhex(j,k) = dlami[isort[j]][k];
	    }

	  break;
	}
      };
    // return ip;    
  }


  void
  TensorProductIntegrationRule :: GetIP_DTrig_DQuad(int i, Mat<2> & dtrig_dquad) const
  {
    int ix = i % nx;
    i /= nx;
    int iy = i;

    switch (type)
      {
      case ET_TRIG:
	{
	  double xi = seg_rule[ix](0);
	  double eta = seg_rule[iy](0);

	  if (useho)
	    {
	      // 	      double lami[3];
	      // 	      lami[2] = eta;
	      // 	      lami[1] = xi * (1-eta);
	      // 	      lami[0] = 1-lami[1]-lami[2];

	      double dlami[3][2];
	      
	      dlami[2][0] = 0;
	      dlami[2][1] = 1;
	      
	      dlami[1][0] = (1-eta);
	      dlami[1][1] = -xi;
	      
	      dlami[0][0] = -(1-eta);
	      dlami[0][1] = -(1-xi);
	      
	      for (int j = 0; j < 2; j++)
		for (int k = 0; k < 2; k++)
		  dtrig_dquad(j,k) = dlami[isort[j]][k];
	    }

	  break;
	}
      };
  }


  void TensorProductIntegrationRule :: GetWeights (FlatArray<double> weights)
  {
    switch (type)
      {
      case ET_TET:
	{
	  for (int iz = 0, i = 0; iz < nx; iz++)
	    for (int iy = 0; iy < nx; iy++)
	      {
		double eta = seg_rule[iy](0);
		double zeta = seg_rule[iz](0);
		
		double fac_yz = 
		  (1-eta)*(1-zeta) *
		  seg_rule[iy].Weight() * (1-zeta) *
		  seg_rule[iz].Weight();
		
		for (int ix = 0; ix < nx; ix++, i++)
		  weights[i] = seg_rule[ix].Weight() * fac_yz;
	      }
	  break;
	}
      case ET_TRIG:
	{
	  for (int iy = 0, i = 0; iy < nx; iy++)
	    {
	      double eta = seg_rule[iy](0);
	      double fac_y = (1-eta)* seg_rule[iy].Weight();
	      for (int ix = 0; ix < nx; ix++, i++)
		weights[i] = seg_rule[ix].Weight() * fac_y;
	    }
	  break;
	}
      }
  }
#endif

  




  template <>
  IntegrationRuleTP<1> :: IntegrationRuleTP (const ElementTransformation & eltrans,
                                             INT<1> vorder, bool compute_duffy, bool compute_points) 
  {
    int order = vorder[0];
    irx = &SelectIntegrationRule (ET_SEGM, order);
    int nip = irx->GetNIP();

    SetSize(nip);
    dxdxi_duffy.SetSize(nip);

    for (int i = 0; i < irx->GetNIP(); i++)
      {
        (*this)[i] = IntegrationPoint ((*irx)[i](0), 0, 0, 
                                       (*irx)[i].Weight());
      }
        
    for (int i = 0; i < nip; i++)
      dxdxi_duffy[i] = 1;
  }


  template <>
  IntegrationRuleTP<2> :: IntegrationRuleTP (const ElementTransformation & eltrans,
                                             INT<2> order, bool compute_duffy, bool compute_points)
  {
    int nip = 0;

    switch (eltrans.GetElementType())
      {
      case ET_TRIG:
        {
          irx = &SelectIntegrationRuleJacobi10 (order[0]);
          iry = &SelectIntegrationRule (ET_SEGM, order[0]);
          
          int sort[3] = { 0, 1, 2 };
          eltrans.GetSort (FlatArray<int> (3, &sort[0]) );
          int isort[3];
          for (int i = 0; i < 3; i++) isort[sort[i]] = i;
          
          nip = irx->GetNIP() * iry->GetNIP();

          SetAllocSize(nip+8);
          SetSize(nip);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++, ii++)
              {
                double 
                  x = (*irx)[i1](0), 
                  y = (*iry)[i2](0);
                
                double lami[] = { x, 
                                  y * (1-x), 
                                  (1-x) * (1-y) };
                
		(*this)[ii] = 
		  IntegrationPoint (lami[isort[0]], 
				    lami[isort[1]], 
				    0,
				    (*irx)[i1].Weight()*(*iry)[i2].Weight()*(1-x));
              }

          // trig permutation transformation
          double dlamdx[3][2] = 
            { { 1, 0 },
              { 0, 1 },
              { -1, -1 }};
        
          Mat<2> trans2;
          for (int k = 0; k < 2; k++)
            for (int l = 0; l < 2; l++)
              trans2(l,k) = dlamdx[sort[l]][k];
        
          // trig permutation plus quad -> trig mapping
          dxdxi_duffy.SetSize(nip);
          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            {
              double x = (*irx)[i1](0);
              double invx = 1.0 / (1-x);
            
              for (int i2 = 0; i2 < iry->GetNIP(); i2++, ii++)
                {
                  double y = (*iry)[i2](0);
                  
                  Mat<2> trans3; // , trans, hmat;
                  
                  // quad -> trig transform
                  trans3(0,0) = 1;
                  trans3(0,1) = 0;
                  
                  trans3(1,0) = y*invx;
                  trans3(1,1) = 1*invx;
                  
                  dxdxi_duffy[ii] = trans3 * trans2;
                }
            }
          break;
        }

      case ET_QUAD:
        {
          irx = &SelectIntegrationRule (ET_SEGM, order[0]);
          iry = &SelectIntegrationRule (ET_SEGM, order[1]);

          if (compute_points)
            {
              nip = irx->GetNIP() * iry->GetNIP();
              SetSize(nip);
              
              int ii = 0;
              for (IntegrationPoint & ipx : *irx)
                for (IntegrationPoint & ipy : *iry)
                  (*this)[ii++] = IntegrationPoint (ipx(0), ipy(0), 0, ipx.Weight()*ipy.Weight());
            }
          
          if (compute_duffy)
            {
              Mat<2> id;
              id = 0;
              id(0,0) = id(1,1) = 1;
              
              dxdxi_duffy.SetSize(nip);
              for (int i = 0; i < nip; i++)
                dxdxi_duffy[i] = id;
            }
          break;
        }
      default:
        {
          throw Exception (ToString("IntegratonRuleTP not available for element type ") +
			   ElementTopology::GetElementName(eltrans.GetElementType()) + "\n");
        }

      }
  }





  template <>
  IntegrationRuleTP<3> :: IntegrationRuleTP (const ElementTransformation & eltrans,
                                             INT<3> vorder, bool compute_duffy, bool compute_points)
  {
    int order = vorder[0];
    int nip = 0;

    switch (eltrans.GetElementType())
      {

      case ET_TET:
        {
          irx = &SelectIntegrationRuleJacobi20 (order);
          iry = &SelectIntegrationRuleJacobi10 (order);
          irz = &SelectIntegrationRule (ET_SEGM, order);

          int sort[4] = { 0, 1, 2, 3}, isort[4];   // compiler warnings
	  eltrans.GetSort (FlatArray<int> (4, sort) );
          for (int i = 0; i < 4; i++) isort[sort[i]] = i;
          
          nip = irx->GetNIP() * iry->GetNIP() * irz->GetNIP();

          SetAllocSize(nip+8);
          SetSize(nip);

          static Timer t ("IntegrationruleTP - compute pts", NoTracing);
          NgProfiler::RegionTimer reg(t);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++)
              for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                {
                  double 
                    x = (*irx)[i1](0), 
                    y = (*iry)[i2](0), 
                    z = (*irz)[i3](0);

                  double lami[] = { x, 
                                    y * (1-x), 
                                    z * (1-x) * (1-y), 
                                    (1-x)*(1-y)*(1-z) };
                
		  (*this)[ii] = IntegrationPoint (lami[isort[0]],
						  lami[isort[1]],
						  lami[isort[2]],
						  (*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight() * 
						  sqr(1-x) * (1-y) );
                }

          // tet permutation transformation
          double dlamdx[4][3] = 
            { { 1, 0, 0 },
              { 0, 1, 0 },
              { 0, 0, 1 },
              { -1, -1, -1 }};
        
          Mat<3> trans2;
          for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
              trans2(l,k) = dlamdx[sort[l]][k];
	  
	  dxdxi_permute = trans2;

          // tet permutation plus hex -> tet mapping
          if (compute_duffy)
            {
              dxdxi_duffy.SetSize(nip);
              for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
                {
                  double x = (*irx)[i1](0);
                  double invx = 1.0 / (1-x);
                  
                  for (int i2 = 0; i2 < iry->GetNIP(); i2++)
                    {
                      double y = (*iry)[i2](0);
                      double invxy = 1.0 / ( (1-x) * (1-y) );
                      
                      for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                        {
                          double z = (*irz)[i3](0);
                          Mat<3> trans3; 
                          
                          // hex -> tet transform
                          trans3(0,0) = 1;
                          trans3(0,1) = 0;
                          trans3(0,2) = 0;
                          
                          trans3(1,0) = y*invx;
                          trans3(1,1) = 1*invx;
                          trans3(1,2) = 0;
                          
                          trans3(2,0) = z*invxy;
                          trans3(2,1) = z*invxy;
                          trans3(2,2) = 1*invxy;
                          
                          dxdxi_duffy[ii] = trans3 * trans2;
                        }
                    }
                }
            }
          break;
        }

      case ET_PRISM:
        {
          irx = &SelectIntegrationRule (ET_SEGM, order);
          iry = &SelectIntegrationRule (ET_SEGM, order+1);
          irz = &SelectIntegrationRule (ET_SEGM, order);

          int sort[6], isort[6];

          eltrans.GetSort (FlatArray<int> (6, sort) );
          for (int i = 0; i < 6; i++) isort[sort[i]] = i;

          nip = irx->GetNIP() * iry->GetNIP() * irz->GetNIP();
	  
	  SetSize(nip);
	  /*
          xi.SetSize(nip);
          weight.SetSize(nip);
	  */
          dxdxi_duffy.SetSize(nip);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++)
              for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                {
                  double 
                    z = (*irx)[i1](0), 
                    x = (*iry)[i2](0), 
                    y = (*irz)[i3](0);
 
                  double lami[] = { x, y *(1-x), (1-x)*(1-y) };
		  
 
		  (*this)[ii] = IntegrationPoint (lami[isort[0]],
						  lami[isort[1]],
						  z,
						  (*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight() * (1-x));
		  /*
                  xi[ii](0) = lami[isort[0]];
                  xi[ii](1) = lami[isort[1]];
                  xi[ii](2) = z;
                  weight[ii] = (*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight() * (1-x);
		  */
                }


          // trig permutation transformation
          double dlamdx[3][2] = 
            { { 1, 0 },
              { 0, 1 },
              { -1, -1 }};
        
          Mat<3> trans2;
          trans2 = 0.0;
          for (int k = 0; k < 2; k++)
            for (int l = 0; l < 2; l++)
              trans2(l,k) = dlamdx[sort[l]][k];
          trans2(2,2) = 1.0;

	  dxdxi_permute = trans2;
        
          // prism permutation plus hex -> tet mapping
          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++)
              for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                {
                  double x = (*iry)[i2](0);
                  double invx = 1.0 / (1-x);
                  double y = (*irz)[i3](0);
                
                  Mat<3> trans3;
                
                  // hex -> prism transform
                  trans3 = 0.0;
                
                  trans3(1,0) = 1;
                  trans3(2,0) = y*invx;
                  trans3(2,1) = 1*invx;
                
                  trans3(0,2) = 1;

                  dxdxi_duffy[ii] = trans3 * trans2;
                }
          break;
        }

      case ET_HEX:
        {
          irx = &SelectIntegrationRule (ET_SEGM, order);
          iry = &SelectIntegrationRule (ET_SEGM, order);
          irz = &SelectIntegrationRule (ET_SEGM, order);

          nip = irx->GetNIP() * iry->GetNIP() * irz->GetNIP();

	  SetSize (nip);
	  /*
          xi.SetSize(nip);
          weight.SetSize(nip);
	  */
          dxdxi_duffy.SetSize(nip);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++)
              for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                {
		  (*this)[ii] = IntegrationPoint ((*irx)[i1](0),
						  (*iry)[i2](0),
						  (*irz)[i3](0),
						  (*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight());

		  /*
		    xi[ii](0) = (*irx)[i1](0);
		    xi[ii](1) = (*iry)[i2](0);
		    xi[ii](2) = (*irz)[i3](0);
		    weight[ii] = (*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight();
		  */
                }
        

	  /*
          Mat<3> id;
          id = 0;
          id(0,0) = id(1,1) = id(2,2) = 1;
	  */
          for (int i = 0; i < nip; i++)
            dxdxi_duffy[i] = Id<3>();
          break;
        }

     case ET_PYRAMID:
        {
          irx = &SelectIntegrationRule (ET_SEGM, order+2); // in z-direction
          irz = &SelectIntegrationRule (ET_SEGM, order);
          iry = &SelectIntegrationRule (ET_SEGM, order);

          nip = irx->GetNIP() * iry->GetNIP() * irz->GetNIP();

	  SetSize (nip);
	  /*
          xi.SetSize(nip);
          weight.SetSize(nip);
	  */
          dxdxi_duffy.SetSize(nip);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++)
              for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                {
		  double z = (*irx)[i1](0);

		  (*this)[ii] = IntegrationPoint ((*iry)[i2](0)*(1-z),
						  (*irz)[i3](0)*(1-z),
						  (*irx)[i1](0),
						  sqr(1-z)*(*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight());
		  /*
                  xi[ii](0) = (*iry)[i2](0)*(1-z);
                  xi[ii](1) = (*irz)[i3](0)*(1-z);
                  xi[ii](2) = (*irx)[i1](0);
                  weight[ii] = sqr(1-z)*(*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight();
		  */
                }
        
          Mat<3> id;
          id = 0;
          id(0,0) = id(1,1) = id(2,2) = 1;
        
          for (int i = 0; i < nip; i++)
            dxdxi_duffy[i] = id;

          break;
        }

      default:
        {
          stringstream str;
          str<< "IntegratonRuleTP not available for element type " 
             << ElementTopology::GetElementName(eltrans.GetElementType()) << endl;
          throw Exception (str.str());
        }

      }

    /*
    if (compute_mapping && !eltrans.Boundary())
      {
        x.SetSize(nip);
        dxdxi.SetSize(nip);
        // eltrans.CalcMultiPointJacobian (*this, x, dxdxi, lh);
	throw Exception ("intruleTP comput_mapping currently not available");
      }
    */
  }














  template <>
  IntegrationRuleTP<1> :: IntegrationRuleTP (ELEMENT_TYPE eltype, FlatArray<int> sort, 
                                             NODE_TYPE nt, int nodenr, int order, LocalHeap & lh)
  { ; }

  template <>
  IntegrationRuleTP<2> :: IntegrationRuleTP (ELEMENT_TYPE eltype, FlatArray<int> sort, 
                                             NODE_TYPE nt, int nodenr, int order, LocalHeap & lh)
  {
    int nip;

    static IntegrationRule ir0, ir1;
    {
      lock_guard<mutex> guard(intruletpfacet_mutex);
      if (ir0.GetNIP() == 0)
	{
	  ir0.AddIntegrationPoint (IntegrationPoint (0.0, 0, 0, 1.0));
	  ir1.AddIntegrationPoint (IntegrationPoint (1.0, 0, 0, 1.0));
	}
    }

    switch (eltype)
      {
      case ET_TRIG:
        {
          if (nt == NT_EDGE)
            {
              // int sort[3];
              // eltrans.GetSort (FlatArray<int> (3, sort) );
              int isort[3];
              for (int i = 0; i < 3; i++) isort[sort[i]] = i;

              const EDGE & edge = ElementTopology::GetEdges (ET_TRIG)[nodenr];
              EDGE sedge;
              sedge[0] = isort[edge[0]];
              sedge[1] = isort[edge[1]];
              if (sedge[0] > sedge[1]) swap (sedge[0], sedge[1]);

              irx = 0; iry = 0;
              if (sedge[0] == 1 && sedge[1] == 2)
                {
                  irx = &ir0;
                  iry = &SelectIntegrationRule (ET_SEGM, order);
                }
                  
              if (sedge[0] == 0 && sedge[1] == 1)
                {
                  irx = &SelectIntegrationRule (ET_SEGM, order);
                  iry = &ir1;
                }

              if (sedge[0] == 0 && sedge[1] == 2)
                {
                  irx = &SelectIntegrationRule (ET_SEGM, order);
                  iry = &ir0;
                }

              nip = irx->GetNIP() * iry->GetNIP();
	      SetSize(nip);

	      /*
              xi.SetSize(nip);
              weight.SetSize(nip);
              */
              for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
                for (int i2 = 0; i2 < iry->GetNIP(); i2++, ii++)
                  {
                    double 
                      x = (*irx)[i1](0), 
                      y = (*iry)[i2](0);
                    
                    double lami[] = { x, 
                                      y * (1-x), 
                                      (1-x) * (1-y) };
                    
		    (*this)[ii] = IntegrationPoint (lami[isort[0]],
						    lami[isort[1]],
						    (*irx)[i1].Weight()*(*iry)[i2].Weight());
		    /*
                    xi[ii](0) = lami[isort[0]];
                    xi[ii](1) = lami[isort[1]];
                    weight[ii] = (*irx)[i1].Weight()*(*iry)[i2].Weight(); 
		    */
                  }
              
              // trig permutation transformation
              double dlamdx[3][2] = 
                { { 1, 0 },
                  { 0, 1 },
                  { -1, -1 }};
              
              Mat<2> trans2;
              for (int k = 0; k < 2; k++)
                for (int l = 0; l < 2; l++)
                  trans2(l,k) = dlamdx[sort[l]][k];
              
              // trig permutation plus quad -> trig mapping
              dxdxi_duffy.SetSize(nip);
              for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
                {
                  double x = (*irx)[i1](0);
                  double invx = 1.0 / (1-x);
                  
                  for (int i2 = 0; i2 < iry->GetNIP(); i2++, ii++)
                    {
                      double y = (*iry)[i2](0);
                      
                      Mat<2> trans3; // , trans, hmat;
                      
                      // quad -> trig transform
                      trans3(0,0) = 1;
                      trans3(0,1) = 0;
                      
                      trans3(1,0) = y*invx;
                      trans3(1,1) = 1*invx;
                      
                      dxdxi_duffy[ii] = trans3 * trans2;
                    }
                }
              break;
            }
          
        }

      default:
        {
          stringstream str;
          str<< "IntegratonRuleTP,2 not available for element type " 
             << ElementTopology::GetElementName(eltype) << endl;
          throw Exception (str.str());
        }
      }
  }













  template <>
  IntegrationRuleTP<3> :: IntegrationRuleTP (ELEMENT_TYPE eltype, FlatArray<int> sort, 
                                             NODE_TYPE nt, int nodenr, int order, LocalHeap & lh)
  {
    int nip;

    static IntegrationRule ir0, ir1;
    {
      lock_guard<mutex> guard(intruletpfacet_mutex);
      if (ir0.GetNIP() == 0)
	{
	  ir0.AddIntegrationPoint (IntegrationPoint (0.0, 0, 0, 1.0));
	  ir1.AddIntegrationPoint (IntegrationPoint (1.0, 0, 0, 1.0));
	}
    }

    switch (eltype)
      {

      case ET_TET:
        {
	  int isort[4];
	  for (int i = 0; i < 4; i++) isort[sort[i]] = i;

	  irx = iry = irz = 0;
	  int powx = 0, powy = 0;


          if (nt == NT_FACE)
	    {
	      const FACE & face = ElementTopology::GetFaces (ET_TET)[nodenr];
	      FACE sface;
	      sface[0] = isort[face[0]];
	      sface[1] = isort[face[1]];
	      sface[2] = isort[face[2]];
	      if (sface[0] > sface[1]) swap (sface[0], sface[1]);
	      if (sface[1] > sface[2]) swap (sface[1], sface[2]);
	      if (sface[0] > sface[1]) swap (sface[0], sface[1]);
	      
	      if (sface[0] == 1 && sface[1] == 2 && sface[2] == 3)
		{
		  powy = 1;
		  irx = &ir0;
		  iry = &SelectIntegrationRule (ET_SEGM, order+1);
		  irz = &SelectIntegrationRule (ET_SEGM, order);
		}
	      if (sface[0] == 0 && sface[1] == 2 && sface[2] == 3)
		{
		  powx = 1;
		  irx = &SelectIntegrationRule (ET_SEGM, order+1);
		  iry = &ir0;
		  irz = &SelectIntegrationRule (ET_SEGM, order);
		}
	      if (sface[0] == 0 && sface[1] == 1 && sface[2] == 3)
		{
		  powx = 1;
		  irx = &SelectIntegrationRule (ET_SEGM, order+1);
		  iry = &SelectIntegrationRule (ET_SEGM, order);
		  irz = &ir0;
		}
	      if (sface[0] == 0 && sface[1] == 1 && sface[2] == 2)
		{
		  powx = 1;
		  irx = &SelectIntegrationRule (ET_SEGM, order+1);
		  iry = &SelectIntegrationRule (ET_SEGM, order);
		  irz = &ir1;
		}
	    }

	  else if (nt == NT_EDGE)
	    {
	      const EDGE & edge = ElementTopology::GetEdges (ET_TET)[nodenr];
	      EDGE sedge;
	      sedge[0] = isort[edge[0]];
	      sedge[1] = isort[edge[1]];
	      if (sedge[0] > sedge[1]) swap (sedge[0], sedge[1]);


	      if (sedge[0] == 0 && sedge[1] == 1)
		{
		  irx = &SelectIntegrationRule (ET_SEGM, order);
		  iry = &ir1;
		  irz = &ir0;
		}
	      if (sedge[0] == 0 && sedge[1] == 2)
		{
		  irx = &SelectIntegrationRule (ET_SEGM, order);
		  iry = &ir0;
		  irz = &ir1;
		}
	      if (sedge[0] == 0 && sedge[1] == 3)
		{
		  irx = &SelectIntegrationRule (ET_SEGM, order);
		  iry = &ir0;
		  irz = &ir0;
		}
	      if (sedge[0] == 1 && sedge[1] == 2)
		{
		  irx = &ir0;
		  iry = &SelectIntegrationRule (ET_SEGM, order);
		  irz = &ir1;
		}
	      if (sedge[0] == 1 && sedge[1] == 3)
		{
		  irx = &ir0;
		  iry = &SelectIntegrationRule (ET_SEGM, order);
		  irz = &ir0;
		}
	      if (sedge[0] == 2 && sedge[1] == 3)
		{
		  irx = &ir0;
		  iry = &ir0;
		  irz = &SelectIntegrationRule (ET_SEGM, order);
		}
	    }




          nip = irx->GetNIP() * iry->GetNIP() * irz->GetNIP();

	  SetSize (nip);

          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            for (int i2 = 0; i2 < iry->GetNIP(); i2++)
              for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                {
                  double 
                    x = (*irx)[i1](0), 
                    y = (*iry)[i2](0), 
                    z = (*irz)[i3](0);

                  double lami[] = { x, 
                                    y * (1-x), 
                                    z * (1-x) * (1-y), 
                                    (1-x)*(1-y)*(1-z) };
                
		  (*this)[ii] = IntegrationPoint (lami[isort[0]],
						  lami[isort[1]],
						  lami[isort[2]],
						  (*irx)[i1].Weight()*(*iry)[i2].Weight()*(*irz)[i3].Weight() * 
						  pow (1-x, powx) * pow(1-y, powy));
                }


          // tet permutation transformation
          double dlamdx[4][3] = 
            { { 1, 0, 0 },
              { 0, 1, 0 },
              { 0, 0, 1 },
              { -1, -1, -1 }};
        
          Mat<3> trans2;
          for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++)
              trans2(l,k) = dlamdx[sort[l]][k];
        
        
          // tet permutation plus hex -> tet mapping
          dxdxi_duffy.SetSize(nip);
          for (int i1 = 0, ii = 0; i1 < irx->GetNIP(); i1++)
            {
              double x = (*irx)[i1](0);
              double invx = 1.0 / (1-x);
            
              for (int i2 = 0; i2 < iry->GetNIP(); i2++)
                {
                  double y = (*iry)[i2](0);
                  double invxy = 1.0 / ( (1-x) * (1-y) );
                
                  for (int i3 = 0; i3 < irz->GetNIP(); i3++, ii++)
                    {
                      double z = (*irz)[i3](0);
                    
                      Mat<3> trans3; // , trans, hmat;
                    
                      // hex -> tet transform
                      trans3(0,0) = 1;
                      trans3(0,1) = 0;
                      trans3(0,2) = 0;
                    
                      trans3(1,0) = y*invx;
                      trans3(1,1) = 1*invx;
                      trans3(1,2) = 0;
                    
                      trans3(2,0) = z*invxy;
                      trans3(2,1) = z*invxy;
                      trans3(2,2) = 1*invxy;
                    
                      dxdxi_duffy[ii] = trans3 * trans2;
                    }
                }
            }
          break;
        }
      default:
        {
          stringstream str;
          str<< "IntegratonRuleTP,2 not available for element type " 
             << ElementTopology::GetElementName(eltype) << endl;
          throw Exception (str.str());
        }
      }
  }









  template class IntegrationRuleTP<1>;
  template class IntegrationRuleTP<2>;
  template class IntegrationRuleTP<3>;






  /** 
      Integration Rules.
      A global class maintaining integration rules. If a rule of specific
      order is requested for the first time, than the rule is generated.
  */
  class IntegrationRules
  {
  public:
    IntegrationRule pointrule;  // 0-dim IR
    Array<IntegrationRule*> segmentrules, segmentrules_inv;
    Array<IntegrationRule*> trigrules;
    Array<IntegrationRule*> quadrules;
    Array<IntegrationRule*> tetrules;
    Array<IntegrationRule*> prismrules;
    Array<IntegrationRule*> pyramidrules;
    Array<IntegrationRule*> hexamidrules;
    Array<IntegrationRule*> hexrules;

    SIMD_IntegrationRule simd_pointrule;
    Array<SIMD_IntegrationRule*> simd_segmentrules, simd_segmentrules_inv;
    Array<SIMD_IntegrationRule*> simd_trigrules;
    Array<SIMD_IntegrationRule*> simd_quadrules;
    Array<SIMD_IntegrationRule*> simd_tetrules;
    Array<SIMD_IntegrationRule*> simd_prismrules;
    Array<SIMD_IntegrationRule*> simd_pyramidrules;
    Array<SIMD_IntegrationRule*> simd_hexamidrules;
    Array<SIMD_IntegrationRule*> simd_hexrules;
    
    Array<IntegrationRule*> jacobirules10;
    Array<IntegrationRule*> jacobirules20;

  public:
    static IntegrationRule intrule0, intrule1;
    static SIMD_IntegrationRule *simd_intrule0, *simd_intrule1;

    ///
    IntegrationRules ();
    ///
    ~IntegrationRules ();
    /// returns the integration rule for the element type and integration order
    const IntegrationRule & SelectIntegrationRule (ELEMENT_TYPE eltyp, int order) const;
    ///
    const IntegrationRule & SelectIntegrationRuleJacobi10 (int order) const;
    ///
    const IntegrationRule & SelectIntegrationRuleJacobi20 (int order) const;
    const SIMD_IntegrationRule & SIMD_SelectIntegrationRule (ELEMENT_TYPE eltyp, int order);
    ///
    const IntegrationRule & GenerateIntegrationRule (ELEMENT_TYPE eltyp, int order);
    const IntegrationRule & GenerateIntegrationRuleJacobi10 (int order);
    const IntegrationRule & GenerateIntegrationRuleJacobi20 (int order);
  };

  IntegrationRule IntegrationRules :: intrule0;
  IntegrationRule IntegrationRules :: intrule1;
  SIMD_IntegrationRule * IntegrationRules :: simd_intrule0;
  SIMD_IntegrationRule * IntegrationRules :: simd_intrule1;




  IntegrationRules :: IntegrationRules ()
  {
    // ************************************
    // ** left/right points
    // ************************************

    intrule0.Append (IntegrationPoint (0,0,0,1));
    intrule1.Append (IntegrationPoint (1,0,0,1));
    simd_intrule0 = new SIMD_IntegrationRule (intrule0);
    simd_intrule1 = new SIMD_IntegrationRule (intrule1);
    
    // ************************************
    // ** Point integration rules
    // ************************************

    pointrule.Append (IntegrationPoint (0, 0, 0, 1));
    // cout << "simd_pointrule: " << simd_pointrule << endl;
    simd_pointrule = SIMD_IntegrationRule (pointrule);
    // cout << "simd_pointrule: " << simd_pointrule << endl;
    
    // ************************************
    // ** Segment integration rules
    // ************************************

    for (int p = 0; p < 20; p++)
      {
        GenerateIntegrationRule (ET_SEGM, p);
        SIMD_SelectIntegrationRule (ET_SEGM, p);
      }

    // ************************************
    // ** Triangle integration rules
    // ************************************


    trigrules.SetSize (7);

    static double qf_trig_order1_points[][3] = 
      {
	{ 1.0/3.0, 1.0/3.0 },
      };
    
    static double qf_trig_order1_weights[] = 
      {  
	0.5
      } ;

    trigrules[0] = new IntegrationRule (1, qf_trig_order1_points, qf_trig_order1_weights);
    trigrules[1] = new IntegrationRule (1, qf_trig_order1_points, qf_trig_order1_weights);


    static double qf_trig_order2_points[][3] = 
      {
	{ 0,   0.5 },
	{ 0.5, 0,  },
	{ 0.5, 0.5 }
      };
   
    static double qf_trig_order2_weights[] = 
      {
	1.0/6.0, 1.0/6.0 , 1.0/6.0
      };

    trigrules[2] = new IntegrationRule (3, qf_trig_order2_points, qf_trig_order2_weights);




    static double qf_trig_order4_points[][3] = 
      {
	{ 0.816847572980459, 0.091576213509771, },
	{ 0.091576213509771, 0.816847572980459, },
	{ 0.091576213509771, 0.091576213509771, },
	{ 0.108103018168070, 0.445948490915965, },
	{ 0.445948490915965, 0.108103018168070, },
	{ 0.445948490915965, 0.445948490915965 }
      };

    static double qf_trig_order4_weights[] = 
      {
	0.054975871827661, 0.054975871827661, 0.054975871827661,
	0.111690794839005, 0.111690794839005, 0.111690794839005
      };

    trigrules[3] = new IntegrationRule (6, qf_trig_order4_points, qf_trig_order4_weights);
    trigrules[4] = new IntegrationRule (6, qf_trig_order4_points, qf_trig_order4_weights);



    static double qf_trig_order6_points[][3] = 
      {
	{ 0.873821971016996, 0.063089014491502, },
	{ 0.063089014491502, 0.873821971016996, },
	{ 0.063089014491502, 0.063089014491502, },
	{ 0.501426509658179, 0.249286745170910, },
	{ 0.249286745170910, 0.501426509658179, },
	{ 0.249286745170910, 0.249286745170910, },

	{ 0.636502499121399, 0.310352451033785, },
	{ 0.310352451033785, 0.053145049844816, },
	{ 0.053145049844816, 0.636502499121399, },
	{ 0.636502499121399, 0.053145049844816, },
	{ 0.310352451033785, 0.636502499121399, },
	{ 0.053145049844816, 0.310352451033785, }
      };

    static double qf_trig_order6_weights[] = 
      {
	0.025422453185103, 0.025422453185103, 0.025422453185103,
	0.058393137863189, 0.058393137863189, 0.058393137863189,

	0.041425537809187, 0.041425537809187, 0.041425537809187,
	0.041425537809187, 0.041425537809187, 0.041425537809187 
      };

    trigrules[5] = new IntegrationRule (12, qf_trig_order6_points, qf_trig_order6_weights);
    trigrules[6] = new IntegrationRule (12, qf_trig_order6_points, qf_trig_order6_weights);


    for (int p = 7; p <= 10; p++)
      GenerateIntegrationRule (ET_TRIG, p);
    for (int p = 0; p <= 10; p++)    
      SIMD_SelectIntegrationRule (ET_TRIG, p);


    // ************************************
    // ** Quadrilateral integration rules
    // ************************************
    
    for (int p = 0; p <= 10; p++)
      {
        GenerateIntegrationRule (ET_QUAD, p);
        SIMD_SelectIntegrationRule (ET_QUAD, p);
      }


    // ************************************
    // ** Tetrahedral integration rules
    // ************************************


    tetrules.SetSize(6);

    static double qf_tetra_order1_points[][3] = 
      { 
	{ 0.25, 0.25, 0.25 },
      };
    
    static double qf_tetra_order1_weights[] = 
      {
	1.0/6.0
      };

    tetrules[0] = new IntegrationRule (1, qf_tetra_order1_points, qf_tetra_order1_weights);
    tetrules[1] = new IntegrationRule (1, qf_tetra_order1_points, qf_tetra_order1_weights);    

    static double qf_tetra_order2_points[][3] = 
      {
	{ 0.585410196624969, 0.138196601125011, 0.138196601125011 },
	{ 0.138196601125011, 0.585410196624969, 0.138196601125011 },
	{ 0.138196601125011, 0.138196601125011, 0.585410196624969 },
	{ 0.138196601125011, 0.138196601125011, 0.138196601125011 }
      };
    
    static double qf_tetra_order2_weights[] = 
      { 1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/24.0 };

    tetrules[2] = new IntegrationRule (4, qf_tetra_order2_points, qf_tetra_order2_weights);    



    static double qf_tetra_order5_points[][3] = 
      {
	{ 0.454496295874350,   0.454496295874350,   0.045503704125650 },
	{ 0.454496295874350,   0.045503704125650,   0.454496295874350 },
	{ 0.454496295874350,   0.045503704125650,   0.045503704125650 },
	{ 0.045503704125650,   0.454496295874350,   0.454496295874350 },
	{ 0.045503704125650,   0.454496295874350,   0.045503704125650 },
	{ 0.045503704125650,   0.045503704125650,   0.454496295874350 },
	{ 0.310885919263301,   0.310885919263301,   0.310885919263301 },   
	{ 0.067342242210098,   0.310885919263301,   0.310885919263301 }, 
	{ 0.310885919263301,   0.067342242210098,   0.310885919263301 },  
	{ 0.310885919263301,   0.310885919263301,   0.067342242210098 },
      
	{ 0.092735250310891,   0.092735250310891,   0.092735250310891 },
	{ 0.721794249067326,   0.092735250310891,   0.092735250310891 },
	{ 0.092735250310891,   0.721794249067326,   0.092735250310891 },
	{ 0.092735250310891,   0.092735250310891,   0.721794249067326 }
      };

    static double qf_tetra_order5_weights[] = 
      {
	0.007091003462847, 0.007091003462847, 0.007091003462847,
	0.007091003462847, 0.007091003462847, 0.007091003462847,
	0.018781320953003, 0.018781320953003, 0.018781320953003, 0.018781320953003,
	0.012248840519394, 0.012248840519394, 0.012248840519394, 0.012248840519394
      };
    
    tetrules[3] = new IntegrationRule (14, qf_tetra_order5_points, qf_tetra_order5_weights);    
    tetrules[4] = new IntegrationRule (14, qf_tetra_order5_points, qf_tetra_order5_weights);    
    tetrules[5] = new IntegrationRule (14, qf_tetra_order5_points, qf_tetra_order5_weights);    


    for (int p = 6; p <= 10; p++)
      GenerateIntegrationRule (ET_TET, p);
    for (int p = 0; p <= 10; p++)
      SIMD_SelectIntegrationRule (ET_TET, p);
    

    // ************************************
    // ** Prismatic integration rules
    // ************************************

    for (int p = 0; p <= 6; p++)
      {
        GenerateIntegrationRule (ET_PRISM, p);
        SIMD_SelectIntegrationRule (ET_PRISM, p);
      }


    // ************************************
    // ** Pyramid integration rules
    // ************************************

    for (int p = 0; p <= 6; p++)
      {
        GenerateIntegrationRule (ET_PYRAMID, p);
        SIMD_SelectIntegrationRule (ET_PYRAMID, p);
      }
    
    // ************************************
    // ** Hexaeder integration rules
    // ************************************

    for (int p = 0; p <= 6; p++)
      {
        GenerateIntegrationRule (ET_HEX, p);
        SIMD_SelectIntegrationRule (ET_HEX, p);

        GenerateIntegrationRule (ET_HEXAMID, p);
        SIMD_SelectIntegrationRule (ET_HEXAMID, p);
      }
    





    // cout << "Check trig intrule:";
    for (int order = 0; order < 6; order++)
      {
	// cout << "order = " << order << endl;

	const IntegrationRule & rule = *trigrules[order];
	const IntegrationRule & rule2 = *trigrules[10];

	for (int ox = 0; ox <= order+2; ox++)
	  for (int oy = 0; ox+oy <= order; oy++)
	    {
	      double sum = 0;
	      for (int j = 0; j < rule.GetNIP(); j++)
		{
		  const IntegrationPoint & ip = rule[j];
		  sum += ip.Weight() * pow (ip(0), ox) * pow (ip(1), oy);
		}
	      double sum2 = 0;
	      for (int j = 0; j < rule2.GetNIP(); j++)
		{
		  const IntegrationPoint & ip = rule2[j];
		  sum2 += ip.Weight() * pow (ip(0), ox) * pow (ip(1), oy);
		}
	      
	      if (fabs (sum - sum2) > 1e-12)
		cout << "ERROR trig rule: \\int x^" << ox << " y^" << oy << " = " << sum << ", diff = " << sum-sum2 << endl;
	    }
      }

    // cout << "Check tet intrule:";
    for (int order = 0; order < 6; order++)
      {
	// cout << "order = " << order << endl;

	const IntegrationRule & rule = *tetrules[order];
	const IntegrationRule & rule2 = *tetrules[10];

	for (int ox = 0; ox <= order+2; ox++)
	  for (int oy = 0; ox+oy <= order+2; oy++)
	    for (int oz = 0; ox+oy+oz <= order; oz++)
	      {
		double sum = 0;
		for (int j = 0; j < rule.GetNIP(); j++)
		  {
		    const IntegrationPoint & ip = rule[j];
		    sum += ip.Weight() * pow (ip(0), ox) * pow (ip(1), oy) * pow(ip(2), oz);
		  }
		double sum2 = 0;
		for (int j = 0; j < rule2.GetNIP(); j++)
		  {
		    const IntegrationPoint & ip = rule2[j];
		    sum2 += ip.Weight() * pow (ip(0), ox) * pow (ip(1), oy) * pow(ip(2), oz);
		  }
		
		if (fabs (sum - sum2) > 1e-12)
		  cout << "ERROR, tet rule: \\int x^" << ox << " y^" << oy << " z^" << oz <<  " = " << sum << ", diff = " << sum-sum2 << endl;
	      }
      }
    



    for (int i = 0; i < 50; i++)
      {
	GenerateIntegrationRuleJacobi10 (i);
	GenerateIntegrationRuleJacobi20 (i);
      }
  }




  IntegrationRules :: ~IntegrationRules ()
  {
    for (int i = 0; i < segmentrules.Size(); i++)
      delete segmentrules[i];

    for (int i = 0; i < trigrules.Size(); i++)
      delete trigrules[i];

    for (int i = 0; i < quadrules.Size(); i++)
      delete quadrules[i];

    for (int i = 0; i < tetrules.Size(); i++)
      delete tetrules[i];

    for (int i = 0; i < prismrules.Size(); i++)
      delete prismrules[i];

    for (int i = 0; i < pyramidrules.Size(); i++)
      delete pyramidrules[i];

    for (int i = 0; i < hexrules.Size(); i++)
      delete hexrules[i];

    for (int i = 0; i < hexamidrules.Size(); i++)
      delete hexamidrules[i];

    for (int i = 0; i < jacobirules10.Size(); i++)
      delete jacobirules10[i];

    for (int i = 0; i < jacobirules20.Size(); i++)
      delete jacobirules20[i];
   
  }



  const IntegrationRule & IntegrationRules :: 
  SelectIntegrationRule (ELEMENT_TYPE eltyp, int order) const
  {
    const Array<IntegrationRule*> * ira;

    switch (eltyp)
      {
      case ET_POINT:
	return pointrule;
      case ET_SEGM:
	ira = &segmentrules; break;
      case ET_TRIG:
	ira = &trigrules; break;
      case ET_QUAD:
	ira = &quadrules; break;
      case ET_TET:
	ira = &tetrules; break;
      case ET_PYRAMID:
	ira = &pyramidrules; break;
      case ET_PRISM:
	ira = &prismrules; break;
      case ET_HEXAMID:
	ira = &hexamidrules; break;
      case ET_HEX:
	ira = &hexrules; break;
      default:
	{
	  stringstream str;
	  str << "no integration rules for element " << int(eltyp) << endl;
	  throw Exception (str.str());
	}
      }

    if (order < 0) 
      { order = 0; }

    if (order >= ira->Size() || (*ira)[order] == 0)
      {
        return const_cast<IntegrationRules&> (*this).
	  GenerateIntegrationRule (eltyp, order);
      }

    return *((*ira)[order]);
  }
 

  const IntegrationRule & IntegrationRules :: SelectIntegrationRuleJacobi10 (int order) const
  {
    const Array<IntegrationRule*> * ira;
  
    ira = &jacobirules10; 

    if (order < 0) { order = 0; }

    if (order >= ira->Size() || (*ira)[order] == 0)
      {
	return const_cast<IntegrationRules&> (*this).
	  GenerateIntegrationRuleJacobi10 (order);
      }

    return *((*ira)[order]);
  }
 

  const IntegrationRule & IntegrationRules :: SelectIntegrationRuleJacobi20 (int order) const
  {
    const Array<IntegrationRule*> * ira;
  
    ira = &jacobirules20; 

    if (order < 0) { order = 0; }

    if (order >= ira->Size() || (*ira)[order] == 0)
      {
	return const_cast<IntegrationRules&> (*this).
	  GenerateIntegrationRuleJacobi20 (order);
      }

    return *((*ira)[order]);
  }
 






  const IntegrationRule & IntegrationRules :: 
  GenerateIntegrationRule (ELEMENT_TYPE eltyp, int order)
  {
    Array<IntegrationRule*> * ira;

    if (eltyp == ET_QUAD || eltyp == ET_TRIG)
      {
        GenerateIntegrationRule (ET_SEGM, order);
      }

    if (eltyp == ET_TET || eltyp == ET_PRISM || eltyp == ET_HEX || eltyp == ET_PYRAMID) 
      {        
        GenerateIntegrationRule (ET_SEGM, order);
        GenerateIntegrationRule (ET_SEGM, order+2);
        GenerateIntegrationRule (ET_TRIG, order);
        GenerateIntegrationRule (ET_QUAD, order);
      }

    {
      lock_guard<mutex> guard(genintrule_mutex);
  
      switch (eltyp)
	{
	case ET_SEGM:
	  ira = &segmentrules; break;
	case ET_TRIG:
	  ira = &trigrules; break;
	case ET_QUAD:
	  ira = &quadrules; break;
	case ET_TET:
	  ira = &tetrules; break;
	case ET_PYRAMID:
	  ira = &pyramidrules; break;
	case ET_PRISM:
	  ira = &prismrules; break;
	case ET_HEXAMID:
	  ira = &hexamidrules; break;
	case ET_HEX:
	  ira = &hexrules; break;
	default:
	    /*
	      {
	      stringstream str;
	      str << "no integration rules for element " << int(eltyp) << endl;
	      throw Exception (str.str());
	      }
	    */
	  throw Exception ("no integration rules for element " + 
			   ToString(int(eltyp)) + "\n"); 
	}

      if (ira -> Size() < order+1)
	{
	  int oldsize = ira -> Size();
	  ira -> SetAllocSize (order+1);
	  for (int i = oldsize; i < order+1; i++)
	    ira -> Append(nullptr);
	}

      if ( (*ira)[order] == 0)
	{
	  switch (eltyp)
	    {
	    case ET_POINT: break;
	    case ET_SEGM:
	      {
		Array<double> xi, wi;
		ComputeGaussRule (order/2+1, xi, wi);
		IntegrationRule * rule = new IntegrationRule; //  (xi.Size());
		double xip[3] = { 0, 0, 0 };
		for (int j = 0; j < xi.Size(); j++)
		  {
		    xip[0] = xi[j];
		    IntegrationPoint ip = IntegrationPoint (xip, wi[j]);
		    ip.SetNr (j);
		    // if (order < 20)
		      // ip.SetGlobNr (segmentpoints.Append (ip)-1);
		    rule->AddIntegrationPoint (ip);
		  }
		segmentrules[order] = rule;
                
                IntegrationRule * rule_inv = new IntegrationRule;
                for (int j = rule->Size()-1; j >= 0; j--)
                  rule_inv->AddIntegrationPoint ((*rule)[j]);
                if (segmentrules_inv.Size() < order+1)
                  segmentrules_inv.SetSize(order+1);
                segmentrules_inv[order] = rule_inv;
		break;
	      }

	    case ET_TRIG:
	      {
		Array<double> xx, wx, xy, wy;
		ComputeGaussJacobiRule (order/2+1, xx, wx, 1, 0);
		// ComputeGaussRule (order/2+2, xx, wx);
		ComputeGaussRule (order/2+1, xy, wy);

		IntegrationRule * trigrule = new IntegrationRule; // (xx.Size()*xy.Size());
	      
		int ii = 0;
                for (int j = 0; j < xy.Size(); j++)
                  for (int i = 0; i < xx.Size(); i++)
		    {
		      IntegrationPoint ip = 
			IntegrationPoint (xx[i], xy[j]*(1-xx[i]), 0,
                                          wx[i]*wy[j]*(1-xx[i]));
		      ip.SetNr (ii); ii++;
		      // if (order <= 10)
			// ip.SetGlobNr (trigpoints.Append (ip)-1);
		      trigrule->AddIntegrationPoint (ip);
		    }
		trigrules[order] = trigrule;
		break;
	      }


	    case ET_QUAD:
	      {
		const IntegrationRule & segmrule = SelectIntegrationRule (ET_SEGM, order);
		IntegrationRule * quadrule = new IntegrationRule; // (segmrule.GetNIP()*segmrule.GetNIP());
	      
		double point[3], weight;

		int ii = 0;
		for (int i = 0; i < segmrule.GetNIP(); i++)
		  for (int j = 0; j < segmrule.GetNIP(); j++)
		    {
		      const IntegrationPoint & ipsegm1 = segmrule[i];
		      const IntegrationPoint & ipsegm2 = segmrule[j];
		      point[0] = ipsegm1.Point()[0];
		      point[1] = ipsegm2.Point()[0];
		      point[2] = 0;
		      weight = ipsegm1.Weight() * ipsegm2.Weight();
		    
		      IntegrationPoint ip = IntegrationPoint (point, weight);
		      ip.SetNr (ii); ii++;
		      // if (order <= 10)
			// ip.SetGlobNr (quadpoints.Append (ip)-1);
		      quadrule->AddIntegrationPoint (ip);
		    }
		quadrules[order] = quadrule;
		break;
	      }
  

	    case ET_TET:
	      {
		// tet-rules by degenerated tensor product rules:

		Array<double> xx, wx, xy, wy, xz, wz;
		ComputeGaussRule (order/2+1, xz, wz);
		ComputeGaussJacobiRule (order/2+1, xy, wy, 1, 0);
		ComputeGaussJacobiRule (order/2+1, xx, wx, 2, 0);

		IntegrationRule * tetrule = new IntegrationRule; // (xx.Size()*xy.Size()*xz.Size());
	      
		int ii = 0;
                /*
		for (int i = 0; i < xx.Size(); i++)
		  for (int j = 0; j < xy.Size(); j++)
		    for (int k = 0; k < xz.Size(); k++)
                */
                // change loop order ... x-pnts are innermost
                for (int k = 0; k < xz.Size(); k++)
		  for (int j = 0; j < xy.Size(); j++)
                    for (int i = 0; i < xx.Size(); i++)
		      {
			IntegrationPoint ip = 
			  IntegrationPoint (xx[i], 
                                            xy[j]*(1-xx[i]),
                                            xz[k]*(1-xx[i])*(1-xy[j]),
                                            wx[i]*wy[j]*wz[k]*sqr(1-xx[i])*(1-xy[j]));
			ip.SetNr (ii); ii++;
			// if (order <= 6)
			  // ip.SetGlobNr (tetpoints.Append (ip)-1);
			tetrule->AddIntegrationPoint (ip);
		      }
		tetrules[order] = tetrule;
		break;
	      }

	    case ET_HEX:
	      {
		const IntegrationRule & segmrule = SelectIntegrationRule (ET_SEGM, order);
                
		IntegrationRule * hexrule = 
		  new IntegrationRule; //(segmrule.GetNIP()*segmrule.GetNIP()*segmrule.GetNIP());
	
		double point[3], weight;
		int ii = 0;
		for (int i = 0; i < segmrule.GetNIP(); i++)
		  for (int j = 0; j < segmrule.GetNIP(); j++)
		    for (int l = 0; l < segmrule.GetNIP(); l++)
		      {
			const IntegrationPoint & ipsegm1 = segmrule[i];
			const IntegrationPoint & ipsegm2 = segmrule[j];
			const IntegrationPoint & ipsegm3 = segmrule[l];
		      
			point[0] = ipsegm1.Point()[0];
			point[1] = ipsegm2.Point()[0];
			point[2] = ipsegm3.Point()[0];
			weight = ipsegm1.Weight() * ipsegm2.Weight() * ipsegm3.Weight();
		      
			IntegrationPoint ip = 
			  IntegrationPoint (point, weight);
		      
			ip.SetNr (ii); ii++;
			hexrule->AddIntegrationPoint (ip);
		      }
		hexrules[order] = hexrule;
		break;
	      }


	    case ET_HEXAMID:
	      {
		const IntegrationRule & segmrule = SelectIntegrationRule (ET_SEGM, order);
                
		IntegrationRule * hexamidrule = 
		  new IntegrationRule; //(segmrule.GetNIP()*segmrule.GetNIP()*segmrule.GetNIP());
	
		double point[3], weight;
		int ii = 0;
		for (int i = 0; i < segmrule.GetNIP(); i++)
		  for (int j = 0; j < segmrule.GetNIP(); j++)
		    for (int l = 0; l < segmrule.GetNIP(); l++)
		      {
			const IntegrationPoint & ipsegm1 = segmrule[i];
			const IntegrationPoint & ipsegm2 = segmrule[j];
			const IntegrationPoint & ipsegm3 = segmrule[l];
                        double x = ipsegm1.Point()[0];
                        double y = ipsegm2.Point()[0];
                        double z = ipsegm3.Point()[0];

                        /*
			point[0] = x*(1-y)*(1-z);
			point[1] = y;
			point[2] = z;
			weight = (1-y)*(1-z)*ipsegm1.Weight() * ipsegm2.Weight() * ipsegm3.Weight();
                        */
			point[0] = x*(1-y*z);
			point[1] = y;
			point[2] = z;
			weight = (1-y*z)*ipsegm1.Weight() * ipsegm2.Weight() * ipsegm3.Weight();
                        
			IntegrationPoint ip = 
			  IntegrationPoint (point, weight);
		      
			ip.SetNr (ii); ii++;
			hexamidrule->AddIntegrationPoint (ip);
		      }
		hexamidrules[order] = hexamidrule;
		break;
	      }


              

	    case ET_PRISM:
	      {
		const IntegrationRule & segmrule = SelectIntegrationRule (ET_SEGM, order);
		const IntegrationRule & trigrule = SelectIntegrationRule (ET_TRIG, order);

		IntegrationRule * prismrule = new IntegrationRule; 
      
		double point[3], weight;
	      
		int ii = 0;
		for (int i = 0; i < segmrule.GetNIP(); i++)
		  for (int j = 0; j < trigrule.GetNIP(); j++)
		    {
		      const IntegrationPoint & ipsegm = segmrule[i]; // .GetIP(i);
		      const IntegrationPoint & iptrig = trigrule[j]; // .GetIP(j);
		      point[0] = iptrig.Point()[0];
		      point[1] = iptrig.Point()[1];
		      point[2] = ipsegm.Point()[0];
		      weight = iptrig.Weight() * ipsegm.Weight();
	      
		      IntegrationPoint ip = 
			IntegrationPoint (point, weight);
		      ip.SetNr (ii); ii++;
		      prismrule->AddIntegrationPoint (ip);
		    }
		prismrules[order] = prismrule;
		break;
	      }


	    case ET_PYRAMID:
	      {
		const IntegrationRule & quadrule = SelectIntegrationRule (ET_QUAD, order);
		const IntegrationRule & segrule = SelectIntegrationRule (ET_SEGM, order+2);

		IntegrationRule * pyramidrule = new IntegrationRule;
	
		double point[3], weight;
	      
		int ii = 0;
		for (int i = 0; i < quadrule.GetNIP(); i++)
		  for (int j = 0; j < segrule.GetNIP(); j++)
		    {
		      const IntegrationPoint & ipquad = quadrule[i]; // .GetIP(i);
		      const IntegrationPoint & ipseg = segrule[j]; // .GetIP(j);
		      point[0] = (1-ipseg(0)) * ipquad(0);
		      point[1] = (1-ipseg(0)) * ipquad(1);
		      point[2] = ipseg(0);
		      weight = ipseg.Weight() * sqr (1-ipseg(0)) * ipquad.Weight();
		    
		      IntegrationPoint ip = IntegrationPoint (point, weight);
		      ip.SetNr (ii); ii++;
		      pyramidrule->AddIntegrationPoint (ip);
		    }
		pyramidrules[order] = pyramidrule;
		break;
	      }
	    }
	}

      if ( (*ira)[order] == 0)
	{
	  stringstream str;
	  str << "could not generate Integration rule of order " << order 
	      << " for element type " 
	      << ElementTopology::GetElementName(eltyp) << endl;
	  throw Exception (str.str());
	}
    }

    return *(*ira)[order];
  }






  const IntegrationRule & IntegrationRules :: GenerateIntegrationRuleJacobi10 (int order)
  {
    Array<IntegrationRule*> * ira;
    ira = &jacobirules10; 

    {
      lock_guard<mutex> guard(genintrule_mutex);
      if (ira -> Size() < order+1)
	{
	  int oldsize = ira -> Size();
	  ira -> SetSize (order+1);
	  for (int i = oldsize; i < order+1; i++)
	    (*ira)[i] = 0;
	}

      if ( (*ira)[order] == 0)
	{
	  Array<double> xi, wi;
	  // ComputeGaussRule (order/2+1, xi, wi);
	  ComputeGaussJacobiRule (order/2+1, xi, wi, 1, 0);
	  IntegrationRule * rule = new IntegrationRule; // (xi.Size());
	  double xip[3] = { 0, 0, 0 };
	  for (int j = 0; j < xi.Size(); j++)
	    {
	      xip[0] = xi[j];
	      IntegrationPoint ip = IntegrationPoint (xip, wi[j]);
	      ip.SetNr (j);
              // 	      if (order < 20)
		// ip.SetGlobNr (segmentpoints.Append (ip)-1);
	      rule->AddIntegrationPoint (ip);
	    }
	  jacobirules10[order] = rule;	      
	}

      if ( (*ira)[order] == 0)
	{
	  stringstream str;
	  str << "could not generate Jacobi-10 integration rule of order " << order 
	      << " for element type " << endl;
	  throw Exception (str.str());
	}
    }
    return *(*ira)[order];
  }






  const IntegrationRule & IntegrationRules :: GenerateIntegrationRuleJacobi20 (int order)
  {
    Array<IntegrationRule*> * ira;
    ira = &jacobirules20; 

    {
      lock_guard<mutex> guard(genintrule_mutex);
      if (ira -> Size() < order+1)
	{
	  int oldsize = ira -> Size();
	  ira -> SetSize (order+1);
	  for (int i = oldsize; i < order+1; i++)
	    (*ira)[i] = 0;
	}

      if ( (*ira)[order] == 0)
	{
	  Array<double> xi, wi;
	  ComputeGaussJacobiRule (order/2+1, xi, wi, 2, 0);
	  IntegrationRule * rule = new IntegrationRule; 
	  double xip[3] = { 0, 0, 0 };
	  for (int j = 0; j < xi.Size(); j++)
	    {
	      xip[0] = xi[j];
	      IntegrationPoint ip = IntegrationPoint (xip, wi[j]);
	      ip.SetNr (j);
	      rule->AddIntegrationPoint (ip);
	    }
	  jacobirules20[order] = rule;	      
	}

      if ( (*ira)[order] == 0)
	{
	  stringstream str;
	  str << "could not generate Jacobi-20 integration rule of order " << order 
	      << " for element type " << endl;
	  throw Exception (str.str());
	}
    }
    return *(*ira)[order];
  }


  const SIMD_IntegrationRule & IntegrationRules :: SIMD_SelectIntegrationRule (ELEMENT_TYPE eltype, int order)
  {
    Array<SIMD_IntegrationRule*> * ira;

    switch (eltype)
      {
      case ET_POINT:
	return simd_pointrule;
      case ET_SEGM:
	ira = &simd_segmentrules; break;
      case ET_TRIG:
	ira = &simd_trigrules; break;
      case ET_QUAD:
	ira = &simd_quadrules; break;
      case ET_TET:
	ira = &simd_tetrules; break;
      case ET_PYRAMID:
	ira = &simd_pyramidrules; break;
      case ET_PRISM:
	ira = &simd_prismrules; break;
      case ET_HEXAMID:
	ira = &simd_hexamidrules; break;
      case ET_HEX:
	ira = &simd_hexrules; break;
      default:
	{
	  stringstream str;
	  str << "no simd-integration rules for element " << int(eltype) << endl;
	  throw Exception (str.str());
	}
      }

    if (order < 0) 
      { order = 0; }

    if (order >= ira->Size() || (*ira)[order] == 0)
      {
        lock_guard<mutex> guard(simd_genintrule_mutex[eltype]);

        if (order >= ira->Size())
          {
            int oldsize = ira->Size();
            ira -> SetAllocSize(order+2);
            for (int i = oldsize; i < order+2; i++)
                ira -> Append(nullptr);
          }

        if ( (*ira)[order] == nullptr)
          {
            IntegrationRule ir(eltype, order);
            auto tmp = new SIMD_IntegrationRule(ir);
            switch (eltype)
              {
              case ET_SEGM:
                {
                  if (simd_segmentrules_inv.Size() < order+2)
                    simd_segmentrules_inv.SetSize(order+2);
                  simd_segmentrules_inv[order] = new SIMD_IntegrationRule(*segmentrules_inv[order]);
                  if (order % 2 == 0)
                    simd_segmentrules_inv[order+1] = simd_segmentrules_inv[order];
                  else
                    simd_segmentrules_inv[order-1] = simd_segmentrules_inv[order];
                  break;
                }
              case ET_QUAD:
                {
                  tmp->SetIRX (&SIMD_SelectIntegrationRule (ET_SEGM, order));
                  tmp->SetIRY (&SIMD_SelectIntegrationRule (ET_SEGM, order));
                  break;
                }
              case ET_TRIG:
                {
                  const IntegrationRule & irx = SelectIntegrationRuleJacobi10(order);
                  tmp->SetIRX (new SIMD_IntegrationRule(irx));
                  tmp->SetIRY (&SIMD_SelectIntegrationRule (ET_SEGM, order));
                  if (tmp->GetNIP() != tmp->GetIRX().GetNIP()*tmp->GetIRY().GetNIP())
                    {
                      tmp->SetIRX(nullptr);
                      tmp->SetIRY(nullptr);
                    }
                  break;
                }
              case ET_TET:
                {
                  const IntegrationRule & irx = SelectIntegrationRuleJacobi20(order);
                  const IntegrationRule & iry = SelectIntegrationRuleJacobi10(order);
                  tmp->SetIRX (new SIMD_IntegrationRule(irx));
                  tmp->SetIRY (new SIMD_IntegrationRule(iry));                  
                  tmp->SetIRZ (&SIMD_SelectIntegrationRule (ET_SEGM, order));
                  if (tmp->GetNIP() != tmp->GetIRX().GetNIP()*tmp->GetIRY().GetNIP()*tmp->GetIRZ().GetNIP())
                    {
                      tmp->SetIRX(nullptr);
                      tmp->SetIRY(nullptr);
                      tmp->SetIRZ(nullptr);
                    }
                  break;
                }
              case ET_HEX:
                {
                  tmp->SetIRX (&SIMD_SelectIntegrationRule (ET_SEGM, order));
                  tmp->SetIRY (&SIMD_SelectIntegrationRule (ET_SEGM, order));
                  tmp->SetIRZ (&SIMD_SelectIntegrationRule (ET_SEGM, order));
                  break;
                }
              default:
                ;
              }
            (*ira)[order] = tmp;
          }
      }

    return *((*ira)[order]);
  }

  SIMD_IntegrationRule::SIMD_IntegrationRule (const IntegrationRule & ir)
    : Array<SIMD<IntegrationPoint>> (0, nullptr)
  {
    nip = ir.GetNIP();
    this->size = (ir.Size()+SIMD<IntegrationPoint>::Size()-1) / SIMD<IntegrationPoint>::Size();
    this -> mem_to_delete = new SIMD<IntegrationPoint>[this->size];
    this->data = this->mem_to_delete;
    
    dimension = ir.Dim();
    for (int i = 0; i < Size(); i++)
      (*this)[i] = [&] (int j) { int nr = i*SIMD<IntegrationPoint>::Size()+j;
                                 bool regularip = nr < ir.Size();
                                 IntegrationPoint ip = ir[regularip ? nr : ir.Size()-1];
                                 if (!regularip) ip.SetWeight(0);
                                 return ip; };
    // return (nr < ir.Size()) ? ir[nr] : IntegrationPoint(0,0,0,0); };
  }

  SIMD_IntegrationRule::SIMD_IntegrationRule (const IntegrationRule & ir, LocalHeap & lh)
    : Array<SIMD<IntegrationPoint>> ( (ir.Size()+SIMD<IntegrationPoint>::Size()-1) / SIMD<IntegrationPoint>::Size(), lh)
  {
    nip = ir.GetNIP();
    dimension = ir.Dim();
    for (int i = 0; i < Size(); i++)
      (*this)[i] = [&] (int j) { int nr = i*SIMD<IntegrationPoint>::Size()+j;
                                 return (nr < ir.Size()) ? ir[nr] : IntegrationPoint(0,0,0,0); };
  }

  SIMD_IntegrationRule::SIMD_IntegrationRule (int _nip, LocalHeap & lh)
    : Array<SIMD<IntegrationPoint>> ( (_nip+SIMD<IntegrationPoint>::Size()-1) / SIMD<IntegrationPoint>::Size(), lh)
  {
    nip = _nip;
  }

  /*
  SIMD_IntegrationRule::~SIMD_IntegrationRule()
  {
    if (mem_to_delete)
      _mm_free(mem_to_delete);
    mem_to_delete = nullptr;
  }
  */
  
  template <int DIM_ELEMENT, int DIM_SPACE>
  SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> :: 
  SIMD_MappedIntegrationRule (const SIMD_IntegrationRule & ir, 
                              const ElementTransformation & aeltrans, 
                              Allocator & lh)
    : SIMD_BaseMappedIntegrationRule (ir, aeltrans), mips(ir.Size(), lh)
  {
    dim_element = DIM_ELEMENT;
    dim_space = DIM_SPACE;
    baseip = (char*)(void*)(SIMD<BaseMappedIntegrationPoint>*)(&mips[0]);
    incr = sizeof (SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>);
    
    FlatArray<SIMD<IntegrationPoint>> hir = ir;
    FlatArray<SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>> hmips = mips;
    for (size_t i = 0; i < hir.Size(); i++)
      new (&hmips[i]) SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>> (hir[i], eltrans, -1);

    /*
    new (&points) BareSliceMatrix<SIMD<double>> (sizeof(SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)/sizeof(SIMD<double>),
                                                 &mips[0].Point()(0),
                                                 DummySize(mips.Size(), DIM_SPACE));
    
    new (&normals) BareSliceMatrix<SIMD<double>> (sizeof(SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)/sizeof(SIMD<double>),
                                                  &mips[0].NV()(0),
                                                  DummySize(mips.Size(), DIM_SPACE));
    */
    new (&points) BareSliceMatrix<SIMD<double>> (mips.Size(), DIM_SPACE,
                                                 sizeof(SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)/sizeof(SIMD<double>),
                                                 &mips[0].Point()(0));
    
    new (&normals) BareSliceMatrix<SIMD<double>> (mips.Size(), DIM_SPACE,
                                                  sizeof(SIMD<MappedIntegrationPoint<DIM_ELEMENT, DIM_SPACE>>)/sizeof(SIMD<double>),
                                                  &mips[0].NV()(0));
    
    eltrans.CalcMultiPointJacobian (ir, *this);

    if (ir.Size())
      if (ir[0].VB() != VOL)
        ComputeNormalsAndMeasure (eltrans.GetElementType(), ir[0].FacetNr());
  }

  template <int DIM_ELEMENT, int DIM_SPACE>
  void SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> :: 
  ComputeNormalsAndMeasure (ELEMENT_TYPE et, int facetnr)
  {   // 
    auto hmips = mips;
    if (hmips.Size() == 0) return;

    // codimension == 0:
    if ((Dim(et) >= 2) && (Dim(et)-int(hmips[0].IP().VB()) == 0))
      {
        for (size_t i = 0; i < hmips.Size(); i++)
          hmips[i].SetMeasure(1);
        return;
      }
    
    if constexpr(DIM_ELEMENT == 3)
      if (hmips[0].IP().VB() == BBND && Dim(et) == 3)
      {
	FlatVector<Vec<3>> points(99,(Vec<3>*)ElementTopology::GetVertices (et));
        auto edges = ElementTopology::GetEdges (et);
        Vec<3> tau_ref = points(edges[facetnr][1]) - points(edges[facetnr][0]);
        for (int i = 0; i < hmips.Size(); i++)
          {
            auto & mip = hmips[i];        
            Vec<3,SIMD<double>> tau = mip.GetJacobian() * tau_ref;
            mip.SetMeasure(L2Norm(tau));
          }
        return;
      }
    
    Vec<DIM_ELEMENT> normal_ref = ElementTopology::GetNormals<DIM_ELEMENT>(et)[facetnr];
    for (int i = 0; i < hmips.Size(); i++)
      {
        auto & mip = hmips[i];
        /*
        Mat<DIM_ELEMENT,DIM_SPACE,SIMD<double>> inv_jac = mip.GetJacobianInverse();
        SIMD<double> det = mip.GetMeasure();
        Vec<DIM_SPACE,SIMD<double>> normal = det * Trans (inv_jac) * normal_ref;
        SIMD<double> len = L2Norm (normal);       // that's the surface measure
        normal *= SIMD<double>(1.0)/len;                           // normal vector on physical element
        */

        Vec<DIM_SPACE,SIMD<double>> normal = mip.GetJacobianCofactor() * normal_ref;
        SIMD<double> len = L2Norm (normal);
        SIMD<double> sign = IfPos(mip.GetJacobiDet(), SIMD<double>(1), SIMD<double>(-1));
        normal *= sign/len;    
        
        // SIMD<double> len = 1.0;
        // mip.SetNV(normal);
        mip.SetMeasure (len);
        
        if constexpr (DIM_ELEMENT == DIM_SPACE)
          {
            mip.SetNV(normal);
            if (DIM_SPACE == 2)
              mip.SetTV(Vec<2,SIMD<double>>(-normal(1),normal(0)));
          }
        else
          {
            if constexpr (DIM_SPACE == 3 && DIM_ELEMENT == 1)
              {
                mip.SetTV(Vec<3,SIMD<double>>(normal));
              }
            else if constexpr ( DIM_SPACE == 3 && DIM_ELEMENT == 2)
              {
                Vec<3,SIMD<double>> tang = Cross(Vec<3,SIMD<double>> (normal), Vec<3,SIMD<double>> (mip.GetNV()));
                mip.SetTV(tang);
              }
            else if constexpr ( DIM_SPACE == 2 && DIM_ELEMENT == 1)
              {
                mip.SetTV(Vec<2,SIMD<double>>(normal));
              }
            else if constexpr (DIM_ELEMENT == 0)
                                ; 
            else
              {
                throw Exception("In ComputeNormalsAndMeasure SIMD: not implemented yet!");
              }
          }
        
      }
  }


  template <int DIM_ELEMENT, int DIM_SPACE>
  void SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> :: Print (ostream & ost) const
  {
    ost << "simd-mir, size = " << mips.Size() << endl;
    for (size_t i = 0; i < mips.Size(); i++)
      mips[i].Print(ost);
  }


  template <int DIM_ELEMENT, int DIM_SPACE>
  void SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> ::
  TransformGradient (BareSliceMatrix<SIMD<double>> grad) const 
  {
    if constexpr(DIM_ELEMENT != DIM_SPACE)
      throw Exception("transformgrad only available for volume mapping");
    else
      for (size_t i = 0; i < mips.Size(); i++)
      {
        Vec<DIM_ELEMENT,SIMD<double>> vref = grad.Col(i);
        // Vec<DIM_SPACE,SIMD<double>> vphys =
        grad.Col(i).Range(0,DIM_SPACE) = Trans(mips[i].GetJacobianInverse()) * vref;
      }
  }
  
  template <int DIM_ELEMENT, int DIM_SPACE>
  void SIMD_MappedIntegrationRule<DIM_ELEMENT,DIM_SPACE> ::
  TransformGradientTrans (BareSliceMatrix<SIMD<double>> grad) const 
  {
    if constexpr (DIM_ELEMENT != DIM_SPACE)
      throw Exception("transformgrad only available for volume mapping");
    else
      for (size_t i = 0; i < mips.Size(); i++)
      {
        Vec<DIM_ELEMENT,SIMD<double>> vref = grad.Col(i);
        // Vec<DIM_SPACE,SIMD<double>> vphys =
        grad.Col(i).Range(DIM_SPACE) = mips[i].GetJacobianInverse() * vref;
      }
  }

  
  
  template class SIMD_MappedIntegrationRule<0,0>;
  template class SIMD_MappedIntegrationRule<0,1>;
  template class SIMD_MappedIntegrationRule<1,1>;
  template class SIMD_MappedIntegrationRule<2,2>;
  template class SIMD_MappedIntegrationRule<3,3>;
  template class SIMD_MappedIntegrationRule<1,2>;
  template class SIMD_MappedIntegrationRule<1,3>;
  template class SIMD_MappedIntegrationRule<2,3>;
  template class SIMD_MappedIntegrationRule<0,2>;
  template class SIMD_MappedIntegrationRule<0,3>;






  
  int Integrator :: common_integration_order = -1;

  static IntegrationRules intrules;
  const IntegrationRules & GetIntegrationRules ()
  {
    return intrules;
  }


  const IntegrationRule & SelectIntegrationRule (ELEMENT_TYPE eltype, int order)
  {
    return GetIntegrationRules ().SelectIntegrationRule (eltype, order);
  }

  const IntegrationRule & SelectIntegrationRuleJacobi20 (int order)
  {
    return GetIntegrationRules ().SelectIntegrationRuleJacobi20 (order);
  }

  const IntegrationRule & SelectIntegrationRuleJacobi10 (int order)
  {
    return GetIntegrationRules ().SelectIntegrationRuleJacobi10 (order);
  }

  const SIMD_IntegrationRule & SIMD_SelectIntegrationRule (ELEMENT_TYPE eltype, int order)
  {
    return const_cast<IntegrationRules&>(GetIntegrationRules()).SIMD_SelectIntegrationRule (eltype, order);
  }








  const SIMD_IntegrationRule & Facet2ElementTrafo :: operator() (int fnr, const SIMD_IntegrationRule & irfacet, LocalHeap & lh)
  {
    if (vb == VOL) return irfacet;
    SIMD_IntegrationRule & irvol = *new (lh) SIMD_IntegrationRule (irfacet.GetNIP(), lh);

    FlatArray<SIMD<IntegrationPoint>> hirfacet = irfacet;
    FlatArray<SIMD<IntegrationPoint>> hirvol = irvol;
    
    // switch (ElementTopology::GetFacetType(eltype, fnr))
    switch (FacetType(fnr))
      {
      case ET_POINT:
        {
          irvol[0](0) = points (fnr)[0];
          irvol[0](1) = points (fnr)[1];
          irvol[0](2) = points (fnr)[2];
          break;
        }
      case ET_SEGM:
        {
          Vec<3> p1 = points (edges[fnr][0]);
          Vec<3> p2 = points (edges[fnr][1]);
          Vec<3> delta = p1-p2;
          for (int i = 0; i < hirfacet.Size(); i++)
            {
              auto ip = hirfacet[i];
              auto & ipvol = hirvol[i];              
              for (int k = 0; k < 3; k++)
                ipvol(k) = p2(k) + delta(k) * ip(0);
            }
          break;
        }

      case ET_TRIG:
        {
          Vec<3> p0 = points(faces[fnr][0]);
          Vec<3> p1 = points(faces[fnr][1]);
          Vec<3> p2 = points(faces[fnr][2]);
          Vec<3> delta1 = p0-p2;
          Vec<3> delta2 = p1-p2;

          /*
          for (int i = 0; i < irfacet.GetNIP(); i++)
            irvol[i] = Vec<3> (p2 + irfacet[i](0) * (p0-p2) + irfacet[i](1)*(p1-p2));
          */
          for (int i = 0; i < hirfacet.Size(); i++)
            {
              auto ip = hirfacet[i];
              auto & ipvol = hirvol[i];              
              for (int k = 0; k < 3; k++)
                ipvol(k) = p2(k) + delta1(k) * ip(0) + delta2(k) * ip(1);
            }
          break;
        }

	case ET_QUAD:
	  {
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][3]);
            Vec<3> delta1 = p1-p0;
            Vec<3> delta2 = p2-p0;

	    // for (int i = 0; i < irfacet.GetNIP(); i++)
            // irvol[i] = Vec<3> (p0 + irfacet[i](0) * (p1-p0) + irfacet[i](1)*(p2-p0));
            for (int i = 0; i < hirfacet.Size(); i++)
              {
                auto ip = hirfacet[i];
                auto & ipvol = hirvol[i];              
                for (int k = 0; k < 3; k++)
                  ipvol(k) = p0(k) + delta1(k) * ip(0) + delta2(k) * ip(1);
              }
	    break;
	  }

      default:
        throw Exception ("undefined facet type in SIMD_Facet2ElementTrafo()\n");
      } 

    for (int i = 0; i < hirfacet.Size(); i++)
      {
        hirvol[i].SetFacetNr(fnr, vb);
        hirvol[i].Weight() = hirfacet[i].Weight();
      }


    // tensor-product
    switch (eltype)
      {
      case ET_TET:
        {
          if (irfacet.IsTP())
          switch (fnr)
            {
            case 0:
              {
                irvol.SetIRX (intrules.simd_intrule0);
                irvol.SetIRY (&irfacet.GetIRX());
                if (!swapped)
                  irvol.SetIRZ (&irfacet.GetIRY());
                else
                  {
                    int order = 2*irfacet.GetIRY().GetNIP()-1;
                    irvol.SetIRZ (intrules.simd_segmentrules_inv[order]);
                  }
                break;
              }
            case 1:
              {
                irvol.SetIRX (&irfacet.GetIRX());
                irvol.SetIRY (intrules.simd_intrule0);
                if (!swapped)
                  irvol.SetIRZ (&irfacet.GetIRY());
                else
                  {
                    int order = 2*irfacet.GetIRY().GetNIP()-1;
                    irvol.SetIRZ (intrules.simd_segmentrules_inv[order]);
                  }
                break;
              }
            case 2:
              {
                irvol.SetIRX (&irfacet.GetIRX());
                irvol.SetIRY (&irfacet.GetIRY());
                irvol.SetIRZ (intrules.simd_intrule0);
                break;
              }
            case 3:
              {
                irvol.SetIRX (&irfacet.GetIRX());
                irvol.SetIRY (&irfacet.GetIRY());
                irvol.SetIRZ (intrules.simd_intrule1);
                break;
              }
            default:
              ;
            }
          break;
        }


      case ET_QUAD:
        {
          FlatVec<3> p0 = points(edges[fnr][0]);
          FlatVec<3> p1 = points(edges[fnr][1]);
          Vec<3> delta1 = p1-p0;
          
          const SIMD_IntegrationRule * ir1d = nullptr;
          for (int dir = 0; dir < 2; dir++)
            {
              if (delta1(dir) == 0)
                {
                  if (p0(dir) < 0.5)
                    ir1d = intrules.simd_intrule0;
                  else
                    ir1d = intrules.simd_intrule1;
                }
              else
                {
                  if (delta1(dir) < 0)
                    ir1d = &irfacet; // .GetIRX();
                  else
                    {
                      int order = 2*irfacet.GetNIP()-1;
                      ir1d = intrules.simd_segmentrules_inv[order];
                    }
                  
                  // ir1d = &irfacet; // .GetIRX();
                }
              switch (dir)
                {
                case 0:
                  irvol.SetIRX (ir1d); break;
                case 1:
                  irvol.SetIRY (ir1d); break;
                }
            }
          /*
            irvol.SetIRX(nullptr);
            irvol.SetIRY(nullptr);
            irvol.SetIRZ(nullptr);
          */
          break;
        }
        
          
      case ET_HEX:
        {
          if (!irfacet.IsTP()) break;
          FlatVec<3> p0 = points(faces[fnr][0]);
          FlatVec<3> p1 = points(faces[fnr][1]);
          FlatVec<3> p2 = points(faces[fnr][3]);
          Vec<3> delta1 = p1-p0;
          Vec<3> delta2 = p2-p0;
          
          const SIMD_IntegrationRule * ir1d = nullptr;
          for (int dir = 0; dir < 3; dir++)
            {
              if (delta1(dir) == 0 && delta2(dir) == 0)
                {
                  if (p0(dir) < 0.5)
                    ir1d = intrules.simd_intrule0;
                  else
                    ir1d = intrules.simd_intrule1;
                }
              else
                {
                  if (delta1(dir) > 0 || delta2(dir) > 0)
                    ir1d = &irfacet.GetIRX();
                  else
                    {
                      int order = 2*irfacet.GetIRX().GetNIP()-1;
                      ir1d = intrules.simd_segmentrules_inv[order];
                    }
                }
              switch (dir)
                {
                case 0:
                  irvol.SetIRX (ir1d); break;
                case 1:
                  irvol.SetIRY (ir1d); break;
                case 2:
                  irvol.SetIRZ (ir1d); break;
                }
            }
          /*
            irvol.SetIRX(nullptr);
            irvol.SetIRY(nullptr);
            irvol.SetIRZ(nullptr);
          */
          break;
        }
      default:
        ;
      }
    
    return irvol;
  }




  SIMD_IntegrationRule & Facet2SurfaceElementTrafo :: operator() (const SIMD_IntegrationRule & irfacet, LocalHeap & lh)
  {
    SIMD_IntegrationRule & irvol = *new (lh) SIMD_IntegrationRule (irfacet.GetNIP(), lh);
    int fnr = 0;
    
    FlatArray<SIMD<IntegrationPoint>> hirfacet = irfacet;
    FlatArray<SIMD<IntegrationPoint>> hirvol = irvol;
    
    switch (eltype) // ElementTopology::GetFacetType(eltype, fnr))
      {
      case ET_POINT:
        {
          irvol[0](0) = points (fnr)[0];
          irvol[0](1) = 0.0;
          irvol[0](2) = 0.0;
          break;
        }
      case ET_SEGM:
        {
          Vec<2> p1 = points (edges[fnr][0]);
          Vec<2> p2 = points (edges[fnr][1]);
          Vec<2> delta = p1-p2;
          for (int i = 0; i < hirfacet.Size(); i++)
            {
              auto ip = hirfacet[i];
              auto & ipvol = hirvol[i];              
              for (int k = 0; k < 2; k++)
                ipvol(k) = p2(k) + delta(k) * ip(0);
              ipvol(2) = 0.0;
            }
          break;
        }

      case ET_TRIG:
        {
          Vec<3> p0 = points(faces[fnr][0]);
          Vec<3> p1 = points(faces[fnr][1]);
          Vec<3> p2 = points(faces[fnr][2]);
          Vec<3> delta1 = p0-p2;
          Vec<3> delta2 = p1-p2;

          /*
          for (int i = 0; i < irfacet.GetNIP(); i++)
            irvol[i] = Vec<3> (p2 + irfacet[i](0) * (p0-p2) + irfacet[i](1)*(p1-p2));
          */
          for (int i = 0; i < hirfacet.Size(); i++)
            {
              auto ip = hirfacet[i];
              auto & ipvol = hirvol[i];              
              for (int k = 0; k < 3; k++)
                ipvol(k) = p2(k) + delta1(k) * ip(0) + delta2(k) * ip(1);
            }
          break;
        }

	case ET_QUAD:
	  {
	    FlatVec<3> p0 = points(faces[fnr][0]);
	    FlatVec<3> p1 = points(faces[fnr][1]);
	    FlatVec<3> p2 = points(faces[fnr][3]);
            Vec<3> delta1 = p1-p0;
            Vec<3> delta2 = p2-p0;

	    // for (int i = 0; i < irfacet.GetNIP(); i++)
            // irvol[i] = Vec<3> (p0 + irfacet[i](0) * (p1-p0) + irfacet[i](1)*(p2-p0));
            for (int i = 0; i < hirfacet.Size(); i++)
              {
                auto ip = hirfacet[i];
                auto & ipvol = hirvol[i];              
                for (int k = 0; k < 3; k++)
                  ipvol(k) = p0(k) + delta1(k) * ip(0) + delta2(k) * ip(1);
              }
	    break;
	  }

      default:
        throw Exception ("undefined facet type in SIMD_Facet2SurfaceElementTrafo()\n");
      } 
    
    for (int i = 0; i < irfacet.Size(); i++)
      {
        irvol[i].SetFacetNr(fnr);
        irvol[i].Weight() = irfacet[i].Weight();
      }

    return irvol;
  }

  
  IntegrationRule & Facet2SurfaceElementTrafo ::
  Inverse (const IntegrationRule & ir, LocalHeap & lh)
  {
    IntegrationRule & irinv = *new (lh) IntegrationRule (ir.GetNIP(), lh);

    int fnr = 0;
    
    FlatArray<IntegrationPoint> hir = ir;
    FlatArray<IntegrationPoint> hirinv = irinv;
    
    switch (eltype)
      {
        /*
      case ET_POINT:
        {
          irvol[0](0) = points (fnr)[0];
          irvol[0](1) = 0.0;
          irvol[0](2) = 0.0;
          break;
        }
        */
      case ET_SEGM:
        {
          Vec<2> p1 = points (edges[fnr][0]);
          Vec<2> p2 = points (edges[fnr][1]);
          double deltainv = p1(0)-p2(0); // +1 or -1, do not need to invert!
          for (int i = 0; i < hir.Size(); i++)
            {
              auto ip = hir[i];
              auto & ipinv = hirinv[i];
              ipinv(0) = deltainv*(ip(0) - p2(0));
              ipinv(1) = 0.0;
              ipinv(2) = 0.0;
            }
          break;
        }
      case ET_TRIG:
        {
          Vec<3> p0 = points(faces[fnr][0]);
          Vec<3> p1 = points(faces[fnr][1]);
          Vec<3> p2 = points(faces[fnr][2]);
          Vec<3> delta1 = p0-p2;
          Vec<3> delta2 = p1-p2;

          Vec<2> p22d(p2(0), p2(1));
          Mat<2> mat;
          mat(0,0) = delta1(0);   mat(0,1) = delta2(0);
          mat(1,0) = delta1(1);   mat(1,1) = delta2(1);
          mat = Inv(mat);

          for (int i = 0; i < hir.Size(); i++)
            {
              auto ip = hir[i];
              auto & ipinv = hirinv[i];
              Vec<2> vip(ip(0), ip(1));
              Vec<2> vipinv = mat * (vip - p22d);
              ipinv(0) = vipinv(0);
              ipinv(1) = vipinv(1);
              ipinv(2) = 0;
              /*
              for (int k = 0; k < 3; k++)
                ipinv(k) = p2(k) + delta1(k) * ip(0) + delta2(k) * ip(1);
              */
            }
          break;
        }

        
      case ET_QUAD:
	{
	  FlatVec<3> p0 = points(faces[fnr][0]);
	  FlatVec<3> p1 = points(faces[fnr][1]);
	  FlatVec<3> p2 = points(faces[fnr][3]);
          Vec<3> delta1 = p1-p0;
          Vec<3> delta2 = p2-p0;

          Vec<2> p02d(p0(0), p0(1));
          Mat<2> mat;
          mat(0,0) = delta1(0);   mat(0,1) = delta2(0);
          mat(1,0) = delta1(1);   mat(1,1) = delta2(1);
          mat = Inv(mat);
          for (int i = 0; i < hir.Size(); i++)
            {
              auto ip = hir[i];
              auto & ipinv = hirinv[i];
              Vec<2> vip(ip(0), ip(1));
              Vec<2> vipinv = mat * (vip - p02d);
              ipinv(0) = vipinv(0);
              ipinv(1) = vipinv(1);
              ipinv(2) = 0;
            }
	  break;
	}
        
      default:
        throw Exception ("undefined facet type in SIMD_Facet2SurfaceElementTrafo()\n");
      } 
    
    for (int i = 0; i < ir.Size(); i++)
      {
        irinv[i].SetFacetNr(fnr);
        irinv[i].SetWeight(ir[i].Weight());
      }


    
    return irinv;
  }
  
  
  SIMD_IntegrationRule & Facet2SurfaceElementTrafo ::
  Inverse (const SIMD_IntegrationRule & ir, LocalHeap & lh)
  {
    SIMD_IntegrationRule & irinv = *new (lh) SIMD_IntegrationRule (ir.GetNIP(), lh);

    int fnr = 0;
    
    FlatArray<SIMD<IntegrationPoint>> hir = ir;
    FlatArray<SIMD<IntegrationPoint>> hirinv = irinv;
    
    switch (eltype)
      {
        /*
      case ET_POINT:
        {
          irvol[0](0) = points (fnr)[0];
          irvol[0](1) = 0.0;
          irvol[0](2) = 0.0;
          break;
        }
        */
      case ET_SEGM:
        {
          Vec<2> p1 = points (edges[fnr][0]);
          Vec<2> p2 = points (edges[fnr][1]);
          double deltainv = p1(0)-p2(0); // +1 or -1, do not need to invert!
          for (int i = 0; i < hir.Size(); i++)
            {
              auto ip = hir[i];
              auto & ipinv = hirinv[i];
              ipinv(0) = deltainv*(ip(0) - p2(0));
              ipinv(1) = 0.0;
              ipinv(2) = 0.0;
            }
          break;
        }
      case ET_TRIG:
        {
          Vec<3> p0 = points(faces[fnr][0]);
          Vec<3> p1 = points(faces[fnr][1]);
          Vec<3> p2 = points(faces[fnr][2]);
          Vec<3> delta1 = p0-p2;
          Vec<3> delta2 = p1-p2;

          Vec<2> p22d(p2(0), p2(1));
          Mat<2> mat;
          mat(0,0) = delta1(0);   mat(0,1) = delta2(0);
          mat(1,0) = delta1(1);   mat(1,1) = delta2(1);
          mat = Inv(mat);

          for (int i = 0; i < hir.Size(); i++)
            {
              auto ip = hir[i];
              auto & ipinv = hirinv[i];
              Vec<2,SIMD<double>> vip(ip(0), ip(1));
              Vec<2,SIMD<double>> vipinv = mat * (vip - p22d);
              ipinv(0) = vipinv(0);
              ipinv(1) = vipinv(1);
              ipinv(2) = 0;
              /*
              for (int k = 0; k < 3; k++)
                ipinv(k) = p2(k) + delta1(k) * ip(0) + delta2(k) * ip(1);
              */
            }
          break;
        }

      case ET_QUAD:
	{
	  FlatVec<3> p0 = points(faces[fnr][0]);
	  FlatVec<3> p1 = points(faces[fnr][1]);
	  FlatVec<3> p2 = points(faces[fnr][3]);
          Vec<3> delta1 = p1-p0;
          Vec<3> delta2 = p2-p0;

          Vec<2> p02d(p0(0), p0(1));
          Mat<2> mat;
          mat(0,0) = delta1(0);   mat(0,1) = delta2(0);
          mat(1,0) = delta1(1);   mat(1,1) = delta2(1);
          mat = Inv(mat);
          for (int i = 0; i < hir.Size(); i++)
            {
              auto ip = hir[i];
              auto & ipinv = hirinv[i];
              Vec<2,SIMD<double>> vip(ip(0), ip(1));
              Vec<2,SIMD<double>> vipinv = mat * (vip - p02d);
              ipinv(0) = vipinv(0);
              ipinv(1) = vipinv(1);
              ipinv(2) = 0;
            }
	  break;
	}
      default:
        throw Exception ("undefined facet type in SIMD_Facet2SurfaceElementTrafo()\n");
      } 
    
    for (int i = 0; i < ir.Size(); i++)
      {
        irinv[i].SetFacetNr(fnr);
        irinv[i].Weight() = ir[i].Weight();
      }


    return irinv;
  }
}


namespace ngcore
{

  FlatVector<SIMD<double>>
  SIMD<ngfem::BaseMappedIntegrationPoint> :: GetPoint() const
  {
    switch (eltrans->SpaceDim())
      {
      case 1: return static_cast<const SIMD<ngfem::DimMappedIntegrationPoint<1>>&> (*this).GetPoint();
      case 2: return static_cast<const SIMD<ngfem::DimMappedIntegrationPoint<2>>&> (*this).GetPoint();
      case 3: return static_cast<const SIMD<ngfem::DimMappedIntegrationPoint<3>>&> (*this).GetPoint();
      }
    throw Exception("BaseMappedIntegrationPoint::GetPoint, illegal dimension");
  }

  int SIMD<ngfem::BaseMappedIntegrationPoint> ::  DimElement() const
  {
    return eltrans->ElementDim();
  }

  int SIMD<ngfem::BaseMappedIntegrationPoint> ::  DimSpace() const
  {
    return eltrans->SpaceDim();
  }

  
  template <int DIMS, int DIMR>
  void SIMD<ngfem::MappedIntegrationPoint<DIMS,DIMR>>::
  CalcHesse (Vec<DIMR,Mat<DIMS,DIMS,SIMD<double>>> & ddx) const
  {
    double eps = 1e-6;
    LocalHeapMem<2500*SIMD<double>::Size()> lh("calchesse");
    ngfem::SIMD_IntegrationRule ir(2*SIMD<double>::Size(), lh);  // number of scalar pnts
    
    for (int dir = 0; dir < DIMS; dir++)
      {
        ir[0] = this->IP();
        ir[0](dir) += eps;
        ir[1] = this->IP();
        ir[1](dir) -= eps;
        ngfem::SIMD_MappedIntegrationRule<DIMS,DIMR> mir(ir, *this->eltrans, lh);
        auto jacr = mir[0].GetJacobian();
        auto jacl = mir[1].GetJacobian();
        for (int k = 0; k < DIMR; k++)
          for (int j = 0; j < DIMS; j++)
            ddx(k)(dir,j) = (jacr(k,j) - jacl(k,j) ) / (2*eps);
      }
  }
  
  template class SIMD<ngfem::MappedIntegrationPoint<1,1>>;
  template class SIMD<ngfem::MappedIntegrationPoint<2,2>>;
  template class SIMD<ngfem::MappedIntegrationPoint<3,3>>;

  template class SIMD<ngfem::MappedIntegrationPoint<1,2>>;
  template class SIMD<ngfem::MappedIntegrationPoint<1,3>>;
  template class SIMD<ngfem::MappedIntegrationPoint<2,3>>;

  
}

/*
namespace ngstd
{
  SIMD<ngfem::BaseMappedIntegrationPoint>::~SIMD<ngfem::BaseMappedIntegrationPoint>()
  {
    // if (owns_trafo) delete eltrans; 
  }
  
}
*/
