/*********************************************************************/
/* File:   meshaccess.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   Access to fe mesh
*/

#include <comp.hpp>
#include "../fem/h1lofe.hpp"
#include <regex>

namespace ngcomp
{
  
  template <int DIMS, int DIMR>
  class Ng_ElementTransformation : public ElementTransformation
  {
    const MeshAccess * mesh;	

  public:
    Ng_ElementTransformation (const MeshAccess * amesh,
                              ELEMENT_TYPE aet, ElementId ei, int elindex) 
      : ElementTransformation(aet, ei, elindex), 
        mesh(amesh) 
    {
      iscurved = true;
    }

    virtual int SpaceDim () const
    {
      return DIMR;
    }

    virtual VorB VB() const
    {
      return VorB(int(DIMR)-int(DIMS));
    }

    virtual bool BelongsToMesh (const void * mesh2) const 
    {
      // return mesh == &(static_cast<const MeshAccess*> (mesh2) -> mesh);
      return mesh == mesh2;
    }

    virtual const void * GetMesh () const { return mesh; }

    virtual void GetSort (FlatArray<int> sort) const
    {
      int vnums[12];

      Ngs_Element nel = mesh -> GetElement<DIMS, (DIMS==DIMR)?VOL:BND> (elnr);
      for (int j = 0; j  < nel.vertices.Size(); j++)
        vnums[j] = nel.vertices[j];

      switch (eltype)
	{
	case ET_TRIG:
	  for (int i = 0; i < 3; i++) sort[i] = i;
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  // vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]]
	  break; 

	case ET_TET:
	  for (int i = 0; i < 4; i++) sort[i] = i;
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  if (vnums[sort[2]] > vnums[sort[3]]) Swap (sort[2], sort[3]);
	  if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
	  if (vnums[sort[1]] > vnums[sort[3]]) Swap (sort[1], sort[3]);
	  if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);

	  // vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]] < vnums[sort[3]]
	  break; 

	case ET_PRISM:
	  for (int i = 0; i < 6; i++) sort[i] = i;

	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
        
	  if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);
	  if (vnums[sort[4]] > vnums[sort[5]]) Swap (sort[4], sort[5]);
	  if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);
	  break;

	default:
	  throw Exception ("undefined eltype in ElementTransformation::GetSort()\n");
	}
      
    }


    virtual void CalcJacobian (const IntegrationPoint & ip,
			       FlatMatrix<> dxdxi) const
    {
      mesh->mesh.ElementTransformation <DIMS,DIMR> (elnr, &ip(0), NULL, &dxdxi(0));
    }
    
    virtual void CalcPoint (const IntegrationPoint & ip,
			    FlatVector<> point) const
    {
      mesh->mesh.ElementTransformation <DIMS,DIMR> (elnr, &ip(0), &point(0), NULL);
    }

    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, FlatMatrix<> dxdxi) const
    {
      mesh->mesh.ElementTransformation <DIMS,DIMR> (elnr, &ip(0), &point(0), &dxdxi(0));
    }

    virtual BaseMappedIntegrationPoint & operator() (const IntegrationPoint & ip, Allocator & lh) const
    {
      return *new (lh) MappedIntegrationPoint<DIMS,DIMR> (ip, *this);
    }

    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, Allocator & lh) const
    {
      return *new (lh) MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }

    virtual SIMD_BaseMappedIntegrationRule & operator() (const SIMD_IntegrationRule & ir, Allocator & lh) const
    {
      return *new (lh) SIMD_MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }

    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & bmir) const
    {
      if (sizeof(IntegrationPoint) % 8 != 0)
        {
          cerr << "Integration must should have 8-byte alignment" << endl;
          exit(1);
        }

      // static Timer t("eltrans::multipointjacobian"); RegionTimer reg(t);
      MappedIntegrationRule<DIMS,DIMR> & mir = 
	static_cast<MappedIntegrationRule<DIMS,DIMR> &> (bmir);
      
      mesh->mesh.MultiElementTransformation <DIMS,DIMR>
        (elnr, ir.Size(),
         &ir[0](0), ir.Size()>1 ? &ir[1](0)-&ir[0](0) : 0,
         &mir[0].Point()(0), ir.Size()>1 ? &mir[1].Point()(0)-&mir[0].Point()(0) : 0, 
         &mir[0].Jacobian()(0,0), ir.Size()>1 ? &mir[1].Jacobian()(0,0)-&mir[0].Jacobian()(0,0) : 0);

      /*
      for (int i = 0; i < ir.Size(); i++)
        mir[i].Compute();
      */
      for (auto & mip : mir)
        mip.Compute();
  }


    virtual void CalcMultiPointJacobian (const SIMD_IntegrationRule & ir,
					 SIMD_BaseMappedIntegrationRule & bmir) const
    {
      // static Timer t("eltrafo - nonconst, calcmultipoint"); RegionTimer reg(t);
      // t.AddFlops (ir.GetNIP());
      if (sizeof(IntegrationPoint) % 8 != 0)
        {
          cerr << "Integration must should have 8-byte alignment" << endl;
          exit(1);
        }

      // static Timer t("eltrans::multipointjacobian"); RegionTimer reg(t);
      SIMD_MappedIntegrationRule<DIMS,DIMR> & mir = 
	static_cast<SIMD_MappedIntegrationRule<DIMS,DIMR> &> (bmir);
      
      mesh->mesh.MultiElementTransformation <DIMS,DIMR>
        (elnr, ir.Size(),
         &ir[0](0).Data(), ir.Size()>1 ? &ir[1](0)-&ir[0](0) : 0,
         &mir[0].Point()(0).Data(), ir.Size()>1 ? &mir[1].Point()(0)-&mir[0].Point()(0) : 0, 
         &mir[0].Jacobian()(0,0).Data(), ir.Size()>1 ? &mir[1].Jacobian()(0,0)-&mir[0].Jacobian()(0,0) : 0);
      
      for (int i = 0; i < ir.Size(); i++)
        mir[i].Compute();
    }

  };
  




  template <int DIMS, int DIMR, typename BASE>
  class ALE_ElementTransformation : public BASE
  {
    const GridFunction * deform;
    const ScalarFiniteElement<DIMS> * fel;
    FlatVector<> elvec;
    FlatMatrix<> elvecs;
  public:
    ALE_ElementTransformation (const MeshAccess * amesh, 
                               ELEMENT_TYPE aet, ElementId ei, int elindex,
                               const GridFunction * adeform,
                               LocalHeap & lh)
      : BASE(amesh, aet, ei, elindex), 
        deform(adeform) 
    {
      fel = dynamic_cast<const ScalarFiniteElement<DIMS>*> (&deform->GetFESpace()->GetFE(ei, lh));

      Array<int> dnums(fel->GetNDof(), lh);
      deform->GetFESpace()->GetDofNrs(ei, dnums);

      elvec.AssignMemory(DIMR*dnums.Size(), lh);
      deform->GetElementVector(dnums, elvec);

      elvecs.AssignMemory(DIMR, dnums.Size(), lh);
      for (int j = 0; j < DIMR; j++)
        elvecs.Row(j) = elvec.Slice(j,DIMR);
    }

    /*
    virtual void SetElement (bool aboundary, int aelnr, int aelindex)
    {
      Ng_ElementTransformation<DIMS,DIMR> :: SetElement (aboundary, aelnr, aelindex);

      ElementId id(aboundary ? BND : VOL, aelnr);
      fel = dynamic_cast<const ScalarFiniteElement<DIMR>*> (&deform->GetFESpace()->GetFE(id, lh));

      Array<int> dnums;
      deform->GetFESpace()->GetDofNrs(id, dnums);

      elvec.AssignMemory(DIMR*dnums.Size(), lh);
      deform->GetElementVector(dnums, elvec);

      elvecs.AssignMemory(DIMR, dnums.Size(), lh);
      for (int j = 0; j < DIMR; j++)
        elvecs.Row(j) = elvec.Slice(j,DIMR);
    }
    */

    virtual ~ALE_ElementTransformation() { ; }



    virtual void CalcJacobian (const IntegrationPoint & ip,
			       FlatMatrix<> dxdxi) const
    {
      /*
      Mat<DIMR,DIMS> tmp, itmp;
      Ng_ElementTransformation<DIMS,DIMR>::CalcJacobian (ip, tmp);
      itmp = Inv (tmp);

      Mat<DIMR,DIMR> def, mdef;
      for (int i = 0; i < DIMR; i++)
        def.Row(i) = fel->EvaluateGrad (ip, elvecs.Row(i));
      
      mdef = def*itmp;
      mdef += Id<DIMR>();
      dxdxi = mdef * tmp;
      */

      Mat<DIMR,DIMS> tmp;
      BASE::CalcJacobian (ip, tmp);

      Mat<DIMR,DIMS> def;
      for (int i = 0; i < DIMR; i++)
        def.Row(i) = fel->EvaluateGrad (ip, elvecs.Row(i));
      dxdxi = def + tmp;
    }
    
    virtual void CalcPoint (const IntegrationPoint & ip,
			    FlatVector<> point) const
    {
      Vec<DIMR> tmp;
      BASE::CalcPoint (ip, tmp);

      Vec<DIMR> def;
      for (int i = 0; i < DIMR; i++)
        def(i) = fel->Evaluate (ip, elvecs.Row(i));
      point = tmp + def;
    }

    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, FlatMatrix<> dxdxi) const
    {
      CalcJacobian (ip, dxdxi);
      CalcPoint (ip, point);
      // this->mesh->mesh.ElementTransformation <DIMS,DIMR> (elnr, &ip(0), &point(0), &dxdxi(0));
    }

    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & bmir) const
    {
      if (sizeof(IntegrationPoint) % 8 != 0)
        {
          cerr << "Integration must should have 8-byte alignment" << endl;
          exit(1);
        }

      MappedIntegrationRule<DIMS,DIMR> & mir = 
	static_cast<MappedIntegrationRule<DIMS,DIMR> &> (bmir);

      for (int i = 0; i < ir.Size(); i++)
        {
          CalcPointJacobian (ir[i], mir[i].Point(), mir[i].Jacobian());
          mir[i].Compute();
        }

      /*
      MappedIntegrationRule<DIMS,DIMR> & mir = 
	static_cast<MappedIntegrationRule<DIMS,DIMR> &> (bmir);
      mesh->mesh.MultiElementTransformation <DIMS,DIMR> (elnr, ir.Size(),
                                                         &ir[0](0), &ir[1](0)-&ir[0](0),
                                                         &mir[0].Point()(0), 
                                                         &mir[1].Point()(0)-&mir[0].Point()(0), 
                                                         &mir[0].Jacobian()(0,0), 
                                                         &mir[1].Jacobian()(0,0)-&mir[0].Jacobian()(0,0));
    
      for (int i = 0; i < ir.Size(); i++)
        mir[i].Compute();
      */
    }
    
    virtual void CalcMultiPointJacobian (const SIMD_IntegrationRule & ir,
					 SIMD_BaseMappedIntegrationRule & bmir) const
    {
      SIMD_MappedIntegrationRule<DIMS,DIMR> & mir = 
	static_cast<SIMD_MappedIntegrationRule<DIMS,DIMR> &> (bmir);
      
      BASE::CalcMultiPointJacobian (ir, bmir);

      // LocalHeapMem<100000> lh("tmp");
      // AFlatVector<double> def(ir.GetNIP(), lh);
      // AFlatMatrix<double> grad(DIMS, ir.GetNIP(), lh);
      STACK_ARRAY(SIMD<double>, mem0, ir.Size());
      FlatVector<SIMD<double>> def(ir.Size(), &mem0[0]);
      STACK_ARRAY(SIMD<double>, mem1, (DIMS*ir.Size()));
      FlatMatrix<SIMD<double>> grad(DIMS, ir.Size(), &mem1[0]);

      for (int i = 0; i < DIMR; i++)
        {
          fel->Evaluate (ir, elvec.Slice(i,DIMR), def);
          fel->EvaluateGrad (ir, elvec.Slice(i,DIMR), grad);
          
          for (size_t k = 0; k < ir.Size(); k++)
            {
              mir[k].Point()(i) += def(k);
              for (int j = 0; j < DIMS; j++)
                mir[k].Jacobian()(i,j) += grad(j,k);
            }
        }
      
      for (int i = 0; i < ir.Size(); i++)
        mir[i].Compute();
    }


  };


  template <int DIMS, int DIMR, typename BASE>
  class PML_ElementTransformation : public BASE
  {
    const PML_Transformation & pml_global_trafo;
  public:
    PML_ElementTransformation (const MeshAccess * amesh, 
                                     ELEMENT_TYPE aet, ElementId ei, int elindex,
                                     const PML_Transformation & _pml_global_trafo)
      : BASE(amesh, aet, ei, elindex), pml_global_trafo(_pml_global_trafo)
    {
      this->is_complex = true;
    }

    virtual ~PML_ElementTransformation() { ; }

    virtual void CalcJacobian (const IntegrationPoint & ip,
			       FlatMatrix<> dxdxi) const
    {
      BASE::CalcJacobian(ip,dxdxi);
    }
    
    virtual void CalcPoint (const IntegrationPoint & ip,
			    FlatVector<> point) const
    {
      BASE::CalcPoint(ip,point);
    }

    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, FlatMatrix<> dxdxi) const
    {
      CalcJacobian (ip, dxdxi);
      CalcPoint (ip, point);
    }


    virtual BaseMappedIntegrationPoint & operator() (const IntegrationPoint & ip, Allocator & lh) const
    {
      auto & mip = *new (lh) MappedIntegrationPoint<DIMS,DIMR,Complex> (ip, *this, -47);
      MappedIntegrationPoint<DIMS,DIMR> hip(ip, *this);
      Vec<DIMR,Complex> point;
      
      Mat<DIMS,DIMR> hjac(hip.Jacobian());
          
      Mat<DIMR,DIMR,Complex> tjac;
      const PML_TransformationDim<DIMR> & dimpml = 
        static_cast<const PML_TransformationDim<DIMR>&> (pml_global_trafo);
      dimpml.MapIntegrationPoint (hip, point, tjac);
                    
      mip.Point() = point; 
      mip.Jacobian() = tjac*hjac;
      mip.Compute();
      return mip;
    }

    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, Allocator & lh) const
    {
      return *new (lh) MappedIntegrationRule<DIMS,DIMR,Complex> (ir, *this, lh);
    }

    virtual SIMD_BaseMappedIntegrationRule & operator() (const SIMD_IntegrationRule & ir, Allocator & lh) const
    {
      throw ExceptionNOSIMD("PML-trafo operator() (SIMD_IntegrationRule & ir) not overloaded");
      // return *new (lh) SIMD_MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }

    
    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & bmir) const
    {
      if (!bmir.IsComplex())
        {
          BASE::CalcMultiPointJacobian (ir, bmir);
          return;
        }
      //LocalHeapMem<1000000> lh("testwise");
      //MappedIntegrationRule<DIMS,DIMR> mir_real(ir, *this, lh);

      auto & mir_complex = dynamic_cast<MappedIntegrationRule<DIMS,DIMR,Complex>&> (bmir);
      const PML_TransformationDim<DIMR> & dimpml = 
            static_cast<const PML_TransformationDim<DIMR>&> (pml_global_trafo);

      Vec<DIMR,Complex> point;
      Mat<DIMR,DIMR,Complex> tjac;
      for (int i = 0; i < ir.Size(); i++)
        {
          MappedIntegrationPoint<DIMS,DIMR> mip(ir[i], *this);
          Mat<DIMS,DIMR> hjac(mip.Jacobian());
          dimpml.MapIntegrationPoint (mip, point, tjac);
                    
          mir_complex[i].Point() = point; 
          mir_complex[i].Jacobian() = tjac*hjac;
          mir_complex[i].Compute();          
        }
    }
    
    virtual void CalcMultiPointJacobian (const SIMD_IntegrationRule & ir,
					 SIMD_BaseMappedIntegrationRule & bmir) const
    {
      throw ExceptionNOSIMD ("pml - trafo cannot calc simd_rule");
    }
  };
  
  template <int DIMS, int DIMR>
  class Ng_ConstElementTransformation : public ElementTransformation
  {
    const MeshAccess * mesh;
    Vec<DIMR> p0;
    Mat<DIMR,DIMS> mat;
  public:
    INLINE Ng_ConstElementTransformation (const MeshAccess * amesh,
                                          ELEMENT_TYPE aet, ElementId ei, int elindex) 
      : ElementTransformation(aet, ei, elindex), 
        mesh(amesh) 
    { 
      iscurved = false;
      if ( (DIMR==3) && (eltype == ET_TET) )
        {
          Ngs_Element nel = mesh -> GetElement<DIMS,VOL> (elnr);
          // p0 = FlatVec<3> (point0[point_delta*nel.Vertices()[0]]);
          p0 = FlatVec<3, const double> (mesh->mesh.GetPoint (nel.Vertices()[3]));
	  for (int j = 0; j < 3; j++)
	    {
	      Vec<3> pj = FlatVec<3, const double>(mesh->mesh.GetPoint(nel.Vertices()[j])) - p0;
	      for (int k = 0; k < 3; k++)
		mat(k,j) = pj(k);
	    }
        }
      
      else if ( (DIMR==2) && (DIMS==2) && (eltype == ET_TRIG) )
        {
          Ngs_Element nel = mesh -> GetElement<DIMS,VOL> (elnr);
          Vec<2> hp0 = FlatVec<2, const double> (mesh->mesh.GetPoint (nel.Vertices()[2]));
          p0 = hp0;
	  for (int j = 0; j < 2; j++)
	    {
	      Vec<2> pj = FlatVec<2, const double>(mesh->mesh.GetPoint(nel.Vertices()[j])) - hp0;
	      for (int k = 0; k < 2; k++)
		mat(k,j) = pj(k);
	    }
        }

      else if ( (DIMR==2) && (DIMS==1) && (eltype == ET_SEGM) )
        {
          Ngs_Element nel = mesh -> GetElement<DIMS,VOL> (elnr);
          p0 = FlatVec<2, const double> (mesh->mesh.GetPoint (nel.Vertices()[1]));
	  for (int j = 0; j < 1; j++)
	    {
	      Vec<2> pj = FlatVec<2, const double>(mesh->mesh.GetPoint(nel.Vertices()[j])) - p0;
	      for (int k = 0; k < 2; k++)
		mat(k,j) = pj(k);
	    }
          //mat.Col(0) = FlatVec<3, const double> (mesh -> GetPoint (nel.Vertices()[0])) - p0;
          //mat.Col(1) = FlatVec<3, const double> (mesh -> GetPoint (nel.Vertices()[1])) - p0;
          //mat.Col(2) = FlatVec<3, const double> (mesh -> GetPoint (nel.Vertices()[2])) - p0;
        }

      else if ( (DIMR==1) && (DIMS==1) && (eltype == ET_SEGM) )
        {
          Ngs_Element nel = mesh -> GetElement<DIMS,VOL> (elnr);
          p0 = FlatVec<1, const double> (mesh->mesh.GetPoint (nel.Vertices()[1]));
	  for (int j = 0; j < 1; j++)
	    {
	      Vec<1> pj = FlatVec<1, const double>(mesh->mesh.GetPoint(nel.Vertices()[j])) - p0;
	      for (int k = 0; k < 1; k++)
		mat(k,j) = pj(k);
	    }
        }

      else
        {
          Vec<DIMS> pref = 0.0;
          mesh->mesh.ElementTransformation <DIMS,DIMR> (elnr, &pref(0), &p0(0), &mat(0));
        }
    }

    /*
    virtual void SetElement (bool aboundary, int aelnr, int aelindex)
    {
      elnr = aelnr;
      elindex = aelindex;
      iscurved = false;

      if ( (DIMR==3) && (eltype == ET_TET) )
        {
          Ngs_Element nel = mesh -> GetElement<DIMS,VOL> (elnr);
          // p0 = FlatVec<3> (point0[point_delta*nel.Vertices()[0]]);
          p0 = FlatVec<3, const double> (mesh->mesh.GetPoint (nel.Vertices()[3]));
	  for (int j = 0; j < 3; j++)
	    {
	      Vec<3> pj = FlatVec<3, const double>(mesh->mesh.GetPoint(nel.Vertices()[j])) - p0;
	      for (int k = 0; k < 3; k++)
		mat(k,j) = pj(k);
	    }
          //mat.Col(0) = FlatVec<3, const double> (mesh -> GetPoint (nel.Vertices()[0])) - p0;
          //mat.Col(1) = FlatVec<3, const double> (mesh -> GetPoint (nel.Vertices()[1])) - p0;
          //mat.Col(2) = FlatVec<3, const double> (mesh -> GetPoint (nel.Vertices()[2])) - p0;
        }
      
      else if ( (DIMR==2) && (DIMS==2) && (eltype == ET_TRIG) )
        {
          Ngs_Element nel = mesh -> GetElement<DIMS,VOL> (elnr);
          p0 = FlatVec<2, const double> (mesh->mesh.GetPoint (nel.Vertices()[2]));
	  for (int j = 0; j < 2; j++)
	    {
	      Vec<2> pj = FlatVec<2, const double>(mesh->mesh.GetPoint(nel.Vertices()[j])) - p0;
	      for (int k = 0; k < 2; k++)
		mat(k,j) = pj(k);
	    }
          //mat.Col(0) = FlatVec<3, const double> (mesh -> GetPoint (nel.Vertices()[0])) - p0;
          //mat.Col(1) = FlatVec<3, const double> (mesh -> GetPoint (nel.Vertices()[1])) - p0;
          //mat.Col(2) = FlatVec<3, const double> (mesh -> GetPoint (nel.Vertices()[2])) - p0;
        }


      else
        {
          Vec<DIMS> pref = 0.0;
          mesh->mesh.ElementTransformation <DIMS,DIMR> (elnr, &pref(0), &p0(0), &mat(0));
        }
    }
    */

    virtual int SpaceDim () const
    {
      return DIMR;
    }

    virtual VorB VB() const
    {
      return VorB(int(DIMR)-int(DIMS));
    }
    
    virtual bool BelongsToMesh (const void * mesh2) const 
    {
      // return mesh == &(static_cast<const MeshAccess*> (mesh2) -> mesh);
      return mesh == mesh2;
    }

    virtual const void * GetMesh () const { return mesh; }

    virtual void GetSort (FlatArray<int> sort) const
    {
      int vnums[12];

      Ngs_Element nel = mesh -> GetElement<DIMS, (DIMS==DIMR)?VOL:BND> (elnr);
      for (int j = 0; j  < nel.vertices.Size(); j++)
        vnums[j] = nel.vertices[j];

      switch (eltype)
	{
	case ET_TRIG:
	  for (int i = 0; i < 3; i++) sort[i] = i;
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  // vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]]
	  break; 

	case ET_TET:
	  for (int i = 0; i < 4; i++) sort[i] = i;
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  if (vnums[sort[2]] > vnums[sort[3]]) Swap (sort[2], sort[3]);
	  if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
	  if (vnums[sort[1]] > vnums[sort[3]]) Swap (sort[1], sort[3]);
	  if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);

	  // vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]] < vnums[sort[3]]
	  break; 

	case ET_PRISM:
	  for (int i = 0; i < 6; i++) sort[i] = i;

	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	  if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
	  if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
        
	  if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);
	  if (vnums[sort[4]] > vnums[sort[5]]) Swap (sort[4], sort[5]);
	  if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);
	  break;

	default:
	  throw Exception ("undefined eltype in ElementTransformation::GetSort()\n");
	}
      
    }


    virtual void CalcJacobian (const IntegrationPoint & ip,
			       FlatMatrix<> dxdxi) const
    {
      // dxdxi = mat;
      FlatMatrixFixWidth<DIMS> (DIMR, &dxdxi(0,0)) = mat;
    }
    
    virtual void CalcPoint (const IntegrationPoint & ip,
			    FlatVector<> point) const
    {
      // point = p0 + mat * FlatVec<DIMS, const double> (&ip(0));
      FlatVec<DIMR> (&point(0)) = p0 + mat * FlatVec<DIMS, const double> (&ip(0));
    }

    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, FlatMatrix<> dxdxi) const
    {
      FlatVec<DIMR> (&point(0)) = p0 + mat * FlatVec<DIMS, const double> (&ip(0));
      FlatMatrixFixWidth<DIMS> (DIMR, &dxdxi(0,0)) = mat;
    }

    virtual BaseMappedIntegrationPoint & operator() (const IntegrationPoint & ip, Allocator & lh) const
    {
      return *new (lh) MappedIntegrationPoint<DIMS,DIMR> (ip, *this);
    }

    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, Allocator & lh) const
    {
      return *new (lh) MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }    

    virtual SIMD_BaseMappedIntegrationRule & operator() (const SIMD_IntegrationRule & ir, Allocator & lh) const
    {
      return *new (lh) SIMD_MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }

    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & bmir) const
    {
      MappedIntegrationRule<DIMS,DIMR> & mir = static_cast<MappedIntegrationRule<DIMS,DIMR> &> (bmir);
      for (int i = 0; i < ir.Size(); i++)
        {
          const IntegrationPoint & ip = ir[i];
          mir[i].Point() = p0 + mat * FlatVec<DIMS, const double> (&ip(0));
          mir[i].Jacobian() = mat;
          mir[i].Compute();
        }
    }

    virtual void CalcMultiPointJacobian (const SIMD_IntegrationRule & ir,
					 SIMD_BaseMappedIntegrationRule & bmir) const
    {
      // static Timer t("eltrafo - const, calcmultipoint"); RegionTimer reg(t);
      // t.AddFlops (ir.GetNIP());

      SIMD_MappedIntegrationRule<DIMS,DIMR> & mir = static_cast<SIMD_MappedIntegrationRule<DIMS,DIMR> &> (bmir);
      FlatArray<SIMD<MappedIntegrationPoint<DIMS,DIMR>>> hmir(mir.Size(), &mir[0]);
      FlatArray<SIMD<IntegrationPoint>> hir (ir);
      
      Vec<DIMR,SIMD<double>> simd_p0(p0);
      Mat<DIMR,DIMS,SIMD<double>> simd_mat(mat);

      for (size_t i = 0; i < hir.Size(); i++)
        {
          hmir[i].Point() = simd_p0 + simd_mat * FlatVec<DIMS, const SIMD<double>> (&hir[i](0));
          hmir[i].Jacobian() = simd_mat;
        }

      // static Timer tcompute("eltrafo - const, compute");
      // RegionTimer r2(tcompute);
      for (int i = 0; i < hir.Size(); i++)
        hmir[i].Compute();
    }
  };
  














  MeshAccess :: MeshAccess (shared_ptr<netgen::Mesh> amesh)
    : mesh(amesh), mesh_comm(ngs_comm)
  {
    // the connection to netgen global variables
    ngstd::testout = netgen::testout;
    ngstd::printmessage_importance = netgen::printmessage_importance;
#ifdef PARALLEL
    // best we can do at the moment to get py-mpi running
    mesh_comm = MPI_COMM_WORLD;
    ngs_comm =  MPI_COMM_WORLD;
#endif
    mesh.SelectMesh();
    mesh.UpdateTopology();  // for netgen/ngsolve stand alone
    UpdateBuffers();
  }


  MeshAccess :: ~MeshAccess ()
  {
    // delete mesh;
    // Ng_LoadGeometry("");
  }


  void MeshAccess :: LoadMesh (const string & filename)
  {
    static Timer t("MeshAccess::LoadMesh"); RegionTimer reg(t);
    mesh.LoadMesh (filename);
    UpdateBuffers();
    if (!mesh.Valid())
      throw Exception ("could not load mesh from '" + filename + "'");
  }

  void MeshAccess :: LoadMesh (istream & str)
  {
    static Timer t("MeshAccess::LoadMesh"); RegionTimer reg(t);    
    mesh.LoadMesh (str);
    UpdateBuffers();
  }

  void MeshAccess :: SaveMesh (ostream & str) const
  {
    mesh.SaveMesh (str);
  }

  void MeshAccess :: ArchiveMesh (Archive & archive) 
  {
    mesh.DoArchive (archive);
    if (archive.Input()) UpdateBuffers();
  }

  void MeshAccess :: SelectMesh() const
  {
    mesh.SelectMesh();
  }


  /*
  void MeshAccess :: LoadMeshFromString(const string & str)
  {
    Ng_LoadMeshFromString(const_cast<char*>(str.c_str()));
    UpdateBuffers();
  }
  */

  void MeshAccess :: GetSElNeighbouringDomains(const int elnr, int & in, int & out) const
  {
    ArrayMem<int, 2> elnums;
    ArrayMem<int, 2> fnums;
    GetElFacets(ElementId(BND,elnr), fnums);
    GetFacetElements ( fnums[0], elnums );
    if (elnums.Size()==1)
    {
      in = GetElIndex(ElementId(VOL,elnums[0]))+1;
      out = 0;
    }
    else
    {
      out = GetElIndex(ElementId(VOL,elnums[0]))+1;
      in = GetElIndex(ElementId(VOL,elnums[1]))+1;
    }
  }


  void MeshAccess :: UpdateBuffers()
  {
    static Timer t("MeshAccess::UpdateBuffers");
    RegionTimer reg(t);

    timestamp = NGS_Object::GetNextTimeStamp();
    
    if (!mesh.Valid())
      {
        for (int i = 0; i < 4; i++)  
          {
            nnodes[i] = 0;
            nelements[i] = 0;
            nnodes_cd[i] = 0;
            nelements_cd[i] = 0;
          }
        nnodes[NT_ELEMENT] = 0;
        nnodes[NT_FACET] = 0;
        dim = -1;
        ne_vb[VOL] = ne_vb[BND] = ne_vb[BBND] = 0;
        return;
      }

    dim = mesh.GetDimension();
    nlevels = mesh.GetNLevels(); 

    if (MyMPI_GetNTasks() > 1 && MyMPI_GetId() == 0)
      {
        for (int i = 0; i < 4; i++)  
          {
            nnodes[i] = 0;
            nelements[i] = 0;
            nnodes_cd[i] = 0;
            nelements_cd[i] = 0;
          }
        ne_vb[VOL] = ne_vb[BND] = ne_vb[BBND] = 0;
      }
    else
      {
	for (int i = 0; i < 4; i++)  
	  {
	    nnodes[i] = mesh.GetNNodes(i);
	    nelements[i] = mesh.GetNElements(i);
	  }
	for (int i = 0; i <= dim; i++)
	  {
	    nnodes_cd[i] = nnodes[dim-i];
	    nelements_cd[i] = nelements[dim-i];
	  }
        ne_vb[VOL] = nelements_cd[0];
        ne_vb[BND] = nelements_cd[1];
	if(dim==1)
	  ne_vb[BBND] = 0;
	else 
	  ne_vb[BBND] = nelements_cd[2];
      }
    nnodes[NT_ELEMENT] = nnodes[StdNodeType (NT_ELEMENT, dim)];
    nnodes[NT_FACET] = nnodes[StdNodeType (NT_FACET, dim)];
    

    ndomains = -1;
    int ne = GetNE(); 
    for (int i = 0; i < ne; i++)
      {
        int elindex = GetElIndex(ElementId(VOL,i));
        if (elindex < 0) throw Exception("mesh with negative element-index");
        ndomains = max2(ndomains, elindex);
      }

    ndomains++;
    ndomains = MyMPI_AllReduce (ndomains, MPI_MAX);
    pml_trafos.SetSize(ndomains);

    nboundaries = -1;
    int nse = GetNSE(); 
    for (int i = 0; i < nse; i++)
      {
        ElementId sei(BND, i);
        int elindex = GetElIndex(sei);
        if (elindex < 0) throw Exception("mesh with negative boundary-condition number");
        nboundaries = max2(nboundaries, elindex);
      }

    nboundaries++;
    nboundaries = MyMPI_AllReduce (nboundaries, MPI_MAX);

    CalcIdentifiedFacets();
    if(mesh.GetDimension() == 1)
      {
        nbboundaries = 0;
      }
    else
      {
        nbboundaries = -1;
        int ncd2e = nelements_cd[2];
        for (int i=0; i< ncd2e; i++)
          {
            ElementId ei(BBND, i);
            int elindex = GetElIndex(ei);
            //if (elindex < 0) throw Exception ("mesh with negative cd2 condition number");
            if (elindex >=0)
              nbboundaries = max2(nbboundaries, elindex);
          }
        nbboundaries++;
        nbboundaries = MyMPI_AllReduce(nbboundaries, MPI_MAX);
      }
    
    // update periodic mappings
    auto nid = mesh.GetMesh()->GetIdentifications().GetMaxNr();
    periodic_node_pairs[NT_VERTEX]->SetSize(nid);
    periodic_node_pairs[NT_EDGE]->SetSize(nid);
    periodic_node_pairs[NT_FACE]->SetSize(nid);
    for (auto idnr : Range(1,nid+1))
      {
        // only if it is periodic
        if (mesh.GetMesh()->GetIdentifications().GetType(idnr)!=2) continue;
        size_t nverts = Ng_GetNPeriodicVertices(idnr); 
        (*periodic_node_pairs[NT_VERTEX])[idnr-1].SetSize(nverts);
        Ng_GetPeriodicVertices(idnr,&(*periodic_node_pairs[NT_VERTEX])[idnr-1][0][0]);
        for (auto i : Range(nverts))
          {
            (*periodic_node_pairs[NT_VERTEX])[idnr-1][i][0]--;
            (*periodic_node_pairs[NT_VERTEX])[idnr-1][i][1]--;
          }
      

        // build vertex map for idnr
        Array<int> vertex_map(GetNV());
        for (auto i : Range(GetNV()))
          vertex_map[i] = i;
        for (auto pair : (*periodic_node_pairs[NT_VERTEX])[idnr-1])
          vertex_map[pair[1]] = pair[0];
        
        // build vertex-pair to edge hashtable:
        HashTable<INT<2>, int> vp2e(GetNEdges());
        
        for (int enr = 0; enr < GetNEdges(); enr++)
          {
            int v1, v2;
            GetEdgePNums (enr, v1, v2);
            if (v1 > v2) Swap (v1, v2);
            vp2e[INT<2>(v1,v2)] = enr;
          }
        int count = 0;
        for (int enr = 0; enr < GetNEdges(); enr++)
          {
            int v1,v2;
            GetEdgePNums(enr,v1,v2);
            int mv1 = vertex_map[v1];
            int mv2 = vertex_map[v2];
            if(mv1 != v1 && mv2 != v2)
              count++;
          }
        (*periodic_node_pairs[NT_EDGE])[idnr-1].SetSize(count);
        count = 0;
        for (int enr = 0; enr < GetNEdges(); enr++)
          {
            int v1, v2;
            GetEdgePNums (enr, v1, v2);
            int mv1 = vertex_map[v1];
            int mv2 = vertex_map[v2];
            if(mv1 != v1 && mv2 != v2)
              {               
                if (mv1 > mv2) Swap(mv1,mv2);
                int menr = vp2e.Get(INT<2>(mv1,mv2));
                (*periodic_node_pairs[NT_EDGE])[idnr-1][count][0] = menr;
                (*periodic_node_pairs[NT_EDGE])[idnr-1][count++][1] = enr;
              }
          }
        // build vertex-triple to face hashtable
        HashTable<INT<3>, int> v2f(GetNFaces());
        Array<int> pnums;
        for (auto fnr : Range(GetNFaces()))
          {
            GetFacePNums (fnr, pnums);
            INT<3> i3(pnums[0], pnums[1], pnums[2]);
            i3.Sort();
            v2f[i3] = fnr;
          }

        count = 0;
        for (auto fnr : Range(GetNFaces()))
          {
            GetFacePNums(fnr,pnums);
            if(vertex_map[pnums[0]] != pnums[0] && vertex_map[pnums[1]] != pnums[1] &&
               vertex_map[pnums[2]] != pnums[2])
              {
                count++;
              }
          }
        (*periodic_node_pairs[NT_FACE])[idnr-1].SetSize(count);
        count = 0;
        for (auto fnr : Range(GetNFaces()))
          {
            GetFacePNums(fnr,pnums);
            INT<3> mv(vertex_map[pnums[0]],vertex_map[pnums[1]],vertex_map[pnums[2]]);
            if(mv[0] != pnums[0] && mv[1] != pnums[1] && mv[2] != pnums[2])
              {
                mv.Sort();
                int mfnr = v2f[mv];
                (*periodic_node_pairs[NT_FACE])[idnr-1][count][0] = mfnr;
                (*periodic_node_pairs[NT_FACE])[idnr-1][count++][1] = fnr;
              }
          }
      }
    
 }


#ifdef ABC
  void MeshAccess::GetTopologicElement (int elnr, TopologicElement & topel) const
  {
    int help[54];
    int nnodes;

    nnodes = Ng_GetElementClosureNodes (dim, elnr, 
                                        NodeSet (NT_VERTEX, NT_EDGE, NT_FACE, NT_CELL), help);
        
    topel.SetElementType (GetElType(elnr));
    topel.Clear();
    for (int i = 0; i < nnodes; i++)
      topel.AddNode (Node (NODE_TYPE (help[2*i]), help[2*i+1]));
    
    /*    
    // 2D case
    ArrayMem<int, 12> nums;
    
    topel.SetElementType (GetElType(elnr));
    topel.Clear();

    GetElVertices (elnr, nums);
    for (int j = 0; j < nums.Size(); j++)
    topel.AddNode (Node (NT_VERTEX, nums[j]));

    GetElEdges (elnr, nums);
    for (int j = 0; j < nums.Size(); j++)
    topel.AddNode (Node (NT_EDGE, nums[j]));

    if (GetDimension() == 3)
    {
    GetElFaces (elnr, nums);
    for (int j = 0; j < nums.Size(); j++)
    topel.AddNode (Node (NT_FACE, nums[j]));
        
    topel.AddNode (Node (NT_CELL, elnr));
    }
    else
    topel.AddNode (Node (NT_FACE, elnr));
    */
  }
#endif


#ifdef ABC
  void MeshAccess :: GetSegmentPNums (int snr, Array<int> & pnums) const
  {
    pnums = GetElement<1> (snr).Points();
    /*
    pnums.SetSize(3);
    int np;
    Ng_GetSegment (snr+1, &pnums[0], &np);
    pnums.SetSize(np);

    for (int i = 0; i < np; i++)
      pnums[i]--;
    */
  }


  int MeshAccess :: GetSegmentIndex (int snr) const
  {
    return mesh -> GetElementIndex<1> (snr);
    //return Ng_GetSegmentIndex(snr+1);
  }
#endif


  void MeshAccess :: 
  GetElEdges (int elnr, Array<int> & ednums, Array<int> & orient) const
  {
    ednums.SetSize (12);
    orient.SetSize (12);
    int ned = 
      Ng_GetElement_Edges (elnr+1, &ednums[0], &orient[0]);
    ednums.SetSize (ned);
    orient.SetSize (ned);
    for (int i = 0; i < ned; i++)
      ednums[i]--;
  }

  void MeshAccess :: 
  GetSElEdges (int selnr, Array<int> & ednums, Array<int> & orient) const
  {
    ednums.SetSize (4);
    orient.SetSize (4);
    int ned = 
      Ng_GetSurfaceElement_Edges (selnr+1, &ednums[0], &orient[0]);
    ednums.SetSize (ned);
    orient.SetSize (ned);

    for (int i = 0; i < ned; i++)
      ednums[i]--;
  }


  void MeshAccess :: GetEdgeElements (int enr, Array<int> & elnums) const
  {
    // static Timer t("getedgeelements"); RegionTimer reg(t);    
    elnums.SetSize0();
    /*
    int p0, p1;
    GetEdgePNums(enr, p0, p1);

    auto velems0 = GetVertexElements(p0);
    auto velems1 = GetVertexElements(p1);
    */
    auto vts = mesh.GetNode<1>(enr).vertices;
    auto velems0 = GetVertexElements(vts[0]);
    auto velems1 = GetVertexElements(vts[1]);
    
    /*
    // n^2 intersection 
    for (auto el : velems0)
      if (velems1.Contains(el))
        elnums.Append(el);
    */
    
    // intersect sorted arrays ...
    int j = 0;
    for (int i = 0; i < velems0.Size(); i++)
      for ( ; j < velems1.Size(); j++)
        {
          if (velems0[i] < velems1[j]) break;
          if (velems0[i] == velems1[j])
            elnums.Append (velems0[i]);
        }
  }

  void MeshAccess :: GetEdgeSurfaceElements (int enr, Array<int> & elnums) const
  {
    // static Timer t("getedgesurfelements"); RegionTimer reg(t);
    elnums.SetSize0();
    // ArrayMem<int,3> pnums;
    int p0, p1;
    // ArrayMem<int,50> velems0, velems1; 
    GetEdgePNums(enr, p0, p1);
    /*
    GetVertexSurfaceElements(p0, velems0);
    GetVertexSurfaceElements(p1, velems1);
    */
    auto velems0 = GetVertexSurfaceElements(p0);
    auto velems1 = GetVertexSurfaceElements(p1);
    // now compare
    for (int i=0; i<velems0.Size(); i++) 
      for (int j=0; j<velems1.Size(); j++) 
	if (velems0[i] == velems1[j]) 
	  {
	    elnums.Append(velems0[i]);
	    continue;
	  }
  }

  

  void MeshAccess :: 
  GetElFaces (int elnr, Array<int> & fnums, Array<int> & orient) const
  {
    fnums.SetSize (6);
    orient.SetSize (6);
    int nfa = 
      Ng_GetElement_Faces (elnr+1, &fnums[0], &orient[0]);
    fnums.SetSize (nfa);
    orient.SetSize (nfa);

    for (int i = 0; i < nfa; i++)
      fnums[i]--;
  }

  int MeshAccess :: 
  GetSElFace (int selnr) const
  {
    return GetElement<2,BND>(selnr).Faces()[0];
    // return Ng_GetSurfaceElement_Face (selnr+1, 0)-1;
  }
  
  void MeshAccess :: 
  GetSElFace (int selnr, int & fnum, int & orient) const
  {
    fnum = Ng_GetSurfaceElement_Face (selnr+1, &orient);
    fnum--;
  }

  void MeshAccess :: GetFacePNums (int fnr, Array<int> & pnums) const
  {
    pnums = ArrayObject (mesh.GetNode<2> (fnr).vertices);
  }

 
  void MeshAccess :: GetFaceEdges (int fnr, Array<int> & edges) const
  {
    edges.SetSize(4);
    int ned = Ng_GetFace_Edges (fnr+1, &edges[0]);
    edges.SetSize(ned);
    for (int i = 0; i < ned; i++) edges[i]--;
  }
 

  void MeshAccess :: GetFaceElements (int fnr, Array<int> & elnums) const
  {
    if (dim == 3)
      {
        elnums.SetSize0();
        auto vnums = ArrayObject(mesh.GetNode<2> (fnr).vertices);
        auto vels = ArrayObject(mesh.GetNode<0> (vnums[0]).elements);
        for (auto el : vels)
          {
            for (int f : ArrayObject(mesh.GetElement<3>(el).faces))
              if (f == fnr)
                elnums.Append (el);
          }
        return;
      }

    // 2D case: as it was before: one volume element, is it still needed ???
    ArrayMem<int, 9> vnums;
    GetFacePNums(fnr, vnums);

    ArrayMem<int, 50> vels;
    GetVertexElements (vnums[0], vels);

    int faces[8];
    elnums.SetSize (0);
    for (int i = 0; i < vels.Size(); i++)
      {
	int nfa = Ng_GetElement_Faces (vels[i]+1, faces, 0);
	for (int j = 0; j < nfa; j++)
	  if (faces[j]-1 == fnr)
	    elnums.Append (vels[i]);
      }
  }


  void MeshAccess :: GetFaceSurfaceElements (int fnr, Array<int> & elnums) const
  {
    /*
    ArrayMem<int, 9> vnums;
    GetFacePNums(fnr, vnums);

    ArrayMem<int, 50> vels;
    GetVertexSurfaceElements (vnums[0], vels);

    elnums.SetSize (0);
    for (int i = 0; i < vels.Size(); i++)
      {
	int sface = Ng_GetSurfaceElement_Face (vels[i]+1)-1;
        if (sface == fnr)
          elnums.Append (vels[i]);
      }
    */
    size_t v0 = GetFacePNums(fnr)[0];
    elnums.SetSize0();    
    for (auto sel : GetVertexSurfaceElements(v0))
      {
        int sface = Ng_GetSurfaceElement_Face (sel+1)-1;
        if (sface == fnr)
          elnums.Append (sel);
      }
  }


  /*
  void MeshAccess :: GetEdgePNums (int enr, int & pn1, int & pn2) const
  {
    auto edge = mesh.GetNode<1>(enr);
    pn1 = edge.vertices[0];
    pn2 = edge.vertices[1];

    // int v2[2];
    // Ng_GetEdge_Vertices (enr+1, v2);
    // pn1 = v2[0]-1;
    // pn2 = v2[1]-1;
  }
  */

  /*
  auto MeshAccess :: GetEdgePNums (int enr) const -> decltype(ArrayObject(INT<2>()))
  {
    int v2[2];
    Ng_GetEdge_Vertices (enr+1, v2);
    return ArrayObject (INT<2> (v2[0]-1, v2[1]-1));
  }
  */

  void MeshAccess :: GetEdgePNums (int enr, Array<int> & pnums) const
  {
    pnums.SetSize(2);
    // Ng_GetEdge_Vertices (enr+1, &pnums[0]);
    // pnums[0] -= 1;
    // pnums[1] -= 1;
    auto edge = mesh.GetNode<1>(enr);
    pnums[0] = edge.vertices[0];
    pnums[1] = edge.vertices[1];
  }

  void MeshAccess :: GetElFacets (ElementId ei, Array<int> & fnums) const
  {
    /*
    if (dim == 1)
      fnums = GetElement(ei).Vertices();
    else
    */
    fnums = GetElement(ei).Facets();
  }

  // some utility for Facets
  void MeshAccess :: GetElFacets (int elnr, Array<int> & fnums) const
  {
    switch (dim)
      {
      case 1: fnums = GetElement<1,VOL> (elnr).Vertices(); break;
      case 2: fnums = GetElement<2,VOL> (elnr).Edges(); break;
      default:
        fnums = GetElement<3,VOL> (elnr).Faces();
      }
  } 
    
  void MeshAccess :: GetSElFacets (int selnr, Array<int> & fnums) const
  {
    switch (dim)
      {
      case 1:
        GetSElVertices(selnr, fnums); break;
      case 2:
        GetSElEdges(selnr, fnums); break;
      default:
        {
          fnums.SetSize(1);
          fnums[0] = GetSElFace(selnr);
        }
      }
  }
  
  void MeshAccess :: GetFacetPNums (int fnr, Array<int> & pnums) const
  {
    switch (dim)
      {
      case 1: pnums.SetSize(1); pnums[0] = fnr; break;
      case 2: GetEdgePNums(fnr, pnums); break;
      case 3: GetFacePNums(fnr, pnums); break;
      }
  }

  ELEMENT_TYPE MeshAccess :: GetFacetType (int fnr) const
  {
    switch (dim)
      {
      case 1: return ET_POINT; 
      case 2: return ET_SEGM;
      default:  // i.e. dim = 3
        /*
	ArrayMem<int, 4> pnums;
	GetFacePNums(fnr, pnums);
	return (pnums.Size() == 3) ? ET_TRIG : ET_QUAD;
        */
        return (mesh.GetNode<2>(fnr).vertices.Size() == 3) ? ET_TRIG : ET_QUAD;
      }
  }

  void MeshAccess::CalcIdentifiedFacets()
  {
    identified_facets.SetSize(nnodes_cd[1]);
    for(auto i : Range(identified_facets.Size()))
      identified_facets[i] = std::tuple<int,int>(i,1);

    auto & idents = mesh.GetMesh()->GetIdentifications();
    int nid = idents.GetMaxNr();
    
    for(auto id : Range(1,nid+1))
      {
	// for periodic identification by now
	if (idents.GetType(id) != 2) continue;
	Array<INT<2>> pairs;
	switch(mesh.GetDimension())
	  {
	  case 1: GetPeriodicVertices(id,pairs); break;
	  case 2: GetPeriodicEdges(id,pairs); break;
	    // default: // here should be a 3D version
	  }
	for(auto pair : pairs)
	  for(auto l : Range(2))
	    identified_facets[pair[l]] = std::tuple<int,int>(pair[1-l],2);
      }
  }

  void MeshAccess::PushStatus (const char * str) const
  { 
    Ng_PushStatus (str);
  }
  void MeshAccess::PopStatus () const
  { 
    Ng_PopStatus ();
  }
  void MeshAccess::SetThreadPercentage (double percent) const
  { 
    Ng_SetThreadPercentage (percent); 
  }
  void MeshAccess :: GetStatus (string & str, double & percent) const
  {
    char * s;
    Ng_GetStatus(&s,percent);
    str = s;
  }

  void MeshAccess :: SetTerminate(void) const
  {
    Ng_SetTerminate();
  }
  void MeshAccess :: UnSetTerminate(void) const
  {
    Ng_UnSetTerminate();
  }
  bool MeshAccess :: ShouldTerminate(void) const
  {
    return (Ng_ShouldTerminate() != 0);
  }
  
  void MeshAccess :: GetVertexElements (int vnr, Array<int> & elnrs) const
  {
    elnrs = ArrayObject(mesh.GetNode<0> (vnr).elements);
    /*
    int nel = Ng_GetNVertexElements (vnr+1);
    elnrs.SetSize (nel);
    Ng_GetVertexElements (vnr+1, &elnrs[0]);
    for (int j = 0; j < nel; j++)
      elnrs[j]--;
    */
  }







  
  template <int DIM>
  ElementTransformation & MeshAccess :: 
  GetTrafoDim (size_t elnr, Allocator & lh) const
  {
    // static Timer t("MeshAccess::GetTrafoDim"); RegionTimer reg(t);
    
    ElementTransformation * eltrans;
    GridFunction * loc_deformation = deformation.get();
    
    Ngs_Element el (mesh.GetElement<DIM> (elnr), ElementId(VOL, elnr));
    
    if (pml_trafos[el.GetIndex()])
      {
        eltrans = new (lh)
          PML_ElementTransformation<DIM, DIM, Ng_ElementTransformation<DIM,DIM>>
          (this, el.GetType(), 
           ElementId(VOL,elnr), el.GetIndex(), *pml_trafos[el.GetIndex()]);
      }
    
    else if (loc_deformation)
      {
        if (el.is_curved)
          eltrans = new (lh)
            ALE_ElementTransformation<DIM, DIM, Ng_ElementTransformation<DIM,DIM>>
            (this, el.GetType(), 
             ElementId(VOL,elnr), el.GetIndex(),
             loc_deformation, 
             dynamic_cast<LocalHeap&> (lh));

        else

          eltrans = new (lh)
            ALE_ElementTransformation<DIM, DIM, Ng_ConstElementTransformation<DIM,DIM>>
            (this, el.GetType(), 
             ElementId(VOL,elnr), el.GetIndex(),
             loc_deformation, 
             dynamic_cast<LocalHeap&> (lh));
      }

    else if ( el.is_curved )

      eltrans = new (lh) Ng_ElementTransformation<DIM,DIM> (this, el.GetType(), 
                                                            ElementId(VOL,elnr), el.GetIndex()); 

    else
      eltrans = new (lh) Ng_ConstElementTransformation<DIM,DIM> (this, el.GetType(), 
                                                                 ElementId(VOL,elnr), el.GetIndex()); 

    /*
    eltrans->SetElementType (el.GetType());
    int elind = el.GetIndex();
    eltrans->SetElement (0, elnr, elind);
    */

    if(higher_integration_order.Size() == GetNE() && higher_integration_order[elnr])
      eltrans->SetHigherIntegrationOrder();
    else
      eltrans->UnSetHigherIntegrationOrder();

    return *eltrans;
  }


  template <int DIM>
  ElementTransformation & MeshAccess :: 
  GetSTrafoDim (size_t elnr, Allocator & lh) const
  {
    // static Timer t("MeshAccess::GetTrafoDim");

    ElementTransformation * eltrans;
    
    Ngs_Element el(mesh.GetElement<DIM-1> (elnr), ElementId(BND, elnr));
    GridFunction * loc_deformation = deformation.get();
    
    if (loc_deformation)

      eltrans = new (lh) ALE_ElementTransformation<DIM-1,DIM, Ng_ElementTransformation<DIM-1,DIM>>
        (this, el.GetType(), 
         ElementId(BND,elnr), el.GetIndex(),
         loc_deformation, 
         dynamic_cast<LocalHeap&> (lh)); 
    
    else if ( el.is_curved )

      eltrans = new (lh) Ng_ElementTransformation<DIM-1,DIM> (this, el.GetType(), 
                                                              ElementId(BND,elnr), el.GetIndex()); 
    
    else
      eltrans = new (lh) Ng_ConstElementTransformation<DIM-1,DIM> (this, el.GetType(), 
                                                                   ElementId(BND,elnr), el.GetIndex()); 

    /*
    eltrans->SetElementType (el.GetType());
    int elind = el.GetIndex();
    eltrans->SetElement (0, elnr, elind);
    */

    if(higher_integration_order.Size() == GetNSE() && higher_integration_order[elnr])
      eltrans->SetHigherIntegrationOrder();
    else
      eltrans->UnSetHigherIntegrationOrder();

    return *eltrans;
  }

  template <int DIM>
  ElementTransformation & MeshAccess ::
  GetCD2TrafoDim (size_t elnr, Allocator & lh) const
  {
    ElementTransformation * eltrans;
    Ngs_Element el(mesh.GetElement<DIM-2>(elnr), ElementId(BBND,elnr));
    GridFunction * loc_deformation = deformation.get();
    if(loc_deformation)
      eltrans = new (lh) ALE_ElementTransformation<DIM-2,DIM, Ng_ElementTransformation<DIM-2,DIM>>
	(this,el.GetType(),
	 ElementId(BBND,elnr), el.GetIndex(),
	 loc_deformation,
	 dynamic_cast<LocalHeap&>(lh));
    else if( el.is_curved )
      eltrans = new (lh) Ng_ElementTransformation<DIM-2,DIM> (this, el.GetType(),
							      ElementId(BBND,elnr), el.GetIndex());

    else
      eltrans = new (lh) Ng_ConstElementTransformation<DIM-2,DIM> (this, el.GetType(),
								   ElementId(BBND,elnr), el.GetIndex());

    if(higher_integration_order.Size() == GetNCD2E() && higher_integration_order[elnr])
      eltrans->SetHigherIntegrationOrder();
    else
      eltrans->UnSetHigherIntegrationOrder();

    return *eltrans;
  }

  template <> ElementTransformation & MeshAccess :: GetTrafo (T_ElementId<VOL,1> ei, Allocator & lh) const;
  template <> ElementTransformation & MeshAccess :: GetTrafo (T_ElementId<VOL,2> ei, Allocator & lh) const;
  template <> ElementTransformation & MeshAccess :: GetTrafo (T_ElementId<VOL,3> ei, Allocator & lh) const;
  template <> ElementTransformation & MeshAccess :: GetTrafo (T_ElementId<BND,1> ei, Allocator & lh) const;
  template <> ElementTransformation & MeshAccess :: GetTrafo (T_ElementId<BND,2> ei, Allocator & lh) const;
  template <> ElementTransformation & MeshAccess :: GetTrafo (T_ElementId<BND,3> ei, Allocator & lh) const;
  template <> ElementTransformation & MeshAccess :: GetTrafo (T_ElementId<BBND,2> ei, Allocator & lh) const;
  template <> ElementTransformation & MeshAccess :: GetTrafo (T_ElementId<BBND,3> ei, Allocator & lh) const;

  ngfem::ElementTransformation & MeshAccess :: GetTrafo (ElementId ei, Allocator & lh) const
  {
    int elnr = ei.Nr();
    VorB vb = ei.VB();
    
    switch(vb)
      {
      case VOL:
        switch (dim)
          {
          case 1: return GetTrafoDim<1> (elnr, lh);
          case 2: return GetTrafoDim<2> (elnr, lh);
          case 3: return GetTrafoDim<3> (elnr, lh);

          default:
            throw Exception ("MeshAccess::GetTrafo, illegal dimension");
          }

      case BND:
        switch (dim)
          {
          case 1: return GetSTrafoDim<1> (elnr, lh);
          case 2: return GetSTrafoDim<2> (elnr, lh);
          case 3: return GetSTrafoDim<3> (elnr, lh);

          default:
            throw Exception ("MeshAccess::GetSTrafo, illegal dimension");
          }

      case BBND:
	switch(dim)
	  {
	  case 2: return GetCD2TrafoDim<2>(elnr,lh);
	  case 3: return GetCD2TrafoDim<3>(elnr,lh);

	  default:
	    throw Exception ("MeshAccess::GetCD2Trafo, illegal dimension");
	  }
      }
  }


  Region :: Region (shared_ptr<MeshAccess> amesh,
                    VorB avb, string pattern)
    : mesh(amesh), vb(avb)
  {
    mask = BitArray(mesh->GetNRegions(vb));
    mask.Clear();
    regex re_pattern(pattern);

    for (int i : Range(mask))
      if (regex_match(mesh->GetMaterial(vb,i), re_pattern))
        mask.Set(i);
    
    /*
    if (vb == VOL)
      {
        mask = BitArray(mesh->GetNDomains());
        mask.Clear();
        regex re_pattern(pattern);
        for (int i : Range(mask))
          if (regex_match(mesh->GetMaterial(VOL,i), re_pattern))
            mask.Set(i);
      }
    else
      if (vb==BND)
      {
        mask = BitArray(mesh->GetNBoundaries());
        mask.Clear();
        regex re_pattern(pattern);
        for (int i : Range(mask))
          if (regex_match(mesh->GetMaterial(BND,i), re_pattern))
            mask.Set(i);
      }
      else
	{
          mask = BitArray(mesh->GetNBBoundaries());
          mask.Clear();
          regex re_pattern(pattern);
	  for(int i : Range(mask))
	    (*testout) << "boundary condition " << i << ": " << mesh->GetMaterial(BBND,i) << endl;
	  for (int i : Range(mask))
	    if (regex_match(mesh->GetMaterial(BBND,i), re_pattern))
	      mask.Set(i);
	  (*testout) << "mask: " << mask << endl;
	}
    */
  }      


  

  double MeshAccess :: ElementVolume (int elnr) const
  {
    static FE_Segm0 segm0;
    static ScalarFE<ET_TRIG,0> trig0;
    static ScalarFE<ET_QUAD,0> quad0;
    static ScalarFE<ET_TET,0> tet0;
    static FE_Prism0 prism0;
    static FE_Pyramid0 pyramid0;
    FE_Hex0 hex0;
    
    const FiniteElement * fe = NULL;
    switch (GetElType (ElementId(VOL, elnr)))
      {
      case ET_SEGM: fe = &segm0; break;
      case ET_TRIG: fe = &trig0; break;
      case ET_QUAD: fe = &quad0; break;
      case ET_TET: fe = &tet0; break;
      case ET_PYRAMID: fe = &pyramid0; break;
      case ET_PRISM: fe = &prism0; break;
	// case ET_HEX: fe = &hex0; break;
      default:
	{
	  cerr << "ElementVolume not implemented for el " << GetElType(ElementId(VOL, elnr)) << endl;
	}
      }
  
    LocalHeapMem<10000> lh("MeshAccess - elementvolume");

    ElementTransformation & trans = GetTrafo (ElementId(VOL, elnr), lh);
    ConstantCoefficientFunction ccf(1);

    if (GetDimension() == 1)
      {
	SourceIntegrator<1> si (&ccf);
	FlatVector<> elvec(fe->GetNDof(), lh);
	si.CalcElementVector (*fe, trans, elvec, lh);
	return elvec(0);
      }
    else if (GetDimension() == 2)
      {
	SourceIntegrator<2> si (&ccf);
	FlatVector<> elvec(fe->GetNDof(), lh);
	si.CalcElementVector (*fe, trans, elvec, lh);
	return elvec(0);
      }
    else
      {
	SourceIntegrator<3> si(&ccf);
	FlatVector<> elvec(fe->GetNDof(), lh);
	si.CalcElementVector (*fe, trans, elvec, lh);
	return elvec(0);
      }
  }


  
  double MeshAccess :: SurfaceElementVolume (int selnr) const
  {
    static ScalarFE<ET_TRIG,0> trig0;
    static ScalarFE<ET_QUAD,0> quad0;
    ElementId sei(BND, selnr);
    const FiniteElement * fe;
    switch (GetElType (sei))
      {
      case ET_TRIG: fe = &trig0; break;
      case ET_QUAD: fe = &quad0; break;
      default:
	{
	  cerr << "SurfaceElementVolume not implemented for el " << GetElType(sei) << endl;
	  return 0;
	}
      }

    LocalHeapMem<10000> lh("MeshAccess - surfaceelementvolume");

    ElementTransformation & trans = GetTrafo (sei, lh);
    ConstantCoefficientFunction ccf(1);

    if (GetDimension() == 2)
      {
	NeumannIntegrator<2> si( &ccf );
	FlatVector<> elvec (fe->GetNDof(), lh);
	si.CalcElementVector (*fe, trans, elvec, lh);
	return elvec(0);
      }
    else
      {
	NeumannIntegrator<3> si( &ccf );
	FlatVector<> elvec (fe->GetNDof(), lh);
	si.CalcElementVector (*fe, trans, elvec, lh);
	return elvec(0);
      }
  }


  void MeshAccess :: SetPointSearchStartElement(const int el) const
  {
    Ng_SetPointSearchStartElement(el+1);
  }

  int MeshAccess :: FindElementOfPoint (FlatVector<double> point,
					IntegrationPoint & ip, 
					bool build_searchtree,
					int index) const
  {
    ArrayMem<int,1> dummy(1);
    dummy[0] = index;
    return FindElementOfPoint(point,ip,build_searchtree,&dummy);
  }

  int MeshAccess :: FindElementOfPoint (FlatVector<double> point,
					IntegrationPoint & ip,
					bool build_searchtree,
					const Array<int> * const indices) const
  {
    static Timer t("FindElementOfPonit");
    RegionTimer reg(t);


    if (indices != NULL)
      {
        switch (dim)
          {
          case 1:
            return mesh.FindElementOfPoint<1> (&point(0), &ip(0), build_searchtree, 
                                               &(*indices)[0],indices->Size());
          case 2:
            return mesh.FindElementOfPoint<2> (&point(0), &ip(0), build_searchtree, 
                                               &(*indices)[0],indices->Size());
          case 3:
            return mesh.FindElementOfPoint<3> (&point(0), &ip(0), build_searchtree,
                                               &(*indices)[0],indices->Size());
          }
      }
    else
      {  
        switch (dim)
          {
          case 1: return mesh.FindElementOfPoint<1> (&point(0), &ip(0), build_searchtree, NULL, 0);
          case 2: return mesh.FindElementOfPoint<2> (&point(0), &ip(0), build_searchtree, NULL, 0);
          case 3: return mesh.FindElementOfPoint<3> (&point(0), &ip(0), build_searchtree, NULL, 0);
          }
      }

    return -1;
  }


  int MeshAccess :: FindSurfaceElementOfPoint (FlatVector<double> point,
					       IntegrationPoint & ip, 
					       bool build_searchtree,
					       int index) const
  {
    ArrayMem<int,1> dummy(1);
    dummy[0] = index;
    return FindSurfaceElementOfPoint(point,ip,build_searchtree,&dummy);
  }

  int MeshAccess :: FindSurfaceElementOfPoint (FlatVector<double> point,
					       IntegrationPoint & ip,
					       bool build_searchtree,
					       const Array<int> * const indices) const
  {
    static Timer t("FindSurfaceElementOfPonit");
    RegionTimer reg(t);
    int elnr;

    if(indices != NULL && indices->Size()>0)
      {
	// int * dummy = new int[indices->Size()];
	// for(int i=0; i<indices->Size(); i++) dummy[i] = (*indices)[i]+1;
	// elnr = Ng_FindSurfaceElementOfPoint (&point(0), lami, build_searchtree,dummy,indices->Size());
	// delete [] dummy;
        elnr = mesh.FindElementOfPoint<2> (&point(0), &ip(0), build_searchtree,
                                           &(*indices)[0],indices->Size());
      }
    else
      {  
	// elnr = Ng_FindSurfaceElementOfPoint (&point(0), lami, build_searchtree);
        elnr = mesh.FindElementOfPoint<2> (&point(0), &ip(0), build_searchtree, NULL, 0);        
      }

    return elnr;
  }


  void NGSolveTaskManager (function<void(int,int)> func)
  {
    // cout << "call ngsolve taskmanager from netgen, tm = " << task_manager << endl;
    if (!task_manager)
      func(0,1);
    else
      task_manager->CreateJob
        ([&](TaskInfo & info)
         {
           func(info.task_nr, info.ntasks);
         }, TasksPerThread(4));
  }

  
  void MeshAccess :: Refine ()
  {
    static Timer t("MeshAccess::Refine"); RegionTimer reg(t);
    nlevels = std::numeric_limits<int>::max();
    mesh.Refine(NG_REFINE_H, &NGSolveTaskManager);
    UpdateBuffers();
  }

  int MeshAccess :: GetNPairsPeriodicVertices () const 
  {
    return Ng_GetNPeriodicVertices(0);
  }

  int MeshAccess :: GetNPairsPeriodicVertices (int idnr) const 
  {
    return Ng_GetNPeriodicVertices(idnr);
  }
 
  void MeshAccess :: GetPeriodicVertices ( Array<INT<2> > & pairs) const
  {
    int npairs;

    npairs = Ng_GetNPeriodicVertices (0);
    pairs.SetSize (npairs);

    Ng_GetPeriodicVertices (0,&pairs[0][0]);
    for (int i = 0; i < pairs.Size(); i++)
      {
	pairs[i][0]--;
	pairs[i][1]--;
      }
  }
 
  void MeshAccess :: GetPeriodicVertices (int idnr, Array<INT<2> > & pairs) const
  {
    int npairs;

    npairs = Ng_GetNPeriodicVertices (idnr);

    pairs.SetSize (npairs);

    Ng_GetPeriodicVertices (idnr,&pairs[0][0]);

    for (int i = 0; i < pairs.Size(); i++)
      {
	pairs[i][0]--;
	pairs[i][1]--;
      }
  }


  const Array<INT<2>> & MeshAccess :: GetPeriodicNodes (NODE_TYPE nt, int idnr) const
  {
    return (*periodic_node_pairs[nt])[idnr-1];
  }


  int MeshAccess :: GetNPairsPeriodicEdges () const 
  {
    return Ng_GetNPeriodicEdges(0);
  }
 
  void MeshAccess :: GetPeriodicEdges ( Array<INT<2> > & pairs) const
  {
    int npairs;

    npairs = Ng_GetNPeriodicEdges (0);
    pairs.SetSize (npairs);

    Ng_GetPeriodicEdges (0,&pairs[0][0]);
    for (int i = 0; i < pairs.Size(); i++)
      {
	pairs[i][0]--;
	pairs[i][1]--;
      }
  }

  int MeshAccess :: GetNPairsPeriodicEdges (int idnr) const 
  {
    return Ng_GetNPeriodicEdges(idnr);
  }
 
  void MeshAccess :: GetPeriodicEdges (int idnr, Array<INT<2> > & pairs) const
  {
    int npairs;

    npairs = Ng_GetNPeriodicEdges (idnr);
    pairs.SetSize (npairs);

    Ng_GetPeriodicEdges (idnr,&pairs[0][0]);
    for (int i = 0; i < pairs.Size(); i++)
      {
	pairs[i][0]--;
	pairs[i][1]--;
      }
  }


///// Added by Roman Stainko ....
void MeshAccess::GetVertexSurfaceElements( int vnr, Array<int>& elems) const
{
    elems = GetVertexSurfaceElements(vnr);
}



  void MeshAccess::SetHigherIntegrationOrder(int elnr)
  {
    if(higher_integration_order.Size() != GetNE())
      {
	higher_integration_order.SetSize(GetNE());
	higher_integration_order = false;
      }
    higher_integration_order[elnr] = true;
  }
  void MeshAccess::UnSetHigherIntegrationOrder(int elnr)
  {
    if(higher_integration_order.Size() != GetNE())
      {
	higher_integration_order.SetSize(GetNE());
	higher_integration_order = false;
      }
    higher_integration_order[elnr] = false;
  }


  void MeshAccess :: InitPointCurve(double red, double green, double blue) const
  {
    Ng_InitPointCurve(red, green, blue);
  }

  void MeshAccess :: AddPointCurvePoint(const Vec<3> & point) const
  {
    Ng_AddPointCurvePoint(&(point(0)));
  }

 

#ifdef PARALLEL
 
  size_t MeshAccess ::GetGlobalNodeNum (Node node) const
  {
    int glob = NgPar_GetGlobalNodeNum (node.GetType(), node.GetNr());
    return glob;
  }
  
  void MeshAccess :: GetDistantProcs (Node node, Array<int> & procs) const
  {
    procs.SetSize( NgPar_GetNDistantNodeNums(node.GetType(), node.GetNr()) );
    NgPar_GetDistantNodeNums ( node.GetType(), node.GetNr(), &procs[0] );
  }


#else

  size_t MeshAccess ::GetGlobalNodeNum (NodeId node) const  
  {
    return -1;
  }

  void MeshAccess :: GetDistantProcs (NodeId node, Array<int> & procs) const
  {
    procs.SetSize (0);
  }
#endif


  function<void()> cleanup_func;
  ProgressOutput :: ProgressOutput (shared_ptr<MeshAccess> ama,
				    string atask, size_t atotal)
    : ma(ama), task(atask), total(atotal)
  {
    is_root = (MyMPI_GetId() == 0);
    prevtime = WallTime();
    size_t glob_total = MyMPI_Reduce (total);
    if (is_root) total = glob_total;

    done_called = false;
    cnt = 0;
    thd_cnt = 0;
    cleanup_func = [this] () {  this->SumUpLocal(); };
    TaskManager::SetCleanupFunction(cleanup_func); 
  }

  ProgressOutput :: ~ProgressOutput ()
  {
    Done();
    TaskManager::SetCleanupFunction();
  }  

  atomic<size_t> ProgressOutput :: cnt;
  thread_local size_t ProgressOutput :: thd_cnt = 0;
  thread_local double ProgressOutput :: thd_prev_time = WallTime();
  
  void ProgressOutput :: Update ()
  {
    thd_cnt++;
    double time = WallTime();
    if (time > thd_prev_time+0.05)
      {
        thd_prev_time = time;
        cnt += thd_cnt;
        thd_cnt = 0;
        Update(cnt);
      }
  }
  
  void ProgressOutput :: SumUpLocal()
  {
    cnt += thd_cnt;
    thd_cnt = 0;
    // Update(cnt);
  }

  void ProgressOutput :: Update (size_t nr)
  {
    static mutex progressupdate_mutex;
    double time = WallTime();
    if (time > prevtime+0.05)
      {
	{
          lock_guard<mutex> guard(progressupdate_mutex);
	  if (is_root)
	    {
	      cout << IM(3) << "\r" << task << " " << nr << "/" << total << flush;
	      ma->SetThreadPercentage ( 100.0*nr / total);
	    }
#ifdef PARALLEL
	  else
	    {
	      static Timer t("dummy - progressreport"); RegionTimer r(t);
	      MPI_Send (&nr, 1, MPI_INT, 0, MPI_TAG_SOLVE, ngs_comm);
              // changed from BSend (VSC-problem)
	    }
#endif
	  prevtime = WallTime();
	}
      }
  }
  


  void ProgressOutput :: Done()
  {
    if (done_called) return;
    done_called = true;

    if (is_root)
      {
#ifdef PARALLEL	  
	int ntasks = MyMPI_GetNTasks();
	if (ntasks > 1)
	  {
	    Array<int> working(ntasks), computed(ntasks);
	    working = 1;
	    computed = 0;
	    while (1)
	      {
		int flag, data, num_working = 0, got_flag = false;
		for (int source = 1; source < ntasks; source++)
		  {
		    if (!working[source]) continue;
		    num_working++;
		    MPI_Iprobe (source, MPI_TAG_SOLVE, ngs_comm, &flag, MPI_STATUS_IGNORE);
		    if (flag)
		      {
			got_flag = true;
			MPI_Recv (&data, 1, MPI_INT, source, MPI_TAG_SOLVE, ngs_comm, MPI_STATUS_IGNORE);
			if (data == -1) 
			  working[source] = 0;
			else
			  computed[source] = data;
		      }
		  }
		int sum = 0;
		for (int j = 1; j < ntasks; j++) 
		  sum += computed[j];
		cout << IM(3) 
		     << "\r" << task << " " << sum << "/" << total
		     << " (" << num_working << " procs working) " << flush;
		ma->SetThreadPercentage ( 100.0*sum / total );
		if (!num_working) break;
		if (!got_flag) usleep (1000);
	      }
	  }
#endif
	cout << IM(3) << "\r" << task << " " << total << "/" << total
	     << "                                 " << endl;
      }
    else
      {
#ifdef PARALLEL
	MPI_Send (&total, 1, MPI_INT, 0, MPI_TAG_SOLVE, ngs_comm);
	int final = -1;
	MPI_Send (&final, 1, MPI_INT, 0, MPI_TAG_SOLVE, ngs_comm);
#endif
      }
  }
  


}




