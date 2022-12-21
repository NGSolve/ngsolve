/*********************************************************************/
/* File:   meshaccess.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   Access to fe mesh
*/
#include <core/register_archive.hpp>
#ifdef NGS_PYTHON
#include <core/python_ngcore.hpp>
#endif // NGS_PYTHON

#include <ngstd.hpp>
#include <nginterface.h>

#include <comp.hpp>
#include "../fem/h1lofe.hpp"
#include <regex>

namespace ngcomp
{


  template <int DIMS, int DIMR, typename BASE> class ALE_ElementTransformation;
  
  
  string Ngs_Element::defaultstring = "default";
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

    virtual int SpaceDim () const override
    {
      return DIMR;
    }

    virtual VorB VB() const override
    {
      return VorB(int(DIMR)-int(DIMS));
    }

    virtual bool BelongsToMesh (const void * mesh2) const override
    {
      // return mesh == &(static_cast<const MeshAccess*> (mesh2) -> mesh);
      return mesh == mesh2;
    }

    virtual const void * GetMesh () const override { return mesh; }

    virtual void GetSort (FlatArray<int> sort) const override
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
			       FlatMatrix<> dxdxi) const override
    {
      mesh->mesh.ElementTransformation <DIMS,DIMR> (elnr, &ip(0), NULL, &dxdxi(0));
    }
    
    virtual void CalcPoint (const IntegrationPoint & ip,
			    FlatVector<> point) const override
    {
      mesh->mesh.ElementTransformation <DIMS,DIMR> (elnr, &ip(0), &point(0), NULL);
    }

    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, FlatMatrix<> dxdxi) const override
    {
      mesh->mesh.ElementTransformation <DIMS,DIMR> (elnr, &ip(0), &point(0), &dxdxi(0));
    }

    virtual BaseMappedIntegrationPoint & operator() (const IntegrationPoint & ip, Allocator & lh) const override
    {
      return *new (lh) MappedIntegrationPoint<DIMS,DIMR> (ip, *this);
    }

    virtual void VCalcHesse (const SIMD<ngfem::IntegrationPoint> & ip, SIMD<double> * hesse) const override
    {
      double eps = 1e-4;
      Mat<4*DIMS,DIMS,SIMD<double>> pnts;
      for (int i = 0; i < DIMS; i++)
        {
          for (int j = 0; j < DIMS; j++)
	    pnts(2*i,j) = pnts(2*i+1, j) = pnts(2*(i+DIMS), j) = pnts(2*(i+DIMS)+1, j)= ip(j); 
          pnts(2*i,i) += eps;
          pnts(2*i+1,i) -= eps;
	  pnts(2*(i+DIMS),i) += 2*eps;
	  pnts(2*(i+DIMS)+1,i) -= 2*eps;
        }

      Mat<4*DIMS,DIMR,SIMD<double>> x;
      Mat<4*DIMS,DIMR*DIMS,SIMD<double>> dx;

      mesh->mesh.MultiElementTransformation <DIMS,DIMR>
        (elnr, 4*DIMS,
         &pnts(0,0), &pnts(1,0)-&pnts(0,0),
         &x(0,0), &x(1,0)-&x(0,0),
         &dx(0,0), &dx(1,0)-&dx(0,0));
      
      for (int i = 0; i < DIMR; i++)
        for (int j = 0; j < DIMS; j++)
          for (int k = 0; k < DIMS; k++)
	    hesse[i*DIMS*DIMS + j*DIMS + k] = (8.0*dx(2*j, i*DIMS+k) - 8.0*dx(2*j+1, i*DIMS+k) - dx(2*(j+DIMS), i*DIMS+k) + dx(2*(j+DIMS)+1, i*DIMS+k)) / (12.0*eps);
      //hesse[i*DIMS*DIMS + j*DIMS + k] = (dx(2*j, i*DIMS+k) - dx(2*j+1, i*DIMS+k) ) / (2*eps);
    }

    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, Allocator & lh) const override
    {
      return *new (lh) MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }

    virtual SIMD_BaseMappedIntegrationRule & operator() (const SIMD_IntegrationRule & ir, Allocator & lh) const override
    {
      return *new (lh) SIMD_MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }

    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & bmir) const override
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

      for (auto & mip : mir)
        mip.Compute();
  }


    virtual void CalcMultiPointJacobian (const SIMD_IntegrationRule & ir,
					 SIMD_BaseMappedIntegrationRule & bmir) const override
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
         &ir[0](0), ir.Size()>1 ? &ir[1](0)-&ir[0](0) : 0,
         &mir[0].Point()(0), ir.Size()>1 ? &mir[1].Point()(0)-&mir[0].Point()(0) : 0,
         &mir[0].Jacobian()(0,0), ir.Size()>1 ? &mir[1].Jacobian()(0,0)-&mir[0].Jacobian()(0,0) : 0);
      
      for (int i = 0; i < ir.Size(); i++)
        mir[i].Compute();
    }

    virtual const ElementTransformation & VAddDeformation (const GridFunction * gf, LocalHeap & lh) const override
    {
      return * new (lh) ALE_ElementTransformation<DIMS,DIMR,Ng_ElementTransformation<DIMS,DIMR>>
        (gf->GetMeshAccess().get(), eltype, GetElementId(), elindex, gf, lh);
    }
    
  };
  




  template <int DIMS, int DIMR, typename BASE>
  class ALE_ElementTransformation : public BASE
  {
    const GridFunction * deform;
    const ScalarFiniteElement<DIMS> * fel;
    // FlatVector<> elvec;
    FlatMatrix<> elvecs;
  public:
    ALE_ElementTransformation (const MeshAccess * amesh, 
                               ELEMENT_TYPE aet, ElementId ei, int elindex,
                               const GridFunction * adeform,
                               Allocator & lh)
      : BASE(amesh, aet, ei, elindex), 
        deform(adeform) 
    {
      // LocalHeap & lh = dynamic_cast<LocalHeap&> (alh);
      this->iscurved = true;

      auto & bfel = deform->GetFESpace()->GetFE(ei, lh);
      auto cfel = dynamic_cast<VectorFiniteElement*> (&bfel);
      if (cfel) // VectorH1
        {
          fel = dynamic_cast<const ScalarFiniteElement<DIMS>*> (&(*cfel)[0]);
          size_t nds = fel->GetNDof();
          ArrayMem<int,100> dnums(cfel->GetNDof());
          deform->GetFESpace()->GetDofNrs(ei, dnums);
          
          // elvec.AssignMemory(dnums.Size(), lh);

          // FlatVector<> helvec(dnums.Size(), lh); // temporary
          VectorMem<100> helvec(dnums.Size()); // temporary
          deform->GetElementVector(dnums, helvec);
          /*
          for (int i = 0; i < DIMR; i++)
            for (int j = 0; j < nds; j++)
              elvec(j*DIMR+i) = helvec(i*nds+j);
          elvecs.AssignMemory(DIMR, nds, lh);
          for (int j = 0; j < DIMR; j++)
            elvecs.Row(j) = elvec.Slice(j,DIMR);
          */
          // elvecs.AssignMemory(DIMR, nds, lh);
          elvecs.AssignMemory(DIMR, nds, new (lh) double[DIMR*nds]);
          for (int j = 0; j < DIMR; j++)
            elvecs.Row(j) = helvec.Range(j*nds, (j+1)*nds);
        }
      else
        {
          fel = dynamic_cast<const ScalarFiniteElement<DIMS>*> (&bfel);
          
          ArrayMem<int,100> dnums(fel->GetNDof());
          deform->GetFESpace()->GetDofNrs(ei, dnums);

          // elvec.AssignMemory(DIMR*dnums.Size(), lh);
          VectorMem<100> elvec(DIMR*dnums.Size());
          deform->GetElementVector(dnums, elvec);
          
          // elvecs.AssignMemory(DIMR, dnums.Size(), lh);
          elvecs.AssignMemory(DIMR, dnums.Size(), new (lh) double[DIMR*dnums.Size()]);          
          for (int j = 0; j < DIMR; j++)
            elvecs.Row(j) = elvec.Slice(j,DIMR);
        }
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
			       FlatMatrix<> dxdxi) const override
    {
      Mat<DIMR,DIMS> tmp;
      BASE::CalcJacobian (ip, tmp);

      Mat<DIMR,DIMS> def;
      for (int i = 0; i < DIMR; i++)
        def.Row(i) = fel->EvaluateGrad (ip, elvecs.Row(i));
      dxdxi = def + tmp;
    }
    
    virtual void CalcPoint (const IntegrationPoint & ip,
			    FlatVector<> point) const override
    {
      Vec<DIMR> tmp;
      BASE::CalcPoint (ip, tmp);

      Vec<DIMR> def;
      for (int i = 0; i < DIMR; i++)
        def(i) = fel->Evaluate (ip, elvecs.Row(i));
      point = tmp + def;
    }

    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, FlatMatrix<> dxdxi) const override
    {
      CalcJacobian (ip, dxdxi);
      CalcPoint (ip, point);
      // this->mesh->mesh.ElementTransformation <DIMS,DIMR> (elnr, &ip(0), &point(0), &dxdxi(0));
    }

    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & bmir) const override
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
					 SIMD_BaseMappedIntegrationRule & bmir) const override
    {
      SIMD_MappedIntegrationRule<DIMS,DIMR> & mir = 
	static_cast<SIMD_MappedIntegrationRule<DIMS,DIMR> &> (bmir);
      
      BASE::CalcMultiPointJacobian (ir, bmir);

      STACK_ARRAY(SIMD<double>, mem0, ir.Size());
      FlatVector<SIMD<double>> def(ir.Size(), &mem0[0]);
      STACK_ARRAY(SIMD<double>, mem1, (DIMS*ir.Size()));
      // FlatMatrix<SIMD<double>> grad(DIMS, ir.Size(), &mem1[0]);
      FlatMatrix<SIMD<double>> grad(DIMS, ir.Size(), mem1);

      for (int i = 0; i < DIMR; i++)
        {
          // fel->Evaluate (ir, elvec.Slice(i,DIMR), def);
          // fel->EvaluateGrad (ir, elvec.Slice(i,DIMR), grad);
          fel->Evaluate (ir, elvecs.Row(i), def);
          fel->EvaluateGrad (ir, elvecs.Row(i), grad);
          
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

    virtual int SpaceDim () const override
    {
      return DIMR;
    }

    virtual VorB VB() const override
    {
      return VorB(int(DIMR)-int(DIMS));
    }

    virtual void VCalcHesse (const SIMD<ngfem::IntegrationPoint> & ip, SIMD<double> * hesse) const override
    {
      for (size_t i = 0; i < DIMS*DIMS*DIMR; i++)
        hesse[i] = SIMD<double>(0);
    }

    virtual bool BelongsToMesh (const void * mesh2) const override 
    {
      // return mesh == &(static_cast<const MeshAccess*> (mesh2) -> mesh);
      return mesh == mesh2;
    }

    virtual const void * GetMesh () const override { return mesh; }

    virtual void GetSort (FlatArray<int> sort) const override
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
			       FlatMatrix<> dxdxi) const override
    {
      NETGEN_CHECK_RANGE(dxdxi.Width(), DIMS, DIMS+1);
      NETGEN_CHECK_RANGE(DIMR, 0, dxdxi.Height()+1);
      // dxdxi = mat;
      FlatMatrixFixWidth<DIMS> (DIMR, &dxdxi(0,0)) = mat;
    }
    
    virtual void CalcPoint (const IntegrationPoint & ip,
			    FlatVector<> point) const override
    {
      NETGEN_CHECK_RANGE(DIMR, 0, point.Size()+1);
      // point = p0 + mat * FlatVec<DIMS, const double> (&ip(0));
      FlatVec<DIMR> (&point(0)) = p0 + mat * FlatVec<DIMS, const double> (&ip(0));
    }

    virtual void CalcPointJacobian (const IntegrationPoint & ip,
				    FlatVector<> point, FlatMatrix<> dxdxi) const override
    {
      NETGEN_CHECK_RANGE(dxdxi.Width(), DIMS, DIMS+1);
      NETGEN_CHECK_RANGE(DIMR, 0, dxdxi.Height()+1);
      NETGEN_CHECK_RANGE(DIMR, 0, point.Size()+1);

      FlatVec<DIMR> (&point(0)) = p0 + mat * FlatVec<DIMS, const double> (&ip(0));
      FlatMatrixFixWidth<DIMS> (DIMR, &dxdxi(0,0)) = mat;
    }

    virtual BaseMappedIntegrationPoint & operator() (const IntegrationPoint & ip, Allocator & lh) const override
    {
      return *new (lh) MappedIntegrationPoint<DIMS,DIMR> (ip, *this);
    }

    virtual BaseMappedIntegrationRule & operator() (const IntegrationRule & ir, Allocator & lh) const override
    {
      return *new (lh) MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }    

    virtual SIMD_BaseMappedIntegrationRule & operator() (const SIMD_IntegrationRule & ir, Allocator & lh) const override
    {
      return *new (lh) SIMD_MappedIntegrationRule<DIMS,DIMR> (ir, *this, lh);
    }

    virtual void CalcMultiPointJacobian (const IntegrationRule & ir,
					 BaseMappedIntegrationRule & bmir) const override
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
					 SIMD_BaseMappedIntegrationRule & bmir) const override
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
          hmir[i].Compute();
        }
    }
    virtual const ElementTransformation & VAddDeformation (const GridFunction * gf, LocalHeap & lh) const override
    {
      return * new (lh) ALE_ElementTransformation<DIMS,DIMR,Ng_ConstElementTransformation<DIMS,DIMR>>
        (gf->GetMeshAccess().get(), eltype, GetElementId(), elindex, gf, lh);
    }

  };
  











  MeshAccess :: MeshAccess ()
    : mesh(shared_ptr<netgen::Mesh>())
  {
    ;
  }


  MeshAccess :: MeshAccess (shared_ptr<netgen::Mesh> amesh)
    : mesh(amesh)
  {
    // the connection to netgen global variables
    // ngstd::printmessage_importance = netgen::printmessage_importance;

    mesh.SelectMesh();

    mesh.UpdateTopology();  
    UpdateBuffers();
  }

  MeshAccess :: MeshAccess (string filename, NgMPI_Comm amesh_comm)
    : mesh(filename, amesh_comm)
  {
    // the connection to netgen global variables
    // ngstd::printmessage_importance = netgen::printmessage_importance;

    UpdateBuffers();
  }
  
  MeshAccess :: ~MeshAccess ()
  {
    ;
  }

  /*
  void MeshAccess :: LoadMesh (const string & filename)
  {
    static Timer t("MeshAccess::LoadMesh"); RegionTimer reg(t);
    mesh.LoadMesh (filename, this->mesh_comm);
    UpdateBuffers();
    if (!mesh.Valid())
      throw Exception ("could not load mesh from '" + filename + "'");
  }

  void MeshAccess :: LoadMesh (istream & str)
  {
    static Timer t("MeshAccess::LoadMesh"); RegionTimer reg(t);    
    mesh.LoadMesh (str, this->mesh_comm);
    UpdateBuffers();
  }
  */
  
  void MeshAccess :: SaveMesh (ostream & str) const
  {
    mesh.SaveMesh (str);
  }

  void MeshAccess :: DoArchive(Archive& ar)
  {
    auto mshptr = mesh.GetMesh();
    ar.Shallow(mshptr);
    if(ar.Input())
      {
        // Create the Ngx_Mesh interface
        mesh = {mshptr};
        UpdateBuffers();
      }
  }

  void MeshAccess :: SelectMesh() const
  {
    mesh.SelectMesh();
  }


  void MeshAccess :: GetSElNeighbouringDomains(const int elnr, int & in, int & out) const
  {
    ArrayMem<int, 2> elnums;
    auto fnums = GetElFacets(ElementId(BND,elnr));
    GetFacetElements ( fnums[0], elnums );
    if (elnums.Size()==0)
      {  // surface mesh only
        in = 0;
        out = 0;
      }
    else if (elnums.Size()==1)
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
        nnodes[NT_GLOBAL] = ( (GetCommunicator().Size() > 1) && (GetCommunicator().Rank() == 0) ? 0 : 1 );
        dim = -1;
        return;
      }

    auto netgen_mesh_timestamp = mesh.GetTimeStamp();
    // cout << "netgen_mesh_timestamp = " << netgen_mesh_timestamp << ", mytime = " << mesh_timestamp << endl;
    if (netgen_mesh_timestamp == mesh_timestamp) return;
    mesh_timestamp = netgen_mesh_timestamp;
    
    timestamp = NGS_Object::GetNextTimeStamp();
    

    dim = mesh.GetDimension();
    nlevels = mesh.GetNLevels(); 

    if (GetCommunicator().Size() > 1 && GetCommunicator().Rank() == 0)
      {
        for (int i = 0; i < 4; i++)  
          {
            nnodes[i] = 0;
            nelements[i] = 0;
            nnodes_cd[i] = 0;
            nelements_cd[i] = 0;
          }
        nnodes[NT_GLOBAL] = 0;
      }
    else
      {
	for (int i = 0; i < 4; i++)  
	  {
	    nnodes[i] = mesh.GetNNodes(i);
	    nelements[i] = mesh.GetNElements(i);
	  }
	for (int i = 0; i < 4; i++)  
	  {
	    nnodes_cd[i] = 0;
	    nelements_cd[i] = 0;
	  }
	for (int i = 0; i <= dim; i++)
	  {
	    nnodes_cd[i] = nnodes[dim-i];
	    nelements_cd[i] = nelements[dim-i];
	  }
        nnodes[NT_GLOBAL] = 1;
      }
    nnodes[NT_ELEMENT] = nnodes[StdNodeType (NT_ELEMENT, dim)];
    nnodes[NT_FACET] = nnodes[StdNodeType (NT_FACET, dim)];

    for (auto & p : trafo_jumptable) p = nullptr;
    Iterate<4> ([&](auto DIM)
                {
                  if (DIM.value == dim)
                    {
                      if constexpr (DIM.value >= 1)
                                     {
                                       trafo_jumptable[VOL] = &MeshAccess::GetTrafoDim<DIM>;
                                       trafo_jumptable[BND] = &MeshAccess::GetSTrafoDim<DIM>;
                                     }
                      if constexpr (DIM.value >= 2)
                                     trafo_jumptable[BBND] = &MeshAccess::GetCD2TrafoDim<DIM>;
                    }
                });
    

    
    int & ndomains = nregions[0];    
    ndomains = -1;
    // int ne = GetNE();
    
    auto minmax =
      ParallelReduce (GetNE(VOL),
                      [&] (size_t i)
                      {
                        auto ind = GetElIndex(ElementId(VOL, i));
                        return make_pair(ind, ind);
                      },
                      [] (pair<int,int> a, pair<int,int> b)
                      {
                        return make_pair(min2(a.first, b.first),
                                         max2(a.second, b.second));
                      },
                      pair<int,int> (std::numeric_limits<int>::max(),
                                     0));

    ndomains = minmax.second;
    if (minmax.first < 0)
      throw Exception("mesh with negative element-index");      
    /*
    for (int i = 0; i < ne; i++)
      {
        int elindex = GetElIndex(ElementId(VOL,i));
        if (elindex < 0) throw Exception("mesh with negative element-index");
        ndomains = max2(ndomains, elindex);
      }
    */

    ndomains++;
    
    // ndomains = MyMPI_AllReduce (ndomains, MPI_MAX, GetCommunicator());
    ndomains = GetCommunicator().AllReduce (ndomains, MPI_MAX);
    
    pml_trafos.SetSize(ndomains);
    
    int nboundaries = -1;
    for (auto el : Elements(BND))
      {
        int elindex = el.GetIndex();
        if (elindex < 0) throw Exception("mesh with negative boundary-condition number");
        nboundaries = max2(nboundaries, elindex);
      }
    nboundaries++;
    nboundaries = GetCommunicator().AllReduce (nboundaries, MPI_MAX);
    nregions[1] = nboundaries;


    int & nbboundaries = nregions[BBND];
    if(mesh.GetDimension() == 1)
      {
        nbboundaries = 0;
      }
    else
      {
        nbboundaries = -1;
        // int ncd2e = nelements_cd[2];
        // for (int i=0; i< ncd2e; i++)
        for (auto el : Elements(BBND))          
          {
            // ElementId ei(BBND, i);
            // int elindex = GetElIndex(ei);
            int elindex = el.GetIndex();            
            //if (elindex < 0) throw Exception ("mesh with negative cd2 condition number");
            if (elindex >=0)
              nbboundaries = max2(nbboundaries, elindex);
          }
        nbboundaries++;
        nbboundaries = GetCommunicator().AllReduce(nbboundaries, MPI_MAX);
      }

    int & nbbboundaries = nregions[BBBND];
    if(mesh.GetDimension() == 1 || mesh.GetDimension() == 2)
      {
        nbbboundaries = 0;
      }
    else
      {
        nbbboundaries = -1;
        for (auto el : Elements(BBBND))          
          {
            int elindex = el.GetIndex();
            if (elindex >=0)
              nbbboundaries = max2(nbbboundaries, elindex);
          }
        nbbboundaries++;
        nbbboundaries = GetCommunicator().AllReduce(nbbboundaries, MPI_MAX);
      }
    
    // update periodic mappings
    auto nid = mesh.GetNIdentifications();
    periodic_node_pairs[NT_VERTEX]->SetSize(0);
    periodic_node_pairs[NT_EDGE]->SetSize(0);
    periodic_node_pairs[NT_FACE]->SetSize(0);

    if (GetCommunicator().Size() > 1 && GetCommunicator().Rank() == 0)
      nid = 0; //hopefully this is enough...

    for (auto idnr : Range(nid))
      {
        // only if it is periodic
        if (mesh.GetIdentificationType(idnr)!=2) continue;
        auto pv_buffer = mesh.GetPeriodicVertices(idnr);
        auto pidnr = periodic_node_pairs[NT_VERTEX]->Append(Array<INT<2,int>>(pv_buffer.Size(), (INT<2,int>*) pv_buffer.Release(), true));

        // build vertex map for idnr
        Array<int> vertex_map(GetNV());
        for (auto i : Range(GetNV()))
          vertex_map[i] = i;
        for (const auto& pair : (*periodic_node_pairs[NT_VERTEX])[pidnr])
          vertex_map[pair[1]] = pair[0];

        // build vertex-pair to edge hashtable:
        HashTable<INT<2>, int> vp2e(GetNEdges());
        
        for (size_t enr = 0; enr < GetNEdges(); enr++)
          {
            INT<2> vts = GetEdgePNums (enr);
            vts.Sort();
            vp2e[vts] = enr;
          }
        size_t count = 0;
        for (size_t enr = 0; enr < GetNEdges(); enr++)
          {
            INT<2> vts = GetEdgePNums(enr);
            size_t mv1 = vertex_map[vts[0]];
            size_t mv2 = vertex_map[vts[1]];
            if(mv1 != vts[0] && mv2 != vts[1])
              count++;
          }
        periodic_node_pairs[NT_EDGE]->Append(Array<INT<2>>(count));
        count = 0;
        for (size_t enr = 0; enr < GetNEdges(); enr++)
          {
            INT<2> vts = GetEdgePNums (enr);
            int mv1 = vertex_map[vts[0]];
            int mv2 = vertex_map[vts[1]];
            if(mv1 != vts[0] && mv2 != vts[1])
              {               
                if (mv1 > mv2) Swap(mv1,mv2);
                int menr = vp2e.Get(INT<2>(mv1,mv2));
                (*periodic_node_pairs[NT_EDGE])[pidnr][count][0] = menr;
                (*periodic_node_pairs[NT_EDGE])[pidnr][count++][1] = enr;
              }
          }
        // build vertex-triple to face hashtable
        HashTable<INT<3>, int> v2f(GetNFaces());
        for (auto fnr : Range(GetNFaces()))
          {
            auto pnums = GetFacePNums (fnr);
            INT<3> i3(pnums[0], pnums[1], pnums[2]);
            i3.Sort();
            v2f[i3] = fnr;
          }

        count = 0;
        for (auto fnr : Range(GetNFaces()))
          {
            auto pnums = GetFacePNums(fnr);
            if(vertex_map[pnums[0]] != pnums[0] && vertex_map[pnums[1]] != pnums[1] &&
               vertex_map[pnums[2]] != pnums[2])
              {
                count++;
              }
          }
        periodic_node_pairs[NT_FACE]->Append(Array<INT<2>>(count));
        count = 0;
        for (auto fnr : Range(GetNFaces()))
          {
            auto pnums = GetFacePNums(fnr);
            INT<3> mv(vertex_map[pnums[0]],vertex_map[pnums[1]],vertex_map[pnums[2]]);
            if(mv[0] != pnums[0] && mv[1] != pnums[1] && mv[2] != pnums[2])
              {
                mv.Sort();
                int mfnr = v2f[mv];
                (*periodic_node_pairs[NT_FACE])[pidnr][count][0] = mfnr;
                (*periodic_node_pairs[NT_FACE])[pidnr][count++][1] = fnr;
              }
          }
      }
    
    CalcIdentifiedFacets();
  }

  void MeshAccess :: BuildNeighbours()
  {
    static Timer timer_neighbours("Build neighbours");
    RegionTimer reg(timer_neighbours);

    const auto& nmesh = *GetNetgenMesh();
    const auto& topology = nmesh.GetTopology();

    for(auto i : Range(4))
      for(auto j : Range(4))
        neighbours[i][j].SetSize(nregions[i]);

    Array<int> edgemap(GetNEdges());
    edgemap = -1;
    for(auto seg : Elements(VorB(dim-1)))
      edgemap[seg.edges[0]] = seg.GetIndex();

    Array<int> vertmap(GetNV());
    vertmap = -1;
    for(auto pel : Elements(VorB(dim)))
      vertmap[pel.Vertices()[0]] = pel.GetIndex();

    if(GetDimension() == 3)
      {
        for(auto sei : Range(nmesh.SurfaceElements()))
          {
            int el1, el2;
            topology.GetSurface2VolumeElement(sei+1, el1, el2);
            const auto& sel = nmesh.SurfaceElements()[sei];
            auto bc = nmesh.GetFaceDescriptor(sel.GetIndex()).BCProperty()-1;
            if(el1 > 0)
              {
                auto index1 = nmesh.VolumeElement(el1).GetIndex()-1;
                neighbours[BND][VOL].AddUnique(bc, index1);
                neighbours[VOL][BND].AddUnique(index1, bc);
              }
            if(el2 > 0)
              {
                auto index2 = nmesh.VolumeElement(el2).GetIndex()-1;
                neighbours[BND][VOL].AddUnique(bc, index2);
                neighbours[VOL][BND].AddUnique(index2, bc);
              }
          }
        for(auto ei : Elements(VOL))
          {
            const auto& el = GetElement(ei);
            auto index = el.GetIndex();
            for (auto edge : el.Edges())
              {
                if(auto eindex = edgemap[edge]; eindex != -1)
                  {
                    neighbours[VOL][BBND].AddUnique(index, eindex);
                    neighbours[BBND][VOL].AddUnique(eindex, index);
                  }
              }
            for(auto v : el.Vertices())
              {
                if(auto vindex = vertmap[v]; vindex != -1)
                  {
                    neighbours[BBBND][VOL].AddUnique(vindex, index);
                    neighbours[VOL][BBBND].AddUnique(index, vindex);
                  }
              }
          }
      }
    auto edge_vb = VorB(GetDimension()-1);
    auto point_vb = VorB(GetDimension());
    if(GetDimension() >= 2)
      {
        auto surf_vb = VorB(GetDimension()-2);
        for(auto sei : Elements(surf_vb))
          {
            const auto& sel = GetElement(sei);
            auto index = sel.GetIndex();
            for(auto edge : sel.Edges())
              {
                if(auto eindex = edgemap[edge]; eindex != -1)
                  {
                    neighbours[surf_vb][edge_vb].AddUnique(index, eindex);
                    neighbours[edge_vb][surf_vb].AddUnique(eindex, index);
                  }
              }
            for(auto v : sel.Vertices())
              {
                if(auto vindex = vertmap[v]; vindex != -1)
                  {
                    neighbours[surf_vb][point_vb].AddUnique(index, vindex);
                    neighbours[point_vb][surf_vb].AddUnique(vindex, index);
                  }
              }
          }
      }
    if(GetDimension() >= 1)
      {
        for(auto ei : Elements(edge_vb))
          {
            const auto& seg = GetElement(ei);
            auto index = seg.GetIndex();
            for(auto v : seg.Vertices())
              {
                if(auto vindex = vertmap[v]; vindex != -1)
                  {
                    neighbours[edge_vb][point_vb].AddUnique(index, vindex);
                    neighbours[point_vb][edge_vb].AddUnique(vindex, index);
                  }
              }
          }
      }

    // same codim neighbours have common codim-1 neighbour
    for(auto i : Range(4))
      if(GetDimension() - i > 0)
        for(auto j : Range(nregions[i]))
          for(auto bnd : neighbours[i][i+1][j])
            for(auto mat : neighbours[i+1][i][bnd])
              if(mat != j)
                neighbours[i][i].AddUnique(j, mat);
  }

  void MeshAccess :: 
  GetElEdges (int elnr, Array<int> & ednums, Array<int> & orient) const
  {
    /*
    ednums.SetSize (12);
    orient.SetSize (12);
    int ned = 
      Ng_GetElement_Edges (elnr+1, &ednums[0], &orient[0]);
    ednums.SetSize (ned);
    orient.SetSize (ned);
    for (int i = 0; i < ned; i++)
      ednums[i]--;
    */
    auto el = GetElement(ElementId(VOL, elnr));
    auto vnums = el.Vertices();
    auto locedges = ElementTopology::GetEdges(el.GetType());
    ednums = el.Edges();
    orient.SetSize(ednums.Size());
    for (auto i : Range(ednums))
      orient[i] = vnums[locedges[i][1]] > vnums[locedges[i][0]] ? 1 : -1;
  }

  void MeshAccess :: 
  GetSElEdges (int selnr, Array<int> & ednums, Array<int> & orient) const
  {
    /*
    ednums.SetSize (4);
    orient.SetSize (4);
    int ned = 
      Ng_GetSurfaceElement_Edges (selnr+1, &ednums[0], &orient[0]);
    ednums.SetSize (ned);
    orient.SetSize (ned);

    for (int i = 0; i < ned; i++)
      ednums[i]--;
    */
    auto el = GetElement(ElementId(BND, selnr));
    auto vnums = el.Vertices();
    auto locedges = ElementTopology::GetEdges(el.GetType());
    ednums = el.Edges();
    orient.SetSize(ednums.Size());
    for (auto i : Range(ednums))
      orient[i] = vnums[locedges[i][1]] > vnums[locedges[i][0]] ? 1 : -1;

  }


  void MeshAccess :: GetEdgeElements (int enr, Array<int> & elnums) const
  {
    // static Timer t("getedgeelements"); RegionTimer reg(t);    
    elnums.SetSize0();

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

    INT<2> pts = GetEdgePNums(enr);

    auto velems0 = GetVertexSurfaceElements(pts[0]);
    auto velems1 = GetVertexSurfaceElements(pts[1]);
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
    /*
    edges.SetSize(4);
    int ned = Ng_GetFace_Edges (fnr+1, &edges[0]);
    edges.SetSize(ned);
    for (int i = 0; i < ned; i++) edges[i]--;
    */
    edges = ArrayObject(mesh.GetFaceEdges(fnr));
  }
  
  void MeshAccess :: GetEdgeFaces (int enr, Array<int> & faces) const
  {
    faces.SetSize0();
    auto [v0,v1] = GetEdgePNums(enr);
    // auto v01 = GetEdgePNums(enr);
    // auto v0 = v01[0], v1 = v01[1];

    for (auto elnr : GetVertexElements(v0))
      {
        auto el = GetElement(ElementId(VOL,elnr));
        if (el.Vertices().Contains(v1))
          for (auto fa : el.Faces())
            {
              auto vnums = GetFacePNums(fa);
              if (vnums.Contains(v0) && vnums.Contains(v1) && !faces.Contains(fa))
                faces.Append(fa);
            }
      }
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
    // ArrayMem<int, 9> vnums;
    auto vnums = GetFacePNums(fnr);

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

  
  void MeshAccess :: GetEdgePNums (int enr, Array<int> & pnums) const
  {
    pnums.SetSize(2);
    auto edge = mesh.GetNode<1>(enr);
    pnums[0] = edge.vertices[0];
    pnums[1] = edge.vertices[1];
  }

  
  void MeshAccess :: GetFacetPNums (int fnr, Array<int> & pnums) const
  {
    switch (dim)
      {
      case 1: pnums.SetSize(1); pnums[0] = fnr; break;
      case 2: pnums = GetEdgePNums(fnr); break;
      case 3: pnums = GetFacePNums(fnr); break;
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
    static Timer t("CalcIdentifiedFacets"); RegionTimer reg(t);
    identified_facets.SetSize(nnodes_cd[1]);
    for(auto i : Range(identified_facets.Size()))
      identified_facets[i] = std::tuple<int,int>(i,1);
 
    // for periodic identification by now
    for(auto id : Range(GetNPeriodicIdentifications()))
      {
        auto dim  = mesh.GetDimension();
	const auto& pairs = GetPeriodicNodes(dim == 3 ? NT_FACE : (dim == 2 ? NT_EDGE : NT_VERTEX),id);

	for(const auto& pair : pairs)
	  for(auto l : Range(2))
	    identified_facets[pair[l]] = std::tuple<int,int>(pair[1-l],2);
      }
  }

  void MeshAccess :: GetVertexElements (size_t vnr, Array<int> & elnrs) const
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

    void MeshAccess :: SetDeformation (shared_ptr<GridFunction> def)
    {
      if (def)
        {
          // if (dim  != def->GetFESpace()->GetDimension())
          // now also VectorH1 is possible !
          if (dim != def->Dimension())
            throw Exception ("Mesh::SetDeformation needs a GridFunction with dim="+ToString(dim));
        }
      deformation = def;
    }
  
    void MeshAccess :: SetPML (const shared_ptr<PML_Transformation> & pml_trafo, int _domnr)
    {
      if (_domnr>=nregions[VOL])
        throw Exception("MeshAccess::SetPML: was not able to set PML, domain index too high!");
      if (pml_trafo->GetDimension()!=dim)
        throw Exception("MeshAccess::SetPML: dimension of PML = "+ToString(pml_trafo->GetDimension())+" does not fit mesh dimension!");
      pml_trafos[_domnr] = pml_trafo; 
    }
    
    void MeshAccess :: UnSetPML (int _domnr)
    {
      if (_domnr>=nregions[VOL])
        throw Exception("MeshAccess::UnSetPML: was not able to unset PML, domain index too high!");
      pml_trafos[_domnr] = nullptr; 
    }
    Array<shared_ptr<PML_Transformation>> & MeshAccess :: GetPMLTrafos()
    { return pml_trafos; }

    shared_ptr<PML_Transformation> MeshAccess :: GetPML(int _domnr)
    {
      if (_domnr>=nregions[VOL])
        throw Exception("MeshAccess::GetPML: was not able to get PML, domain index too high!");
      return pml_trafos[_domnr];
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
             loc_deformation, lh); 
        // dynamic_cast<LocalHeap&> (lh));

        else

          eltrans = new (lh)
            ALE_ElementTransformation<DIM, DIM, Ng_ConstElementTransformation<DIM,DIM>>
            (this, el.GetType(), 
             ElementId(VOL,elnr), el.GetIndex(),
             loc_deformation, lh);
        //dynamic_cast<LocalHeap&> (lh));
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
         loc_deformation, lh);
         // dynamic_cast<LocalHeap&> (lh)); 
    
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
	 loc_deformation, lh);
    // dynamic_cast<LocalHeap&>(lh));
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

  ngfem::ElementTransformation & MeshAccess :: GetTrafoOld (ElementId ei, Allocator & lh) const
  {
    auto elnr = ei.Nr();
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

      case BBBND: default:
        {
          ElementTransformation * eltrans;
          Ngs_Element el(mesh.GetElement<0>(elnr), ElementId(BBBND,elnr));
          GridFunction * loc_deformation = deformation.get();
          if(loc_deformation)
            eltrans = new (lh) ALE_ElementTransformation<0,3, Ng_ElementTransformation<0,3>>
              (this,el.GetType(),
               ElementId(BBBND,elnr), el.GetIndex(),
               loc_deformation, lh);
          // dynamic_cast<LocalHeap&>(lh));
          else 
            eltrans = new (lh) Ng_ConstElementTransformation<0,3> (this, el.GetType(),
                                                                   ElementId(BBBND,elnr), el.GetIndex());
          
          return *eltrans;
        }
      }
  }


  Region :: Region (const shared_ptr<MeshAccess> & amesh,
                    VorB avb, string pattern)
    : mesh(amesh), vb(avb)
  {
    mask = make_shared<BitArray>(mesh->GetNRegions(vb));
    mask->Clear();
    regex re_pattern(pattern);

    for (int i : Range(*mask))
      if (regex_match(mesh->GetMaterial(vb,i), re_pattern))
        mask->SetBit(i);
  }      

  Region :: Region (const shared_ptr<MeshAccess> & amesh, VorB avb, const BitArray & amask)
    : mesh(amesh), vb(avb), mask(make_shared<BitArray>(amask))
  { ; }

  Region :: Region (const shared_ptr<MeshAccess> & amesh, VorB avb, bool all)
    : mesh(amesh), vb(avb),
      mask(make_shared<BitArray>(amesh->GetNRegions(avb)))
  {
    if(all)
      mask->Set();
    else
      mask->Clear();
  }

  Region Region :: GetNeighbours(VorB other_vb)
  {
    if(mesh->neighbours[vb][other_vb].Size() == 0)
      {
        static mutex calc_neighbour_mutex;
        lock_guard<mutex> guard(calc_neighbour_mutex);
        if(mesh->neighbours[vb][other_vb].Size() == 0)
          mesh->BuildNeighbours();
      }
    Region neighbours(mesh, other_vb);
    for(auto i : Range(mask->Size()))
      if(mask->Test(i))
        for(auto j : mesh->neighbours[vb][other_vb][i])
          neighbours.Mask().SetBit(j);
    return neighbours;
  }

  bool Region :: operator==(const Region& other) const
  {
    if(mesh.get() != other.Mesh().get())
      return false;
    if(vb != other.VB())
      return false;
    return *mask == other.Mask();
  }

  size_t Region :: Hash() const
  {
    HashArchive ar;
    ar & *mask;
    return ar.GetHash();
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
    RegionTracer reg(TaskManager::GetThreadId(), t);


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
    static Timer t("FindSurfaceElementOfPoint");
    RegionTimer reg(t);
    int elnr = -1;

    if(indices != NULL && indices->Size()>0)
      {
	// int * dummy = new int[indices->Size()];
	// for(int i=0; i<indices->Size(); i++) dummy[i] = (*indices)[i]+1;
	// elnr = Ng_FindSurfaceElementOfPoint (&point(0), lami, build_searchtree,dummy,indices->Size());
	// delete [] dummy;
        switch (dim)
          {
          case 1:
            throw Exception("FindSurfaceElementOfPoint for mesh-dim = 1 not implemented yet!");
            break;
          case 2:
            elnr = mesh.FindElementOfPoint<1> (&point(0), &ip(0), build_searchtree,
                                               &(*indices)[0],indices->Size());
            break;
          case 3:
            elnr = mesh.FindElementOfPoint<2> (&point(0), &ip(0), build_searchtree,
                                               &(*indices)[0],indices->Size());
            break;
          }
      }
    else
      {  
	// elnr = Ng_FindSurfaceElementOfPoint (&point(0), lami, build_searchtree);
        switch (dim)
          {
          case 1:
            throw Exception("FindSurfaceElementOfPoint for mesh-dim = 1 not implemented yet!");
            break;
          case 2:
            elnr = mesh.FindElementOfPoint<1> (&point(0), &ip(0), build_searchtree, NULL, 0);
            break;
          case 3:
            elnr = mesh.FindElementOfPoint<2> (&point(0), &ip(0), build_searchtree, NULL, 0);
           break;
          }
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

  map<string, unique_ptr<Timer<>>> ngtimers;
  void NGSolveTracer (string name, bool stop)
  {
    // cout << "************* tracer: " << name << ", stop = " << stop << endl;
    int count = ngtimers.count(name);
    if (count == 0)
      ngtimers[name] = make_unique<Timer<>> (name);
    Timer<> * timer = ngtimers[name].get();
    if (!stop)
      timer->Start();
    else
      timer->Stop();
  }
  
  
  void MeshAccess :: Refine (bool onlyonce)
  {
    static Timer t("MeshAccess::Refine"); RegionTimer reg(t);
    nlevels = std::numeric_limits<int>::max();
    mesh.Refine(NG_REFINE_H, onlyonce, &NGSolveTaskManager, &NGSolveTracer);
    UpdateBuffers();
    updateSignal.Emit();
  }

  void MeshAccess :: Curve (int order)
  {
    mesh.Curve(order);
  } 
  
  int MeshAccess :: GetCurveOrder ()
  {
    return mesh.GetCurveOrder();
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

  size_t MeshAccess :: GetNPeriodicNodes(NODE_TYPE nt) const
  {
    size_t npairs = 0;
    for(const auto& ar : (*periodic_node_pairs[nt]))
      npairs += ar.Size();
    return npairs;
  }

  void MeshAccess :: GetPeriodicNodes(NODE_TYPE nt, Array<INT<2>>& pairs) const
  {
    pairs.SetSize(0);
    pairs.SetAllocSize(GetNPeriodicNodes(nt));
    for(const auto& ar : (*periodic_node_pairs[nt]))
      for(const auto& val : ar)
        pairs.Append(val);
  }

  const Array<INT<2>> & MeshAccess :: GetPeriodicNodes (NODE_TYPE nt, int idnr) const
  {
    if(idnr > periodic_node_pairs[nt]->Size())
      throw Exception("Not enough periodic identifications, GetPeriodicNodes is 0 based and has only periodic Identifications!");
    return (*periodic_node_pairs[nt])[idnr];
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

  void MeshAccess::GetVertexSurfaceElements (size_t vnr, Array<int> & elems) const
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


  const Table<size_t> & MeshAccess :: GetElementsOfClass()
  {
    if (int(elements_of_class_timestamp) >= mesh.GetTimeStamp())
      return elements_of_class;

    Array<short> classnr(GetNE());
    LocalHeap lh(10000);  // should use ParallelForRange, then we don't need the lh
    IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         classnr[el.Nr()] = 
           SwitchET<ET_TRIG,ET_TET>
           (el.GetType(),
            [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
       });
    
    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);

    elements_of_class = creator.MoveTable();
    elements_of_class_timestamp = mesh.GetTimeStamp();
    return elements_of_class;
  }

  shared_ptr<CoefficientFunction> MeshAccess ::
  RegionCF(VorB vb, shared_ptr<CoefficientFunction> default_value, const Array<pair<variant<string, Region>, shared_ptr<CoefficientFunction>>>& region_values)
  {
    Array<shared_ptr<CoefficientFunction>> cfs(GetNRegions(vb));
    shared_ptr<MeshAccess> spm(this, NOOP_Deleter);
    for(const auto& val : region_values)
      {
        const auto& key = val.first;
        Region region;
        if(auto skey = get_if<string>(&key); skey)
          region = Region(spm, vb, *skey);
        else
          {
            region = *get_if<Region>(&key);
            if(region.VB() != vb)
              throw Exception("Region with false vb given to region wise CF!");
          }
        for(auto i : Range(region.Mask().Size()))
          {
            if(region.Mask()[i])
              cfs[i] = val.second;
          }
      }
    for(auto i : Range(cfs))
      if(cfs[i] == nullptr)
        cfs[i] = default_value;
    return MakeDomainWiseCoefficientFunction(move(cfs));
  }

  class BoundaryFromVolumeCoefficientFunction :
    public T_CoefficientFunction<BoundaryFromVolumeCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> vol_cf;
    typedef T_CoefficientFunction<BoundaryFromVolumeCoefficientFunction> BASE;

  public:
    BoundaryFromVolumeCoefficientFunction() = default;
    BoundaryFromVolumeCoefficientFunction (shared_ptr<CoefficientFunction> avol_cf)
      : BASE(avol_cf->Dimension(), avol_cf->IsComplex()), vol_cf(avol_cf)
    {
      SetDimensions(avol_cf->Dimensions());
    }

    ~BoundaryFromVolumeCoefficientFunction() { ; }
    
    void DoArchive(Archive& ar) override
    {
      BASE::DoArchive(ar);
      ar.Shallow(vol_cf);
    }

    virtual bool DefinedOn (const ElementTransformation & trafo) override
    {
      auto ei = trafo.GetElementId();
      if (ei.IsVolume())
        return vol_cf->DefinedOn(trafo);

      const MeshAccess & ma = *static_cast<const MeshAccess*> (trafo.GetMesh());
      int facet = ma.GetElFacets(ei)[0];
      ArrayMem<int,2> elnums;
      ma.GetFacetElements (facet, elnums);
      LocalHeapMem<1000> lh("definedonlh");

      for (auto elnr : elnums)
        {
          ElementId vei(VOL, elnr);
          ElementTransformation & vol_trafo = ma.GetTrafo (vei, lh);
          if (vol_cf->DefinedOn(vol_trafo))
            return true;
        }

      return false;
    }

    using BASE::Evaluate;
    template <typename MIR, typename T, ORDERING ORD>            
    void T_Evaluate (const MIR & ir, BareSliceMatrix<T,ORD> values) const
    {
      // static Timer t("BFV - evaluate T-eval", 1);
      // RegionTimer reg(t); // , TaskManager::GetThreadId());
      
      LocalHeapMem<100000> lh("lhbfv");
      const ElementTransformation & trafo = ir.GetTransformation();
      auto ei = trafo.GetElementId();
      if (ei.IsVolume())
        {
          vol_cf->Evaluate(ir, values);
          return;
        }
      const MeshAccess & ma = *static_cast<const MeshAccess*> (trafo.GetMesh());
      int facet = ma.GetElFacets(ei)[0];
      ArrayMem<int,2> elnums;
      ma.GetFacetElements (facet, elnums);

      for (auto elnr : elnums)
        {
          ElementId vei(VOL, elnr);
          int locfacnr = ma.GetElFacets(vei).Pos(facet);

          ElementTransformation & vol_trafo = ma.GetTrafo (vei, lh);
          if (!vol_cf->DefinedOn(vol_trafo)) continue;

          Facet2ElementTrafo f2el(vol_trafo.GetElementType(), ma.GetElVertices(vei));
          Array<int> surfvnums { ma.GetElVertices(ei) };
          Facet2SurfaceElementTrafo f2sel(trafo.GetElementType(), surfvnums);
          auto & ir_ref = f2sel.Inverse(ir.IR(), lh);
          auto & ir_vol = f2el(locfacnr, ir_ref, lh);
          auto & mir_vol = vol_trafo(ir_vol, lh);
          mir_vol.ComputeNormalsAndMeasure (vol_trafo.GetElementType(), locfacnr);
          vol_cf -> Evaluate (mir_vol, values);
          return;
        }

    }
    
    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    {
      throw Exception ("Surface-from-Volume inout not implemented");
    }
  };

  shared_ptr<CoefficientFunction>
  MakeBoundaryFromVolumeCoefficientFunction  (shared_ptr<CoefficientFunction> avol_cf)
  {
    return make_shared<BoundaryFromVolumeCoefficientFunction> (avol_cf);
  }

  class TrafoCF : public CoefficientFunctionNoDerivative
  {
    shared_ptr<CoefficientFunction> func, trafo;
    Region region;
    unique_ptr<netgen::BoxTree<3>> searchtree;
  public:
    TrafoCF(shared_ptr<CoefficientFunction> _func,
            shared_ptr<CoefficientFunction> _trafo,
            Region _region)
      : CoefficientFunctionNoDerivative(_func->Dimension()),
        func(_func), trafo(_trafo), region(_region)
    {
      if(region.VB() == BND)
        {
          netgen::Box<3> box(netgen::Box<3>::EMPTY_BOX);
          for(const auto& el : region.GetElements())
            for(auto v : el.Points())
              {
                auto p = region.Mesh()->GetPoint<3>(v);
                box.Add({ p[0], p[1], p[2] });
              }

          box.Scale(1.4);
          searchtree = make_unique<netgen::BoxTree<3>>(box);
          for(const auto& el : region.GetElements())
            {
              netgen::Box<3> box (netgen::Box<3>::EMPTY_BOX);
              for(auto v : el.Points())
                {
                  auto p = region.Mesh()->GetPoint<3>(v);
                  box.Add({ p[0], p[1], p[2] });
                }
              box.Scale(1.4);
              searchtree->Insert(box, el.Nr());
            }
        }
    }

    using CoefficientFunctionNoDerivative::Evaluate;
    double Evaluate(const BaseMappedIntegrationPoint& mip) const override
    {
      Vec<3> mapped_point;
      trafo->Evaluate(mip, mapped_point);
      IntegrationPoint ip;
      int elnr;

      if(region.VB() == BND)
        {
      netgen::Point<3> p = {mapped_point[0], mapped_point[1], mapped_point[2] };
      searchtree->GetFirstIntersecting(p, p, [&](int si)
      {
        if(region.Mesh()->GetNetgenMesh()
           ->PointContainedIn2DElement({mapped_point[0],
               mapped_point[1],
               mapped_point[2]},
             &ip(0),
             si+1, false))
          {
            elnr = si;
            return true;
          }
        return false;
      });
      if(region.Mesh()->GetNetgenMesh()->SurfaceElement(elnr+1).GetType() == netgen::TRIG)
        {
          double lam0 = 1-ip(0)-ip(1);
          ip(1) = ip(0);
          ip(0) = lam0;
        }
        }
      else if(region.VB() == VOL)
        {
          Array<int> indices;
          for(auto i : Range(region.Mask().Size()))
            if(region.Mask().Test(i))
              indices.Append(i);
          elnr = region.Mesh()->FindElementOfPoint(mapped_point, ip, true, &indices);
        }
      else
        throw Exception("Only VOL and BND implemented yet!");
      auto& eltrafo = region.Mesh()->GetTrafo(ElementId(region.VB(), elnr), global_alloc);
      auto& mapped_mip = eltrafo(ip, global_alloc);
      return func->Evaluate(mapped_mip);
    }
  };

  shared_ptr<CoefficientFunction> MakeTrafoCF(shared_ptr<CoefficientFunction> func,
                                              shared_ptr<CoefficientFunction> trafo,
                                              Region region)
  {
    return make_shared<TrafoCF>(func, trafo, region);
  }


  

#ifdef PARALLEL
 
  size_t MeshAccess ::GetGlobalNodeNum (NodeId node) const
  {
    /** There is only one global node **/
    if (node.GetType() == NT_GLOBAL)
      { return 0; }
    int glob = NgPar_GetGlobalNodeNum (node.GetType(), node.GetNr());
    return glob;
  }

  size_t MeshAccess :: GetGlobalVertexNum (int locnum) const
  {
    // return NgPar_GetGlobalNodeNum (NT_VERTEX, locnum);
    return mesh.GetGlobalVertexNum (locnum);
  }
  
  void MeshAccess :: GetDistantProcs (NodeId node, Array<int> & procs) const
  {
    if (GetCommunicator().Size() > 1) {
      if (node.GetType() == NT_GLOBAL) {
	/** The global node is shared by everyone except rank 0 **/
	if (GetCommunicator().Rank() == 0) {
	  procs.SetSize0();
	  return;
	}
	else {
	  procs.SetSize(GetCommunicator().Size()-2);
	  for (int k = 1; k < GetCommunicator().Rank(); k++)
	    { procs[k-1] = k; }
	  for (int k = GetCommunicator().Rank() + 1; k < GetCommunicator().Size(); k++)
	    { procs[k-2] = k; }
	}
      }
      else { // node.GetType() == NT_GLOBAL
	auto dps = mesh.GetDistantProcs(node.GetType(), node.GetNr());
	procs.SetSize(dps.Size());
	for (auto k : Range(procs.Size()))
          procs[k] = dps[k]; 
      }
    }
    else // GetCommunicator().Size() > 1
      { procs.SetSize0(); }
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


  function<void()> cleanup_func = ProgressOutput :: SumUpLocal;
  ProgressOutput :: ProgressOutput (shared_ptr<MeshAccess> ama,
				    string atask, size_t atotal)
    : ma(ama), comm(ama->GetCommunicator()), task(atask), total(atotal)
  {
    use_mpi = false;   // hardcoded for the moment 
    is_root = comm.Rank() == 0;
    prevtime = WallTime();

    if (use_mpi)
      {
        size_t glob_total = comm.Reduce (total, MPI_SUM);
        if (is_root) total = glob_total;
      }
    
    done_called = false;
    cnt = 0;
    thd_cnt = 0;
    // cleanup_func = [this] () {  this->SumUpLocal(); };
    TaskManager::SetCleanupFunction(cleanup_func);
  }

  ProgressOutput :: ~ProgressOutput ()
  {
    Done();
    TaskManager::SetCleanupFunction();
  }  

  atomic<size_t> ProgressOutput :: cnt;
  thread_local size_t ProgressOutput :: thd_cnt = 0;
  // thread_local double ProgressOutput :: thd_prev_time = WallTime();
  thread_local size_t ProgressOutput :: thd_prev_time = GetTimeCounter();
  size_t tsc_wait = 0.05*(1.0/seconds_per_tick);
  void ProgressOutput :: Update ()
  {
    thd_cnt++;
    // double time = WallTime();
    size_t time = GetTimeCounter();
    // if (time > thd_prev_time+0.05)
    if (time > thd_prev_time+tsc_wait)
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
	  else if (use_mpi)
	    {
	      static Timer t("dummy - progressreport"); RegionTimer r(t);
	      comm.Send (nr, 0, MPI_TAG_SOLVE);
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
	int ntasks = comm.Size();
	if (ntasks > 1 && use_mpi)
	  {
	    Array<int> working(ntasks), computed(ntasks);
	    working = 1;
	    computed = 0;
	    while (1)
	      {
		int flag, num_working = 0, got_flag = false;
		size_t data;
		for (int source = 1; source < ntasks; source++)
		  {
		    if (!working[source]) continue;
		    num_working++;
		    MPI_Iprobe (source, MPI_TAG_SOLVE, comm, &flag, MPI_STATUS_IGNORE);
		    if (flag)
		      {
			got_flag = true;
			comm.Recv (data, source, MPI_TAG_SOLVE);
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
    else if (use_mpi)
      {
	comm.Send (total, 0, MPI_TAG_SOLVE);
	size_t final = -1;
	comm.Send (final, 0, MPI_TAG_SOLVE);
      }
  }
  
static RegisterClassForArchive<BoundaryFromVolumeCoefficientFunction, CoefficientFunction> regbfvcf;
static RegisterClassForArchive<Region> regregion;

}




