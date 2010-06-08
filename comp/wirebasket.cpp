#include <comp.hpp>
#include <solve.hpp>


namespace ngcomp
{
  template <class SCAL> class ElementByElementMatrix;

  template <class SCAL>
  ElementByElement_BilinearForm<SCAL> :: 
  ElementByElement_BilinearForm (const FESpace & afespace, const string & aname,
				 const Flags & flags)
    : S_BilinearForm<SCAL> (afespace, aname, flags)
  { ; }

  template <class SCAL>
  ElementByElement_BilinearForm<SCAL> :: ~ElementByElement_BilinearForm ()
  { ; }



  
  template <class SCAL>
  void ElementByElement_BilinearForm<SCAL> :: AllocateMatrix ()
  {
    cout << "alloc matrix" << endl;
    const FESpace & fespace = this->fespace;
    this->mats.Append (new ElementByElementMatrix<SCAL> (fespace.GetNDof(), this->ma.GetNE()+this->ma.GetNSE() ));
  }


  template<class SCAL>
  BaseVector * ElementByElement_BilinearForm<SCAL> :: CreateVector() const
  {
    return new VVector<SCAL> (this->fespace.GetNDof());
  }

  template<class SCAL>
  void ElementByElement_BilinearForm<SCAL> :: AddElementMatrix (const Array<int> & dnums1,
								const Array<int> & dnums2,
								const FlatMatrix<SCAL> & elmat,
								bool inner_element, int elnr,
								LocalHeap & lh)
  {
    cout << "ebe bilinearform ???" << endl;
    /*
      (*testout) << "inner_element = " << inner_element << endl;
      (*testout) << "elnr = " << elnr << endl;
      (*testout) << "elmat = " << endl << elmat << endl;
      (*testout) << "dnums1 = " << endl << dnums1 << endl;
    */
    
    dynamic_cast<ElementByElementMatrix<SCAL>&> (this->GetMatrix()) . AddElementMatrix (elnr, dnums1, dnums2, elmat);
  }
  


  template class ElementByElement_BilinearForm<double>;
  template class ElementByElement_BilinearForm<Complex>;




  
  template <class SCAL>
  ElementByElementMatrix<SCAL> :: ElementByElementMatrix (int h, int ane) 
  {
    height = h; 
    ne = ane; 
    elmats.SetSize(ne);
    dnums.SetSize(ne);
    for (int i = 0; i < ne; i++)
      {
        elmats[i].AssignMemory (0, 0, NULL);
        dnums[i] = FlatArray<int> (0, NULL);
      }
  }
  
  
  template <class SCAL>
  void ElementByElementMatrix<SCAL> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static int timer = NgProfiler::CreateTimer ("EBE-matrix::MultAdd");
    NgProfiler::RegionTimer reg (timer);

    int maxs = 0;
    for (int i = 0; i < dnums.Size(); i++)
      maxs = max2 (maxs, dnums[i].Size());

    ArrayMem<SCAL, 100> mem1(maxs), mem2(maxs);
      
    FlatVector<SCAL> vx = dynamic_cast<const S_BaseVector<SCAL> & >(x).FVScal();
    FlatVector<SCAL> vy = dynamic_cast<S_BaseVector<SCAL> & >(y).FVScal();

    for (int i = 0; i < dnums.Size(); i++)
      {
        FlatArray<int> di (dnums[i]);
        FlatVector<SCAL> hv1(di.Size(), &mem1[0]);
        FlatVector<SCAL> hv2(di.Size(), &mem2[0]);
	  
        for (int j = 0; j < di.Size(); j++)
          hv1(j) = vx (di[j]);

        hv2 = elmats[i] * hv1;
        hv2 *= s;

        for (int j = 0; j < dnums[i].Size(); j++)
          vy (di[j]) += hv2[j];
      }
  }

  template <class SCAL>
  BaseMatrix *  ElementByElementMatrix<SCAL> :: InverseMatrix ( BitArray * subset ) const
  {
    cout << "wird das tatsaechlich verwendet ???" << endl;
    ElementByElementMatrix<SCAL> * invmat = new ElementByElementMatrix<SCAL> (height, ne);

    int maxs = 0;
    for (int i = 0; i < dnums.Size(); i++)
      maxs = max2 (maxs, dnums[i].Size());

    LocalHeap lh (maxs*maxs*sizeof(SCAL)+100);

    for ( int i = 0; i < dnums.Size(); i++ )
      {
        int nd = dnums[i] . Size();
        FlatMatrix<SCAL> mat(nd, nd, lh);
        Array<int> dnumsarray(nd);
        for ( int j = 0; j < nd; j++ )
          dnumsarray[j] = dnums[i][j];
        mat = elmats[i];

        LapackInverse(mat);

        invmat -> AddElementMatrix(i, dnumsarray, dnumsarray, mat);
        lh.CleanUp();
      }
  
    return invmat;
  }
  
  template <class SCAL>
  void ElementByElementMatrix<SCAL> :: AddElementMatrix (int elnr,
                                                         const Array<int> & dnums1,
                                                         const Array<int> & dnums2,
                                                         const FlatMatrix<SCAL> & elmat)
  {
    ArrayMem<int,50> used;
    for (int i = 0; i < dnums1.Size(); i++)
      if (dnums1[i] >= 0) used.Append(i);

    int s = used.Size();

    FlatMatrix<SCAL> mat (s, new SCAL[s*s]);
    for (int i = 0; i < s; i++)
      for (int j = 0; j < s; j++)
        mat(i,j) = elmat(used[i], used[j]);

    FlatArray<int> dn(s, new int[s]);
    for (int i = 0; i < s; i++)
      dn[i] = dnums1[used[i]];

    if (elnr < elmats.Size())
      {
        dnums[elnr] = dn;
        elmats[elnr].AssignMemory (s, s, &mat(0,0));
      }
  }


  template <class SCAL>
  BaseBlockJacobiPrecond * ElementByElementMatrix<SCAL> :: 
  CreateBlockJacobiPrecond (Table<int> & blocks,
                            const BaseVector * constraint, int * paralleloptions) const
  { 
    return 0;//new helmholtz_exp_cpp::BlockJacobiPrecond_ElByEl<SCAL> (*this, blocks );
  }

  template class ElementByElementMatrix<double>;
  template class ElementByElementMatrix<Complex>;

  template <class SCAL>
  HO_BilinearForm<SCAL> :: HO_BilinearForm (const FESpace & afespace, const string & aname,
                                            const Flags & flags)
    : S_BilinearForm<SCAL> (afespace, aname, flags)
  { ; }

  template <class SCAL>
  HO_BilinearForm<SCAL> :: ~HO_BilinearForm ()
  { ; }



  
  template <class SCAL>
  void HO_BilinearForm<SCAL> :: AllocateMatrix ()
  {
    cout << "HO_Biforlm::Allocatematrix" << endl;
    const FESpace & fespace = this->fespace;

    this->mats.Append (new ElementByElementMatrix<SCAL> (fespace.GetNDof(), this->ma.GetNE()+this->ma.GetNSE()));

    // cout << "generate inexact schur matrix" << endl;
    // generate inexact schur matrix

    MatrixGraph * graph = 0;

    if (this->ma.GetDimension() == 3)
      {
        int ne = this->ma.GetNE();
        int nfa = this->ma.GetNFaces();
        Array<int> vnums, ednums, fanums, elnums, dnums;
        Array<int> cnt(nfa+ne);
        cnt = 0;

	
        for (int i = 0; i < nfa; i++)
          {
            this->ma.GetFacePNums (i, vnums);

            this->ma.GetFaceEdges (i, ednums);
	    this->ma.GetFaceElements (i, elnums);

            /*
            // only for shells
            if ( this->ma.GetElType(elnums[0]) == ET_PRISM 
            && this->ma.GetElType(elnums[1]) == ET_PRISM )
            if (vnums.Size() == 3) continue;
            */

            for (int j = 0; j < vnums.Size(); j++)
              {
                fespace.GetVertexDofNrs (vnums[j], dnums);
                for (int k = 0; k < dnums.Size(); k++)
                  if (dnums[k] != -1) cnt[i]++;
              }

            for (int j = 0; j < ednums.Size(); j++)
              {
                fespace.GetEdgeDofNrs (ednums[j], dnums);
                for (int k = 0; k < dnums.Size(); k++)
                  if (dnums[k] != -1) cnt[i]++;
              }

            fespace.GetFaceDofNrs (i, dnums);
            for (int k = 0; k < dnums.Size(); k++)
              if (dnums[k] != -1) cnt[i]++;


            for (int j = 0; j < elnums.Size(); j++)
              {
                fespace.GetWireBasketDofNrs (elnums[j], dnums);
                for (int k = 0; k < dnums.Size(); k++)
                  if (dnums[k] != -1) cnt[i]++;
              }
          }

        /*
          for (int i = 0; i < ne; i++)
          {
          fespace.GetExternalDofNrs (i, dnums);
          for (int j = 0; j < dnums.Size(); j++)
          if (dnums[j] != -1)
          cnt[nfa+i]++;
          }
        */

        Table<int> fa2dof(cnt);

        for (int i = 0; i < nfa; i++)
          {
            this->ma.GetFacePNums (i, vnums);
            this->ma.GetFaceElements (i, elnums);

            /*
            // only for shells
            if ( this->ma.GetElType(elnums[0]) == ET_PRISM 
            && this->ma.GetElType(elnums[1]) == ET_PRISM )
            if (vnums.Size() == 3) continue;
            */

            this->ma.GetFaceEdges (i, ednums);
            int ii = 0;
            for (int j = 0; j < vnums.Size(); j++)
              {
                fespace.GetVertexDofNrs (vnums[j], dnums);
                for (int k = 0; k < dnums.Size(); k++)
                  if (dnums[k] != -1) fa2dof[i][ii++] = dnums[k];
              }
            for (int j = 0; j < ednums.Size(); j++)
              {
                fespace.GetEdgeDofNrs (ednums[j], dnums);
                for (int k = 0; k < dnums.Size(); k++)
                  if (dnums[k] != -1) fa2dof[i][ii++] = dnums[k];
              }

            fespace.GetFaceDofNrs (i, dnums);
            for (int k = 0; k < dnums.Size(); k++)
              if (dnums[k] != -1) fa2dof[i][ii++] = dnums[k];


            for (int j = 0; j < elnums.Size(); j++)
              {
                fespace.GetWireBasketDofNrs (elnums[j], dnums);
                for (int k = 0; k < dnums.Size(); k++)
                  if (dnums[k] != -1) fa2dof[i][ii++] = dnums[k];
              }
          }

        graph = new MatrixGraph (fespace.GetNDof(), fa2dof, this->symmetric);
      }
    
    else
      
      {
        // int nel = this->ma.GetNE();
        int ned = this->ma.GetNEdges();
        // int nsel = this->ma.GetNSE();
        Array<int> vnums, ednums, elnums, dnums;
        Array<int> cnt(2*ned);//+nsel);
        cnt = 0;

        for (int i = 0; i < ned; i++)
          {
            int vi[2];
            this->ma.GetEdgePNums (i, vi[0], vi[1]);

            for (int j = 0; j < 2; j++)
              {
                fespace.GetVertexDofNrs (vi[j], dnums);
                for (int k = 0; k < dnums.Size(); k++)
                  if (dnums[k] != -1) cnt[i]++;
              }

            fespace.GetEdgeDofNrs (i, dnums);
            for (int k = 0; k < dnums.Size(); k++)
              if (dnums[k] != -1) cnt[i]++;

            this->ma.GetEdgeElements (i, elnums);
            for (int j = 0; j < elnums.Size(); j++)
              {
                fespace.GetWireBasketDofNrs (elnums[j], dnums);
                for (int k = 0; k < dnums.Size(); k++)
                  if (dnums[k] != -1) cnt[i]++;
              }
          }


        Table<int> fa2dof(cnt);

        for (int i = 0; i < ned; i++)
          {
            int ii = 0;
            int vi[2];
            this->ma.GetEdgePNums (i, vi[0], vi[1]);

            for (int j = 0; j < 2; j++)
              {
                fespace.GetVertexDofNrs (vi[j], dnums);
                for (int k = 0; k < dnums.Size(); k++)
                  if (dnums[k] != -1)  fa2dof[i][ii++] = dnums[k];
              }

            fespace.GetEdgeDofNrs (i, dnums);
            for (int k = 0; k < dnums.Size(); k++)
              if (dnums[k] != -1) fa2dof[i][ii++] = dnums[k];

            this->ma.GetEdgeElements (i, elnums);
            for (int j = 0; j < elnums.Size(); j++)
              {
                fespace.GetWireBasketDofNrs (elnums[j], dnums);
                for (int k = 0; k < dnums.Size(); k++)
                  if (dnums[k] != -1) fa2dof[i][ii++] = dnums[k];
              }
          }


        // (*testout) << "face graph table: " << endl << fa2dof << endl;
        /*    
              for (int i = 0; i < ne; i++)
              {
              fespace.GetExternalDofNrs (i, dnums);
              int ii = 0;
              for (int j = 0; j < dnums.Size(); j++)
              if (dnums[j] != -1)
              fa2dof[nfa+i][ii++] = dnums[j];
              }
        */
        // (*testout) << "fa2dof = " << endl << fa2dof << endl;
    

        graph = new MatrixGraph (fespace.GetNDof(), fa2dof, this->symmetric);
      }


    if (this->symmetric)
      inexact_schur = new SparseMatrixSymmetric<SCAL> (*graph, 1);
    else
      inexact_schur = new SparseMatrix<SCAL> (*graph, 1);

    delete graph;

    inexact_schur -> AsVector() = 0.0;
    
    cout << "matrix complete" << endl;
  }


  template<class SCAL>
  BaseVector * HO_BilinearForm<SCAL> :: CreateVector() const
  {
    return new VVector<SCAL> (this->fespace.GetNDof());
  }

  template<class SCAL>
  void HO_BilinearForm<SCAL> :: AddElementMatrix (const Array<int> & dnums1,
                                                  const Array<int> & dnums2,
                                                  const FlatMatrix<SCAL> & elmat,
                                                  bool inner_element, int elnr,
                                                  LocalHeap & lh)
  {
    {
      /*
        (*testout) << "inner_element = " << inner_element << endl;
        (*testout) << "elnr = " << elnr << endl;
        (*testout) << "elmat = " << endl << elmat << endl;
        (*testout) << "dnums1 = " << endl << dnums1 << endl;
      */
      
      int nr = elnr;
      if (!inner_element) nr += this->ma.GetNE();
      dynamic_cast<ElementByElementMatrix<SCAL>&> (this->GetMatrix()) . AddElementMatrix (nr, dnums1, dnums2, elmat);
    
      
      return;

      // rest needed for Wirebasket
      


      if (inner_element)
        {
          static int timer = NgProfiler::CreateTimer ("Wirebasket matrix processing");
          NgProfiler::RegionTimer reg (timer);

          // build Schur complements for facets

        

          Array<int> usedi, unusedi, wirebasketdofs;
        
#ifdef WB_EXTRA
          this->fespace.GetWireBasketDofNrs (elnr, wirebasketdofs);
          
	  
          for (int j = 0; j < dnums1.Size(); j++)
            if (dnums1[j] != -1)
              {
                bool has = wirebasketdofs.Contains(dnums1[j]);
                
                if (has)  usedi.Append (j);
                if (!has) unusedi.Append (j);
              }
          
          int su = usedi.Size();
          int sunu = unusedi.Size();
          
          FlatMatrix<SCAL> a(su, su, lh);
          FlatMatrix<SCAL> b(su, sunu, lh);
          FlatMatrix<SCAL> c(su, sunu, lh);
          FlatMatrix<SCAL> idc(su, sunu, lh);
          FlatMatrix<SCAL> d(sunu, sunu, lh);
          
          for (int k = 0; k < su; k++)
            for (int l = 0; l < su; l++)
              a(k,l) = elmat(usedi[k], usedi[l]);
          
          for (int k = 0; k < su; k++)
            for (int l = 0; l < sunu; l++)
              {
                b(k,l) = elmat(usedi[k], unusedi[l]);
                c(k,l) = elmat(unusedi[l], usedi[k]);
              }
          
          for (int k = 0; k < sunu; k++)
            for (int l = 0; l < sunu; l++)
              d(k,l) = elmat(unusedi[k], unusedi[l]);
          
#ifdef LAPACK
          if ( sunu > 0 )
            {
              LapackInverse (d);
              LapackMultABt (c, d, idc);
              LapackMultAddABt (b, idc, -1, a);
            }
#else
          FlatMatrix<SCAL> invd(sunu, sunu, lh);
          CalcInverse (d, invd);
          d = invd;
          idc = c * Trans (invd);
          a -= b * Trans (idc);
#endif


#pragma omp critical (addinexact_schur)
          if (this->symmetric)
            {
              for (int k = 0; k < usedi.Size(); k++)
                for (int l = 0; l < usedi.Size(); l++)
                  if (dnums1[usedi[k]] >= dnums1[usedi[l]])
                    (*inexact_schur)(dnums1[usedi[k]], dnums1[usedi[l]]) += a(k,l);
            }
          else
            {
              for (int k = 0; k < usedi.Size(); k++)
                for (int l = 0; l < usedi.Size(); l++)
                  (*inexact_schur)(dnums1[usedi[k]], dnums1[usedi[l]]) += a(k,l);
            }
#endif

          if (this->ma.GetDimension() == 3)
            {

              Array<int> useddofs;
              Array<int> vnums, ednums, fanums, dnums;	

              this->ma.GetElPNums (elnr, vnums);
              this->ma.GetElEdges (elnr, ednums);
              this->ma.GetElFaces (elnr, fanums);

              for (int i = 0; i < fanums.Size(); i++)
                {
                  useddofs.SetSize (0);
                  usedi.SetSize (0);
                  unusedi.SetSize (0);

                  this->ma.GetFacePNums (fanums[i], vnums);
		
                  /*
                  // only for shells
                  if ( this->ma.GetElType(elnr) == ET_PRISM )
                  if (vnums.Size() == 3) continue;
                  */

                  // this->ma.GetElPNums (elnr, vnums);   // all vertices

                  this->fespace.GetWireBasketDofNrs (elnr, useddofs);

                  this->ma.GetFaceEdges (fanums[i], ednums);
	    
                  for (int j = 0; j < vnums.Size(); j++)
                    {
                      this->fespace.GetVertexDofNrs (vnums[j], dnums);
                      for (int k = 0; k < dnums.Size(); k++)
                        if (dnums[k] != -1) useddofs.Append (dnums[k]);
                    }
                  for (int j = 0; j < ednums.Size(); j++)
                    {
                      this->fespace.GetEdgeDofNrs (ednums[j], dnums);
                      for (int k = 0; k < dnums.Size(); k++)
                        if (dnums[k] != -1) useddofs.Append (dnums[k]);
                    }
                  this->fespace.GetFaceDofNrs (fanums[i], dnums);
                  for (int k = 0; k < dnums.Size(); k++)
                    if (dnums[k] != -1) useddofs.Append (dnums[k]);


                  for (int j = 0; j < dnums1.Size(); j++)
                    if (dnums1[j] != -1)
                      {
                        bool has = useddofs.Contains (dnums1[j]);
                        /*
                          bool wb = wirebasketdofs.Contains (dnums1[j]);
                          if (has && !wb) usedi.Append (j);
                          if (!has && !wb) unusedi.Append (j);
                        */
                        if (has) usedi.Append (j);
                        if (!has) unusedi.Append (j);
                      }

                  int su = usedi.Size();
                  int sunu = unusedi.Size();


                  FlatMatrix<SCAL> a(su, su, lh);
                  FlatMatrix<SCAL> b(su, sunu, lh);
                  FlatMatrix<SCAL> c(su, sunu, lh);
                  FlatMatrix<SCAL> d(sunu, sunu, lh);
			    
                  for (int k = 0; k < su; k++)
                    for (int l = 0; l < su; l++)
                      a(k,l) = elmat(usedi[k], usedi[l]);
			    
                  for (int k = 0; k < su; k++)
                    for (int l = 0; l < sunu; l++)
                      {
                        b(k,l) = elmat(usedi[k], unusedi[l]);
                        c(k,l) = elmat(unusedi[l], usedi[k]);
                      }
	    
                  for (int k = 0; k < sunu; k++)
                    for (int l = 0; l < sunu; l++)
                      d(k,l) = elmat(unusedi[k], unusedi[l]);

#ifdef LAPACK
                  if ( sunu > 0 )
                    {
                      LapackAInvBt (d, b);
                      LapackMultAddABt (b, c, -1, a);
                    }
#else
                  FlatMatrix<SCAL> invd(sunu, sunu, lh);
                  FlatMatrix<SCAL> idc(su, sunu, lh);
                  CalcInverse (d, invd);
                  d = invd;
                  idc = c * Trans (invd);
                  a -= b * Trans (idc);
#endif
	
#pragma omp critical (addinexact_schur)
                  {
                    if (this->symmetric)
                      {
                        for (int k = 0; k < usedi.Size(); k++)
                          for (int l = 0; l < usedi.Size(); l++)
                            if (dnums1[usedi[k]] >= dnums1[usedi[l]])
                              (*inexact_schur)(dnums1[usedi[k]], dnums1[usedi[l]]) += a(k,l);
                      }
                    else
                      {
                        for (int k = 0; k < usedi.Size(); k++)
                          for (int l = 0; l < usedi.Size(); l++)
                            (*inexact_schur)(dnums1[usedi[k]], dnums1[usedi[l]]) += a(k,l);
                      }
                  }
                }
            }
          else
            {
              Array<int> useddofs;
              Array<int> dnums, ednums;

              this->ma.GetElEdges (elnr, ednums);
              Matrix<SCAL> helmat(elmat.Height());
              helmat = elmat;

                
              for (int i = 0; i < ednums.Size(); i++)
                {
                  useddofs.SetSize (0);
                  usedi.SetSize (0);
                  unusedi.SetSize (0);

                  this->fespace.GetWireBasketDofNrs (elnr, useddofs);
                    
                  int vi[2];
                  this->ma.GetEdgePNums (ednums[i], vi[0], vi[1]);

                  for (int j = 0; j < 2; j++)
                    {
                      this->fespace.GetVertexDofNrs (vi[j], dnums);
                      for (int k = 0; k < dnums.Size(); k++)
                        if (dnums[k] != -1) useddofs.Append (dnums[k]);
                    }
                    
                  this->fespace.GetEdgeDofNrs (ednums[i], dnums);
                  for (int k = 0; k < dnums.Size(); k++)
                    if (dnums[k] != -1) useddofs.Append (dnums[k]);
                    
                  for (int j = 0; j < dnums1.Size(); j++)
                    if (dnums1[j] != -1)
                      {
                        bool has = useddofs.Contains(dnums1[j]);
                        /*
                          bool wb = wirebasketdofs.Contains(dnums1[j]);
                          if (has && !wb) usedi.Append (j);
                          if (!has && !wb) unusedi.Append (j);
                        */

                        if (has) usedi.Append (j);
                        if (!has) unusedi.Append (j);
                      }
                    
                  int su = usedi.Size();
                  int sunu = unusedi.Size();
                    
                  FlatMatrix<SCAL> a(su, su, lh);
                  FlatMatrix<SCAL> b(su, sunu, lh);
                  FlatMatrix<SCAL> c(su, sunu, lh);
                  FlatMatrix<SCAL> idc(su, sunu, lh);
                  FlatMatrix<SCAL> d(sunu, sunu, lh);
                    
                  for (int k = 0; k < su; k++)
                    for (int l = 0; l < su; l++)
                      a(k,l) = helmat(usedi[k], usedi[l]);
                    
                  for (int k = 0; k < su; k++)
                    for (int l = 0; l < sunu; l++)
                      {
                        b(k,l) = helmat(usedi[k], unusedi[l]);
                        c(k,l) = helmat(unusedi[l], usedi[k]);
                      }
                    
                  for (int k = 0; k < sunu; k++)
                    for (int l = 0; l < sunu; l++)
                      d(k,l) = helmat(unusedi[k], unusedi[l]);
                    

                  /*
                    (*testout) << "usedi = " << endl << usedi << endl;
                    (*testout) << "unusedi = " << endl << unusedi << endl;
                    (*testout) << "a = " << endl << a << endl;
                    (*testout) << "b = " << endl << b << endl;
                    (*testout) << "c = " << endl << c << endl;
                    (*testout) << "d = " << endl << d << endl;
                  */

                    
#ifdef LAPACK
                  if ( sunu > 0 )
                    {
                      LapackInverse (d);
                      LapackMultABt (c, d, idc);
                      LapackMultAddABt (b, idc, -1, a);
                    }
#else
                  FlatMatrix<SCAL> invd(sunu, sunu, lh);
                  CalcInverse (d, invd);
                  d = invd;
                  idc = c * Trans (invd);
                  a -= b * Trans (idc);
#endif

                  // *testout << "schur = " << endl << a << endl;

                  /*
                    for (int k = 0; k < su; k++)
                    for (int l = 0; l < su; l++)
                    helmat(usedi[k], usedi[l]) -= 0.99999 * a(k,l);
                  */


#pragma omp critical (addinexact_schur)
                  if (this->symmetric)
                    {
                      for (int k = 0; k < usedi.Size(); k++)
                        for (int l = 0; l < usedi.Size(); l++)
                          if (dnums1[usedi[k]] >= dnums1[usedi[l]])
                            (*inexact_schur)(dnums1[usedi[k]], dnums1[usedi[l]]) += a(k,l);
                    }
                  else
                    {
                      for (int k = 0; k < usedi.Size(); k++)
                        for (int l = 0; l < usedi.Size(); l++)
                          (*inexact_schur)(dnums1[usedi[k]], dnums1[usedi[l]]) += a(k,l);
                    }
                }
              // *testout << "helmat = " << endl << helmat << endl;
            }
        }
      else
        {

#pragma omp critical (addinexact_schur)
          {
	  
            if (this -> symmetric)
              {
                for (int i = 0; i < dnums1.Size(); i++)
                  for (int j = 0; j < dnums2.Size(); j++)
                    if (dnums1[i] != -1 && dnums2[j] != -1)
                      if (dnums1[i] >= dnums1[j])
                        (*inexact_schur)(dnums1[i], dnums2[j]) += elmat(i, j);
              }
            else
              {
                for (int i = 0; i < dnums1.Size(); i++)
                  for (int j = 0; j < dnums2.Size(); j++)
                    if (dnums1[i] != -1 && dnums2[j] != -1)
                      (*inexact_schur)(dnums1[i], dnums2[j]) += elmat(i, j);
              }
          }
        }
    }
  }





  template class HO_BilinearForm<double>;
  template class HO_BilinearForm<Complex>;





  template <class SCAL>
  WireBasketPreconditioner<SCAL> ::
  WireBasketPreconditioner (const PDE * pde, const Flags & aflags, const string aname)
    : Preconditioner (pde, aflags, aname)
  {
    cout << "\ncreate wire-basket precond" << endl;
    bfa = dynamic_cast<const HO_BilinearForm<SCAL>*>(pde->GetBilinearForm (aflags.GetStringFlag ("bilinearform", NULL)));
    inversetype = flags.GetStringFlag("inverse", "sparsecholesky");
    smoothingtype = int(flags.GetNumFlag("blocktype", -1));
  }


  template <class SCAL>
  void WireBasketPreconditioner<SCAL> ::
  Update ()
  {
    cout << "update wirebasket" << endl;

    const SparseMatrix<SCAL> & mat = bfa -> InexactSchur();

    //      const SparseMatrix<SCAL> & mat = 
    //	dynamic_cast<const SparseMatrix<SCAL>& > (bfa -> GetMatrix());

    //       (*testout) << "type = " << typeid(mat).name() << endl;
    //       (*testout) << "mat = " << mat << endl;
      
    mat.SetInverseType(inversetype);
    cout << "inversetype = " << inversetype << endl;
    if ( smoothingtype == -1 )
      {
        pre = mat.InverseMatrix(); // const_cast<SparseMatrix<SCAL>*> (&mat); // 
      }
    else
      {
        Array<int> cnt(ma.GetNE());
        // const ElementByElementMatrix<SCAL> & elbyelmat =
        // dynamic_cast<const ElementByElementMatrix<SCAL>&> (bfa -> GetMatrix() );
        // 	  const helmholtz_exp_cpp::HybridHelmholtzFESpace & hhspace = 
        // 	    dynamic_cast<const helmholtz_exp_cpp::HybridHelmholtzFESpace & > 
        // 	    (bfa -> GetFESpace() );

        Flags flags;
        flags.SetFlag ("blocktype", smoothingtype);
	if (bfa->UsesEliminateInternal())
	  flags.SetFlag("eliminate_internal");	
        Table<int> * blocktable = bfa->GetFESpace().CreateSmoothingBlocks(flags);
        pre = mat.  CreateBlockJacobiPrecond (*blocktable);
      }
    if (test) Test();
  }


  template class WireBasketPreconditioner<double>;
  template class WireBasketPreconditioner<Complex>;



  class BDDCMatrix : public BaseMatrix
  {
    const HO_BilinearForm<double> & bfa;
    Array<int> restrict;
    Array<int> multiple;
    BaseMatrix * inv;
    string inversetype;
    BitArray * free_dofs;
    const MeshAccess & ma;
  public:
    BDDCMatrix (const HO_BilinearForm<double> & abfa, const string & inversetype)
      : ma(abfa.GetFESpace().GetMeshAccess()),bfa(abfa) 
    {
      // LocalHeap lh(10000000);
      const FESpace & fes = bfa.GetFESpace();
      const MeshAccess & ma = fes.GetMeshAccess();

      int ne = ma.GetNE();
      Array<int> cnt(ne);

      Array<int> wbdofs(fes.GetNDof()), lwbdofs, dnums;
      wbdofs = 0;

      for (int i = 0; i < ne; i++)
        {
          fes.GetExternalDofNrs (i, dnums);
          cnt[i] = dnums.Size();

          fes.GetWireBasketDofNrs (i, lwbdofs);
          for (int j = 0; j < lwbdofs.Size(); j++)
            wbdofs[lwbdofs[j]] = 1;
	  *testout << "dnums = " << endl << dnums << endl;
	  *testout << "lwbdofs = " << endl << lwbdofs << endl;
	  
        }

      Table<int> dcdofs(cnt);   // discontinuous dofs

      

      int firstdcdof = fes.GetNDof();
      // for (int i = 0; i < wbdofs.Size(); i++)
      // if (wbdofs[i]) firstdcdof = i+1;
      int ndcdof = fes.GetNDof();

      for (int i = 0; i < ne; i++)
        {
          fes.GetExternalDofNrs (i, dnums);

          for (int j = 0; j < dnums.Size(); j++)
            {
              if (dnums[j] == -1) continue;
              if (wbdofs[dnums[j]]) continue;
              dnums[j] = ndcdof;
              firstdcdof++;
	      ndcdof++;
            }

          for (int j = 0; j < dnums.Size(); j++)
            dcdofs[i][j] = dnums[j];
        }

      // *testout << "dcdofs = " << endl << dcdofs << endl;

      restrict.SetSize(ndcdof);
      restrict = -1;
      for (int i = 0; i < ne; i++)
        {
          fes.GetExternalDofNrs (i, dnums);

          for (int j = 0; j < dnums.Size(); j++)
            {
              if (dnums[j] != -1)
                restrict[dcdofs[i][j]] = dnums[j];
            }
        }

      // *testout << "restrict = " << endl << restrict << endl;

      multiple.SetSize (fes.GetNDof());
      multiple = 0;
      for (int i = 0; i < restrict.Size(); i++)
        if (restrict[i] != -1)
          multiple[restrict[i]]++;


      cout << "now build graph" << endl;

      MatrixGraph graph(firstdcdof, dcdofs, bfa.IsSymmetric());

      cout << "now allocate matrix" << endl;

      SparseMatrix<double> * pdcmat;
      if (bfa.IsSymmetric()){
	pdcmat = new SparseMatrixSymmetric<double>(graph,1);
	cout << "symmetric" << endl;
      }
      else{
	pdcmat = new SparseMatrix<double>(graph,1);
	cout << "nonsymmetric" << endl;
      }
      SparseMatrix<double>& dcmat=*pdcmat;

      dcmat.SetInverseType (inversetype);
      
      cout << "have matrix" << endl;

      for (int i = 0; i < ne; i++)
        {
          FlatMatrix<> elmat = 
            dynamic_cast<const ElementByElementMatrix<double>&> (bfa.GetMatrix()) . GetElementMatrix (i);

          FlatArray<int> dofs = dcdofs[i];
          for (int k = 0; k < dofs.Size(); k++)
            for (int l = 0; l < dofs.Size(); l++)
              if (dofs[k] != -1 && dofs[l] != -1 && (!bfa.IsSymmetric() || dofs[k] >= dofs[l]))
                dcmat(dofs[k], dofs[l]) += elmat(k,l);
        }
      
      cout << "matrix filled" << endl;

      for (int i = 0; i < restrict.Size(); i++)
	if (restrict[i] == -1)
	  dcmat(i,i) = 1;

      // *testout << "dcmat = " << endl << dcmat << endl;
      *testout << "freedofs = " << *fes.GetFreeDofs() << endl;


      free_dofs = new BitArray (ndcdof);
      free_dofs->Clear();
      
      for (int i = 0; i < ndcdof; i++)
	{
	  if (restrict[i] != -1)
	    if (fes.GetFreeDofs()->Test(restrict[i]))
	      free_dofs->Set(i);
	}
      

      // *testout << "dcmat = " << endl << dcmat << endl;
      // *testout << "restrict = " << endl << restrict << endl;
      cout << "call inverse" << endl;
      inv = dcmat.InverseMatrix(free_dofs);
      cout << "has inverse" << endl;
      *testout << "inverse = " << (*inv) << endl;
      delete pdcmat;
    }

    ~BDDCMatrix(){delete inv;}
    
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      FlatVector<> fx = x.FVDouble();
      FlatVector<> fy = y.FVDouble();

      VVector<> lx(restrict.Size());
      VVector<> ly(restrict.Size());

      for (int i = 0; i < restrict.Size(); i++)
        if (restrict[i] != -1)
          lx(i) = fx(restrict[i]) / multiple[restrict[i]];

      ly = (*inv) * lx;
            
      fy = 0.0;
      for (int i = 0; i < restrict.Size(); i++)
        if (restrict[i] != -1)
          fy(restrict[i]) += ly(i) / multiple[restrict[i]];
    }

  };











  class BDDCMatrixRefElement : public BaseMatrix
  {
    const S_BilinearForm<double> & bfa;
    Array<int> restrict;
    Array<int> multiple;
    BaseMatrix * inv;
    
    Table<int> *wbdofs, *dcdofs, *internaldofs, *externaldofs;
    Table<int> *globwbdofs, *globdcdofs, *globinternaldofs, *globexternaldofs;

    Array<int> elclassnr;
    Array<Matrix<>*> invdc_ref, extwb_ref, schurwb_ref;
    Array<Matrix<>*> inv_int_ref, extension_int_ref, extension_int_ref_trans;
    int ne, ndof;
    
    string inversetype;
    bool print;

  public:
    BDDCMatrixRefElement (const S_BilinearForm<double> & abfa, const string & ainversetype)
      : bfa(abfa), inversetype (ainversetype)
    {
      cout << "BDDC MatrixRefElement" << endl;
      print = true;
      LocalHeap lh(50000000);

      if (print)
	*testout << "BDDC - Constructor" << endl;

      const FESpace & fes = bfa.GetFESpace();
      const MeshAccess & ma = fes.GetMeshAccess();

      ne = ma.GetNE();
      ndof = fes.GetNDof();
      Array<int> cnt(ne), cntwb(ne), cntdc(ne), cntinternal(ne), cntexternal(ne);

      Array<int> gwbdofs(ndof);
      Array<int> wbdnums, dnums, extdnums, dcdnums, vnums;
      gwbdofs = 0;

      for (int i = 0; i < ne; i++)
        {
          fes.GetDofNrs (i, dnums);
          fes.GetExternalDofNrs (i, extdnums);
          fes.GetWireBasketDofNrs (i, wbdnums);

          cnt[i] = dnums.Size();
          cntwb[i] = wbdnums.Size();
          cntdc[i] = extdnums.Size()-wbdnums.Size();
          cntinternal[i] = dnums.Size()-extdnums.Size();
          cntexternal[i] = extdnums.Size();

          for (int j = 0; j < wbdnums.Size(); j++)
            gwbdofs[wbdnums[j]] = 1;
        }

      if (print)
	{
	  *testout << "Glocal wirebasket - dofs:" << endl;
	  *testout << "ndof = " << ndof << endl;
	  *testout << gwbdofs << endl;
	}

      // Table<int> dcdofnrs(cnt);   // discontinuous dofs
      Table<int> dcdofnrs(cntexternal);   // discontinuous dofs

      wbdofs = new Table<int>(cntwb);
      dcdofs = new Table<int>(cntdc);
      internaldofs = new Table<int>(cntinternal);
      externaldofs = new Table<int>(cntexternal);

      globwbdofs = new Table<int>(cntwb);
      globdcdofs = new Table<int>(cntdc);
      globinternaldofs = new Table<int>(cntinternal);
      globexternaldofs = new Table<int>(cntexternal);      

      int totdofs  = fes.GetNDof();
      for (int i = 0; i < ne; i++)
        {
	  *testout << "element " << i << endl;

          fes.GetDofNrs (i, dnums);
          fes.GetExternalDofNrs (i, extdnums);

	  *testout << "dnums = " << endl << dnums << endl;
	  *testout << "ext dnums = " << endl << extdnums << endl;

	  int nint = 0, next = 0;

	  for (int j = 0; j < dnums.Size(); j++)
	    {
	      if (extdnums.Contains (dnums[j]))
		{
		  (*externaldofs)[i][next] = j;
		  (*globexternaldofs)[i][next] = dnums[j];
		  next++;
		}
	      else
		{
		  (*internaldofs)[i][nint] = j;
		  (*globinternaldofs)[i][nint] = dnums[j];
		  nint++;
		}
	    }
	  
          dcdnums.SetSize (extdnums.Size());

          for (int j = 0; j < extdnums.Size(); j++)
            {
	      if (extdnums[j] == -1)
		dcdnums[j] = -1;
	      else
		{
		  if (gwbdofs[extdnums[j]])
		    dcdnums[j] = extdnums[j];
		  else
		    {
		      dcdnums[j] = totdofs;
		      totdofs++;
		    }
		}
	      /*
              dcdnums[j] = extdnums[j];
              if (dcdnums[j] == -1) continue;
              if (gwbdofs[dcdnums[j]]) continue;
              dcdnums[j] = totdofs;
              totdofs++;
	      */
            }

	  *testout << "dcdnums = " << endl << dcdnums << endl;

          int wbi = 0, dci = 0;
          for (int j = 0; j < extdnums.Size(); j++)
            {
              if (extdnums[j] == -1) continue;
              if (gwbdofs[extdnums[j]])
                {
                  (*wbdofs)[i][wbi] = j;
                  (*globwbdofs)[i][wbi] = extdnums[j];
                  wbi++;
                }
              else
                {
                  (*dcdofs)[i][dci] = j;
                  (*globdcdofs)[i][dci] = dcdnums[j];
                  dci++;
                }
            }
          
          for (int j = 0; j < extdnums.Size(); j++)
            dcdofnrs[i][j] = dcdnums[j];
        }


      if (print)
	{
	  *testout << "dcdofnrs = " << endl;
	  *testout << dcdofnrs << endl;
	}

      int firstdcdof = fes.GetNDof();      
      
      MatrixGraph graph(firstdcdof, *globwbdofs, 1);
      SparseMatrixSymmetric<double> dcmat(graph, 1);
      dcmat = 0.0;
      // dcmat.SetInverseType ("sparsecholesky");
      dcmat.SetInverseType (inversetype);
      *testout << "dcmat, init = " << dcmat << endl;

      elclassnr.SetSize(ne);
      int maxclass = 32;

      invdc_ref.SetSize(maxclass);
      extwb_ref.SetSize(maxclass);
      schurwb_ref.SetSize(maxclass);
      
      inv_int_ref.SetSize(maxclass);
      extension_int_ref.SetSize(maxclass);
      extension_int_ref_trans.SetSize(maxclass);

      invdc_ref = NULL;
      extwb_ref = NULL;
      schurwb_ref = NULL;
      inv_int_ref = NULL;
      extension_int_ref = NULL;

      for (int i = 0; i < ne; i++)
        {
          HeapReset hr(lh);
          ma.GetElVertices(i, vnums);

          int classnr = 0;

          int sort[4] = { 0, 1, 2, 3 };
          if (vnums.Size() == 4)
            {
              if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 1; }
              if (vnums[sort[2]] > vnums[sort[3]]) { Swap (sort[2], sort[3]); classnr += 2; }
              if (vnums[sort[0]] > vnums[sort[2]]) { Swap (sort[0], sort[2]); classnr += 4; }
              if (vnums[sort[1]] > vnums[sort[3]]) { Swap (sort[1], sort[3]); classnr += 8; }
              if (vnums[sort[1]] > vnums[sort[2]]) { Swap (sort[1], sort[2]); classnr += 16; }
            }
          else
            {
              if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 1; }
              if (vnums[sort[1]] > vnums[sort[2]]) { Swap (sort[1], sort[2]); classnr += 2; }
              if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 4; }
            }
	  //  classnr = i;

          elclassnr[i] = classnr;
          if (invdc_ref[classnr]) continue;

          cout << "assemble class " << classnr << endl;

	  if (print)
	    *testout << "assemble class " << classnr << endl;

          fes.GetDofNrs (i, dnums);
          fes.GetExternalDofNrs (i, extdnums);
	  
          Matrix<> elmat(dnums.Size());
          Matrix<> partelmat(dnums.Size());
          elmat = 0.0;

          ElementTransformation eltrans;
          ma.GetElementTransformation (i, eltrans, lh);
          const FiniteElement & fel = fes.GetFE(i, lh);

          for (int j = 0; j < bfa.NumIntegrators(); j++)
            {
	      if (bfa.GetIntegrator(j) -> BoundaryForm()) continue;
              bfa.GetIntegrator(j)->AssembleElementMatrix (fel, eltrans, partelmat, lh);
              elmat += partelmat;
            }

	  

	  if (print)
	    *testout << "elmat = " << endl << elmat << endl;
	  
          int nwb = (*wbdofs)[i].Size();
          int ndc = (*dcdofs)[i].Size();
          int nint = (*internaldofs)[i].Size();
          int next = (*externaldofs)[i].Size();

	  
	  // external/internal extension


	  *testout << "nint = " << nint << ", next = " << next << ", nwb = " << nwb << ", ndc = " << ndc << endl;
	  *testout << "wbdofs = " << (*wbdofs)[i] << endl;

          Matrix<> matee(next);
          Matrix<> matie(nint, next);
          Matrix<> extie(nint, next);
          Matrix<> matii(nint);

          for (int j = 0; j < next; j++)
            for (int k = 0; k < next; k++)
              matee(j,k) = elmat((*externaldofs)[i][j], (*externaldofs)[i][k]);
          for (int j = 0; j < nint; j++)
            for (int k = 0; k < next; k++)
              matie(j,k) = elmat((*internaldofs)[i][j], (*externaldofs)[i][k]);
          for (int j = 0; j < nint; j++)
            for (int k = 0; k < nint; k++)
              matii(j,k) = elmat((*internaldofs)[i][j], (*internaldofs)[i][k]);


	  *testout << "matee = " << endl << matee << endl;
	  *testout << "matie = " << endl << matie << endl;
	  *testout << "matii = " << endl << matii << endl;

          LapackInverse (matii);
          extie = matii * matie;
          matee -= Trans (matie) * extie;
	  

	  *testout << "schur-ext = " << endl << matee << endl;

          Matrix<> matwbwb(nwb);
          Matrix<> matdcwb(ndc, nwb);
          Matrix<> extdcwb(ndc, nwb);
          Matrix<> matdcdc(ndc);

          for (int j = 0; j < nwb; j++)
            for (int k = 0; k < nwb; k++)
              matwbwb(j,k) = matee((*wbdofs)[i][j], (*wbdofs)[i][k]);
          for (int j = 0; j < ndc; j++)
            for (int k = 0; k < nwb; k++)
              matdcwb(j,k) = matee((*dcdofs)[i][j], (*wbdofs)[i][k]);
          for (int j = 0; j < ndc; j++)
            for (int k = 0; k < ndc; k++)
              matdcdc(j,k) = matee((*dcdofs)[i][j], (*dcdofs)[i][k]);

          LapackInverse (matdcdc);
          extdcwb = matdcdc * matdcwb;
          matwbwb -= Trans (matdcwb) * extdcwb;

	  *testout << "schur-wb = " << endl << matwbwb << endl;
          
	  invdc_ref[classnr] = new Matrix<> (ndc, ndc);
	  extwb_ref[classnr] = new Matrix<> (ndc, nwb);
	  schurwb_ref[classnr] = new Matrix<> (nwb, nwb);
	  
	  *(invdc_ref[classnr]) = matdcdc;
	  *(extwb_ref[classnr]) = extdcwb;
	  *(schurwb_ref[classnr]) = matwbwb;

	  inv_int_ref[classnr] = new Matrix<> (nint, nint);
	  *(inv_int_ref[classnr]) = matii;
	  extension_int_ref[classnr] = new Matrix<> (nint, next);
	  *(extension_int_ref[classnr]) = extie;

	  extension_int_ref_trans[classnr] = new Matrix<> (next, nint);
	  *(extension_int_ref_trans[classnr]) = Trans(extie);
        }


      for (int i = 0; i < ne; i++)
        {
          FlatArray<int> dofs = (*globwbdofs)[i]; 
          FlatMatrix<> matwbwb = *schurwb_ref[elclassnr[i]];
          for (int k = 0; k < dofs.Size(); k++)
            for (int l = 0; l < dofs.Size(); l++)
              if (dofs[k] != -1 && dofs[l] != -1 && dofs[k] >= dofs[l])
                dcmat(dofs[k], dofs[l]) += matwbwb(k,l);
        }

      if (print)
	*testout << "disconnected matrix: " << endl 
		 << dcmat << endl;


      restrict.SetSize(totdofs);
      restrict = -1;
      for (int i = 0; i < ne; i++)
        {
          // fes.GetDofNrs (i, dnums);
          fes.GetExternalDofNrs (i, extdnums);

          for (int j = 0; j < extdnums.Size(); j++)
            {
              if (extdnums[j] != -1)
                restrict[dcdofnrs[i][j]] = extdnums[j];
            }
        }

      multiple.SetSize (fes.GetNDof());
      multiple = 0;
      for (int i = 0; i < restrict.Size(); i++)
        if (restrict[i] != -1)
          multiple[restrict[i]]++;

      for (int i = 0; i < firstdcdof; i++)
        if (gwbdofs[i] == 0)
          dcmat(i,i) = 1;


      const BitArray & free_dofs = *fes.GetFreeDofs();
      for (int i = 0; i < restrict.Size(); i++)
	if (restrict[i] != -1 && !free_dofs.Test(restrict[i]))
	  restrict[i] = -1;


      cout << "call inverse, type = " << inversetype << endl;
      cout << "dim = " << dcmat.Height() << endl;

      
      inv = dcmat.InverseMatrix(&free_dofs);

#define sometests
#ifdef sometests
      Array<int> & directblocks = *new Array<int> (free_dofs.Size());
      directblocks = 0;
      // int nv = ma.GetNV();
      int ned = ma.GetNEdges();
      // int nfa = ma.GetNFaces();
      for (int i = 0; i < directblocks.Size(); i++)
	if (free_dofs.Test(i))
	  directblocks[i] = 1;

      for (int i = 0; i < ned; i++)
	{
	  fes.GetEdgeDofNrs (i, dnums);
	  for (int j = 1; j < dnums.Size(); j++)
	    if (free_dofs.Test(dnums[j]))
	      directblocks[dnums[j]] = 2+i;
	}


      /*
      for (int i = 0; i < nfa; i++)
	{
	  fes.GetFaceDofNrs (i, dnums);
	  for (int j = 0; j < dnums.Size(); j++)
	    if (free_dofs.Test(dnums[j]))
	      directblocks[dnums[j]] = 1;
	}
      */

      *testout << "directblocks = " << endl << directblocks << endl;


      inv = dcmat.InverseMatrix(&directblocks);
      // inv = dcmat.InverseMatrix(&free_dofs);
      cout << "complete" << endl;

      EigenSystem eigen(dcmat, *inv);
      eigen.Calc();
      eigen.PrintEigenValues (cout);
      if (print)
	{
	  *testout << "restrict = " << endl  << restrict << endl;
	  *testout << "dcmat = " << endl << dcmat << endl;
	  *testout << "freedofs = " << endl << *fes.GetFreeDofs() << endl;
	}
#endif
    }



    /*
      virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
      {
      static int timer = NgProfiler::CreateTimer ("Apply BDDC preconditioner");
      static int timerrest = NgProfiler::CreateTimer ("Apply BDDC preconditioner, restrict");
      static int timerdc = NgProfiler::CreateTimer ("Apply BDDC preconditioner, decoupled");
      static int timerext = NgProfiler::CreateTimer ("Apply BDDC preconditioner, extend");
      NgProfiler::RegionTimer reg (timer);

      LocalHeap lh(100000);
      FlatVector<> fx = x.FVDouble();
      FlatVector<> fy = y.FVDouble();

      VVector<> lx(restrict.Size());
      VVector<> ly(restrict.Size());

      lx = 0.0;
      for (int i = 0; i < restrict.Size(); i++)
      if (restrict[i] != -1)
      lx(i) = fx(restrict[i]) / multiple[restrict[i]];


      // restrict
      NgProfiler::StartTimer (timerrest);
      #pragma omp parallel
      {
      LocalHeap lh(100000);
      #pragma omp for
      for (int i = 0; i < ne; i++)
      {
      HeapReset hr(lh);
      FlatArray<int> wbd = (*globwbdofs)[i];
      FlatArray<int> dcd = (*globdcdofs)[i];

      FlatVector<> wb(wbd.Size(), lh);
      FlatVector<> dc(dcd.Size(), lh);

      for (int j = 0; j < dcd.Size(); j++)
      dc(j) = lx(dcd[j]);
          
      // wb = Trans(*extwb[i]) * dc;
      wb = Trans(*extwb_ref[elclassnr[i]]) * dc;
          
      for (int j = 0; j < wbd.Size(); j++)
      lx(wbd[j]) -= wb(j);

      NgProfiler::AddFlops (timer, 2*wbd.Size()*dcd.Size());
      NgProfiler::AddFlops (timerrest, wbd.Size()*dcd.Size());
      }
      }
      NgProfiler::StopTimer (timerrest);

      ly = 0.0;
      ly = (*inv) * lx;


      // solve decoupled
      NgProfiler::StartTimer (timerdc);
      #pragma omp parallel
      {
      LocalHeap lh(100000);
      #pragma omp for
      for (int i = 0; i < ne; i++)
      {
      HeapReset hr(lh);
      FlatArray<int> dcd = (*globdcdofs)[i];
            
      FlatVector<> dc1(dcd.Size(), lh);
      FlatVector<> dc2(dcd.Size(), lh);
            
      for (int j = 0; j < dcd.Size(); j++)
      dc1(j) = lx(dcd[j]);
            
      // dc2 = Trans(*invdc[i]) * dc1;
      dc2 = (*invdc_ref[elclassnr[i]]) * dc1;
            
      for (int j = 0; j < dcd.Size(); j++)
      ly(dcd[j]) = dc2(j);
            
      NgProfiler::AddFlops (timer, sqr(dcd.Size()));
      NgProfiler::AddFlops (timerdc, sqr(dcd.Size()));
      }
      }
      NgProfiler::StopTimer (timerdc);

      // extend
      NgProfiler::StartTimer (timerext);
      #pragma omp parallel
      {
      LocalHeap lh(100000);
      #pragma omp for
      for (int i = 0; i < ne; i++)
      {
      HeapReset hr(lh);
      FlatArray<int> wbd = (*globwbdofs)[i];
      FlatArray<int> dcd = (*globdcdofs)[i];

      FlatVector<> wb(wbd.Size(), lh);
      FlatVector<> dc(dcd.Size(), lh);
      for (int j = 0; j < wbd.Size(); j++)
      wb(j) = ly(wbd[j]);

      // dc = (*extwb[i]) * wb;
      dc = (*extwb_ref[elclassnr[i]]) * wb;

      for (int j = 0; j < dcd.Size(); j++)
      ly(dcd[j]) -= dc(j);
      NgProfiler::AddFlops (timerext, wbd.Size()*dcd.Size());
      }
      }
      NgProfiler::StopTimer (timerext);

      fy = 0.0;
      for (int i = 0; i < restrict.Size(); i++)
      if (restrict[i] != -1)
      fy(restrict[i]) += ly(i) / multiple[restrict[i]];
      }
    */



    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      static int timer = NgProfiler::CreateTimer ("Apply BDDC preconditioner");
      static int timertransx = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transx");
      static int timertransy = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transy");
      static int timerrest = NgProfiler::CreateTimer ("Apply BDDC preconditioner, restrict");
      static int timerdc = NgProfiler::CreateTimer ("Apply BDDC preconditioner, decoupled");
      static int timerext = NgProfiler::CreateTimer ("Apply BDDC preconditioner, extend");
      static int timerinner = NgProfiler::CreateTimer ("Apply BDDC preconditioner, inner");
      static int timersolve = NgProfiler::CreateTimer ("Apply BDDC preconditioner, solve");
      NgProfiler::RegionTimer reg (timer);



      LocalHeap lh(100000);
      FlatVector<> fx = x.FVDouble();
      FlatVector<> fy = y.FVDouble();

      VVector<> transx(ndof);
      VVector<> transy(ndof);

      VVector<> lx(restrict.Size());
      VVector<> ly(restrict.Size());
      
      {
	NgProfiler::RegionTimer reg(timertransx);
	transx = 1.0 * x;

	/*

#pragma omp parallel
	{
	  LocalHeap slh = lh.Split();
#pragma omp for
	  for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	  for (int i = 0; i < ne; i++)
	    {
	      if (elclassnr[i] != cl) continue;

	      HeapReset hr(slh);
	      FlatArray<int> dint = (*globinternaldofs)[i];
	      FlatArray<int> dext = (*globexternaldofs)[i];
	      
	      FlatVector<> vi(dint.Size(), slh);
	      FlatVector<> ve(dext.Size(), slh);
	      
	      for (int j = 0; j < dint.Size(); j++)
		vi(j) = transx(dint[j]);
	      
	      ve = (*extension_int_ref_trans[elclassnr[i]]) * vi;

#pragma omp critical (bddc_transx)
	      {
		for (int j = 0; j < ve.Size(); j++)
		  transx(dext[j]) -= ve(j);
		NgProfiler::AddFlops (timertransx, dint.Size()*dext.Size());
	      }
	    }
	}
	*/

	for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	  {
	    int cnt = 0;
	    for (int i = 0; i < ne; i++)
	      if (elclassnr[i] == cl) cnt++;
	    if (cnt == 0) continue;

	    const Matrix<> & ext = (*extension_int_ref_trans[cl]);
	    Matrix<> vim(cnt, ext.Width());
	    Matrix<> vem(cnt, ext.Height());

	    NgProfiler::AddFlops (timertransx, ext.Width()*ext.Height()*cnt);
		
	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dint = (*globinternaldofs)[i];
		
		for (int j = 0; j < dint.Size(); j++)
		  vim(cnt, j) = transx(dint[j]);
		cnt++;
	      }

	    // vem = vim * Trans(ext);
	    LapackMultABt (vim, ext, vem);

	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dext = (*globexternaldofs)[i];

		for (int j = 0; j < dext.Size(); j++)
		  transx(dext[j]) -= vem(cnt,j);
		cnt++;
	      }
	  }
      }

      lx = 0.0;
      for (int i = 0; i < restrict.Size(); i++)
        if (restrict[i] != -1)
          lx(i) = transx(restrict[i]) / multiple[restrict[i]];


      // restrict
      {
	NgProfiler::RegionTimer reg(timerrest);

	// #pragma omp parallel
	{
	  LocalHeap slh = lh.Split();
	  // #pragma omp for
	  for (int i = 0; i < ne; i++)
	    {
	      HeapReset hr(slh);
	      FlatArray<int> wbd = (*globwbdofs)[i];
	      FlatArray<int> dcd = (*globdcdofs)[i];
	      
	      FlatVector<> wb(wbd.Size(), slh);
	      FlatVector<> dc(dcd.Size(), slh);
	    
	      for (int j = 0; j < dcd.Size(); j++)
		dc(j) = lx(dcd[j]);
	      
	      wb = Trans(*extwb_ref[elclassnr[i]]) * dc;
	      
	      // #pragma omp critical (bddc_restrict)
	      {
		for (int j = 0; j < wbd.Size(); j++)
		  lx(wbd[j]) -= wb(j);
	      }
	    }
	}
      }
      
      {
	NgProfiler::RegionTimer reg (timersolve);
	ly = 0.0;
	ly = (*inv) * lx;
      }
	

      // solve decoupled
      {
	NgProfiler::RegionTimer reg(timerdc);

	/*


#pragma omp parallel
	{
	  LocalHeap slh = lh.Split();
#pragma omp for
	  for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	  for (int i = 0; i < ne; i++)
	    {
	      if (elclassnr[i] != cl) continue;
	      HeapReset hr(slh);
	      FlatArray<int> dcd = (*globdcdofs)[i];
	      
	      FlatVector<> dc1(dcd.Size(), slh);
	      FlatVector<> dc2(dcd.Size(), slh);
	      
	      for (int j = 0; j < dcd.Size(); j++)
		dc1(j) = lx(dcd[j]);
	      
	      // dc2 = Trans(*invdc[i]) * dc1;
	      dc2 = (*invdc_ref[elclassnr[i]]) * dc1;
	      
	      for (int j = 0; j < dcd.Size(); j++)
		ly(dcd[j]) += s * dc2(j);
	      
	      NgProfiler::AddFlops (timerdc, sqr (dcd.Size()));
	    }
	   }
	*/

	for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	  {
	    int cnt = 0;
	    for (int i = 0; i < ne; i++)
	      if (elclassnr[i] == cl) cnt++;
	    if (cnt == 0) continue;

	    const Matrix<> & invdc = (*invdc_ref[cl]);
	    Matrix<> dc1(cnt, invdc.Height());
	    Matrix<> dc2(cnt, invdc.Height());

	    NgProfiler::AddFlops (timerdc, invdc.Width()*invdc.Height()*cnt);
		
	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dcd = (*globdcdofs)[i];

		for (int j = 0; j < dcd.Size(); j++)
		  dc1(cnt, j) = lx(dcd[j]);
		cnt++;
	      }

	    // vem = vim * Trans(ext);
	    LapackMultAB (dc1, invdc, dc2);
	    // dc2 = dc1 * invdc;

	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dcd = (*globdcdofs)[i];
		
		for (int j = 0; j < dcd.Size(); j++)
		  ly(dcd[j]) += s*dc2(cnt, j);
		cnt++;
	      }

	  }




      }

      // extend
      {
	NgProfiler::RegionTimer reg(timerext);

	// #pragma omp parallel
	{
	  LocalHeap slh = lh.Split();
	  // #pragma omp for
	  for (int i = 0; i < ne; i++)
	    {
	      HeapReset hr(slh);
	      FlatArray<int> wbd = (*globwbdofs)[i];
	      FlatArray<int> dcd = (*globdcdofs)[i];
	      
	      FlatVector<> wb(wbd.Size(), slh);
	      FlatVector<> dc(dcd.Size(), slh);
	      for (int j = 0; j < wbd.Size(); j++)
		wb(j) = ly(wbd[j]);
	      
	      // dc = (*extwb[i]) * wb;
	      dc = (*extwb_ref[elclassnr[i]]) * wb;
	      
	      for (int j = 0; j < dcd.Size(); j++)
		ly(dcd[j]) -= dc(j);
	    }
	}
      }

      transy = 0.0;
      for (int i = 0; i < restrict.Size(); i++)
	if (restrict[i] != -1)
	  transy(restrict[i]) += s * ly(i) / multiple[restrict[i]];

      
      {
	NgProfiler::RegionTimer reg(timertransy);
	/*
#pragma omp parallel
	{
	  LocalHeap slh = lh.Split();
#pragma omp for
	  for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	    for (int i = 0; i < ne; i++)
	    {
	      if (elclassnr[i] != cl) continue;
	      HeapReset hr(slh);
	      FlatArray<int> dint = (*globinternaldofs)[i];
	      FlatArray<int> dext = (*globexternaldofs)[i];
	      
	      FlatVector<> vi(dint.Size(), slh);
	      FlatVector<> ve(dext.Size(), slh);
	      
	      for (int j = 0; j < dext.Size(); j++)
		ve(j) = transy(dext[j]);
	      
	      vi = (*extension_int_ref[elclassnr[i]]) * ve;
	      
	      for (int j = 0; j < dint.Size(); j++)
		transy(dint[j]) -= vi(j);
	    }
	}
	*/



	for (int cl = 0; cl < inv_int_ref.Size(); cl++)
	  {
	    int cnt = 0;
	    for (int i = 0; i < ne; i++)
	      if (elclassnr[i] == cl) cnt++;
	    if (cnt == 0) continue;

	    const Matrix<> & ext = (*extension_int_ref_trans[cl]);
	    Matrix<> vim(cnt, ext.Width());
	    Matrix<> vem(cnt, ext.Height());

	    NgProfiler::AddFlops (timertransy, ext.Width()*ext.Height()*cnt);
		
	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dext = (*globexternaldofs)[i];

		for (int j = 0; j < dext.Size(); j++)
		  vem(cnt, j) = transy(dext[j]);
		cnt++;
	      }

	    // vem = vim * Trans(ext);
	    LapackMultAB (vem, ext, vim);

	    cnt = 0;	
	    for (int i = 0; i < ne; i++)
	      {
		if (elclassnr[i] != cl) continue;
		FlatArray<int> dint = (*globinternaldofs)[i];
		
		for (int j = 0; j < dint.Size(); j++)
		  transy(dint[j]) -= vim(cnt, j);
		cnt++;
	      }

	  }




      }
      
      y += s * transy;

      {
	NgProfiler::RegionTimer reg(timerinner);
	// #pragma omp parallel
	{
	  LocalHeap slh = lh.Split();
	  // #pragma omp for
	  for (int i = 0; i < ne; i++)
	    {
	      HeapReset hr(slh);
	      FlatArray<int> dint = (*globinternaldofs)[i];
	      
	      FlatVector<> v1(dint.Size(), slh);
	      FlatVector<> v2(dint.Size(), slh);
            
	      for (int j = 0; j < dint.Size(); j++)
		v1(j) = fx(dint[j]);
            
	      v2 = (*inv_int_ref[elclassnr[i]]) * v1;
            
	      for (int j = 0; j < dint.Size(); j++)
		fy(dint[j]) += s * v2(j);
	    }
	}
      }
    }

  };













  template <class SCAL>
  BDDCPreconditioner<SCAL> ::
  BDDCPreconditioner (const PDE * pde, const Flags & aflags, const string aname)
    : Preconditioner (pde, aflags, aname)
  {
    bfa = dynamic_cast<const S_BilinearForm<SCAL>*>(pde->GetBilinearForm (aflags.GetStringFlag ("bilinearform", NULL)));
    inversetype = flags.GetStringFlag("inverse", "sparsecholesky");
    refelement = flags.GetDefineFlag("refelement");
  }


  template <class SCAL>
  void BDDCPreconditioner<SCAL> ::
  Update ()
  {
    cout << "update bddc" << endl;
    if (refelement)
      pre = new BDDCMatrixRefElement(*bfa, inversetype);
    else
      pre = new BDDCMatrix(dynamic_cast<const HO_BilinearForm<double>&> (*bfa), inversetype);

    if (test) Test();
    
  }  



  template class BDDCPreconditioner<double>;
    
  namespace wirebasket_cpp
  {
   
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetPreconditionerClasses().AddPreconditioner ("wirebasket", WireBasketPreconditioner<double>::Create);
      GetPreconditionerClasses().AddPreconditioner ("wirebasket_complex", WireBasketPreconditioner<Complex>::Create);

      GetPreconditionerClasses().AddPreconditioner ("bddc", BDDCPreconditioner<double>::Create);
    }
    
    Init init;
  }


}
