#include <comp.hpp>
#include <solve.hpp>

#include <bddc.hpp>


#include <parallelngs.hpp>



namespace ngcomp
{

  class BDDCMatrix : public BaseMatrix
  {
    // const ElementByElement_BilinearForm<double> & bfa;
    const BilinearForm & bfa;
    Array<Matrix<double>*> & elmats;
    Array<int> restrict;
    Array<int> multiple;
    BaseMatrix * inv;
    BaseMatrix * inv_coarse;
    string inversetype;
    BitArray * free_dofs;
    const MeshAccess & ma;
    BaseMatrix * subassembled_harmonicext;
    BaseMatrix * subassembled_harmonicexttrans;
    BaseMatrix * subassembled_innersolve;    
    ///number of global wirebasket degrees of freedom    
    int nglobalwbdof;
    bool block;
    bool ebe;
    SparseMatrix<double> * pwbmat;    
    BaseVector * tmp;
    BaseVector * tmp2;
  public:
    
    void NEBEConstructor()  
    {
      *testout << "NEBE Constructor" << endl;

      if (id == 0) 
	cout << "BDDC-marix NEBE const" << endl;
      const FESpace & fes = bfa.GetFESpace();
      const MeshAccess & ma = fes.GetMeshAccess();
      int ne = ma.GetNE();
      int ndof = fes.GetNDof();      
      
      if (!bfa.UsesEliminateInternal()) throw Exception("please use eliminate_internal for the bilinearform");
      
      Array<int> wbdcnt(ne); //count number of wirebasket dofs on each element
      Array<int> ifcnt(ne); //count number of interface dofs on each element

      Array<int> wbdofs(fes.GetNDof()); //assignment: dof -> global number of the wirebasket dof (-1 for non-wirebasket-dofs)
      Array<int> ifdofs(fes.GetNDof()); //assignment: dof -> global number of the interface dof (-1 for non-wirebasket-dofs)
      wbdofs = -1;
      ifdofs = -1;
      Array<int> lwbdofs, lifdofs, dnums; //local (on one element) wirebasket/interface/all dofs
      
      for (int i = 0; i < ne; i++)
        {
          fes.GetDofNrs (i, lifdofs, INTERFACE_DOF);
	  int lifcnt = 0;
	  for (int j = 0; j < lifdofs.Size(); j++)
	    if(fes.GetFreeDofs()->Test(lifdofs[j])){ lifcnt++;}
          ifcnt[i] = lifcnt;
	  fes.GetDofNrs(i,lwbdofs,WIREBASKET_DOF);
	  wbdcnt[i] = lwbdofs.Size();	
        }

      Table<int> el2wbdofs(wbdcnt);   // wirebasket dofs on each element
      Table<int> el2ifdofs(ifcnt);   // interface dofs on each element
      
      BitArray ifdof(ndof);
      ifdof.Clear();
      BitArray wbdof(ndof);
      wbdof.Clear();
      
      int nifdofs = 0;
      int nwbdofs = 0;

      multiple.SetSize (fes.GetNDof());
      multiple = 0;

      *testout << "freedofs pointer = " << fes.GetFreeDofs() << endl;
      
      //Check how often each interfacedof/wirebasketdof appears 
      //TODO: (CL): Wann/Wo kommt die Gewichtung rein?
      for (int i = 0; i < ne; i++)
        {
	  fes.GetDofNrs (i, lifdofs, INTERFACE_DOF);
	  int lifcnt = 0;
	  for (int j = 0; j< lifdofs.Size();j++)
	  {
	    if (!ifdof.Test(lifdofs[j])&&fes.GetFreeDofs()->Test(lifdofs[j])){
	      ifdof.Set(lifdofs[j]);
	      nifdofs++;
	      multiple[lifdofs[j]]++;
	    }else{
	      multiple[lifdofs[j]]++;
	    }
	    if (fes.GetFreeDofs()->Test(lifdofs[j]))
	      el2ifdofs[i][lifcnt++] = lifdofs[j];
	  }

	  fes.GetDofNrs(i,lwbdofs,WIREBASKET_DOF);
	  for (int j = 0; j< lwbdofs.Size();j++)
	  {
	    el2wbdofs[i][j] = lwbdofs[j];
	    if (!wbdof.Test(lwbdofs[j]))
	    {
	      wbdof.Set(lwbdofs[j]);
	      nwbdofs++;
	    }	  
	  }
        }


#ifdef PARALLEL
      // accumulate multiple
      ParallelVVector<double> pv_multiple (fes.GetNDof(), &fes.GetParallelDofs(), DISTRIBUTED);

      if (ntasks > 1 && id == 0)
	{
	  pv_multiple.Cumulate();  
	}
      else
	{
	  for (int i = 0; i < multiple.Size(); i++)
	    pv_multiple(i) = multiple[i];
	  pv_multiple.Cumulate();  
	  for (int i = 0; i < multiple.Size(); i++)
	    multiple[i] = int(pv_multiple(i));
	}
#endif


      *testout << " number of external dofs: " << nifdofs+nwbdofs << endl;
      *testout << " number of interface dofs: " << nifdofs << endl;
      *testout << " number of wirebasket dofs: " << nwbdofs << endl;
      *testout << " el2wbdofs = " << endl << el2wbdofs << endl;
      *testout << " el2ifdofs = " << endl << el2ifdofs << endl;
       
      
      MatrixGraph graph_harmonicext(ndof, el2ifdofs, el2wbdofs, false);
      MatrixGraph graph_innersolve(ndof, el2ifdofs, el2ifdofs, bfa.IsSymmetric());
      MatrixGraph graph_wbschur(ndof, el2wbdofs, el2wbdofs, bfa.IsSymmetric());
      
      
      if (!bfa.IsSymmetric()){
	MatrixGraph graph_harmonicexttrans(ndof, el2wbdofs, el2ifdofs, false);
	subassembled_harmonicexttrans = new SparseMatrix<double>(graph_harmonicexttrans, 1);
	subassembled_innersolve = new SparseMatrix<double>(graph_innersolve, 1);
      }
      else
	subassembled_innersolve = new SparseMatrixSymmetric<double>(graph_innersolve, 1);
      
      subassembled_harmonicext = new SparseMatrix<double>(graph_harmonicext, 1);


      if (bfa.IsSymmetric()){
	pwbmat = new SparseMatrixSymmetric<double>(graph_wbschur,1);
	if (id == 0)
	  cout << "symmetric" << endl << endl;
      }
      else{
	pwbmat = new SparseMatrix<double>(graph_wbschur,1);
	if (id == 0)
	  cout << "nonsymmetric" << endl << endl;
      }
      
      SparseMatrix<double>& wbmat=*pwbmat;


      wbmat.SetInverseType (inversetype);
      
      // cout << "have matrix" << endl << endl;


      if (ntasks == 1 || id > 0)
      for (int i = 0; i < ne; i++)
        {
          FlatMatrix<> elmat = *elmats[i]; 
	  // dynamic_cast<const ElementByElementMatrix<double>&> (bfa.GetMatrix()) . GetElementMatrix (i);
	  
	  Array<int> interfacedofs; interfacedofs = el2ifdofs[i];
	  Array<int> wirebasketdofs; wirebasketdofs = el2wbdofs[i];
	  Array<int> localwbdofs; localwbdofs.SetSize(0); //local dofs
	  Array<int> localintdofs; localintdofs.SetSize(0); //local dofs 
	  Array<int> ldnums;
	  fes.GetDofNrs (i, ldnums, EXTERNAL_DOF);
	  
          for (int k = 0; k < ldnums.Size(); k++){
	    if (wbdof.Test(ldnums[k])){
	      localwbdofs.Append(k);
	    }
	    else if (ifdof.Test(ldnums[k])){
	      localintdofs.Append(k);
	    }
	  }
	  
	  int sizew = localwbdofs.Size();
	  int sizei = localintdofs.Size();

	  Matrix<double> a(sizew, sizew);
	  Matrix<double> b(sizew, sizei);
	  Matrix<double> c(sizei, sizew);
	  Matrix<double> d(sizei, sizei);
	  for (int k = 0; k < sizew; k++)
	    for (int l = 0; l < sizew; l++)
	      a(k,l) = elmat(localwbdofs[k], localwbdofs[l]);
	    
	  if (el2ifdofs[i].Size())
	  {      
	    for (int k = 0; k < sizew; k++)
	      for (int l = 0; l < sizei; l++)
		{
		  b(k,l) = elmat(localwbdofs[k], localintdofs[l]);
		  c(l,k) = elmat(localintdofs[l], localwbdofs[k]);
		}

	    for (int k = 0; k < sizei; k++)
	      for (int l = 0; l < sizei; l++)
		d(k,l) = elmat(localintdofs[k], localintdofs[l]);
	      
#ifdef LAPACK
	    LapackInverse (d);
#else
	    Matrix<double> invd(sizei);
	    CalcInverse (d, invd);	  
	    d = invd;
#endif //LAPACK
	    Matrix<double> he (sizei, sizew);
#ifdef LAPACK
	    he = 0.0;
	    LapackMultAddAB(d,c,-1.0,he);
#else	  
	    he = -1.0 * d * c;
#endif	 

#ifdef LAPACK	  
	    LapackMultAddAB (b, he, 1.0, a);
#else
	    a += b*he;
#endif //LAPACK

	    //R * E
	    for (int k = 0; k < sizei; k++)
	      for (int l = 0; l < sizew; l++)
		if ((fes.GetFreeDofs()) && (!fes.GetFreeDofs()->Test(el2ifdofs[i][k]))){ //TEST:NIE
		  cout << "((fes.GetFreeDofs()) && (!fes.GetFreeDofs()->Test(el2ifdofs[i][k]))) == true" << endl;
		  he(k,l) = 0.0;
		}
		else
		  he(k,l) /= multiple[ el2ifdofs[i][k] ];
		
	    dynamic_cast<SparseMatrix<double>*>(subassembled_harmonicext)->AddElementMatrix(interfacedofs,wirebasketdofs,he);
	    
	    if (!bfa.IsSymmetric()){
	      Matrix<double> het (sizew, sizei);
#ifdef LAPACK	    
	      het = 0.0;
	      LapackMultAddAB(b,d,-1.0,het);
#else	    
	      het = -1.0 * b * d;
#endif	    
	      //E * R^T
	      for (int k = 0; k < sizew; k++)
		for (int l = 0; l < sizei; l++)
		  if ((fes.GetFreeDofs()) && (!fes.GetFreeDofs()->Test(el2ifdofs[i][l]))) { //TEST:NIE
		    cout << "((fes.GetFreeDofs()) && (!fes.GetFreeDofs()->Test(el2ifdofs[i][l]))) == true" << endl;
		    het(k,l) = 0.0;
		  }
		  else
		    het(k,l) /= multiple[ el2ifdofs[i][l] ];
	      
	      dynamic_cast<SparseMatrix<double>*>(subassembled_harmonicexttrans)->AddElementMatrix(wirebasketdofs,interfacedofs,het);
	    }
	    //R * A_ii^(-1) * R^T
	    for (int k = 0; k < sizei; k++)
	      for (int l = 0; l < sizei; l++)
		if ((fes.GetFreeDofs()) && ((!fes.GetFreeDofs()->Test(el2ifdofs[i][l])) || (!fes.GetFreeDofs()->Test(el2ifdofs[i][k])))) { //TEST:NIE
		  cout << "((fes.GetFreeDofs()) && ((!fes.GetFreeDofs()->Test(el2ifdofs[i][l])) || (!fes.GetFreeDofs()->Test(el2ifdofs[i][k])))) == true" << endl;
		  d(k,l) = 0.0;
		}
		else
		  d(k,l) /= multiple[ el2ifdofs[i][k] ] * multiple[ el2ifdofs[i][l] ];	    
	    if (bfa.IsSymmetric())
	      dynamic_cast<SparseMatrixSymmetric<double>*>(subassembled_innersolve)->AddElementMatrix(interfacedofs,d);
	    else
	      dynamic_cast<SparseMatrix<double>*>(subassembled_innersolve)->AddElementMatrix(interfacedofs,interfacedofs,d);
	  }
	  
	  if (bfa.IsSymmetric())
	    dynamic_cast<SparseMatrixSymmetric<double>&>(wbmat).AddElementMatrix(wirebasketdofs,a);
	  else
	    dynamic_cast<SparseMatrix<double>&>(wbmat).AddElementMatrix(wirebasketdofs,wirebasketdofs,a);

        }
      

      // cout << "matrix filed" << endl;


      free_dofs = new BitArray (ndof);
      free_dofs->Clear();

      *free_dofs = wbdof;
      if (fes.GetFreeDofs())
	free_dofs -> And (*fes.GetFreeDofs());
      int cntfreedofs=0;
      for (int i = 0; i < free_dofs->Size(); i++)
	if (free_dofs->Test(i)) cntfreedofs++;


      *testout << "fes.GetFreeDofs() = " << fes.GetFreeDofs() << endl;
      *testout << "BDDC - precond, orig freedofs = " << *fes.GetFreeDofs() << endl;
      *testout << "local freedofs = " << *free_dofs << endl;
      *testout << "wbdof = " << wbdof << endl;


      if (block){
	//Smoothing Blocks
	Flags flags;
	flags.SetFlag("eliminate_internal");
	flags.SetFlag("subassembled");
	cout << "call Create Smoothing Blocks of " << bfa.GetFESpace().GetName() << endl;
	Table<int> & blocks = *(bfa.GetFESpace().CreateSmoothingBlocks(flags));
	cout << "has blocks" << endl << endl;
 	*testout << " blocks = \n " << blocks << endl;

	cout << "call block-jacobi inverse" << endl;
	inv = wbmat.CreateBlockJacobiPrecond(blocks, 0, 0, 0);      
	cout << "has inverse" << endl << endl;

	//Coarse Grid of Wirebasket
	cout << "call directsolverclusters inverse" << endl;
	Array<int> & clusters = *(bfa.GetFESpace().CreateDirectSolverClusters(flags));
	cout << "has clusters" << endl << endl;
 	*testout << " clusters = \n " << clusters << endl;

	cout << "call coarse wirebasket grid inverse" << endl;
	inv_coarse = wbmat.InverseMatrix(&clusters);
	cout << "has inverse" << endl << endl;
	
	tmp = new VVector<>(ndof);
	tmp2 = new VVector<>(ndof);
      }
      else
      {
#ifdef PARALLEL
	if (ntasks > 1)
	  {
	    // *testout << "bddc, local wbmatrix = " << endl << wbmat << endl;
	    inv = new MasterInverse<double> (wbmat, free_dofs, &bfa.GetFESpace().GetParallelDofs());
	    tmp = new ParallelVVector<>(ndof, &bfa.GetFESpace().GetParallelDofs());
	    subassembled_innersolve = new ParallelMatrix (*subassembled_innersolve, bfa.GetFESpace().GetParallelDofs());
	    subassembled_harmonicext = new ParallelMatrix (*subassembled_harmonicext, bfa.GetFESpace().GetParallelDofs());
	    subassembled_harmonicexttrans = new ParallelMatrix (*subassembled_harmonicexttrans, bfa.GetFESpace().GetParallelDofs());
	  }
	else
#endif
	  {
	    cout << "call wirebasket inverse ( with " << cntfreedofs << " free dofs )" << endl;
	    *testout << "wbmat = " << endl << wbmat << endl;
	    inv = wbmat.InverseMatrix(free_dofs);
	    cout << "has inverse" << endl;
	    tmp = new VVector<>(ndof);
	  }
      }
    };

    
   void EBEConstructor()
   {
      const FESpace & fes = bfa.GetFESpace();
      const MeshAccess & ma = fes.GetMeshAccess();
      int ne = ma.GetNE();
      
      if (!bfa.UsesEliminateInternal()) throw Exception("please use eliminate_internal for the bilinearform");
      
      Array<int> cnt(ne); //count number of (external) dofs on each element
      Array<int> wbdcnt(ne); //count number of wirebasket dofs on each element
      Array<int> wifcnt(ne); //count number of wirebasket dofs on each element

      Array<int> wbdofs(fes.GetNDof()); //assignment: dof -> global number of the wirebasket dof (-1 for non-wirebasket-dofs)
      Array<int> ifdofs(fes.GetNDof()); //assignment: dof -> global number of the interface dof (-1 for non-wirebasket-dofs)
      wbdofs = -1;
      ifdofs = -1;
      Array<int> lwbdofs, lifdofs, dnums; //local (on one element) wirebasket (all) dofs
      nglobalwbdof = 0; 
      
      for (int i = 0; i < ne; i++)
        {
          fes.GetDofNrs (i, dnums, EXTERNAL_DOF);	  
          cnt[i] = dnums.Size();

	  fes.GetDofNrs(i,lwbdofs,WIREBASKET_DOF);
	  wbdcnt[i] = lwbdofs.Size();	  
	  fes.GetDofNrs(i,lifdofs,INTERFACE_DOF);
	  wifcnt[i] = lifdofs.Size();	  
          for (int j = 0; j < lwbdofs.Size(); j++){
	    if (wbdofs[lwbdofs[j]] == -1){
	      wbdofs[lwbdofs[j]] = nglobalwbdof++;
	    }
	  }
// 	  *testout << "dnums = " << endl << dnums << endl;
// 	  *testout << "lwbdofs = " << endl << lwbdofs << endl;
        }

      Table<int> dcdofs(cnt);   // discontinuous dofs on each element
      Table<int> el2wbdofs(wbdcnt);   // wirebasket dofs on each element
      Table<int> el2ifdofs(wifcnt);   // interface dofs on each element

      int ndcdof = nglobalwbdof;

      for (int i = 0; i < ne; i++)
      {
	fes.GetDofNrs (i, dnums, EXTERNAL_DOF);	  
	int wbd = 0;
	for (int j = 0; j < dnums.Size(); j++)
	  {
	    if (dnums[j] == -1) continue;
	    if (wbdofs[dnums[j]]!=-1) {
	      el2wbdofs[i][wbd++] = wbdofs[dnums[j]];
	      dnums[j] = wbdofs[dnums[j]];
	      continue;
	    }
	    dnums[j] = ndcdof;
	    ndcdof++;
	  }

	for (int j = 0; j < dnums.Size(); j++){
	  dcdofs[i][j] = dnums[j];
	}
      }

       *testout << "dcdofs = " << endl << dcdofs << endl;
       *testout << "el2wbdofs = " << endl << el2wbdofs << endl;

//        *testout << "test " << endl;
      restrict.SetSize(ndcdof);
      restrict = -1;
      for (int i = 0; i < ne; i++)
        {
          fes.GetDofNrs (i, dnums, EXTERNAL_DOF);	  

          for (int j = 0; j < dnums.Size(); j++)
            {
              if (dnums[j] != -1){
                restrict[dcdofs[i][j]] = dnums[j];
	      }
            }
        }

//       *testout << "restrict = " << endl << restrict << endl;

      multiple.SetSize (fes.GetNDof());
      multiple = 0;
      for (int i = 0; i < restrict.Size(); i++)
        if (restrict[i] != -1)
          multiple[restrict[i]]++;

//       *testout << "multiple = " << endl << multiple << endl;
	cout << "now build graph" << endl;
	MatrixGraph graph(nglobalwbdof, el2wbdofs, el2wbdofs, bfa.IsSymmetric());
	cout << "now allocate matrix" << endl;
	if (bfa.IsSymmetric()){
	  pwbmat = new SparseMatrixSymmetric<double>(graph,1);
	  cout << "symmetric" << endl;
	}
	else{
	  pwbmat = new SparseMatrix<double>(graph,1);
	  cout << "nonsymmetric" << endl;
	}

      SparseMatrix<double>& wbmat=*pwbmat;

      wbmat.SetInverseType (inversetype);
      
      cout << "have matrix" << endl;
      
      subassembled_harmonicext = new ElementByElementMatrix<double>(ndcdof, ne);
      if (!bfa.IsSymmetric())
	subassembled_harmonicexttrans = new ElementByElementMatrix<double>(ndcdof, ne);
      subassembled_innersolve = new ElementByElementMatrix<double>(ndcdof, ne, true);


      for (int i = 0; i < ne; i++)
        {
          FlatMatrix<> elmat = *elmats[i];
	  // dynamic_cast<const ElementByElementMatrix<double>&> (bfa.GetMatrix()) . GetElementMatrix (i);

	  FlatArray<int> dofs2 = dcdofs[i];
	  Array<int> idnums; idnums.SetSize(0); //global dofs
	  Array<int> wdnums; wdnums.SetSize(0); //global dofs
	  Array<int> localwbdofs; localwbdofs.SetSize(0); //local dofs
	  Array<int> localintdofs; localintdofs.SetSize(0); //local dofs  
          for (int k = 0; k < dofs2.Size(); k++){
	    if (dofs2[k]<nglobalwbdof){
	      localwbdofs.Append(k);
	      wdnums.Append(dofs2[k]);
	    }
	    else{
	      localintdofs.Append(k);
	      idnums.Append(dofs2[k]);
	    }
	  }

	  int sizew = localwbdofs.Size();
	  int sizei = localintdofs.Size();
	  Matrix<double> a(sizew, sizew);
	  Matrix<double> b(sizew, sizei);
	  Matrix<double> c(sizei, sizew);
	  Matrix<double> d(sizei, sizei);
	  for (int k = 0; k < sizew; k++)
	    for (int l = 0; l < sizew; l++)
	      a(k,l) = elmat(localwbdofs[k], localwbdofs[l]);
	    
	  if (idnums.Size())
	  {      
	    for (int k = 0; k < sizew; k++)
	      for (int l = 0; l < sizei; l++)
		{
		  b(k,l) = elmat(localwbdofs[k], localintdofs[l]);
		  c(l,k) = elmat(localintdofs[l], localwbdofs[k]);
		}

	    for (int k = 0; k < sizei; k++)
	      for (int l = 0; l < sizei; l++)
		d(k,l) = elmat(localintdofs[k], localintdofs[l]);
	      
#ifdef LAPACK
	    LapackInverse (d);
#else
	    Matrix<double> invd(size);
	    CalcInverse (d, invd);	  
	    d = invd;
#endif //LAPACK
	    Matrix<double> he (sizei, sizew);
#ifdef LAPACK	  
	    he = 0.0;
	    LapackMultAddAB(d,c,-1.0,he);
#else	  
	    he = -1.0 * d * c;
#endif	 
	    static_cast<ElementByElementMatrix<double>*>(subassembled_harmonicext)->AddElementMatrix(i,idnums,wdnums,he);
	    
	    if (!bfa.IsSymmetric()){
	      Matrix<double> het (sizew, sizei);
#ifdef LAPACK	    
	      het = 0.0;
	      LapackMultAddAB(b,d,-1.0,het);
#else	    
	      het = -1.0 * b * d;
#endif	    
	      static_cast<ElementByElementMatrix<double>*>(subassembled_harmonicexttrans)->AddElementMatrix(i,wdnums,idnums,het);
	    }
	    static_cast<ElementByElementMatrix<double>*>(subassembled_innersolve)->AddElementMatrix(i,idnums,idnums,d);
#ifdef LAPACK	  
	    LapackMultAddAB (b, he, 1.0, a);
#else
	    a += b*he;
#endif //LAPACK
	  }
	  for (int k = 0; k < wdnums.Size(); k++)
            for (int l = 0; l < wdnums.Size(); l++)
              if (wdnums[k] != -1 && wdnums[l] != -1 && (!bfa.IsSymmetric() || wdnums[k] >= wdnums[l]))
                wbmat(wdnums[k], wdnums[l]) += a(k,l);	  

        }
      
      cout << "matrix filled" << endl;
//       *testout << "wbmat = " << endl << wbmat << endl;
      
      free_dofs = new BitArray (nglobalwbdof);
      free_dofs->Clear();
      
      for (int i = 0; i < nglobalwbdof; i++)
	{
	  if (restrict[i] != -1)
	    if ((fes.GetFreeDofs()==NULL) || fes.GetFreeDofs()->Test(restrict[i]))
	      free_dofs->Set(i);
	}
      *testout << "wb_free_dofs = " << *free_dofs << endl;
      if (block){
	Flags flags;
	flags.SetFlag("eliminate_internal");
	flags.SetFlag("subassembled");
	Table<int> & blocks = *(bfa.GetFESpace().CreateSmoothingBlocks(flags));
	for (int i = 0; i < blocks.Size(); i++){
	  for (int j = 0; j < blocks[i].Size(); j++)
	    blocks[i][j] = wbdofs[blocks[i][j]];
	}
	cout << "call block-jacobi inverse" << endl;
	inv = wbmat.CreateBlockJacobiPrecond(blocks, 0, 0, 0);      
	cout << "has inverse" << endl;
	cout << "call directsolverclusters inverse" << endl;
	Array<int> & clusters = *(bfa.GetFESpace().CreateDirectSolverClusters(flags));
 	*testout << " clusters = \n " << clusters << endl;
	Array<int> & condclusters = *new Array<int>(nglobalwbdof);
	for (int i=0; i< clusters.Size(); i++)
	  condclusters[wbdofs[i]] = clusters[i];
 	*testout << " condclusters = \n " << condclusters << endl;
	
	inv_coarse = wbmat.InverseMatrix(&condclusters);
	tmp = new VVector<>(nglobalwbdof);
      }
      else
      {
	inv = wbmat.InverseMatrix(free_dofs);
      }
      cout << "has inverse" << endl;
//       *testout << "inverse2 = " << (*inv) << endl;
    }
      
          
    
    
    BDDCMatrix (const BilinearForm & abfa, Array<Matrix<double>*> & aelmats, const string & ainversetype, bool ablock, bool aebe)
      : bfa(abfa), elmats(aelmats), inversetype(ainversetype), ma(abfa.GetFESpace().GetMeshAccess()),block(ablock), ebe(aebe)
   {
      pwbmat = NULL;
      inv = NULL;
      inv_coarse = NULL;
      tmp = NULL;
      tmp2 = NULL;
      // LocalHeap lh(10000000);     
     
      if(ebe) 
	EBEConstructor();
      else
	NEBEConstructor();
   }

      
      

    ~BDDCMatrix()
    {
      if (inv) delete inv;
      if (pwbmat) delete pwbmat;
      if (inv_coarse) delete inv_coarse;
      if (tmp) delete tmp;
      if (tmp2) delete tmp2;
    }
    

    
    virtual void MultAddNEBE (double s, const BaseVector & x, BaseVector & y) const
    {
      static Timer timer ("Apply BDDC preconditioner");
      static Timer timerifs ("Apply BDDC preconditioner - apply ifs");
      static Timer timerwb ("Apply BDDC preconditioner - wb solve");
      static Timer timeretc ("Apply BDDC preconditioner - etc");
      static Timer timerharmonicext ("Apply BDDC preconditioner - harmonic extension");
      static Timer timerharmonicexttrans ("Apply BDDC preconditioner - harmonic extension trans");
      
      NgProfiler::RegionTimer reg (timer);
      
      x.Cumulate();
    
      y = x;

      timerharmonicexttrans.Start();

      if (bfa.IsSymmetric())
	y += Transpose(*subassembled_harmonicext) * x; 
      else
	y += *subassembled_harmonicexttrans * x;
      
      timerharmonicexttrans.Stop();

      timerwb.Start();
      *tmp = 0;
      if (block){
	if (true) //GS
	  {
	    dynamic_cast<BaseBlockJacobiPrecond*>(inv)->GSSmoothResiduum (*tmp, y, *tmp2 ,1);
	    if (inv_coarse)
		  *tmp += (*inv_coarse) * *tmp2; 
	    dynamic_cast<BaseBlockJacobiPrecond*>(inv)->GSSmoothBack (*tmp, y);
	  }else{ //jacobi only (old)
	  *tmp = (*inv) * y;
	  *tmp += (*inv_coarse) * y; 
	}
      }
      else
	{
	  *tmp = (*inv) * y;
	}
      timerwb.Stop();


      timerifs.Start();
      *tmp += *subassembled_innersolve * x;
      timerifs.Stop();

      
      timerharmonicext.Start();
      
      y = *tmp;
      y += *subassembled_harmonicext * *tmp;
    
      timerharmonicext.Stop();

      y.Cumulate();
    }    
    
    
    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      if (!ebe) {MultAddNEBE(s,x,y);return;}
      static int timer = NgProfiler::CreateTimer ("Apply BDDC preconditioner");
      static int timerifs = NgProfiler::CreateTimer ("Apply BDDC preconditioner - apply ifs");
      static int timerwb = NgProfiler::CreateTimer ("Apply BDDC preconditioner - wb solve");
      static int timeretc = NgProfiler::CreateTimer ("Apply BDDC preconditioner - etc");
      NgProfiler::RegionTimer reg (timer);

      const FESpace & fes = bfa.GetFESpace();      
      FlatVector<> fx = x.FVDouble();
      FlatVector<> fy = y.FVDouble();

      NgProfiler::StartTimer (timeretc);

      VVector<> lx2(restrict.Size());
      VVector<> ly2(restrict.Size());

      for (int i = 0; i < restrict.Size(); i++)
//         if (restrict[i] != -1 && fes.GetFreeDofs()->Test(restrict[i]))
	if (restrict[i] != -1 && ( (!fes.GetFreeDofs()) || (fes.GetFreeDofs()->Test(restrict[i])) ))	  
          lx2(i) = fx(restrict[i]) / multiple[restrict[i]];
	else
	   lx2(i) = 0.0;


      NgProfiler::StopTimer (timeretc);
	
      if (bfa.IsSymmetric())
	lx2 += Transpose(*subassembled_harmonicext) * lx2; 
      else
	lx2 += *subassembled_harmonicexttrans * lx2;

      NgProfiler::StartTimer (timerwb);
      
      BaseVector & subx = *(lx2.Range(0,nglobalwbdof));
      BaseVector & suby = *(ly2.Range(0,nglobalwbdof));
      BaseVector & res = *tmp;
       ly2 = 0.0;
      if (block){
	if (false) //GS
	{
	  dynamic_cast<BaseBlockJacobiPrecond*>(inv)->GSSmoothResiduum (suby, subx, res,1);
	  if (inv_coarse)
	    suby += (*inv_coarse) * res; 
	  dynamic_cast<BaseBlockJacobiPrecond*>(inv)->GSSmoothBack (suby, subx);
	}else{ //jacobi only (old)
	  suby = (*inv) * subx;
	  suby += (*inv_coarse) * subx; 
	}
      }
      else
      {
	suby = (*inv) * subx;
      }

      NgProfiler::StopTimer (timerwb);

      NgProfiler::StartTimer (timerifs);
      ly2 += *subassembled_innersolve * lx2;
      NgProfiler::StopTimer (timerifs);

      ly2 += *subassembled_harmonicext * ly2;


      NgProfiler::StartTimer (timeretc);

      for (int i = 0; i < restrict.Size(); i++)
	{
	  if (restrict[i] != -1 && ( (fes.GetFreeDofs()) && (!fes.GetFreeDofs()->Test(restrict[i])) ))
	    ly2(i) = 0.0;
	}

      fy = 0.0;
      for (int i = 0; i < restrict.Size(); i++){
        if (restrict[i] != -1)
          fy(restrict[i]) += s * ly2(i) / multiple[restrict[i]];
      }

      NgProfiler::StopTimer (timeretc);
    }

  };












  class BDDCMatrixRefElement : public BaseMatrix
  {
    const S_BilinearForm<double> & bfa;
    Array<int> restrict;
    Array<int> multiple;
    BaseMatrix * inv, *inv_coarse;
    
    Table<int> *wbdofs, *dcdofs, *internaldofs, *externaldofs;
    Table<int> *globwbdofs, *globdcdofs, *globinternaldofs, *globexternaldofs;

    Table<int> *interface_dofs;
    Table<double> *weighting;

    Array<int> elclassnr;
    Table<int> * els_of_class;
    Array<Matrix<>*> invdc_ref, extwb_ref, schurwb_ref;
    Array<Matrix<>*> inv_int_ref, extension_int_ref, extension_int_ref_trans;
    int ne, ndof;
    
    string inversetype;
    bool print;

    bool block;
    BitArray * free_dofs_wb;

  public:
    BDDCMatrixRefElement (const S_BilinearForm<double> & abfa, const string & ainversetype)
      : bfa(abfa), inversetype (ainversetype)
    {
      cout << "BDDC MatrixRefElement" << endl;
      print = true;
      block = true;
      LocalHeap lh(500000000, "BDDC - MatrixRefEleemnt");

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
	  fes.GetDofNrs (i, extdnums, EXTERNAL_DOF);	  
	  fes.GetDofNrs (i, wbdnums, WIREBASKET_DOF);	  

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
          fes.GetDofNrs (i, extdnums, EXTERNAL_DOF);


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
      
      MatrixGraph graph(firstdcdof, *globwbdofs,  *globwbdofs, 1);
      SparseMatrixSymmetric<double> & dcmat 
	= *new SparseMatrixSymmetric<double> (graph, 1);

      dcmat = 0.0;
      // dcmat.SetInverseType ("sparsecholesky");
      dcmat.SetInverseType (inversetype);
      *testout << "dcmat, init = " << dcmat << endl;

      elclassnr.SetSize(ne);
      int maxclass = ne; // 32;

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
	  // classnr = i;

          elclassnr[i] = classnr;
          if (invdc_ref[classnr]) continue;

          cout << "assemble class " << classnr << " " << flush;

	  if (print)
	    *testout << "assemble class " << classnr << endl;

          fes.GetDofNrs (i, dnums);
          fes.GetDofNrs (i, extdnums, EXTERNAL_DOF);
	  
          Matrix<> elmat(dnums.Size());
          Matrix<> partelmat(dnums.Size());
          elmat = 0.0;

          ElementTransformation eltrans;
          ma.GetElementTransformation (i, eltrans, lh);
          const FiniteElement & fel = fes.GetFE(i, lh);

          for (int j = 0; j < bfa.NumIntegrators(); j++)
            {
	      if (bfa.GetIntegrator(j) -> BoundaryForm()) continue;
              bfa.GetIntegrator(j)->CalcElementMatrix (fel, eltrans, partelmat, lh);
              elmat += partelmat;
            }

	  

	  if (print)
	    *testout << "elmat = " << endl << elmat << endl;
	  
          int nwb = (*wbdofs)[i].Size();
          int ndc = (*dcdofs)[i].Size();
          int nint = (*internaldofs)[i].Size();
          int next = (*externaldofs)[i].Size();

	  cout << " wb = " << nwb << " nifs = " << ndc << " loc = " << nint << endl;
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

      TableCreator<int> creator;
      for ( ; !creator.Done(); creator++)
	{
	  for (int i = 0; i < ne; i++)
	    creator.Add (elclassnr[i], i);
	}
      els_of_class = creator.GetTable();

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
          fes.GetDofNrs (i, extdnums, EXTERNAL_DOF);

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

      
      // Table<int> *interface_dofs;
      // Table<double> *weighting;
      Array<int> cnt_interface(ne);
      cnt_interface = 0;

      for (int i = 0; i < ne; i++)
	cnt_interface[i] = (*globdcdofs)[i].Size();
      
      interface_dofs = new Table<int> (cnt_interface);
      weighting = new Table<double> (cnt_interface);

      for (int i = 0; i < ne; i++)
	{
	  FlatArray<int> dc = (*globdcdofs)[i];
	  for (int j = 0; j < dc.Size(); j++)
	    {
	      int rest = restrict[dc[j]];

	      if (rest == -1) cout << "rest == -1" << endl;
	      if (multiple[rest] <= 0) cout << "multiplie = " << multiple[rest] << endl;

	      (*interface_dofs)[i][j] = rest;

	      if (free_dofs.Test(rest))
		(*weighting)[i][j] = 1.0 / multiple[rest];
	      else
		(*weighting)[i][j] = 0.0;
	    }
	}

      cout << "tables filled" << endl;




      for (int i = 0; i < restrict.Size(); i++)
	if (restrict[i] != -1 && !free_dofs.Test(restrict[i]))
	  restrict[i] = -1;


      free_dofs_wb = new BitArray (free_dofs.Size());

      free_dofs_wb -> Clear();
      for (int i = 0; i < firstdcdof; i++)
        if (gwbdofs[i] != 0)
	  free_dofs_wb->Set(i);
      free_dofs_wb -> And (free_dofs);

      




      cout << "call inverse, type = " << inversetype << endl;
      cout << "dim = " << dcmat.Height() << endl;



      if (block)
	{
	  Flags flags;
	  flags.SetFlag("eliminate_internal");
	  flags.SetFlag("subassembled");
	  Table<int> & blocks = *(bfa.GetFESpace().CreateSmoothingBlocks(flags));
	  /*
	  for (int i = 0; i < blocks.Size(); i++){
	    for (int j = 0; j < blocks[i].Size(); j++)
	      blocks[i][j] = wbdofs[blocks[i][j]];
	  }
	  */
	  cout << "call block-jacobi inverse" << endl;
	  inv = dcmat.CreateBlockJacobiPrecond(blocks, 0, 0, 0);      

	  cout << "dc.height = " << dcmat.Height() << endl;
	  cout << "inv.height = " << inv->Height() << endl;

	  cout << "has inverse" << endl;
	  cout << "call directsolverclusters inverse" << endl;
	  Array<int> & clusters = *(bfa.GetFESpace().CreateDirectSolverClusters(flags));
	  // 	*testout << " clusters = \n " << clusters << endl;
	  /*
	  Array<int> & condclusters = *new Array<int>(nglobalwbdof);
	  for (int i=0; i< clusters.Size(); i++)
	    condclusters[wbdofs[i]] = clusters[i];
	  inv_coarse = dcmat.InverseMatrix(&condclusters);
	  */
	  // 	*testout << " condclusters = \n " << condclusters << endl;
	  inv_coarse = dcmat.InverseMatrix(&clusters);
	  // tmp = new VVector<>(nglobalwbdof);

	  for (int i = 0; i < clusters.Size(); i++)
	    if (clusters[i] && !free_dofs_wb->Test(i))
	      cout << "wrong dirichet bc" << endl;
	}
      else
	{
	  cout << "call inverse" << endl;
	  inv = dcmat.InverseMatrix(free_dofs_wb);
	}

      // inv = dcmat.InverseMatrix(&free_dofs);
      return;




#define sometestsxx
#ifdef sometests

      Array<int> & directblocks = *new Array<int> (free_dofs.Size());
      directblocks = 0;
      // int nv = ma.GetNV();
      int ned = ma.GetNEdges();
      int nfa = ma.GetNFaces();
      for (int i = 0; i < directblocks.Size(); i++)
	if (free_dofs.Test(i))
	  directblocks[i] = 1;

      
      directblocks = 0;

      /*
      for (int i = 0; i < ne; i++)
	{
	  fes.GetDofNrs (i, dnums);
	  if (!dnums.Size()) continue;

	  if (free_dofs.Test(dnums[0]))
	    directblocks[dnums[0]] = 1;

	  for (int j = 0; j < dnums.Size(); j++)
	    if (free_dofs.Test(dnums[j]))
	      directblocks[dnums[j]] = 1;
	}
      */

      for (int i = 0; i < -nfa; i++)
	{
	  fes.GetFaceDofNrs (i, dnums);
	  if (!dnums.Size()) continue;

	  if (free_dofs.Test(dnums[0]))
	    directblocks[dnums[0]] = 1;

	  /*
	  for (int j = 1; j < dnums.Size(); j++)
	    if (free_dofs.Test(dnums[j]))
	      directblocks[dnums[j]] = 1;
	  */
	}



      for (int i = 0; i < ned; i++)
	{
	  fes.GetEdgeDofNrs (i, dnums);
	  if (!dnums.Size()) continue;

	  if (free_dofs.Test(dnums[0]))
	    directblocks[dnums[0]] = 1;
	  for (int j = 1; j < dnums.Size(); j++)
	    if (free_dofs.Test(dnums[j]))
	      directblocks[dnums[j]] = i+2;
	}
      
      cout << "directblocks.Size() = " << directblocks.Size() << endl;
      cout << "freedofs.Size() = " << free_dofs.Size() << endl;
      cout << "mat.size = " << dcmat.Height() << endl;

      /*
      *testout << "directblocks,1 = " << endl << directblocks << endl;
      
      directblocks = 0;
      for (int i = 0; i < directblocks.Size(); i++)
	if (free_dofs.Test(i)) 
	  directblocks[i] = 1;

      *testout << "dcmat = " << endl << dcmat << endl;
      *testout << "directblocks = " << endl << directblocks << endl;
      *testout << "freedofs = " << endl << free_dofs << endl;
      */

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




    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const
    {
      // cout << "Multadd, |x| = " << L2Norm(x) << ", |y| = " << L2Norm(y) << endl;
      static int timer = NgProfiler::CreateTimer ("Apply BDDC preconditioner");
      static int timertransx = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transx");
      static int timertransx1a = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transx - lapack");
      static int timertransx1b = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transx - read");
      static int timertransx1c = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transx - write");
      // static int timertransx1d = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transx - indices");

      static int timertransy = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transy");
      static int timertransyr = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transy - read");
      static int timertransyw = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transy - write");
      // static int timertransyi = NgProfiler::CreateTimer ("Apply BDDC preconditioner, transy - ind");


      static int timerrest = NgProfiler::CreateTimer ("Apply BDDC preconditioner, restrict");
      static int timerdc = NgProfiler::CreateTimer ("Apply BDDC preconditioner, decoupled");
      static int timerext = NgProfiler::CreateTimer ("Apply BDDC preconditioner, extend");
      static int timerinner = NgProfiler::CreateTimer ("Apply BDDC preconditioner, inner");
      static int timersolve = NgProfiler::CreateTimer ("Apply BDDC preconditioner, solve");
      // static int timeretc = NgProfiler::CreateTimer ("Apply BDDC preconditioner, etc");
      NgProfiler::RegionTimer reg (timer);


      FlatVector<> fx = x.FVDouble();
      FlatVector<> fy = y.FVDouble();

      VVector<> transx(ndof);
      VVector<> transy(ndof);
      
      int maxmem = 0;
      for (int cl = 0; cl < els_of_class->Size(); cl++)
	{
	  FlatArray<int> els = (*els_of_class)[cl];
	  if (!els.Size()) continue;
	  
	  const Matrix<> & ext = (*extension_int_ref_trans[cl]);
	  int mem = els.Size() * max (ext.Height(), ext.Width());

	  maxmem = max(maxmem,mem);
	}
      
      Vector<> hmem1(maxmem);
      Vector<> hmem2(maxmem);
      

      {
	NgProfiler::RegionTimer reg(timertransx);
	transx.FV() = fx;  



	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & ext = (*extension_int_ref_trans[cl]);

	    FlatMatrix<> vim(els.Size(), ext.Width(), &hmem1(0));
	    FlatMatrix<> vem(els.Size(), ext.Height(), &hmem2(0));

	    NgProfiler::AddFlops (timertransx, ext.Width()*ext.Height()*els.Size());
	    NgProfiler::AddFlops (timertransx1a, ext.Width()*ext.Height()*els.Size());

	    NgProfiler::StartTimer(timertransx1b);

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> dint = (*globinternaldofs)[els[ii]];
		FlatVector<> vim_row = vim.Row(ii);

		for (int j = 0; j < dint.Size(); j++)
		  vim_row(j) = transx(dint[j]);
	      }

	    NgProfiler::StopTimer(timertransx1b);


	    NgProfiler::AddLoads (timertransx1b, els.Size()*vim.Width()*12);
	    NgProfiler::AddStores (timertransx1b, els.Size()*vim.Width()*8);



	    NgProfiler::StartTimer(timertransx1a);

	    LapackMultABt (vim, ext, vem);

	    NgProfiler::StopTimer(timertransx1a);

	    NgProfiler::StartTimer(timertransx1c);
	    for (int i = 0, ii = 0; i < els.Size(); i++)
	      {
		FlatArray<int> dext = (*globexternaldofs)[els[i]];

		for (int j = 0; j < dext.Size(); j++, ii++)
		  transx(dext[j]) -= vem(ii);
	      }

	    NgProfiler::StopTimer(timertransx1c);









	    NgProfiler::AddLoads (timertransx1c, els.Size()*vem.Width()*20);
	    NgProfiler::AddStores (timertransx1c, els.Size()*vem.Width()*8);
	  }
      }


      // restrict
      { 
	NgProfiler::RegionTimer reg(timerrest);

	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & ext = (*extwb_ref[cl]);
	    Matrix<> v1(els.Size(), ext.Height());
	    Matrix<> v2(els.Size(), ext.Width());
	    
	    NgProfiler::AddFlops (timerrest, v1.Width()*v2.Width()*els.Size());
		
	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> ifs = (*interface_dofs)[els[ii]];
		FlatArray<double> w = (*weighting)[els[ii]];

		for (int j = 0; j < ifs.Size(); j++)
		  v1(ii,j) = w[j] * transx(ifs[j]);
	      }

	    LapackMultAB (v1, ext, v2);

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> wb = (*globwbdofs)[els[ii]];
		
		for (int j = 0; j < wb.Size(); j++)
		  transx(wb[j]) -= v2(ii, j);
	      }
	}
      }


      transy = 0.0;

  

      {
	// cout << "|lx| = " << L2Norm(lx) << flush;
	NgProfiler::RegionTimer reg (timersolve);

	BaseVector & subx = *(transx.Range(0, inv->Height()));
	BaseVector & suby = *(transy.Range(0, inv->Height()));
	VVector<> res(inv->Height());

	if (block)
	  {
	    if (false) //GS
	      {
		// cout << "start GS" << endl;
		dynamic_cast<BaseBlockJacobiPrecond*>(inv)->GSSmoothResiduum (suby, subx, res,1);
	    // cout << "GS done" << endl;

		if (inv_coarse)
		  suby += (*inv_coarse) * res; 

		// cout << "inv_coarse done" << endl;
		dynamic_cast<BaseBlockJacobiPrecond*>(inv)->GSSmoothBack (suby, subx);
	      }
	    else{ //jacobi only (old)
	      suby = (*inv) * subx;
	      suby += (*inv_coarse) * subx; 
	  }
	}
	else
	  {
	    suby = (*inv) * subx;
	  }


	// cout << ", |ly| = " << L2Norm(ly) << endl;
      }
      

      // extend
      {
	NgProfiler::RegionTimer reg(timerext);

	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & ext = (*extwb_ref[cl]);
	    Matrix<> v1(els.Size(), ext.Width());
	    Matrix<> v2(els.Size(), ext.Height());
	    
	    NgProfiler::AddFlops (timerext, v1.Width()*v2.Width()*els.Size());
		
	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> wb = (*globwbdofs)[els[ii]];
		for (int j = 0; j < wb.Size(); j++)
		  v1(ii, j) = transy(wb[j]);
	      }

	    LapackMultABt (v1, ext, v2);

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> ifs = (*interface_dofs)[els[ii]];
		FlatArray<double> w = (*weighting)[els[ii]];

		for (int j = 0; j < ifs.Size(); j++)
		  transy(ifs[j]) -= w[j] * v2(ii,j);
	      }
	}
      }

 


    // solve decoupled
      {
	NgProfiler::RegionTimer reg(timerdc);

	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & invdc = (*invdc_ref[cl]);
	    Matrix<> dc1(els.Size(), invdc.Height());
	    Matrix<> dc2(els.Size(), invdc.Height());

	    NgProfiler::AddFlops (timerdc, invdc.Width()*invdc.Height()*els.Size());
		
	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> ifs = (*interface_dofs)[els[ii]];
		FlatArray<double> w = (*weighting)[els[ii]];

		for (int j = 0; j < ifs.Size(); j++)
		  dc1(ii,j) = w[j] * transx(ifs[j]);
	      }

	    LapackMultAB (dc1, invdc, dc2);

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> ifs = (*interface_dofs)[els[ii]];
		FlatArray<double> w = (*weighting)[els[ii]];

		for (int j = 0; j < ifs.Size(); j++)
		  transy(ifs[j]) += w[j] * dc2(ii,j);
	      }
	  }
      }





      {
	NgProfiler::RegionTimer reg(timertransy);

	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & ext = (*extension_int_ref_trans[cl]);
	    Matrix<> vim(els.Size(), ext.Width());
	    Matrix<> vem(els.Size(), ext.Height());

	    NgProfiler::AddFlops (timertransy, ext.Width()*ext.Height()*els.Size());
		


	    NgProfiler::StartTimer (timertransyr);

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> dext = (*globexternaldofs)[els[ii]];
		for (int j = 0; j < dext.Size(); j++)
		  vem(ii, j) = transy(dext[j]);
	      }

	    NgProfiler::StopTimer (timertransyr);


	    LapackMultAB (vem, ext, vim);


	    NgProfiler::StartTimer (timertransyw);

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> dint = (*globinternaldofs)[els[ii]];
		
		for (int j = 0; j < dint.Size(); j++)
		  transy(dint[j]) -= vim(ii, j);
	      }

	    NgProfiler::StopTimer (timertransyw);

	  }
      }
      

      {
	NgProfiler::RegionTimer reg(timerinner);

	for (int cl = 0; cl < els_of_class->Size(); cl++)
	  {
	    FlatArray<int> els = (*els_of_class)[cl];
	    if (!els.Size()) continue;

	    const Matrix<> & inv_int = (*inv_int_ref[cl]);
	    Matrix<> v1(els.Size(), inv_int.Width());
	    Matrix<> v2(els.Size(), inv_int.Width());

	    NgProfiler::AddFlops (timerinner, v1.Width()*v2.Width()*els.Size());
		
	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> dint = (*globinternaldofs)[els[ii]];
		for (int j = 0; j < dint.Size(); j++)
		  v1(ii, j) = transx(dint[j]);
	      }

	    LapackMultAB (v1, inv_int, v2);

	    for (int ii = 0; ii < els.Size(); ii++)
	      {
		FlatArray<int> dint = (*globinternaldofs)[els[ii]];
		for (int j = 0; j < dint.Size(); j++)
		  transy(dint[j]) += v2(ii, j);
	      }

	  }
      }

      y += s * transy;
    }

  };













  template <class SCAL>
  BDDCPreconditioner<SCAL> ::
  BDDCPreconditioner (const PDE & pde, const Flags & aflags) // , const string aname)
    : Preconditioner (&pde, aflags)
  {
    bfa = dynamic_cast<const S_BilinearForm<SCAL>*>(pde.GetBilinearForm (aflags.GetStringFlag ("bilinearform", NULL)));
    const_cast<S_BilinearForm<SCAL>*> (bfa) -> SetPreconditioner (this);
    inversetype = flags.GetStringFlag("inverse", "sparsecholesky");
    refelement = flags.GetDefineFlag("refelement");
    block = flags.GetDefineFlag("block");
    ebe = flags.GetDefineFlag("ebe");
  }


  template <class SCAL>
  void BDDCPreconditioner<SCAL> ::
  AddElementMatrix (const Array<int> & dnums,
		    const FlatMatrix<> & elmat,
		    bool inner_element, int elnr,
		    LocalHeap & lh)
  {
#pragma omp critical (addbddcelmat)
    {
      if (elmats.Size() < ma.GetNE())
	elmats.SetSize (ma.GetNE());
    }
    
    int used = 0;
    for (int i = 0; i < dnums.Size(); i++)
      if (dnums[i] != -1) used++;
     
    elmats[elnr] = new Matrix<double> (used);

    for (int i = 0, ii = 0; i < dnums.Size(); i++)
      if (dnums[i] != -1) 
	{
	  for (int j = 0, jj = 0; j < dnums.Size(); j++)
	    if (dnums[j] != -1) 
	      {
		(*elmats[elnr])(ii,jj) = elmat(i,j);
		jj++;
	      }
	  ii++;
	}
  }
  


  template <class SCAL>
  void BDDCPreconditioner<SCAL> ::
  Update ()
  {
    if (id == 0)
      cout << "update bddc, inversetype = " << inversetype << endl;

    if (refelement)
      pre = new BDDCMatrixRefElement(*bfa, inversetype);
    else
      pre = new BDDCMatrix(*bfa, elmats, inversetype, block,ebe);

    if (test) Test();
  }  

  template class BDDCPreconditioner<double>;
    


  static RegisterPreconditioner<BDDCPreconditioner<double> > initpre ("bddc");
}
