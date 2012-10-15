#include <comp.hpp>
#include <multigrid.hpp>

#include <parallelngs.hpp>
#include <stdlib.h>

namespace ngcomp
{
  using namespace ngmg;


  GridFunction :: GridFunction (const FESpace & afespace, const string & name,
				const Flags & flags)
    : NGS_Object (afespace.GetMeshAccess(), name), fespace(afespace)
  { 
    nested = flags.GetDefineFlag ("nested");
    visual = !flags.GetDefineFlag ("novisual");
    multidim = int (flags.GetNumFlag ("multidim", 1));
    level_updated = -1;
    cacheblocksize = 1;
  }

  GridFunction :: ~GridFunction()
  {
    for (int i = 0; i < vec.Size(); i++)
      delete vec[i];
  }


  void GridFunction :: Update ()
  {
    throw Exception("GridFunction::Update not overloaded");
  }

  bool GridFunction :: IsUpdated () const
  {
    int ndof = GetFESpace().GetNDof();
    for (int i = 0; i < multidim; i++)
      {
	if (!vec[i]) return false;
	if (ndof != vec[i]->Size()) return false;
      }
    return true;
  }

  void GridFunction :: PrintReport (ostream & ost)
  {
    ost << "on space " << GetFESpace().GetName() << endl
	<< "nested = " << nested << endl;
  }

  void GridFunction :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  {
    //if (&const_cast<GridFunction&> (*this).GetVector())
    if (&(*this).GetVector())
      {
	int olds = mu.Size();
	//const_cast<GridFunction&> (*this).GetVector().MemoryUsage (mu);
	(*this).GetVector().MemoryUsage (mu);
	for (int i = olds; i < mu.Size(); i++)
	  mu[i]->AddName (string(" gf ")+GetName());
      }
  }





  void GridFunction :: Visualize(const string & given_name)
  {
    if (!visual) return;

    /*
    for (int i = 0; i < compgfs.Size(); i++)
      {
	stringstream child_name;
	child_name << given_name << "_" << i+1;
	compgfs[i] -> Visualize (child_name.str());
      }
    */

    const BilinearFormIntegrator * bfi2d = 0, * bfi3d = 0;

    if (ma.GetDimension() == 2)
      {
	bfi2d = fespace.GetIntegrator();
      }
    else
      {
	bfi3d = fespace.GetIntegrator();
	bfi2d = fespace.GetBoundaryIntegrator();
      }

    if (bfi2d || bfi3d)
      {
	// if (ma.GetNLevels() > 1) return;
        // if (!vis)

	netgen::SolutionData * vis;
	if (!fespace.IsComplex())
	  vis = new VisualizeGridFunction<double> (ma, this, bfi2d, bfi3d, 0);
	else
	  vis = new VisualizeGridFunction<Complex> (ma, this, bfi2d, bfi3d, 0);

	Ng_SolutionData soldata;
	Ng_InitSolutionData (&soldata);
	
	soldata.name = given_name.c_str();
	soldata.data = 0;
	soldata.components = vis->GetComponents();
	soldata.iscomplex = vis->IsComplex();
	soldata.draw_surface = bfi2d != 0;
	soldata.draw_volume  = bfi3d != 0;
	soldata.dist = soldata.components;
	soldata.soltype = NG_SOLUTION_VIRTUAL_FUNCTION;
	soldata.solclass = vis;
	Ng_SetSolutionData (&soldata);    
      }
  }



  GridFunction * GridFunction :: 
  GetComponent (int compound_comp) const
  {
    if (!compgfs.Size() || compgfs.Size() < compound_comp)
      throw Exception("GetComponent: compound_comp does not exist!");
    else
      return compgfs[compound_comp];
  }







template <int N>
bool Less( Vec<N, int>& a, Vec<N, int>& b)
{
	bool ret = false;
	for( int i = 0; i < N; i++)
	{
		if( a[i] < b[i])
		{
			ret = true;
			break;
		}
		else if(a[i] > b[i])
			break;
	}
	
	return ret;  
}



template <int N>
void Quicksort(  Array<int>& pos, Array<Vec<N, int> >& vals,int  first, int last);

template <int N>
int Divide(  Array<int>& pos, Array<Vec<N, int> >& vals, int first,int last);



template <int N>
void Quicksort( Array<int>& pos, Array<Vec<N, int> >& vals,int first,int last)
{
  if( first < last)
  {
    int mid = Divide<N> ( pos, vals, first, last);
    Quicksort<N>(  pos,  vals, first, mid-1);
    Quicksort<N>(  pos,  vals, mid+1, last);
  }
}



template <int N>
int Divide( Array<int>& pos, Array<Vec<N, int> >& vals,int first,int last)
{
  int i = first;
  int j = last -1;
  
  Vec<N, int> pivot = vals[  pos[last]  ];
  
  while(i < j)
  {
    while( (i < last) && (  !Less<N>( pivot, vals[ pos[i] ])   )) //pivot >= vals[ pos[i] ]
			i++;
    
    while((j > first) && ( !Less<N>(vals[ pos[j] ], pivot)   ) )  // pivot <= vals[ pos[j] ]
			j--;
    
    if( i < j)
    {
      int temp = pos[i];
      pos[i] = pos[j];
      pos[j] = temp;
    }
    
  }
  
  
  if(  Less<N>(pivot,vals[ pos[i] ]) )
  {
    int temp = pos[i];
    pos[i] = pos[last];
    pos[last] = temp;
  }
  
  return i;
  
}




  
  template <class SCAL>
  void S_GridFunction<SCAL> :: Load (istream & ist)
  {
    int ntasks = MyMPI_GetNTasks();
  
    if (ntasks==1)
      { 
	const FESpace & fes = GetFESpace();
	Array<int> dnums;
	
	for (NODE_TYPE nt = NT_VERTEX; nt <= NT_CELL; nt++)
	  {
	    int nnodes = ma.GetNNodes (nt);
	    for( int i = 0; i < nnodes; i++)
	      {
		fes.GetNodeDofNrs (nt, i,  dnums); 
		Vector<SCAL> elvec(dnums.Size());
		
		for (int i=0; i< elvec.Size(); i++)
		  if (ist.good())
		    LoadBin<SCAL>(ist, elvec(i));			
	  
		SetElementVector (dnums, elvec);
	      }
	  }
      }
    else
      {  
#ifdef PARALLEL	    
	GetVector().SetParallelStatus (DISTRIBUTED);    
	LoadNodeType<1,NT_VERTEX> (ist);
	LoadNodeType<2,NT_EDGE> (ist);
	LoadNodeType<4,NT_FACE> (ist);
	LoadNodeType<8,NT_CELL> (ist);
	GetVector().Cumulate();
      }
#endif
  }

  template <class SCAL>  template <int N, NODE_TYPE NTYPE>
  void  S_GridFunction<SCAL> :: LoadNodeType (istream & ist) 
  {
#ifdef PARALLEL
    int id = MyMPI_GetId();
    int ntasks = MyMPI_GetNTasks();
    
    const FESpace & fes = GetFESpace();
    ParallelDofs & par = fes.GetParallelDofs ();
    
    if(id > 0)
      {
	int nnodes = ma.GetNNodes (NTYPE);
	
	Array<Vec<N+1, int> > nodenums(0);
	Array<int> master_nodes(0);
	Array<int> dnums;
	
	for( int i = 0; i < nnodes; i++)
	  {
	    fes.GetNodeDofNrs (NTYPE, i,  dnums);
	    
	    if( dnums.Size() == 0)
	      continue;
	    
	    if( !par.IsMasterDof ( dnums[0] ))
	      continue;
	    
	    master_nodes.Append(i);
	    
	    Array<int> pnums;
	    if (NTYPE == NT_VERTEX)
	      {
		pnums.SetSize(1);
		pnums[0]=i;
	      }
	    else if( NTYPE == NT_EDGE)
	      ma.GetEdgePNums (i, pnums);
	    else if (NTYPE == NT_FACE)
	      ma.GetFacePNums (i, pnums);
	    else if (NTYPE == NT_CELL)
	      ma.GetElPNums (i, pnums);
	    else
	      cout<<"Error in SaveSolution Node Type not known"<<endl;
	    
	    Vec<N+1, int> points;
	    points = -1;
	    for( int j = 0; j < pnums.Size(); j++)
	      points[j] = ma.GetGlobalNodeNum (Node(NT_VERTEX, pnums[j]));
	    points[N] = dnums.Size();
	    
	    nodenums.Append(points);	
	  }
	
	MyMPI_Send(nodenums,0,12);
	
	Array<SCAL> loc_data;
	MyMPI_Recv(loc_data,0,13);
	
	int cnt=0;
	for( int i = 0; i < master_nodes.Size(); i++)
	  {
	    fes.GetNodeDofNrs (NTYPE, master_nodes[i],  dnums); 
	    Vector<SCAL> elvec(dnums.Size());
	    
	    for (int i=0; i< elvec.Size(); i++)
	      elvec(i)=loc_data[cnt++];
	    
	    SetElementVector (dnums, elvec);
	  }
      }
    else
      {
	Array<Vec<N, int> > points(0);
	Array<int>  nrdofs_per_node(0);
	
	Array<Array<int>* > nodenums_of_procs (ntasks-1);
	
	int actual=0;
	for( int proc = 1; proc < ntasks; proc++)
	  {
	    Array<Vec<N+1, int> > nodenums_proc;
	    MyMPI_Recv(nodenums_proc,proc,12);
	    nodenums_of_procs[proc-1] = new Array<int> (nodenums_proc.Size());
	    
	    Vec<N, int>  temp;
	    
	    for (int j=0; j<nodenums_proc.Size(); j++)
	      {
		for( int k = 0; k < N; k++)
		  temp[k] = nodenums_proc[j][k];
		
		points.Append( temp );
		
		nrdofs_per_node.Append(nodenums_proc[j][N]);
		
		(*(nodenums_of_procs[proc-1]))[j]=actual++;
	      }
	    /*cout << "proc " << proc << endl;
	      cout << "nodenums " << *(nodenums_of_procs[proc-1]) << endl;
	      cout << "points " << points << endl;*/
	  }
	
	Array<int> index(points.Size());
	for( int i = 0; i < index.Size(); i++)
	  index[i] = i;
	
	Quicksort<N>(index,points,0, points.Size()-1);
	
	Array<int> inverse_index(index.Size());
	for (int i = 0; i < index.Size(); i++ ) 	
	  inverse_index[index[i]] = i;
	
	
	/*cout << "points " << points << endl;      
	  cout << "index " << index << endl;      
	  cout << "inverse index " << inverse_index << endl;*/
	
	int nnodes=points.Size();
	int nrdofs_all_nodes=0;
	Array<int> first_node_dof (nnodes+1);
	
	for (int i=0; i<nnodes; i++)
	  {
	    first_node_dof[i]=nrdofs_all_nodes;
	    nrdofs_all_nodes+=nrdofs_per_node[index[i]];
	  }
	first_node_dof[nnodes]=nrdofs_all_nodes;
	
      //cout << "first_node_dof " << first_node_dof << endl;
	
	Array<SCAL> node_data (nrdofs_all_nodes);     
	for (int i=0; i<nrdofs_all_nodes; i++)
	  if (ist.good())
	    LoadBin<SCAL>(ist, node_data[i]);
	//cout << "node_data " << node_data << endl;
	
	for( int proc = 1; proc < ntasks; proc++)
	  {
	    Array<SCAL> loc_data (0);
	    int nr_local_nodes= (*(nodenums_of_procs[proc-1])).Size();	
	    //cout << "proc " << proc << " nr_local_nodes " << nr_local_nodes << endl;
	    for (int i=0; i<nr_local_nodes; i++)
	      {
		int node=inverse_index[(*(nodenums_of_procs[proc-1]))[i]];
		//cout << "i " << i << " inverse i " << inverse_index[i] << " node " << node << endl;
		int first = first_node_dof[node];
		int last = first_node_dof[node+1];
		for (int j=first; j<last; j++)
		  loc_data.Append(node_data[j]);
	      }
	    MyMPI_Send(loc_data,proc,13); 		
	    //	cout << "proc " << proc << " loc_data " << loc_data << endl;
	  }
      }
#endif	
  }
  


  template <class SCAL>
  void S_GridFunction<SCAL> :: Save (ostream & ost) const
  {
    
    int ntasks = MyMPI_GetNTasks();
    const FESpace & fes = GetFESpace();
  
  
    if (ntasks==1)
      {
	Array<int> dnums;	    
	for (NODE_TYPE nt = NT_VERTEX; nt <= NT_CELL; nt++)
	  {
	    int nnodes = ma.GetNNodes (nt);
	    for( int i = 0; i < nnodes; i++)
	      {
		fes.GetNodeDofNrs (nt, i,  dnums); 
		Vector<SCAL> elvec(dnums.Size());
		GetElementVector (dnums, elvec);
		
		for (int i=0; i< elvec.Size(); i++)
		  SaveBin<SCAL>(ost, elvec(i));			
	      }
	  }
      }
    else
      {  
#ifdef PARALLEL	 
	GetVector().Cumulate();        
        
	SaveNodeType<1,NT_VERTEX>(ost);
	SaveNodeType<2,NT_EDGE>  (ost);
	SaveNodeType<4,NT_FACE>  (ost);
	SaveNodeType<8,NT_CELL>  (ost);
#endif
      }
  }

  template <class SCAL>  template <int N, NODE_TYPE NTYPE>
  void  S_GridFunction<SCAL> :: SaveNodeType (ostream & ost) const
  {
#ifdef PARALLEL
    int id = MyMPI_GetId();
    int ntasks = MyMPI_GetNTasks();
    
    const FESpace & fes = GetFESpace();
    ParallelDofs & par = fes.GetParallelDofs ();
    
    if(id > 0)
      { 
	int nnodes = ma.GetNNodes (NTYPE);
	
	Array<Vec<N+1, int> > nodenums(0);
	Array<SCAL> data(0);
    
	Array<int> dnums;
    
	for( int i = 0; i < nnodes; i++)
	  {
	    fes.GetNodeDofNrs (NTYPE, i,  dnums);
	    
	    if( dnums.Size() == 0)
	      continue;
	    
	    if( !par.IsMasterDof ( dnums[0] ))
	      continue;      
	    
	    Array<int> pnums;
	    if (NTYPE == NT_VERTEX)
	      {
		pnums.SetSize(1);
		pnums[0]=i;
	      }
	    else if( NTYPE == NT_EDGE)
	      ma.GetEdgePNums (i, pnums);
	    else if (NTYPE == NT_FACE)
	      ma.GetFacePNums (i, pnums);
	    else if (NTYPE == NT_CELL)
	      ma.GetElPNums (i, pnums);
	    else
	      cout<<"Error in SaveSolution Node Type not known"<<endl;
	    
	    Vec<N+1, int> points;
	    points = -1;
	    for( int j = 0; j < pnums.Size(); j++)
	      points[j] = ma.GetGlobalNodeNum (Node(NT_VERTEX, pnums[j]));
	    points[N] = dnums.Size();
	    
	    nodenums.Append(points);
	    
	    Vector<SCAL> elvec(dnums.Size());
	    GetElementVector (dnums, elvec);
	    
	    for( int j = 0; j < dnums.Size(); j++)
	      data.Append(elvec(j));
	  }    
	
	MyMPI_Send(nodenums,0,22);
	MyMPI_Send(data,0,23);
      }
  
    if( id == 0 )
      {
	Array<Vec<N, int> > points(0);
	Array<Vec<2, int> > positions(0);
	Array<SCAL> data(0);
	
	int size = 0;
	for( int proc = 1; proc < ntasks; proc++)
	  {
	    Array<Vec<N+1, int> > locpoints;
	    Array<SCAL> locdata;
	    MyMPI_Recv(locpoints,proc, 22);
	    MyMPI_Recv(locdata,proc, 23);
	    
	    
	    Vec<N, int>  temp;
	    
	    int actual = 0;
	    for( int j = 0; j < locpoints.Size(); j++ )
	      {
		int nodesize = locpoints[j][N];
		for( int k = 0; k < nodesize; k++)
		  data.Append(locdata[actual++]);
		
		positions.Append(  Vec<2, int> (  size, nodesize) );
		
		for( int k = 0; k < N; k++)
		  temp[k] = locpoints[j][k];
		
		points.Append( temp );
		
		size += nodesize;
	      }
	  }    
	
	Array<int> index(points.Size());
	for( int i = 0; i < index.Size(); i++)
	  index[i] = i;
	
	Quicksort<N>(index,points,0, points.Size()-1);
	
	for( int i = 0; i < points.Size(); i++)
	  {
	    int start = positions[index[i]][0];
	    int end = positions[index[i]][1];
	    
	    for( int j = 0; j < end; j++)
	      SaveBin<SCAL>(ost, data[start++]);
	  }
	
      }
#endif	   
  }







  template <class SCAL>
  S_ComponentGridFunction<SCAL> :: 
  S_ComponentGridFunction (const S_GridFunction<SCAL> & agf_parent, int acomp)
    : S_GridFunction<SCAL> (*dynamic_cast<const CompoundFESpace&> (agf_parent.GetFESpace())[acomp], 
			    agf_parent.GetName()+"."+ToString (acomp+1), Flags()), 
      gf_parent(agf_parent), comp(acomp)
  { 
    const CompoundFESpace * cfe = dynamic_cast<const CompoundFESpace *>(&this->GetFESpace());
    if (cfe)
      {
	int nsp = cfe->GetNSpaces();
	this->compgfs.SetSize(nsp);
	for (int i = 0; i < nsp; i++)
	  this->compgfs[i] = new S_ComponentGridFunction<SCAL> (*this, i);
      }    
    this->Visualize (this->name);
  }

  template <class SCAL>
  S_ComponentGridFunction<SCAL> :: 
  ~S_ComponentGridFunction ()
  {
    this -> vec = NULL;  // base-class desctructor must not delete the vector
  }


  template <class SCAL>
  void S_ComponentGridFunction<SCAL> :: Update()
  {
    const CompoundFESpace & cfes = dynamic_cast<const CompoundFESpace&> (gf_parent.GetFESpace());

    this -> vec.SetSize (gf_parent.GetMultiDim());
    for (int i = 0; i < gf_parent.GetMultiDim(); i++)
      (this->vec)[i] = gf_parent.GetVector(i).Range (cfes.GetRange(comp));
  
    this -> level_updated = this -> ma.GetNLevels();

    for (int i = 0; i < this->compgfs.Size(); i++)
      this->compgfs[i]->Update();
  }



  template <class TV>
  T_GridFunction<TV> ::
  T_GridFunction (const FESpace & afespace, const string & aname, const Flags & flags)
    : S_GridFunction<TSCAL> (afespace, aname, flags)
  {
    vec.SetSize (this->multidim);
    vec = 0;

    const CompoundFESpace * cfe = dynamic_cast<const CompoundFESpace *>(&this->GetFESpace());
    if (cfe)
      {
	int nsp = cfe->GetNSpaces();
	compgfs.SetSize(nsp);
	for (int i = 0; i < nsp; i++)
	  compgfs[i] = new S_ComponentGridFunction<TSCAL> (*this, i);
      }    

    this->Visualize (this->name);
  }

  template <class TV>
  T_GridFunction<TV> :: ~T_GridFunction()
  {
    ;
  }


  template <class TV>
  void T_GridFunction<TV> :: Update () 
  {
    try
      {
	int ndof = this->GetFESpace().GetNDof();

	for (int i = 0; i < this->multidim; i++)
	  {
	    if (vec[i] && ndof == vec[i]->Size())
	      break;
	    
	    BaseVector * ovec = vec[i];
	
#ifdef PARALLEL
	    if ( & this->GetFESpace().GetParallelDofs() )
	      vec[i] = new ParallelVVector<TV> (ndof, &this->GetFESpace().GetParallelDofs(), CUMULATED);
	    else
#endif
 	      vec[i] = new VVector<TV> (ndof);
	    

	    (*vec[i]) = TSCAL(0);

	    if (this->nested && ovec && this->GetFESpace().GetProlongation())
	      {
		*vec[i]->Range (0, ovec->Size()) += (*ovec);

		const_cast<ngmg::Prolongation&> (*this->GetFESpace().GetProlongation()).Update();
		
		this->GetFESpace().GetProlongation()->ProlongateInline
		  (this->GetMeshAccess().GetNLevels()-1, *vec[i]);
	      }

	    //	    if (i == 0)
            // cout << "visualize" << endl;
            // Visualize (this->name);
	    
	    delete ovec;
	  }
	
	this -> level_updated = this -> ma.GetNLevels();

	// const CompoundFESpace * cfe = dynamic_cast<const CompoundFESpace *>(&GridFunction :: GetFESpace());
	// if (cfe)
	  for (int i = 0; i < compgfs.Size(); i++)
	    compgfs[i]->Update();
      }
    catch (exception & e)
      {
	Exception e2 (e.what());
	e2.Append ("\nIn GridFunction::Update()\n");
	throw e2;
      }
    catch (Exception & e)
      {
	e.Append ("In GridFunction::Update()\n");
	throw e;
      }    
  }



  GridFunction * CreateGridFunction (const FESpace * space,
				     const string & name, const Flags & flags)
  {
    GridFunction * gf = 
      CreateVecObject  <T_GridFunction, GridFunction, const FESpace, const string, const Flags>
      (space->GetDimension() * int(flags.GetNumFlag("cacheblocksize",1)), 
       space->IsComplex(), *space, name, flags);
  
    gf->SetCacheBlockSize(int(flags.GetNumFlag("cacheblocksize",1)));

    return gf;
  }












  GridFunctionCoefficientFunction :: 
  GridFunctionCoefficientFunction (GridFunction & agf, int acomp)
    : gf(agf), diffop (NULL), comp (acomp) 
  {
    diffop = gf.GetFESpace().GetEvaluator();
  }

  GridFunctionCoefficientFunction :: 
  GridFunctionCoefficientFunction (GridFunction & agf, 
				   DifferentialOperator * adiffop, int acomp)
    : gf(agf), diffop (adiffop), comp (acomp) 
  {
    ;
  }

  GridFunctionCoefficientFunction :: 
  ~GridFunctionCoefficientFunction ()
  {
    ;
  }

  int GridFunctionCoefficientFunction::Dimension() const
  { 
    if (diffop) return diffop->Dim();
    return gf.GetFESpace().GetIntegrator()->DimFlux();
  }

  bool GridFunctionCoefficientFunction::IsComplex() const
  { 
    return gf.GetFESpace().IsComplex(); 
  }


  double GridFunctionCoefficientFunction :: 
  Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    VectorMem<10> flux(Dimension());
    Evaluate (ip, flux);
    return flux(0);
  }

  void GridFunctionCoefficientFunction :: 
  Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> result) const
  {
    LocalHeapMem<100000> lh2 ("GridFunctionCoefficientFunction, Eval 2");
    static Timer timer ("GFCoeffFunc::Eval-scal");
    RegionTimer reg (timer);

    
    const int elnr = ip.GetTransformation().GetElementNr();
    bool boundary = ip.GetTransformation().Boundary();

    const FESpace & fes = gf.GetFESpace();

    if (!boundary)    
      if (!fes.DefinedOn (ip.GetTransformation().GetElementIndex()))
	{ result = 0.0; return;};
    
    const FiniteElement & fel = fes.GetFE (elnr, boundary, lh2);
    int dim = fes.GetDimension();
    
    ArrayMem<int, 50> dnums;
    fes.GetDofNrs (elnr, boundary, dnums);
    
    VectorMem<50> elu(dnums.Size()*dim);

    gf.GetElementVector (comp, dnums, elu);
    fes.TransformVec (elnr, boundary, elu, TRANSFORM_SOL);

    if (diffop)
      diffop->Apply (fel, ip, elu, result, lh2);
    else
      fes.GetIntegrator(boundary) -> CalcFlux (fel, ip, elu, result, false, lh2);
  }

  void GridFunctionCoefficientFunction :: 
  Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<Complex> result) const
  {
    LocalHeapMem<100000> lh2 ("GridFunctionCoefficientFunction, Eval complex");
    static Timer timer ("GFCoeffFunc::Eval-scal");
    RegionTimer reg (timer);

    
    const int elnr = ip.GetTransformation().GetElementNr();
    bool boundary = ip.GetTransformation().Boundary();

    const FESpace & fes = gf.GetFESpace();

    if (!boundary)    
      if (!fes.DefinedOn (ip.GetTransformation().GetElementIndex()))
	{ result = 0.0; return;};
    
    const FiniteElement & fel = fes.GetFE (elnr, boundary, lh2);
    int dim = fes.GetDimension();
    
    ArrayMem<int, 50> dnums;
    fes.GetDofNrs (elnr, boundary, dnums);
    
    VectorMem<50, Complex> elu(dnums.Size()*dim);

    gf.GetElementVector (comp, dnums, elu);
    fes.TransformVec (elnr, boundary, elu, TRANSFORM_SOL);

    if (diffop)
      diffop->Apply (fel, ip, elu, result, lh2);
    else
      fes.GetIntegrator(boundary) -> CalcFlux (fel, ip, elu, result, false, lh2);
  }


  void GridFunctionCoefficientFunction :: 
  Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<double> values) const
  {
    LocalHeapMem<100000> lh2("GridFunctionCoefficientFunction - Evalute 3");
    static Timer timer ("GFCoeffFunc::Eval-vec");
    RegionTimer reg (timer);

    
    const int elnr = ir.GetTransformation().GetElementNr();
    bool boundary = ir.GetTransformation().Boundary();

    const FESpace & fes = gf.GetFESpace();
    
    if (!boundary)
      if (!fes.DefinedOn(ir.GetTransformation().GetElementIndex())) 
	{ values = 0.0; return;};
    
    const FiniteElement & fel = fes.GetFE (elnr, boundary, lh2);
    int dim = fes.GetDimension();

    ArrayMem<int, 50> dnums;
    fes.GetDofNrs (elnr, boundary, dnums);
     
    VectorMem<50> elu(dnums.Size()*dim);

    gf.GetElementVector (comp, dnums, elu);
    fes.TransformVec (elnr, boundary, elu, TRANSFORM_SOL);
    
    if (diffop)
      diffop->Apply (fel, ir, elu, values, lh2);
    else
      fes.GetIntegrator(boundary) ->CalcFlux (fel, ir, elu, values, false, lh2);
  }





  template <class SCAL>
  VisualizeGridFunction<SCAL> ::
  VisualizeGridFunction (const MeshAccess & ama,
			 const GridFunction * agf,
			 const BilinearFormIntegrator * abfi2d,
			 const BilinearFormIntegrator * abfi3d,
			 bool aapplyd)

    : SolutionData (agf->GetName(), -1, agf->GetFESpace().IsComplex()),
      ma(ama), gf(dynamic_cast<const S_GridFunction<SCAL>*> (agf)), 
      applyd(aapplyd), cache_elnr(-1), lh(10000013, "VisualizedGridFunction 2"), fel(NULL)
  { 
    if(abfi2d)
      bfi2d.Append(abfi2d);
    if(abfi3d)
      bfi3d.Append(abfi3d);

    if (abfi2d) components = abfi2d->DimFlux();
    if (abfi3d) components = abfi3d->DimFlux();
    if (iscomplex) components *= 2;
    multidimcomponent = 0;
  }

  template <class SCAL>
  VisualizeGridFunction<SCAL> ::
  VisualizeGridFunction (const MeshAccess & ama,
			 const GridFunction * agf,
			 const Array<BilinearFormIntegrator *> & abfi2d,
			 const Array<BilinearFormIntegrator *> & abfi3d,
			 bool aapplyd)

    : SolutionData (agf->GetName(), -1, agf->GetFESpace().IsComplex()),
      ma(ama), gf(dynamic_cast<const S_GridFunction<SCAL>*> (agf)), 
      applyd(aapplyd), cache_elnr(-1), lh(10000002, "VisualizeGridFunction"), fel(NULL)
  { 
    for(int i=0; i<abfi2d.Size(); i++)
      bfi2d.Append(abfi2d[i]);
    for(int i=0; i<abfi3d.Size(); i++)
      bfi3d.Append(abfi3d[i]);
    

    if (bfi2d.Size()) components = bfi2d[0]->DimFlux();
    if (bfi3d.Size()) components = bfi3d[0]->DimFlux();

    if (iscomplex) components *= 2;
    multidimcomponent = 0;
  }

  template <class SCAL>
  VisualizeGridFunction<SCAL> :: ~VisualizeGridFunction ()
  {}
  

  template <class SCAL>
  bool VisualizeGridFunction<SCAL> :: GetValue (int elnr, 
						double lam1, double lam2, double lam3,
						double * values) 
  { 
    // cout << "VisGF::GetValue" << endl;
    if (!bfi3d.Size()) return 0;
    if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;

    const FESpace & fes = gf->GetFESpace();

    int dim     = fes.GetDimension();

    if ( !fes.DefinedOn(ma.GetElIndex(elnr)) ) 
      return false;

    cache_elnr = -1;
    if (cache_elnr != elnr || cache_bound)
      {
	lh.CleanUp();

	fes.GetDofNrs (elnr, dnums);
	elu.AssignMemory (dnums.Size() * dim, lh);

	if(gf->GetCacheBlockSize() == 1)
	  {
	    gf->GetElementVector (multidimcomponent, dnums, elu);
	  }
	else
	  {
	    FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
	    //gf->GetElementVector (multidimcomponent, dnums, elu2);
	    gf->GetElementVector (0, dnums, elu2);
	    for(int i=0; i<elu.Size(); i++)
	      elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
	  }

	fes.TransformVec (elnr, 0, elu, TRANSFORM_SOL);

	cache_elnr = elnr;
	cache_bound = 0;
      }

    HeapReset hr(lh);
    ElementTransformation & eltrans = ma.GetTrafo (elnr, false, lh);
    const FiniteElement & fel = fes.GetFE (elnr, lh);

    IntegrationPoint ip(lam1, lam2, lam3, 0);
    MappedIntegrationPoint<3,3> sip (ip, eltrans);

    for(int j = 0; j < bfi3d.Size(); j++)
      {
	HeapReset hr(lh);
	FlatVector<SCAL> flux(bfi3d[j] -> DimFlux(), lh);
	bfi3d[j]->CalcFlux (fel, sip, elu, flux, applyd, lh);

	for (int i = 0; i < components; i++)
	  {
	    if(j == 0) values[i] = 0;
	    values[i] += ((double*)(void*)&flux(0))[i];
	  }
      }

    return true; 
  }


  template <class SCAL>
  bool VisualizeGridFunction<SCAL> :: GetValue (int elnr, 
						const double xref[], const double x[], const double dxdxref[],
						double * values) 
  { 
    // cout << "VisGF::GetValue2" << endl;
    if (!bfi3d.Size()) return 0;
    if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;


    const FESpace & fes = gf->GetFESpace();

    int dim     = fes.GetDimension();

    if ( !fes.DefinedOn(ma.GetElIndex(elnr)) ) return 0;

    if (cache_elnr != elnr || cache_bound)
      {
	lh.CleanUp();
	fes.GetDofNrs (elnr, dnums);
    
	elu.AssignMemory (dnums.Size() * dim, lh);
	if(gf->GetCacheBlockSize() == 1)
	  {
	    gf->GetElementVector (multidimcomponent, dnums, elu);
	  }
	else
	  {
	    FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
	    //gf->GetElementVector (multidimcomponent, dnums, elu2);
	    gf->GetElementVector (0, dnums, elu2);
	    for(int i=0; i<elu.Size(); i++)
	      elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
	  }
	fes.TransformVec (elnr, 0, elu, TRANSFORM_SOL);

	fel = &fes.GetFE (elnr, lh);

	cache_elnr = elnr;
	cache_bound = 0;
      }

    HeapReset hr(lh);

    Vec<3> vx;
    Mat<3,3> mdxdxref;
    for (int i = 0; i < 3; i++)
      {
	vx(i) = x[i];
	for (int j = 0; j < 3; j++)
	  mdxdxref(i,j) = dxdxref[3*i+j];
      }

    ElementTransformation & eltrans = ma.GetTrafo (elnr, false, lh);
    IntegrationPoint ip(xref[0], xref[1], xref[2], 0);
    MappedIntegrationPoint<3,3> sip (ip, eltrans, vx, mdxdxref);

    for(int j = 0; j < bfi3d.Size(); j++)
      {
	FlatVector<SCAL> flux (bfi3d[j]->DimFlux(), lh);
	bfi3d[j]->CalcFlux (*fel, sip, elu, flux, applyd, lh);

	for (int i = 0; i < components; i++)
	  {
	    if(j == 0) values[i] = 0;
	    values[i] += ((double*)(void*)&flux(0))[i];
	  }
      }

    return true; 
  }



  template <class SCAL>
  bool VisualizeGridFunction<SCAL> ::
  GetMultiValue (int elnr, int npts,
		 const double * xref, int sxref,
		 const double * x, int sx,
		 const double * dxdxref, int sdxdxref,
		 double * values, int svalues)
  {
    // cout << "VisGF::GetMultiValue" << endl;
    try
      {
        if (!bfi3d.Size()) return 0;
        if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;

        const FESpace & fes = gf->GetFESpace();
        int dim = fes.GetDimension();
	
        HeapReset hr(lh);
	ElementTransformation & eltrans = ma.GetTrafo (elnr, false, lh);

	fes.GetDofNrs (elnr, dnums);
	fel = &fes.GetFE (elnr, lh);
        
        elu.AssignMemory (dnums.Size() * dim, lh);

        if(gf->GetCacheBlockSize() == 1)
          {
            gf->GetElementVector (multidimcomponent, dnums, elu);
          }
        else
          {
            FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
            gf->GetElementVector (0, dnums, elu2);
            for(int i=0; i<elu.Size(); i++)
              elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
          }
        
        fes.TransformVec (elnr, false, elu, TRANSFORM_SOL);
        
	if (!fes.DefinedOn(eltrans.GetElementIndex()))return 0;

	IntegrationRule ir; 
	ir.SetAllocSize(npts);
	for (int i = 0; i < npts; i++)
	  ir.Append (IntegrationPoint (xref[i*sxref], xref[i*sxref+1], xref[i*sxref+2]));

	MappedIntegrationRule<3,3> mir(ir, eltrans, 1, lh);

	for (int k = 0; k < npts; k++)
	  {
	    Mat<3,3> & mdxdxref = *new((double*)(dxdxref+k*sdxdxref)) Mat<3,3>;
	    FlatVec<3> vx((double*)x + k*sx);
	    mir[k] = MappedIntegrationPoint<3,3> (ir[k], eltrans, vx, mdxdxref);
	  }

	for (int k = 0; k < npts; k++)
	  for (int i = 0; i < components; i++)
	    values[k*svalues+i] = 0.0;

	for(int j = 0; j < bfi3d.Size(); j++)
	  {
	    HeapReset hr(lh);

	    FlatMatrix<SCAL> flux(npts, bfi3d[j]->DimFlux(), lh);
	    bfi3d[j]->CalcFlux (*fel, mir, elu, flux, applyd, lh);
	    
	    for (int k = 0; k < npts; k++)
	      for (int i = 0; i < components; i++)
		values[k*svalues+i] += ((double*)(void*)&flux(k,0))[i];
          }

        return true;
      }
    catch (Exception & e)
      {
        cout << "GetMultiValue caught exception" << endl
             << e.What();
        return 0;
      }
  }
































  template <class SCAL>
  bool VisualizeGridFunction<SCAL> :: GetSurfValue (int elnr, int facetnr, 
                                                    double lam1, double lam2, 
                                                    double * values) 
  { 
    // cout << "VisGF::GetSurfValue" << endl;
    if (!bfi2d.Size()) return 0;
    if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;

    bool bound = (ma.GetDimension() == 3);
    const FESpace & fes = gf->GetFESpace();


    int dim = fes.GetDimension();

    /*
      if ( bound ? 
      !fes.DefinedOnBoundary(ma.GetSElIndex(elnr)) :
      !fes.DefinedOn(ma.GetElIndex(elnr)) ) return 0;
    */

    if (cache_elnr != elnr || !cache_bound)
      {
	lh.CleanUp();

	fel = &fes.GetFE (elnr, bound, lh);
	fes.GetDofNrs (elnr, bound, dnums);

	elu.AssignMemory (dnums.Size() * dim, lh);

	if(gf->GetCacheBlockSize() == 1)
	  {
	    gf->GetElementVector (multidimcomponent, dnums, elu);
	  }
	else
	  {
	    FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
	    //gf->GetElementVector (multidimcomponent, dnums, elu2);
	    gf->GetElementVector (0, dnums, elu2);
	    for(int i=0; i<elu.Size(); i++)
	      elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
	  }

	fes.TransformVec (elnr, bound, elu, TRANSFORM_SOL);

	cache_elnr = elnr;
	cache_bound = 1;
      }

    HeapReset hr(lh);
    ElementTransformation & eltrans = ma.GetTrafo (elnr, bound, lh);

    if ( bound ? 
	 !fes.DefinedOnBoundary(eltrans.GetElementIndex()) : 
	 !fes.DefinedOn(eltrans.GetElementIndex()) ) return false;

    IntegrationPoint ip(lam1, lam2, 0, 0);
    ip.FacetNr() = facetnr;

    BaseMappedIntegrationPoint & mip = eltrans(ip, lh);
    for(int j = 0; j < bfi2d.Size(); j++)
      {
	FlatVector<SCAL> flux(bfi2d[j]->DimFlux(), lh);
	bfi2d[j]->CalcFlux (*fel, mip, elu, flux, applyd, lh);
	
	for (int i = 0; i < components; i++)
	  {
	    if(j == 0) values[i] = 0;
	    values[i] += ((double*)(void*)&flux(0))[i];
	  }
      }

    return true;
  }




  template <class SCAL>
  bool VisualizeGridFunction<SCAL> :: 
  GetSurfValue (int elnr, int facetnr, 
		const double xref[], const double x[], const double dxdxref[],
		double * values) 
  { 
    // cout << "VisGF::GetSurfValue2" << endl;
    try
      {
        if (!bfi2d.Size()) return 0;
        if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;

        bool bound = (ma.GetDimension() == 3);
        const FESpace & fes = gf->GetFESpace();

        int dim     = fes.GetDimension();

        if (cache_elnr != elnr || !cache_bound)
          {
            lh.CleanUp();

	    fes.GetDofNrs (elnr, bound, dnums);
	    fel = &fes.GetFE (elnr, bound, lh);
            elu.AssignMemory (dnums.Size() * dim, lh);

            if(gf->GetCacheBlockSize() == 1)
              {
                gf->GetElementVector (multidimcomponent, dnums, elu);
              }
            else
              {
                FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
                //gf->GetElementVector (multidimcomponent, dnums, elu2);
                gf->GetElementVector (0, dnums, elu2);
                for(int i=0; i<elu.Size(); i++)
                  elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
              }

            fes.TransformVec (elnr, bound, elu, TRANSFORM_SOL);

            cache_elnr = elnr;
            cache_bound = 1;
          }

	HeapReset hr(lh);
	ElementTransformation & eltrans = ma.GetTrafo (elnr, bound, lh);

        if ( bound ? 
             !fes.DefinedOnBoundary(eltrans.GetElementIndex()) : 
             !fes.DefinedOn(eltrans.GetElementIndex()) ) return 0;

        IntegrationPoint ip(xref[0], xref[1], 0, 0);
	ip.FacetNr() = facetnr;
        if (bound)
          {
            // Vec<3> vx;
            Mat<3,2> mdxdxref;
            for (int i = 0; i < 3; i++)
              {
                // vx(i) = x[i];
                for (int j = 0; j < 2; j++)
                  mdxdxref(i,j) = dxdxref[2*i+j];
              }
            MappedIntegrationPoint<2,3> mip (ip, eltrans, (double*)x, mdxdxref); 
            for (int i = 0; i < components; i++)
              values[i] = 0.0;
            for(int j = 0; j<bfi2d.Size(); j++)
              {
		FlatVector<SCAL> flux(bfi2d[j]->DimFlux(), lh);
                bfi2d[j]->CalcFlux (*fel, mip, elu, flux, applyd, lh);
                for (int i = 0; i < components; i++)
                  values[i] += ((double*)(void*)&flux(0))[i];
              }
          }
        else
          {
            // Vec<2> vx;
            Mat<2,2> mdxdxref;
            for (int i = 0; i < 2; i++)
              {
                // vx(i) = x[i];
                for (int j = 0; j < 2; j++)
                  mdxdxref(i,j) = dxdxref[2*i+j];
              }
            MappedIntegrationPoint<2,2> mip (ip, eltrans, (double*)x, mdxdxref); 

            for (int i = 0; i < components; i++)
              values[i] = 0.0;
            for(int j = 0; j<bfi2d.Size(); j++)
              {
                FlatVector<SCAL> flux(bfi2d[j]->DimFlux(), lh);
                bfi2d[j]->CalcFlux (*fel, mip, elu, flux, applyd, lh);
                for (int i = 0; i < components; i++)
                  values[i] += ((double*)(void*)&flux(0))[i];
              }
          }

        return 1; 
      }
    catch (Exception & e)
      {
        cout << "GetSurfaceValue2 caught exception" << endl
             << e.What();

        return 0;
      }
      
  }






  template <class SCAL>
  bool VisualizeGridFunction<SCAL> ::
  GetMultiSurfValue (int elnr, int facetnr, int npts,
                     const double * xref, int sxref,
                     const double * x, int sx,
                     const double * dxdxref, int sdxdxref,
                     double * values, int svalues)
  {
    // cout << "VisGF::GetMultiSurfValue" << endl;
    try
      {
        if (!bfi2d.Size()) return 0;
        if (gf -> GetLevelUpdated() < ma.GetNLevels()) return 0;

        bool bound = (ma.GetDimension() == 3);
	
        const FESpace & fes = gf->GetFESpace();
        int dim = fes.GetDimension();

        
        HeapReset hr(lh);

	ElementTransformation & eltrans = ma.GetTrafo (elnr, bound, lh);

	fes.GetDofNrs (elnr, bound, dnums);
	fel = &fes.GetFE (elnr, bound, lh);
	/*
        if (bound)
          {
            fes.GetSDofNrs (elnr, dnums);
            fel = &fes.GetSFE (elnr, lh);
          }
        else
          {
            fes.GetDofNrs (elnr, dnums);
            fel = &fes.GetFE (elnr, lh);
          }
	*/
        elu.AssignMemory (dnums.Size() * dim, lh);

        if(gf->GetCacheBlockSize() == 1)
          {
            gf->GetElementVector (multidimcomponent, dnums, elu);
          }
        else
          {
            FlatVector<SCAL> elu2(dnums.Size()*dim*gf->GetCacheBlockSize(),lh);
            gf->GetElementVector (0, dnums, elu2);
            for(int i=0; i<elu.Size(); i++)
              elu[i] = elu2[i*gf->GetCacheBlockSize()+multidimcomponent];
          }
        
        fes.TransformVec (elnr, bound, elu, TRANSFORM_SOL);

        
        if ( bound ? 
             !fes.DefinedOnBoundary(eltrans.GetElementIndex()) : 
             !fes.DefinedOn(eltrans.GetElementIndex()) ) return 0;

	for (int k = 0; k < npts; k++)
	  for (int i = 0; i < components; i++)
	    values[k*svalues+i] = 0.0;

        if (bound)
          {
            for (int k = 0; k < npts; k++)
              {
                HeapReset hr(lh);
                
                IntegrationPoint ip(xref[k*sxref], xref[k*sxref+1], 0, 0);
		FlatVec<3> vx( (double*)x + k*sx);
		Mat<3,2> & mdxdxref = *new((double*)(dxdxref+k*sdxdxref)) Mat<3,2>;

		MappedIntegrationPoint<2,3> mip (ip, eltrans, vx, mdxdxref); 
                
                for(int j = 0; j<bfi2d.Size(); j++)
                  {
		    FlatVector<SCAL> flux (bfi2d[j]->DimFlux(), lh);
                    bfi2d[j]->CalcFlux (*fel, mip, elu, flux, applyd, lh);
                    for (int i = 0; i < components; i++)
                      values[k*svalues+i] += ((double*)(void*)&flux(0))[i];
                  }
              }
          }
        else
          {
	    IntegrationRule ir(npts, lh);
	    for (int i = 0; i < npts; i++)
	      {
		ir[i] = IntegrationPoint (xref[i*sxref], xref[i*sxref+1]);
		ir[i].FacetNr() = facetnr;
	      }
	    MappedIntegrationRule<2,2> mir(ir, eltrans, 1, lh);

	    for (int k = 0; k < npts; k++)
	      {
		Mat<2,2> & mdxdxref = *new((double*)(dxdxref+k*sdxdxref)) Mat<2,2>;
		FlatVec<2> vx( (double*)x + k*sx);
		mir[k] = MappedIntegrationPoint<2,2> (ir[k], eltrans, vx, mdxdxref);
	      }

	    for(int j = 0; j < bfi2d.Size(); j++)
	      {
		FlatMatrix<SCAL> flux(npts, bfi2d[j]->DimFlux(), lh);
		bfi2d[j]->CalcFlux (*fel, mir, elu, flux, applyd, lh);

		for (int k = 0; k < npts; k++)
		  for (int i = 0; i < components; i++)
		    values[k*svalues+i] += ((double*)(void*)&flux(k,0))[i];
	      }
          }

        return true; 
      }
    catch (Exception & e)
      {
        cout << "GetMultiSurfaceValue caught exception" << endl
             << e.What();

        return 0;
      }
  }














  
  template <class SCAL>
  void VisualizeGridFunction<SCAL> :: 
  Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages, int component)
  {
    cout << "VisGF::Analyze" << endl;
    int ndomains = 0;

    if (bfi3d.Size()) 
      ndomains = ma.GetNDomains();
    else if(bfi2d.Size()) 
      ndomains = ma.GetNBoundaries();

    Array<double> volumes(ndomains);

    Analyze(minima,maxima,averages,volumes,component);
    
    for(int i=0; i<ndomains; i++)
      {
	if(component == -1)
	  for(int j=0; j<components; j++)
	    averages[i*components+j] /= volumes[i];
	else
	  averages[i] /= volumes[i];
      }
  }
  

  template <class SCAL>
  void VisualizeGridFunction<SCAL> :: 
  Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages_times_volumes, Array<double> & volumes, int component)
  {
    cout << "VisGF::Analyze2" << endl;
    const FESpace & fes = gf->GetFESpace();

    int domain;
    double *val;
    int pos;
    double vol;

    /*
    int ndomains;
    if(bfi3d.Size()) ndomains = ma.GetNDomains();
    else if(bfi2d.Size()) ndomains = ma.GetNBoundaries();
    */

    Array<double> posx;
    Array<double> posy;
    Array<double> posz;
    ELEMENT_TYPE cache_type = ET_SEGM;
	
    LocalHeapMem<10000> lh2("Gridfunction - Analyze");
	
    val = new double[components];
			

    for(int i=0; i<minima.Size(); i++)
      {
	minima[i] = 1e100;
	maxima[i] = -1e100;
	averages_times_volumes[i] = 0;
      }

    for(int i=0; i<volumes.Size(); i++) volumes[i] = 0;
    
    void * heapp = lh2.GetPointer();
    if(bfi3d.Size())
      {
	for(int i=0; i<ma.GetNE(); i++)
	  {
	    const FiniteElement & fel = fes.GetFE(i,lh2);
	    
	    domain = ma.GetElIndex(i);
	    
	    vol = ma.ElementVolume(i);
	    
	    volumes[domain] += vol;
	    
	    if(fel.ElementType() != cache_type)
	      {
		posx.DeleteAll(); posy.DeleteAll(); posz.DeleteAll(); 
		switch(fel.ElementType())
		  {
		  case ET_TET:
		    posx.Append(0.25); posy.Append(0.25); posz.Append(0.25); 
		    posx.Append(0); posy.Append(0); posz.Append(0); 
		    posx.Append(1); posy.Append(0); posz.Append(0);
		    posx.Append(0); posy.Append(1); posz.Append(0);
		    posx.Append(0); posy.Append(0); posz.Append(1);
		    break;
		  case ET_HEX:
		    posx.Append(0.5); posy.Append(0.5); posz.Append(0.5);  
		    posx.Append(0); posy.Append(0); posz.Append(0); 
		    posx.Append(0); posy.Append(0); posz.Append(1);
		    posx.Append(0); posy.Append(1); posz.Append(0);
		    posx.Append(0); posy.Append(1); posz.Append(1);  
		    posx.Append(1); posy.Append(0); posz.Append(0); 
		    posx.Append(1); posy.Append(0); posz.Append(1);
		    posx.Append(1); posy.Append(1); posz.Append(0);
		    posx.Append(1); posy.Append(1); posz.Append(1);
		    break;
                  default:
                    {
                      bool firsttime = 1;
                      if (firsttime)
                        {
                          cerr << "WARNING::VisGridFunction::Analyze: unsupported element "
                               << ElementTopology::GetElementName(fel.ElementType()) << endl;
                          firsttime = 0;
                        }
                      break;
                    } 
		  }
		cache_type = fel.ElementType();
	      }
	    
	    
	    for(int k=0; k<posx.Size(); k++)
	      {
		GetValue(i,posx[k],posy[k],posz[k],val);


		if(component == -1)
		  {
		    for(int j=0; j<components; j++)
		      {
			pos = domain*components+j;
			if(val[j] > maxima[pos]) maxima[pos] = val[j];
			if(val[j] < minima[pos]) minima[pos] = val[j];
			if(k==0) averages_times_volumes[pos] += val[j]*vol;
		      }
		  }
		else
		  {
		    pos = domain;
		    if(val[component] > maxima[pos]) maxima[pos] = val[component];
		    if(val[component] < minima[pos]) minima[pos] = val[component];
		    if(k==0) averages_times_volumes[pos] += val[component]*vol;
		  }
	      }
	    lh2.CleanUp(heapp);
	  }
      }
    else if (bfi2d.Size())
      {
	for(int i=0; i<ma.GetNSE(); i++)
	  {
	    const FiniteElement & fel = fes.GetSFE(i,lh2);

	    domain = ma.GetSElIndex(i);

	    vol = ma.SurfaceElementVolume(i);

	    volumes[domain] += vol;

	    if(fel.ElementType() != cache_type)
	      {
		posx.DeleteAll(); posy.DeleteAll(); posz.DeleteAll(); 
		switch(fel.ElementType())
		  {
		  case ET_POINT:  // ??
		  case ET_SEGM: 
		  case ET_TET: case ET_HEX: case ET_PRISM: case ET_PYRAMID:
		    break;
		    
		  case ET_TRIG:
		    posx.Append(0.33); posy.Append(0.33);
		    posx.Append(0); posy.Append(0);
		    posx.Append(0); posy.Append(1);
		    posx.Append(1); posy.Append(0);
		    break;
		  case ET_QUAD:
		    posx.Append(0.5); posy.Append(0.5);
		    posx.Append(0); posy.Append(0);
		    posx.Append(0); posy.Append(1); 
		    posx.Append(1); posy.Append(0);
		    posx.Append(1); posy.Append(1);
		    break;
		  }
		cache_type = fel.ElementType();
	      }
	    for(int k=0; k<posx.Size(); k++)
	      {
		GetSurfValue(i,-1, posx[k],posy[k],val);
		if(component == -1)
		  {
		    for(int j=0; j<components; j++)
		      {
			pos = domain*components+j;
			if(val[j] > maxima[pos]) maxima[pos] = val[j];
			if(val[j] < minima[pos]) minima[pos] = val[j];
			if(k==0) averages_times_volumes[pos] += val[j]*vol;
		      }
		  }
		else
		  {
		    pos = domain;
		    if(val[component] > maxima[pos]) maxima[pos] = val[component];
		    if(val[component] < minima[pos]) minima[pos] = val[component];
		    if(k==0) averages_times_volumes[pos] += val[component]*vol;
		  }
	      }
	    lh2.CleanUp(heapp);
	  }
      }

    delete [] val;
  }







  VisualizeCoefficientFunction :: 
  VisualizeCoefficientFunction (const MeshAccess & ama,
				const CoefficientFunction * acf)
    : SolutionData ("coef", -1, false /* complex */),
      ma(ama), cf(acf), lh(100000, "viscoef-lh")
  { ; }
  
  VisualizeCoefficientFunction :: ~VisualizeCoefficientFunction ()
  {
    ;
  }
  
  bool VisualizeCoefficientFunction :: GetValue (int elnr, 
						 double lam1, double lam2, double lam3,
						 double * values) 
  {
    cout << "visualizecoef, getvalue (lam1,lam2,alm3) not implemented" << endl;
    return false;
  }
  
  bool VisualizeCoefficientFunction :: 
  GetValue (int elnr, 
	    const double xref[], const double x[], const double dxdxref[],
	    double * values) 
  {
    cout << "visualizecoef, getvalue (xref) not implemented" << endl;
    return false;
  }

  bool VisualizeCoefficientFunction :: 
  GetMultiValue (int elnr, int npts,
		 const double * xref, int sxref,
		 const double * x, int sx,
		 const double * dxdxref, int sdxdxref,
		 double * values, int svalues)
  {
    cout << "visualizecoef, GetMultiValue not implemented" << endl;
    return false;
  }
  
  bool VisualizeCoefficientFunction ::  
  GetSurfValue (int elnr, int facetnr,
		double lam1, double lam2, 
		double * values) 
  {
    HeapReset hr(lh);
    IntegrationPoint ip(lam1, lam2);
    ip.FacetNr() = facetnr;
    bool bound = ma.GetDimension() == 3;
    ElementTransformation & trafo = ma.GetTrafo (elnr, bound, lh);
    BaseMappedIntegrationPoint & mip = trafo(ip, lh);
    values[0] = cf -> Evaluate (mip);
    return true; 
  }

  bool VisualizeCoefficientFunction ::  GetSurfValue (int selnr, int facetnr, 
			       const double xref[], const double x[], const double dxdxref[],
			       double * values)
  {
    cout << "visualizecoef, getsurfvalue (xref) not implemented" << endl;
    return false;
  }

  bool VisualizeCoefficientFunction ::  
  GetMultiSurfValue (int selnr, int facetnr, int npts,
		     const double * xref, int sxref,
		     const double * x, int sx,
		     const double * dxdxref, int sdxdxref,
		     double * values, int svalues)
  {
    for (int i = 0; i < npts; i++)
      GetSurfValue (selnr, facetnr, xref[i*sxref], xref[i*sxref+1], &values[i*svalues]);
    return true;
  }


  void VisualizeCoefficientFunction ::  
  Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages, 
	  int component)
  {
    cout << "visualizecoef, analyze1 not implemented" << endl;
  }

  void VisualizeCoefficientFunction :: 
  Analyze(Array<double> & minima, Array<double> & maxima, Array<double> & averages_times_volumes, 
	  Array<double> & volumes, int component)
  {
    cout << "visualizecoef, analyzed2 not implemented" << endl;
  }










  template class T_GridFunction<double>;
  template class T_GridFunction<Vec<2> >;
  template class T_GridFunction<Vec<3> >;
  template class T_GridFunction<Vec<4> >;
  

  template class  VisualizeGridFunction<double>;
  template class  VisualizeGridFunction<Complex>;

}
