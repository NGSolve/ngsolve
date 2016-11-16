#include "Distribution.hpp"
#include "polynomials.hpp"
#include "omp.h"

void MatVec(AFlatMatrix<double> m, FlatVector<double> x, FlatVector<double> y);

Timer timern2h("nodal2hermite");
Timer timern2htrans("nodal2hermite trans");

Timer timerh2p("hermite2polar");
Timer timerph2p("hermite2polar - 2");
Timer timerh2ph("hermite2polar - 1");
Timer timerh2ptrans("hermite2polar trans");
Timer timerph2ptrans("hermite2polar trans - 2");
Timer timerh2phtrans("hermite2polar trans - 1");


template<int dim>
Distribution<dim>::Distribution(int aorder,bool trafos)
{
    order = aorder;
    n_nodes = order+1;
    n_dof_nodal = pow(order+1,dim);
    if(dim == 3)
        n_dof_hierarchical = (order+1)*(order+2)*(order+3)/6;
    if(dim == 2)
        n_dof_hierarchical = (order+1)*(order+2)/2;
    if(dim == 1)
        n_dof_hierarchical = n_dof_nodal;
    ndof = n_dof_nodal;
    ansatztemp = 1.0;
    ansatzv.SetSize(dim);
    ansatzv = 0.0;
    Array<double> nodes,weights;
    ComputeHermiteRule (n_nodes, nodes, weights);
    // Precalculate the denomenator of the Lagrange polynomials
    fac.SetSize(n_nodes);
    for (int i = 0; i < n_nodes; i++)
    {
        double prod = 1;
        for (int j = 0; j < n_nodes; j++)
            if (j != i) 
            {
                prod *= (nodes[i]-nodes[j]);
            }
        fac(i) = 1.0 / prod;
    }
    // The inverse of the mass matrix in nodal representation
    // In Polar or Hermite representation the mass matrix is the unit matrix
    cout << "Calculating Inverse Mass matrix...";
    diagmassinv.SetSize(GetNDof<NODAL>());
    diagmass.SetSize(GetNDof<NODAL>());
    if(dim == 3)
        for(int i=0,ii=0;i<order+1;i++)
            for(int j=0;j<order+1;j++)
                for(int k=0;k<order+1;k++,ii++)
                {
                    diagmassinv(ii)=1.0/(weights[i]*weights[j]*weights[k]);
                    diagmass(ii)=weights[i]*weights[j]*weights[k];
                }
    else if(dim == 2)
        for(int i=0,ii=0;i<order+1;i++)
            for(int j=0;j<order+1;j++,ii++)
            {
                diagmassinv(ii)=1.0/(weights[i]*weights[j]);
                diagmass(ii)=weights[i]*weights[j];
            }
    else if(dim == 1)
        for(int i=0;i<order+1;i++)
        {
            diagmassinv(i)=1.0/(weights[i]);
            diagmass(i)=weights[i];
        }
    cout << "done"<<endl;
    cout << "Calculate nodal2hermite trafo... ";
    // The transformation of a distribution from its nodal to its hermite representation 
    // and the inverse 
    nodal2hermite.SetSize(order+1);
    hermite2nodal.SetSize(order+1);
    scale.SetSize(GetNDof<NODAL>());
    for (int i = 0; i <= order; i++)
    {
        HermiteFunction (order, nodes[i], nodal2hermite.Row(i));
        nodal2hermite.Row(i)*= (weights[i] *exp(nodes[i]*nodes[i]) );
        scale(i) = exp ( -sqr(nodes[i])/2 );
    }
    hermite2nodal = nodal2hermite;
    CalcInverse(hermite2nodal);
    cout << "done"<<endl;
    if(trafos == false)
        return;
    // The Polar trafo is split into 2 transformations. 
    // The first trafo is block diagonal. To have the second trafo
    // also block diagonal a rearrangement of the dofs in the intermediate basis is necessarry
    phtoph2.SetSize(GetNDof<POLAR>());
    Tensor<3> polartensor(order+1,order+1,2*order+1);
    for(int n = 0,ii = 0;n < order+1;n++)
    {
        for (int k = 0; k <= n; k++)
        {
            for(int s=0;s<=k/2;s++)
                polartensor(n,k,2*s+k%2) = ii++;     // cos( (2*s+k%2) * phi)*...
            for (int s = 0; s <= k/2; s++)
            {
                if (2*s+k%2 != 0)
                    polartensor(n,k,n+2*s+k%2) = ii++; // sin( (2*s+k%2) * phi)*...
            }
        }
    }
    for(int n = 0,ii = 0; n < order+1;n++)
        for(int s = 0;s <= n/2; s++)
        {
            for(int k=2*s;k<=n;k+=2)
                phtoph2[ii++] = (int)polartensor(n,k,2*s + k%2);    // cos( (2*s+k%2) * phi)*...
            for(int k=2*s;k<=n;k+=2)
                if (s != 0)
                    phtoph2[ii++] = (int)polartensor(n,k,n+2*s + k%2);  // sin( (2*s+k%2) * phi)*...  
            for(int k=2*s+1;k<=n;k+=2)
                phtoph2[ii++] = (int)polartensor(n,k,2*s + k%2);        // cos( (2*s+k%2) * phi)*...
            for(int k=2*s+1;k<=n;k+=2)
                phtoph2[ii++] = (int)polartensor(n,k,n+2*s + k%2);      // sin( (2*s+k%2) * phi)*...
        }
    // Calcualation of the trafo from hermite to Polar hermite representation. (x,y) components are
    // transformed, z component stays the same.
    cout << "\rCalculating hermite to polar hermite trafo... "<<flush;
    hermite2polarhermites.SetSize(order+1);
    for(int k=0;k<order+1;k++)
        hermite2polarhermites[k] = new Matrix<>(k+1,k+1);
    Vector<> shapehermite1((order+1)*(order+2)/2);
    Vector<> shapepolarhermite1((order+1)*(order+2)/2);
    Matrix<> ordern((order+1)*(order+2)/2);
    ordern = 0.0;
    for(int i=0;i<order+1;i++)
        for(int j=0;j<order+1;j++)
            for(int k=0;k<order+1;k++)
            {
                CalcShapeOrder<HERMITE>(IntegrationPoint(nodes[i],nodes[j],nodes[k]),order,shapehermite1);
                CalcShapeOrder<POLARHERMITE>(IntegrationPoint(nodes[i],nodes[j],nodes[k]),order,shapepolarhermite1);
                ordern += weights[i]*weights[j]*weights[k]*shapehermite1*Trans(shapepolarhermite1);
            }
    for(int k=0,ii=0;k<order+1;k++)
    {
        *hermite2polarhermites[k] = ordern.Rows(ii,ii+k+1).Cols(ii,ii+k+1);
        ii += k+1;
    }
    cout << "done"<<endl;
    // Finally the trafo from polar hermite to polar representation
    nops_ph2p = 0;
    Array<Array< Matrix<> * >* > ph2ptemp;
    ph2ptemp.SetSize(order+1);
    #pragma omp parallel for
    for(int k=0;k<order+1;k++)
    {
        (ph2ptemp[k]) = new Array<Matrix<> *>();
        int iiii = 0;
        //if(omp_get_thread_num() == 0) TODO
            cout << "\rCalculating polar hermite to polar trafo of order "<<k<<"/"<<order<<"... "<<flush;
        Vector<> shapepolarhermite((k+1)*(k+2)/2);
        Vector<> shapepolar((k+1)*(k+2)/2);
        Matrix<> temp((k+1)*(k+2)/2);
        temp = 0.0;
        Array<double> nodestemp,weightstemp;
        ComputeHermiteRule(k+1,nodestemp,weightstemp);
        for(int i=0;i<k+1;i++)
            for(int j=0;j<k+1;j++)
                for(int kk=0;kk<k+1;kk++)
                {
                    CalcShapeOrder<POLAR>(IntegrationPoint(nodestemp[i],nodestemp[j],nodestemp[kk]),k,shapepolar);
                    CalcShapeOrder<POLARHERMITE2>(IntegrationPoint(nodestemp[i],nodestemp[j],nodestemp[kk]),k,shapepolarhermite);
                    temp += weightstemp[i]*weightstemp[j]*weightstemp[kk]*shapepolarhermite*Trans(shapepolar);
                }
        if(k==0)
            temp = 1.0;
        int blstart = 0;
        int blsize = k/2+1;
        nops_ph2p += blsize*blsize;
        (*ph2ptemp[k]).Append(new Matrix<>(blsize));
        (*(*ph2ptemp[k])[iiii++]) = Trans(temp).Rows(blstart,blstart+blsize).Cols(blstart,blstart+blsize);
        blstart += blsize;
        if(k%2 == 1)
            for(int bl = 0;bl<2;bl ++)
            {
                int blsize = k/2+1;
                nops_ph2p += blsize*blsize;
                (*ph2ptemp[k]).Append(new Matrix<>(blsize));
                (*(*ph2ptemp[k])[iiii++]) = Trans(temp).Rows(blstart,blstart+blsize).Cols(blstart,blstart+blsize);
                blstart += blsize;
            }
        for(int blsize = k/2;blsize>=1;blsize--)
        {
            for(int bl = 0;bl<4;bl ++)
            {
                nops_ph2p += blsize*blsize;
                (*ph2ptemp[k]).Append(new Matrix<>(blsize));
                (*(*ph2ptemp[k])[iiii++]) = Trans(temp).Rows(blstart,blstart+blsize).Cols(blstart,blstart+blsize);
                blstart += blsize;
            }
        }
    }
    for(int i=0,ii=0;i<order+1;i++)
    {
        for(int j=0;j<(*ph2ptemp[i]).Size();j++,ii++)
        {
            polarhermites2polar.Append(new Matrix<>( (*(*ph2ptemp[i])[j]).Height(),(*(*ph2ptemp[i])[j]).Width() ));
            *polarhermites2polar[ii] = (*(*ph2ptemp[i])[j]);
            delete (*ph2ptemp[i])[j];
        }
        delete ph2ptemp[i];
    }
    cout << "done"<<endl;
    CalculateDofTable();
    //CalculateMultr();
}

template<int dim>
void Distribution<dim>::CalculateDofTable()
{
    polardoftable.SetSize((2*order+1)*(order+1));
    for(int u=0;u<2*order+1;u++)
    {
        for(int v=0;v<order+1;v++)
        {
            polardoftable(u*(order+1)+v) = new Array<Vec<3> >(0);
            polardoftable(u*(order+1)+v)->Append(Vec<3>(u,0,0));
            polardoftable(u*(order+1)+v)->Append(Vec<3>(v,0,0));
            int nr = 0;
            for(int n = 0; n<=order;n++)
            {
                for(int j=0;j<=n;j++)
                {
                    for(int k=j/2+(1-n%2)*j%2;k<=n/2;k++)
                    {
                        if(j == 0)
                        {
                            if(u == 0 && 2*k+n%2 == v)
                                polardoftable(u*(order+1)+v)->Append(Vec<3>(nr,2*k+n%2,(n-n%2-2*k)/2));
                            //shape(nr++)   = sqrt(2.0)*radial(2*k+n%2,(n-n%2-2*k)/2)*angular(0,2*k+n%2);
                        }
                        else
                        {
                            if(u == 2*j-1 && 2*k+n%2 == v)
                                polardoftable(u*(order+1)+v)->Append(Vec<3>(nr,2*k+n%2,(n-n%2-2*k)/2));
                        //        cout << u<< " " <<v<< " " << n <<" "<< j << " "<<k<<endl;
                        }
                        nr++;
                            //shape(nr++) = 2.0*radial(2*k+n%2,(n-n%2-2*k)/2)*angular(2*j-1,2*k+n%2);
                    }
                    for(int k=j/2+(1-n%2)*j%2;k<=n/2;k++)
                        if(j!=0)
                        {
                            if(u == 2*j && 2*k+n%2 == v)
                                polardoftable(u*(order+1)+v)->Append(Vec<3>(nr,2*k+n%2,(n-n%2-2*k)/2));
                            //shape(nr++) = 2.0*radial(2*k+n%2,(n-n%2-2*k)/2)*angular(2*j,2*k+n%2);
                            nr++;
                        }
                }
            }
        }
    }
    //for(int i=0;i<polardoftable.Size();i++)
    //    cout << *polardoftable(i)<<endl;
}

template<int dim>
void Distribution<dim>::Project(const FlatVector<> f_in, FlatVector<> f_out, double T_in, double T_out, Vec<dim> v_in, Vec<dim> v_out,LocalHeap & lh) const
{
    FlatMatrix<> f_inmat((int)f_in.Size()/GetNDof<NODAL>(), GetNDof<NODAL>(), &f_in(0));
    Array<double> nodes, weights;
    HermiteRule(nodes,weights,2*order);
    FlatMatrix<> shapes(nodes.Size(),GetNDof<NODAL>(),lh);
    FlatMatrix<> shapesmod(nodes.Size(),GetNDof<NODAL>(),lh);
    FlatMatrix<> masshm(GetNDof<NODAL>(),lh);
    for(int ip=0;ip<nodes.Size();ip++)
    {
        CalcShape1D(nodes[ip],shapes.Row(ip));
        double scaled_weight = weights[ip]*exp(nodes[ip]*nodes[ip]);
        if(use_hm_funcs)
            shapes.Row(ip)*=  scaled_weight;
        else
            shapes.Row(ip)*=  scaled_weight * exp(-nodes[ip]*nodes[ip]/2.0 +  ( sqrt(T_in) * nodes[ip] + v_in(0)-v_out(0) )*( sqrt(T_in) * nodes[ip] + v_in(0)-v_out(0) )/( 2.0*T_out ) );
        CalcShape1D( ( sqrt(T_in) * nodes[ip] + v_in(0)-v_out(0) )/sqrt(T_out),shapesmod.Row(ip));
    }
    masshm = Trans(shapes)*shapesmod;
    FlatMatrix<> f_outmat(f_inmat.Height(),f_inmat.Width(),&f_out(0));
    f_outmat = f_inmat*masshm;
    for(int i=0;i<f_outmat.Width();i++)
        f_outmat.Col(i)*=diagmassinv(i);    
    f_out*=sqrt(T_in/T_out);
}


template<int dim>
void Distribution<dim>::ProjectTrans(const FlatVector<> f_in, FlatVector<> f_out, double T_in, double T_out, Vec<dim> v_in, Vec<dim> v_out,LocalHeap & lh) const
{
    FlatMatrix<> f_inmat((int)f_in.Size()/GetNDof<NODAL>(), GetNDof<NODAL>(), &f_in(0));
    Array<double> nodes, weights;
    HermiteRule(nodes,weights,2*order);
    FlatMatrix<> shapes(nodes.Size(),GetNDof<NODAL>(),lh);
    FlatMatrix<> shapesmod(nodes.Size(),GetNDof<NODAL>(),lh);
    FlatMatrix<> masshm(GetNDof<NODAL>(),lh);
    for(int ip=0;ip<nodes.Size();ip++)
    {
        CalcShape1D(nodes[ip],shapes.Row(ip));
        double scaled_weight = weights[ip]*exp(nodes[ip]*nodes[ip]);
        if(use_hm_funcs)
            shapes.Row(ip)*=  scaled_weight;
        else
            shapes.Row(ip)*=  scaled_weight * exp(-nodes[ip]*nodes[ip]/2.0 +  ( sqrt(T_in) * nodes[ip] + v_in(0)-v_out(0) )*( sqrt(T_in) * nodes[ip] + v_in(0)-v_out(0) )/( 2.0*T_out ) );
        CalcShape1D( ( sqrt(T_in) * nodes[ip] + v_in(0)-v_out(0) )/sqrt(T_out),shapesmod.Row(ip));
    }
    masshm = Trans(shapes)*shapesmod;
    FlatMatrix<> f_outmat(f_inmat.Height(),f_inmat.Width(),&f_out(0));
    for(int i=0;i<f_inmat.Width();i++)
        f_inmat.Col(i)*=diagmassinv(i);
    f_outmat = f_inmat*masshm;
}


template<int dim>
void Distribution<dim>::GetDofNrFixedAngularDof(int m,int n,Array<int> & dofs) const
{
    int nr = (m*(order+1)+n);
    dofs.SetSize(polardoftable(nr)->Size()-2);
    for(int i=2;i<polardoftable(nr)->Size();i++)
        dofs[i-2] = (*polardoftable(nr))[i](0);
}

template<int dim>
void Distribution<dim>::GetRadialIndexFixedAngularDof(int m,int n,Array<Vec<2> > & indices) const
{
    int nr = (m*(order+1)+n);
    indices.SetSize(polardoftable(nr)->Size()-2);
    for(int i=2;i<polardoftable(nr)->Size();i++)
        indices[i-2] = Vec<2>( (*polardoftable(nr))[i](1),(*polardoftable(nr))[i](2) );
}

template<int dim>
void Distribution<dim>::CalculateMultr()
{
    multrs.SetSize(polardoftable.Size());
    for(int i=0;i<multrs.Size();i++)
    {
        if(polardoftable[i]->Size()<=2)
            continue;
        multrs[i] = new Matrix<>(polardoftable[i]->Size()-2);
        (*multrs[i]) = 0.0;
        double alpha = (*polardoftable(i))[2](1);
        Array<double> nodeslag, weightslag;        
        int ordermax = (*polardoftable(i))[ (*polardoftable(i)).Size()-1 ](2);
        ComputeLaguerreRule(ordermax+1,alpha+0.5 + 0.5,nodeslag,weightslag);
        cout << *polardoftable(i) << endl;
        for(int ip=0;ip<nodeslag.Size();ip++)
        {
            double r = nodeslag[ip];
            Vector<> radial(ordermax+1);
            // Calculate Laguerre polynomials times r^k
            LaguerrePolynomial ( ordermax, r, alpha+0.5, radial);
            for (int j = 0; j <= ordermax; j++)
                radial(j) *= exp(-0.5*(my_gammln(alpha+j+1.5)-my_gammln(j+1)));
            (*multrs[i]) += weightslag[ip]*radial*Trans(radial);
        }
    }
}


template<int dim>
inline void Distribution<dim>::Nodes(Array<double> & nodes,int aorder) const
{
    Array<double> weights;
    if(aorder==-1)
        aorder=n_nodes;
    ComputeHermiteRule(aorder,nodes,weights);
}

template<int dim>
inline void Distribution<dim>::Weights(Array<double> & weights,int aorder) const
{
    Array<double> nodes;
    if(aorder==-1)
        aorder=n_nodes;
    ComputeHermiteRule(aorder,nodes,weights);
}

template<int dim>
inline void Distribution<dim>::HermiteRule(Array<double> & nodes,Array<double> & weights,int aorder) const
{
    if(aorder==-1)
        aorder=n_nodes;
    ComputeHermiteRule(n_nodes,nodes,weights);
}

template<int dim>
void Distribution<dim>::CalcShape1D(const double x,BareSliceVector<> shape) const
{
    VectorMem<50> prod1x(n_nodes), prod2x(n_nodes);
    double px = 1;
    Array<double> nodes,weights;
    HermiteRule(nodes,weights);
    for (int i = 0; i < order; i++)
    {
        prod1x[i] = px;
        px *= x-nodes[i];
    }
    prod1x[order] = px; 
    px = 1; 
    for (int i = order; i > 0; i--)
    {
        prod2x[i] = px;
        px *= x-nodes[i];
    } 
    prod2x[0] = px;
    for (int i = 0; i < n_nodes; i++)
        shape(i) = exp(-0.5*x*x)*prod1x(i) * prod2x(i) * fac(i);
}

template<int dim>
template<Representation REP>
void Distribution<dim>::CalcShape1(const IntegrationPoint & ip, BareSliceVector<> shape) const
{
    double x = ip(0);
    double y = ip(1);
    //shape = 0.0;
    if(REP == NODAL)
    {
        Array<double> nodes,weights;
        HermiteRule(nodes,weights);
        Vector<> shapex(n_nodes), shapey(n_nodes);
        if(dim == 1)
            CalcShape1D(ip(0),shape);
        if(dim == 2)
        {
            for (int i = 0; i < n_nodes; i++)
            {
                double prodx = 1;
                double prody = 1;
                for (int j = 0; j < n_nodes; j++)
                    if (j != i) 
                    {
                        prodx *= (x-nodes[j]);
                        prody *= (y-nodes[j]);
                    }
                shapex(i) = prodx*fac(i);
                shapey(i) = prody*fac(i);
            }
            for (int i = 0, ii = 0; i < n_nodes; i++)
                for (int j = 0; j < n_nodes; j++,ii++)
                    shape(ii) = shapex(i) * shapey(j);// * exp( -0.5*(nodes[i]*nodes[i] + nodes[j]*nodes[j]) );
        }
        if(dim == 3)
        {
            double z = ip(2);
            Vector<> shapez(n_nodes);
            for (int i = 0; i < n_nodes; i++)
            {
                double prodx = 1;
                double prody = 1;
                double prodz = 1;
                for (int j = 0; j < n_nodes; j++)
                    if (j != i) 
                    {
                        prodx *= (x-nodes[j]);
                        prody *= (y-nodes[j]);
                        prodz *= (z-nodes[j]);
                    }
                shapex(i) = prodx*fac(i);
                shapey(i) = prody*fac(i);
                shapez(i) = prodz*fac(i);
            }
            for (int i = 0, ii = 0; i < n_nodes; i++)
                for (int j = 0; j < n_nodes; j++)
                    for (int k = 0; k < n_nodes; k++, ii++)
                        shape(ii) = shapex(i) * shapey(j) * shapez(k);// * exp( -0.5*(nodes[i]*nodes[i] + nodes[j]*nodes[j] + nodes[k]*nodes[k]) );
        }
        return;
    }
    if(REP == HERMITE)
    {
        Vector<> valuesx(order+1),valuesy(order+1);
        ScaledHermitePolynomial(order,x, valuesx);
        ScaledHermitePolynomial(order,y, valuesy);
        if(dim == 1)
            for(int l=0;l<=order;l++)
                    shape(l) = valuesx(l);
        if(dim == 2)
            for(int l=0,ii=0;l<=order;l++)
                for(int m=0;m<=l;m++,ii++)
                    shape(ii) = valuesx(m)*valuesy(l-m);
        if(dim == 3)
        {
            double z = ip(2);
            Vector<> valuesz(order+1);
            ScaledHermitePolynomial(order,z, valuesz);
            for(int l=0,ii=0;l<=order;l++)
                for(int m=0;m<=l;m++)
                    for(int n=0;n<=m;n++,ii++)
                        shape(ii) = valuesx(n)*valuesy(m-n)*valuesz(l-m);
        }
        return;
    }
    if(REP == POLAR)
    {
        double z = ip(2);
        double phi = atan2(y, x);
        Matrix<> radial(order+1);
        if(dim == 2)
        {
            double r = sqrt(x*x+y*y);
            for (int al = 0; al <= order; al++)
            {
                LaguerrePolynomial ( (order-al)/2, r*r, al, radial.Row(al));
                for (int j = 0; j <= (order-al)/2; j++)
                    if(al!=0)
                        radial(al,j) *= exp(al*log(r)-0.5*(my_gammln(al+j+1)-my_gammln(j+1)));
            }
            for (int k = 0, ii = 0; k <= order; k++)
                for (int j = k%2; j <= k; j+=2)
                {
                    if (j != 0)
                        shape(ii++) = sin(j*phi) * radial(j, (k-j)/2);
                    shape(ii++) = cos(j*phi) * radial(j, (k-j)/2);
                }
        }
        if(dim == 3)
        {
            double r = sqrt(x*x+y*y+z*z);
            double z = ip(2);
            double theta;
            theta = acos(z/r);
            Matrix<> angular(2*order+1,order+1);
            // Calculate Laguerre polynomials times r^k
            for (int k = 0; k <= order; k++)
            {
                LaguerrePolynomial ( (order-k)/2, r*r, k+0.5, radial.Row(k));
                for (int j = 0; j <= (order-k)/2; j++)
                    radial(k,j) *= exp(k*log(r)-0.5*(my_gammln(k+j+1.5)-my_gammln(j+1)));
            }
            // Calculate the Spherical harmonics (real version)
            SphericalHarmonicReal(order,theta,phi,angular);
            //cout << angular<<endl;
            int nr = 0;
            for(int n = 0; n<=order;n++)
            {
                for(int j=0;j<=n;j++)
                {
                    for(int k=j/2+(1-n%2)*j%2;k<=n/2;k++)
                        if(j == 0)
                            shape(nr++)   = sqrt(2.0)*radial(2*k+n%2,(n-n%2-2*k)/2)*angular(0,2*k+n%2);
                        else
                            shape(nr++) = 2.0*radial(2*k+n%2,(n-n%2-2*k)/2)*angular(2*j-1,2*k+n%2);
                    for(int k=j/2+(1-n%2)*j%2;k<=n/2;k++)
                        if(j!=0)
                            shape(nr++) = 2.0*radial(2*k+n%2,(n-n%2-2*k)/2)*angular(2*j,2*k+n%2);
                }
            }
        }
        return;
    }
    if( REP == POLARHERMITE)
    {
        double z = ip(2);
        Vector<> hermitepoly(order+1);
        ScaledHermitePolynomial(order,z,hermitepoly);
        double phi = atan2(y, x);
        double r = sqrt(x*x+y*y);
        Matrix<> radial(order+1);
        for (int al = 0; al <= order; al++)
        {
            LaguerrePolynomial ( (order-al)/2, r*r, al, radial.Row(al));
            for (int j = 0; j <= (order-al)/2; j++)
                if(al!=0)
                    radial(al,j) *= exp(al*log(r)-0.5*(my_gammln(al+j+1)-my_gammln(j+1)));
        }
        for(int n = 0,ii = 0;n < order+1;n++)
            for (int k = 0; k <= n; k++)
            {
                for(int s=0;s<=k/2;s++)
                {
                    int j = 2*s+k%2;
                    shape(ii++) = cos(j*phi) * radial(j, (k-j)/2) * hermitepoly(n-k) * sqrt(2.0/M_PI);
                    if( j == 0)
                        shape(ii-1) *=  sqrt(1.0/2.0);
                }
                for (int s = 0; s <= k/2; s++)
                {
                    int j = 2*s+k%2;
                    if (j != 0)
                        shape(ii++) = sin(j*phi) * radial(j, (k-j)/2) * hermitepoly(n-k) * sqrt(2.0/M_PI);
                }
            }
        return;
    }
    if( REP == POLARHERMITE2)
    {
        double z = ip(2);
        Vector<> hermitepoly(order+1);
        ScaledHermitePolynomial(order,z,hermitepoly);
        double phi = atan2(y, x);
        Matrix<> radial(order+1);
        double r = sqrt(x*x+y*y);
        for (int al = 0; al <= order; al++)
        {
            LaguerrePolynomial ( (order-al)/2, r*r, al, radial.Row(al));
            for (int j = 0; j <= (order-al)/2; j++)
                if(al!=0)
            radial(al,j) *= exp(al*log(r)-0.5*(my_gammln(al+j+1)-my_gammln(j+1)));
        }
        for(int n = 0,ii = 0; n < order+1;n++)
            for(int s = 0;s <= n/2; s++)
            {
                for(int k=2*s;k<=n;k+=2)
                    if( s == 0)
                        shape(ii++) = cos(2*s*phi) * radial(2*s, (k-2*s)/2) * hermitepoly(n-k) * sqrt(1.0/M_PI);
                    else
                        shape(ii++) = cos(2*s*phi) * radial(2*s, (k-2*s)/2) * hermitepoly(n-k) * sqrt(2.0/M_PI);
                for(int k=2*s;k<=n;k+=2)
                    if (s != 0)
                        shape(ii++) = sin(2*s*phi) * radial(2*s, (k-2*s)/2) * hermitepoly(n-k) * sqrt(2.0/M_PI);
                for(int k=2*s+1;k<=n;k+=2)
                    shape(ii++) = cos((2*s+1)*phi) * radial(2*s+1, (k-(2*s+1))/2) * hermitepoly(n-k) * sqrt(2.0/M_PI);
                for(int k=2*s+1;k<=n;k+=2)
                    shape(ii++) = sin((2*s+1)*phi) * radial(2*s+1, (k-(2*s+1))/2) * hermitepoly(n-k) * sqrt(2.0/M_PI);
            }
            return;
    }
    cout << "CalcShape with dim = "<<dim<<", REP = "<<REP<<" not yet implemented"<<endl;
}

template<int dim>
template<Representation REP>
void Distribution<dim>::CalcShapeOrder(const IntegrationPoint & ip, int k, SliceVector<> shape) const
{
    double x = ip(0);
    double y = ip(1);
    shape = 0.0;
    if(REP == HERMITE)
    {
        Vector<> valuesx(k+1),valuesy(k+1);
        ScaledHermitePolynomial(k,x,valuesx);
        ScaledHermitePolynomial(k,y,valuesy);
        if(dim == 2)            
            for(int m=0;m<=k;m++)
                shape(m) = valuesx(m)*valuesy(k-m);
        if(dim == 3)
        {
            double z = ip(2);
            Vector<> valuesz(k+1);
            ScaledHermitePolynomial(k,z,valuesz);
            for(int m=0,ii=0;m<=k;m++)
                for(int n=0;n<=m;n++,ii++)
                    shape(ii) = valuesx(n)*valuesy(m-n)*valuesz(k-m);
        }
        return;
    }
    if(REP == POLAR)
    {
        double z = ip(2);
        double phi = atan2(y, x);
        Matrix<> radial(k+1);
        radial = 0.0;
        if(dim == 2)
        {
            double r = sqrt(x*x+y*y);
            for (int al = 0; al <= k; al++)
            {
                LaguerrePolynomial ( (k-al)/2, r*r, al, radial.Row(al));
                for (int j = 0; j <= (k-al)/2; j++)
                    if(al!=0)
                        radial(al,j) *= exp(al*log(r)-0.5*(my_gammln(al+j+1)-my_gammln(j+1)));
            }            
            for (int j = k%2, ii=0; j <= k; j+=2)
            {
                if (j != 0)
                    shape(ii++) = sin(j*phi) * radial(j, (k-j)/2);
                shape(ii++) = cos(j*phi) * radial(j, (k-j)/2);
            }
        }
        if(dim == 3)
        {
            double r = sqrt(x*x+y*y+z*z);
            double z = ip(2);
            double theta;
            theta = acos(z/r);
            Matrix<> angular(2*k+1,k+1);
            angular = 0.0;            
            // Calculate Laguerre polynomials times r^k
            for (int i = 0; i <= k; i++)
            {
                LaguerrePolynomial ( (k-i)/2, r*r, i+0.5, radial.Row(i));
                for (int j = 0; j <= (k-i)/2; j++)
                {
                    if(r!=0)
                        radial(i,j) *= exp(i*log(r)-0.5*(my_gammln(i+j+1.5)-my_gammln(j+1)));
                    else
                        radial(i,j) *= exp(-0.5*(my_gammln(i+j+1.5)-my_gammln(j+1)));
                }
            }
            // Calculate the Spherical harmonics (real version)
            SphericalHarmonicReal(k,theta,phi,angular);
            int nr = 0;
            for(int j=0;j<=k;j++)
            {
                for(int i=j/2+(1-k%2)*j%2;i<=k/2;i++)
                    if(j == 0)
                        shape(nr++)   = sqrt(2.0)*radial(2*i+k%2,(k-k%2-2*i)/2)*angular(0,2*i+k%2);
                    else
                        shape(nr++) =         2.0*radial(2*i+k%2,(k-k%2-2*i)/2)*angular(2*j-1,2*i+k%2);
                for(int i=j/2+(1-k%2)*j%2;i<=k/2;i++)
                    if(j!=0)
                        shape(nr++) = 2.0*radial(2*i+k%2,(k-k%2-2*i)/2)*angular(2*j,2*i+k%2);
            }            
            
        }
        return;
    }
    if( REP == POLARHERMITE )
    {
        double z = ip(2);
        Vector<> hermitepoly(order+1);
        ScaledHermitePolynomial(order,z,hermitepoly);
        double phi = atan2(y, x);
        double r = sqrt(x*x+y*y);
        Matrix<> radial(order+1);
        for (int al = 0; al <= order; al++)
        {
            LaguerrePolynomial ( (order-al)/2, r*r, al, radial.Row(al));
            for (int j = 0; j <= (order-al)/2; j++)
                if(al!=0)
                    radial(al,j) *= exp(al*log(r)-0.5*(my_gammln(al+j+1)-my_gammln(j+1)));
        }
        for (int i = 0,ii=0; i <= k; i++)
        {
            for(int s=0;s<=i/2;s++)
            {
                int j = 2*s+i%2;
                shape(ii++) = cos(j*phi) * radial(j, (i-j)/2) * hermitepoly(k-i) * sqrt(2.0/M_PI);
                if( j == 0)
                    shape(ii-1) *=  sqrt(1.0/2.0);
            }
            for (int s = 0; s <= i/2; s++)
            {
                int j = 2*s+i%2;
                if (j != 0)
                    shape(ii++) = sin(j*phi) * radial(j, (i-j)/2) * hermitepoly(k-i) * sqrt(2.0/M_PI);
            }
        }
        return;        
    }
    if( REP == POLARHERMITE2)
    {
        double z = ip(2);
        Vector<> hermitepoly(order+1);
        ScaledHermitePolynomial(order,z,hermitepoly);
        double phi = atan2(y, x);
        Matrix<> radial(order+1);
        double r = sqrt(x*x+y*y);
        for (int al = 0; al <= order; al++)
        {
            LaguerrePolynomial ( (order-al)/2, r*r, al, radial.Row(al));
            for (int j = 0; j <= (order-al)/2; j++)
                if(al!=0)
            radial(al,j) *= exp(al*log(r)-0.5*(my_gammln(al+j+1)-my_gammln(j+1)));
        }
        for(int s = 0,ii=0;s <= k/2; s++)
        {
            for(int i=2*s;i<=k;i+=2)
                if( s == 0)
                    shape(ii++) = cos(2*s*phi) * radial(2*s, (i-2*s)/2) * hermitepoly(k-i) * sqrt(1.0/M_PI);
                else
                    shape(ii++) = cos(2*s*phi) * radial(2*s, (i-2*s)/2) * hermitepoly(k-i) * sqrt(2.0/M_PI);
            for(int i=2*s;i<=k;i+=2)
                if (s != 0)
                    shape(ii++) = sin(2*s*phi) * radial(2*s, (i-2*s)/2) * hermitepoly(k-i) * sqrt(2.0/M_PI);
            for(int i=2*s+1;i<=k;i+=2)
                shape(ii++) = cos((2*s+1)*phi) * radial(2*s+1, (i-(2*s+1))/2) * hermitepoly(k-i) * sqrt(2.0/M_PI);
            for(int i=2*s+1;i<=k;i+=2)
                shape(ii++) = sin((2*s+1)*phi) * radial(2*s+1, (i-(2*s+1))/2) * hermitepoly(k-i) * sqrt(2.0/M_PI);
        }
        return;
    }    
}

template<int dim>
void Distribution<dim>::Set(const std::function<double(IntegrationPoint &)> function, FlatVector<> coeffs) const
{
    Array<double> nodes,weights;
    HermiteRule(nodes,weights);
    int ii=0;
    for(int i=0;i<=order;i++)
        for(int j=0;j<=order;j++)
        {
            if(dim == 3)
                for(int k=0;k<=order;k++)
                {
                    //cout << "\rSetting "<<i<<" / "<<order<<flush;
                    IntegrationPoint ip(nodes[i],nodes[j],nodes[k]);
                    IntegrationPoint iptrans(sqrt(ansatztemp)*nodes[i]+ansatzv(0),sqrt(ansatztemp)*nodes[j]+ansatzv(1),
                                             sqrt(ansatztemp)*nodes[k]+ansatzv(2));
                    coeffs(ii++) = function(iptrans)*exp(L2Norm2(FlatVector<>(dim,&ip(0))));
                }
            else if(dim==2)
            {
                IntegrationPoint ip(nodes[i],nodes[j]);
                IntegrationPoint iptrans(sqrt(ansatztemp)*nodes[i]+ansatzv(0),sqrt(ansatztemp)*nodes[j]+ansatzv(1));
                coeffs(ii++) = function(iptrans)*exp(L2Norm2(FlatVector<>(dim,&ip(0))));
            }
        }
}






template<>
INLINE void Distribution<3>::Tensor2Hierarchical(FlatTensor<3> tensorrep, FlatVector<> vectorrep) const
{
    for(int l=0,ii=0;l<=order;l++)
        for(int m=0;m<=l;m++)
            for(int n=0;n<=m;n++,ii++)
                vectorrep(ii) = tensorrep(n,m-n,l-m);
}

template<>
INLINE void Distribution<2>::Tensor2Hierarchical(FlatTensor<2> tensorrep, FlatVector<> vectorrep) const
{
    for(int l=0,ii=0;l<=order;l++)
        for(int m=0;m<=l;m++,ii++)
            vectorrep(ii) = tensorrep(m,l-m);
}
template<>
INLINE void Distribution<1>::Tensor2Hierarchical(FlatTensor<1> tensorrep, FlatVector<> vectorrep) const
{
    for(int i=0;i<=order;i++)
            vectorrep(i) = tensorrep(i);
}

template<>
INLINE void Distribution<3>::Hierarchical2Tensor(FlatTensor<3> tensorrep, FlatVector<> vectorrep) const
{
    tensorrep = 0.0;
    for(int l=0,ii=0;l<=order;l++)
        for(int m=0;m<=l;m++)
            for(int n=0;n<=m;n++,ii++)
                tensorrep(n,m-n,l-m) = vectorrep(ii);
}

template<>
INLINE void Distribution<2>::Hierarchical2Tensor(FlatTensor<2> tensorrep, FlatVector<> vectorrep) const
{
    tensorrep = 0.0;
    for(int l=0,ii=0;l<=order;l++)
        for(int m=0;m<=l;m++,ii++)
            tensorrep(m,l-m) = vectorrep(ii);
}

template<>
INLINE void Distribution<1>::Hierarchical2Tensor(FlatTensor<1> tensorrep, FlatVector<> vectorrep) const
{
    for(int i=0;i<=order;i++)
        tensorrep(i) = vectorrep(i);
}

template<>
void Distribution<3>::N2H(FlatTensor<3> nodal, FlatTensor<3> hermite, LocalHeap & lh) const
{
    timern2h.Start();
    for (int i = 0; i <= order; i++)        
        for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
                nodal(i,j,k) *= scale(i)*scale(j)*scale(k);

    FlatTensor<3> temp(lh,order+1,order+1,order+1);
    FlatMatrix<> tempmat(order+1,order+1,lh);
    for(int l=0;l<= order;l++)
    {
        tempmat = nodal(l,STAR,STAR)*nodal2hermite | Lapack;
        temp(l,STAR,STAR) = Trans(nodal2hermite)*tempmat | Lapack;
    }
    for(int l=0;l<= order;l++)
        hermite(STAR,l,STAR) = Trans(nodal2hermite)*temp(STAR,l,STAR) | Lapack;
    timern2h.Stop();
    timern2h.AddFlops(3*pow(order+1,4) + pow(order+1,3));
}

template<>
void Distribution<3>::N2H(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const
{
    FlatTensor<3> nodalt(order+1,order+1,order+1,&nodal(0));
    Tensor<3> hermitet(order+1,order+1,order+1);
    N2H(nodalt,hermitet,lh);
    Tensor2Hierarchical(hermitet,hermite);
}

template<>
void Distribution<2>::N2H(FlatTensor<2> nodal, FlatTensor<2> hermite, LocalHeap & lh) const
{
    for (int i = 0, ii = 0; i <= order; i++)
        for (int j = 0; j <= order; j++, ii++)
            nodal(i,j) *= scale(i)*scale(j);
    FlatMatrix<> temp(order+1,order+1,lh);
    temp = nodal(STAR,STAR)*nodal2hermite | Lapack;
    hermite(STAR,STAR) = Trans(nodal2hermite)*temp | Lapack;
}

template<>
void Distribution<2>::N2H(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const
{
    FlatTensor<2> nodalt(order+1,order+1,&nodal(0));
    Tensor<2> hermitet(order+1,order+1);
    N2H(nodalt,hermitet,lh);
    Tensor2Hierarchical(hermitet,hermite);
}

template<>
void Distribution<1>::N2H(FlatTensor<1> nodal, FlatTensor<1> hermite, LocalHeap & lh) const
{
    for (int i = 0; i <= order; i++)
            nodal(i) *= scale(i);
    FlatVector<> temp(order+1,lh);
    hermite(STAR) = nodal(STAR)*nodal2hermite;
}

template<>
void Distribution<1>::N2H(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const
{
    hermite = Trans(nodal2hermite)*nodal;
}


template<>
void Distribution<3>::H2N(FlatTensor<3> nodal, FlatTensor<3> hermite, LocalHeap & lh) const
{
    FlatTensor<3> temp(lh,order+1,order+1,order+1);
    FlatMatrix<> tempmat(order+1,order+1,lh);
    for(int l=0;l<= order;l++)
        temp(STAR,l,STAR) = Trans(hermite2nodal)*hermite(STAR,l,STAR) | Lapack;
    for(int l=0;l<= order;l++)
    {
        tempmat = Trans(hermite2nodal)*temp(l,STAR,STAR) | Lapack;
        nodal(l,STAR,STAR) = tempmat*hermite2nodal | Lapack;
    }
    for (int i = 0; i <= order; i++)
        for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
                nodal(i,j,k) /= (scale(i)*scale(j)*scale(k));
}

template<>
void Distribution<3>::H2N(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const
{
    Tensor<3> hermitet(order+1,order+1,order+1);
    FlatTensor<3> nodalt(order+1,order+1,order+1,&nodal(0));
    hermitet = 0.0;
    Hierarchical2Tensor(hermitet,hermite);
    H2N(nodalt,hermitet,lh);
}

template<>
void Distribution<2>::H2N(FlatTensor<2> nodal, FlatTensor<2> hermite, LocalHeap & lh) const
{
    FlatTensor<2> temp(lh,order+1,order+1);
    FlatMatrix<> tempmat(order+1,order+1,lh);
    tempmat = Trans(hermite2nodal)*temp(STAR,STAR) | Lapack;
    nodal(STAR,STAR) = tempmat*hermite2nodal | Lapack;
    for (int i = 0; i <= order; i++)
        for (int j = 0; j <= order; j++)
                nodal(i,j) /= (scale(i)*scale(j));
}

template<>
void Distribution<2>::H2N(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const
{
    Tensor<2> hermitet(order+1,order+1);
    FlatTensor<2> nodalt(order+1,order+1,&nodal(0));
    hermitet = 0.0;
    Hierarchical2Tensor(hermitet,hermite);
    H2N(nodalt,hermitet,lh);
}

template<>
void Distribution<1>::H2N(FlatTensor<1> nodal, FlatTensor<1> hermite, LocalHeap & lh) const
{
    nodal(STAR) = Trans(hermite2nodal)*hermite(STAR);// | Lapack;
    for (int i = 0; i <= order; i++)
                nodal(i) /= scale(i);
}

template<>
void Distribution<1>::H2N(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const
{
    Tensor<1> hermitet(order+1);
    FlatTensor<1> nodalt(order+1,&nodal(0));
    hermitet = 0.0;
    Hierarchical2Tensor(hermitet,hermite);
    H2N(nodalt,hermitet,lh);
}

template<>
void Distribution<3>::N2HTrans(FlatTensor<3> nodal, FlatTensor<3> hermite, LocalHeap & lh) const
{
    timern2htrans.Start();
    FlatTensor<3> temp(lh,order+1,order+1,order+1);
    FlatMatrix<> tempmat(order+1,order+1,lh);
    for(int l=0;l<= order;l++)
    {
        tempmat = hermite(l,STAR,STAR)*Trans(nodal2hermite) | Lapack;
        temp(l,STAR,STAR) = nodal2hermite*tempmat | Lapack;
    }
    for(int l=0;l<= order;l++)
        nodal(STAR,l,STAR) = nodal2hermite*temp(STAR,l,STAR) | Lapack;

    for (int i = 0; i <= order; i++)
        for (int j = 0; j <= order; j++)
            for (int k = 0; k <= order; k++)
                nodal(i,j,k) *= scale(i)*scale(j)*scale(k);
    timern2htrans.Stop();
    timern2htrans.AddFlops(3*pow(order+1,4) + pow(order+1,3));
}

template<>
void Distribution<3>::N2HTrans(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const
{
    Tensor<3> hermitet(order+1,order+1,order+1);
    hermitet = 0.0;
    FlatTensor<3> nodalt(order+1,order+1,order+1,&nodal(0));
    Hierarchical2Tensor(hermitet,hermite);
    N2HTrans(nodalt,hermitet,lh);
}

template<>
void Distribution<2>::N2HTrans(FlatTensor<2> nodal, FlatTensor<2> hermite, LocalHeap & lh) const
{
    timern2htrans.Start();
    FlatTensor<2> temp(lh,order+1,order+1);
    FlatMatrix<> tempmat(order+1,order+1,lh);
        tempmat = hermite(STAR,STAR)*Trans(nodal2hermite) | Lapack;
        temp(STAR,STAR) = nodal2hermite*tempmat | Lapack;
    
    for (int i = 0; i <= order; i++)
        for (int j = 0; j <= order; j++)
                nodal(i,j) *= scale(i)*scale(j);
    timern2htrans.Stop();
}

template<>
void Distribution<2>::N2HTrans(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const
{
    Tensor<2> hermitet(order+1,order+1);
    hermitet = 0.0;
    FlatTensor<2> nodalt(order+1,order+1,&nodal(0));
    Hierarchical2Tensor(hermitet,hermite);
    N2HTrans(nodalt,hermitet,lh);
}

template<>
void Distribution<1>::N2HTrans(FlatTensor<1> nodal, FlatTensor<1> hermite, LocalHeap & lh) const
{
    timern2htrans.Start();
    FlatTensor<1> temp(lh,order+1);
    FlatMatrix<> tempmat(order+1,lh);
        nodal(STAR) = nodal2hermite*hermite(STAR);
    for (int i = 0; i <= order; i++)
        nodal(i) *= scale(i);
    timern2htrans.Stop();
}

template<>
void Distribution<1>::N2HTrans(FlatVector<> nodal, FlatVector<> hermite, LocalHeap & lh) const
{
    Tensor<1> hermitet(order+1);
    hermitet = 0.0;
    FlatTensor<1> nodalt(order+1,&nodal(0));
    Hierarchical2Tensor(hermitet,hermite);
    N2HTrans(nodalt,hermitet,lh);
}

template<int dim>
void Distribution<dim>::H2PH(FlatVector<> hermite,FlatVector<> polarhermite,LocalHeap & lh) const
{
    timerh2ph.Start();
    // FlatMatrix<> hermitem((order+1)*(order+2)/2,order/2+1,lh);
    // FlatMatrix<> polarm((order+1)*(order+2)/2,order/2+1,lh);
    FlatMatrix<> hermitem((order+1)*(order+2)/2,order+1,lh);
    FlatMatrix<> polarm((order+1)*(order+2)/2,order+1,lh);
    hermitem = 0.0;
    polarm = 0.0;
    for(int k=0,ii=0;k<order+1;k++)
    {
        //if( k % 2 == 0 )
        //    hermitem.Col(k/2).Rows(0,(k+1)*(k+2)/2) = hermite.Range(ii,ii+(k+1)*(k+2)/2);
        hermitem.Col(k).Rows(0,(k+1)*(k+2)/2) = hermite.Range(ii,ii+(k+1)*(k+2)/2);
        ii += (k+1)*(k+2)/2;
    }
    // 0-th order and first order trafos
    polarm.Row(0) = hermitem.Row(0);
    polarm.Row(2) = hermitem.Row(1);
    polarm.Row(1) = hermitem.Row(2);
    // higher order trafos:
    for(int k=2,ii=3;k<order+1;k++)
    {
        polarm.Rows(ii,ii+k+1) = Trans(*hermite2polarhermites[k]) * hermitem.Rows(ii,ii+k+1) | Lapack;
        ii += k+1;
    }
    for(int k=0,ii=0;k<order+1;k++)
    {
        //if( k % 2 == 0)
        //    polarhermite.Range(ii,ii+(k+1)*(k+2)/2) = polarm.Col(k/2).Rows(0,(k+1)*(k+2)/2);
        polarhermite.Range(ii,ii+(k+1)*(k+2)/2) = polarm.Col(k).Rows(0,(k+1)*(k+2)/2);
        ii += (k+1)*(k+2)/2;
    }    
    timerh2ph.Stop();
    timerh2ph.AddFlops(pow(order+1,2)*(order+2)*(2*order+3)/6);
}

template<int dim>
void Distribution<dim>::PH2P(FlatVector<> polarhermite,FlatVector<> polar,LocalHeap & lh) const
{
    timerph2p.Start();
    FlatVector<> polartemp(GetNDof<POLAR>(),lh);
    for(int i=0;i<GetNDof<POLAR>();i++)
        polartemp(i) = polarhermite( phtoph2[i] );
    int iii = 0;
    polar = 0.0;
    for(int k=0,ii=0;k<order+1;k++)
    {
        FlatVector<> coeff_k((k+1)*(k+2)/2,&polartemp(ii));
        int blstart = 0;
        int blsize = k/2+1;
        // The first trafo is of size k/2+1:
        polar.Range(ii+blstart,ii+blstart+blsize) = (*polarhermites2polar[iii])*coeff_k.Range(blstart,blstart+blsize);
        iii++;
        blstart += blsize;
        // For odd k, there two additional trafos of size k/2+1
        if(k%2 == 1) 
        {
            blsize = k/2+1;
            FlatMatrix<double> tempin(2,blsize,&coeff_k(blstart));
            FlatMatrix<double> tempout(2,blsize,&polar(ii+blstart));
            tempout = tempin * Trans( *polarhermites2polar[iii] );
            iii+=2;
            blstart += 2*blsize;
        }
        for(int blsize = k/2;blsize>=1;blsize--)
        {
            for(int bl = 0;bl<2;bl++)
            {
                FlatMatrix<double> tempin(2,blsize,&coeff_k(blstart));
                FlatMatrix<double> tempout(2,blsize,&polar(ii+blstart));
                if( k % 2 == 0 )
                    tempout = tempin * Trans( *polarhermites2polar[iii] );
                iii+=2;
                blstart += 2*blsize;
            }
        }
        ii += (k+1)*(k+2)/2;// + (k+2)*(k+3)/2;
    }
    timerph2p.Stop();
    timerph2p.AddFlops(nops_ph2p);
}

template<int dim>
inline void Distribution<dim>::H2P(FlatVector<> hermite, FlatVector<> polar,LocalHeap & lh) const
{
    timerh2p.Start();
    FlatVector<> temp(GetNDof<POLAR>(),lh);
    H2PH(hermite,temp,lh);
    PH2P(temp,polar,lh);
    timerh2p.Stop();
}

template<int dim>
void Distribution<dim>::H2PHTrans(FlatVector<> hermite,FlatVector<> polarhermite,LocalHeap & lh) const
{
    timerh2phtrans.Start();
    //FlatMatrix<> hermitem((order+1)*(order+2)/2,order/2+1,lh);
    //FlatMatrix<> polarhermitem((order+1)*(order+2)/2,order/2+1,lh);
    FlatMatrix<> hermitem((order+1)*(order+2)/2,order+1,lh);
    FlatMatrix<> polarhermitem((order+1)*(order+2)/2,order+1,lh);
    hermitem = 0.0;
    polarhermitem = 0.0;
    for(int k=0,ii=0;k<order+1;k++)
    {
        // if( k % 2 == 0 )
            // polarhermitem.Col(k/2).Rows(0,(k+1)*(k+2)/2) = polarhermite.Range(ii,ii+(k+1)*(k+2)/2);
        if( k % 2 == 0 )
            polarhermitem.Col(k).Rows(0,(k+1)*(k+2)/2) = polarhermite.Range(ii,ii+(k+1)*(k+2)/2);
        ii += (k+1)*(k+2)/2;
    }
    // 0-th order and first order trafos
    hermitem.Row(0) = polarhermitem.Row(0);
    hermitem.Row(2) = polarhermitem.Row(1);
    hermitem.Row(1) = polarhermitem.Row(2);
    // higher order trafos:
    for(int k=2,ii=3;k<order+1;k++)
    {
        hermitem.Rows(ii,ii+k+1) = (*hermite2polarhermites[k]) * polarhermitem.Rows(ii,ii+k+1) | Lapack;
        ii += k+1;
    }
    for(int k=0,ii=0;k<order+1;k++)
    {
        // if( k % 2 == 0 )
            // hermite.Range(ii,ii+(k+1)*(k+2)/2) = hermitem.Col(k/2).Rows(0,(k+1)*(k+2)/2);
        hermite.Range(ii,ii+(k+1)*(k+2)/2) = hermitem.Col(k).Rows(0,(k+1)*(k+2)/2);
            
        ii += (k+1)*(k+2)/2;
    }    
    timerh2phtrans.Stop();
    timerh2phtrans.AddFlops(pow(order+1,2)*(order+2)*(2*order+3)/6);
}

template<int dim>
void Distribution<dim>::PH2PTrans(FlatVector<> polarhermite,FlatVector<> polar,LocalHeap & lh) const
{
    timerph2ptrans.Start();
    FlatVector<> polarhermitetemp(GetNDof<POLAR>(),lh);
    polarhermitetemp = 0.0;
    int iii = 0;
    for(int k=0,ii=0;k<order+1;k++)
    {
        FlatVector<> coeff_k((k+1)*(k+2)/2,&polar(ii));
        int blstart = 0;
        for(int bl = 0;bl<1;bl ++)
        {
            int blsize = k/2+1;
            // The first trafo is of size k/2+1:
            polarhermitetemp.Range(ii+blstart,ii+blstart+blsize) = Trans(*polarhermites2polar[iii])*coeff_k.Range(blstart,blstart+blsize);
            iii++;
            blstart += blsize;
        }
        // For odd k, there two additional trafos of size k/2+1
        if(k%2 == 1) 
        {
            int blsize = k/2+1;
            FlatMatrix<double> tempin(2,blsize,&coeff_k(blstart));
            FlatMatrix<double> tempout(2,blsize,&polarhermitetemp(ii+blstart));
            tempout = tempin * ( *polarhermites2polar[iii] );
            iii+=2;
            blstart += 2*blsize;
        }
        for(int blsize = k/2;blsize>=1;blsize--)
        {
            for(int bl = 0;bl<2;bl++)
            {
                FlatMatrix<double> tempin(2,blsize,&coeff_k(blstart));
                FlatMatrix<double> tempout(2,blsize,&polarhermitetemp(ii+blstart));
                if( k % 2 == 0 )
                    tempout = tempin * ( *polarhermites2polar[iii] );
                iii+=2;
                blstart += 2*blsize;
            }
        }
        ii += (k+1)*(k+2)/2;// + (k+2)*(k+3)/2;
    }
    for(int i=0;i<GetNDof<POLAR>();i++)
        polarhermite( phtoph2[i] ) = polarhermitetemp( i );
    timerph2ptrans.Stop();
    timerph2ptrans.AddFlops(nops_ph2p);
}

template<int dim>
inline void Distribution<dim>::H2PTrans(FlatVector<> hermite, FlatVector<> polar,LocalHeap & lh) const
{
    timerh2ptrans.Start();
    FlatVector<> temp(GetNDof<POLAR>(),lh);
    PH2PTrans(temp,polar,lh);
    H2PHTrans(hermite,temp,lh);
    timerh2ptrans.Stop();
}

template<int dim>
int Distribution<dim>::GetDimension() const
{
    return dim;
}

template<int dim>
inline void Distribution<dim>::SolveM(FlatVector<> vecin)
{
    for(int i=0;i<GetNDof<NODAL>();i++)
        vecin(i) *= exp(log(GetAnsatzTemp())*dim/2.0)*diagmassinv(i);
}

template<int dim>
template<Representation REP>
void Distribution<dim>::TestMass(FlatMatrix<> mass)
{
    Vector<> shape(GetNDof<REP>());
    // mass matrix:
    Array<double> nodes,weights;
    HermiteRule(nodes,weights);
    mass = 0.0;
    if(dim == 3)
        for(int i = 0;i<nodes.Size();i++)
            for(int j = 0;j<nodes.Size();j++)
                for(int k = 0;k<nodes.Size();k++)
                {
                    CalcShape<REP>(IntegrationPoint(nodes[i],nodes[j],nodes[k]),shape);                
                    mass += weights[i]*weights[j]*weights[k]*shape*Trans(shape);
                }
    else if (dim == 2)
        for(int i = 0;i<nodes.Size();i++)
            for(int j = 0;j<nodes.Size();j++)
            {
                CalcShape<REP>(IntegrationPoint(nodes[i],nodes[j]),shape);
                mass += weights[i]*weights[j]*shape*Trans(shape);
            }
        
    for(int i=0;i<GetNDof<REP>();i++)
        for(int j=0;j<GetNDof<REP>();j++)
            if(fabs(mass(i,j))<1e-12)
                mass(i,j)=0.0;
}

template<int dim>
template<Representation REP>
void Distribution<dim>::Evaluate(FlatVector<> coeffs, IntegrationRule & pts, FlatVector<> values,double Tref,FlatVector<> Vref) const
{
    Matrix<> shapes(pts.Size(),GetNDof<REP>());
    IntegrationRule ir;
    for(int i=0;i<pts.Size();i++)
    {
        IntegrationPoint temp(1.0/sqrt(Tref)*( (pts[i])(0) - Vref(0) ),1.0/sqrt(Tref)*( (pts[i])(1) - Vref(1) ),
                              1.0/sqrt(Tref)*( (pts[i])(2) - Vref(2) ));
        CalcShape1<REP>(temp,shapes.Row(i));
        ir.Append(temp);
    }
    values = shapes*coeffs;
    if(use_hm_funcs)
      return;   
    for(int i=0;i<pts.Size();i++)
      values(i)*=exp(-0.5*L2Norm2(FlatVector<>(dim,&ir[i](0))));
    return;
}

template<int dim>
template<Representation REP>
double Distribution<dim>::Evaluate(FlatVector<> coeffs, IntegrationPoint pt)
{
    Vector<> shape(GetNDof<REP>());
    IntegrationPoint temp(1.0/sqrt(ansatztemp)*(pt(0)-ansatzv(0)),1.0/sqrt(ansatztemp)*(pt(1)-ansatzv(1)),
                          1.0/sqrt(ansatztemp)*(pt(2)-ansatzv(2)));
    CalcShape1<REP>(temp,shape);
    shape*=exp(-L2Norm2(FlatVector<>(dim,&temp(0))));
    return InnerProduct(shape,coeffs);
}

template<>
template<Representation REP>
void Distribution<3>::Macroscopics(FlatVector<> coeff,FlatVector<> result,double ansatztemp, FlatVector<> ansatzv,LocalHeap & lh) const
{
    FlatVector<> mycoeffs(5,lh);
    if(REP == NODAL)
    {
        FlatVector<> coeff1(coeff.Size(),lh);
        coeff1 = coeff;
        FlatTensor<3> nodalt(order+1,order+1,order+1,&coeff1(0));
        FlatTensor<3> hermitet(lh,order+1,order+1,order+1);
        FlatVector<> hermitev(GetNDof<HERMITE>(),lh);
        N2H(nodalt,hermitet,lh);
        Tensor2Hierarchical(hermitet,hermitev);
        mycoeffs(0) = hermitev(0);
        mycoeffs(1) = hermitev(3);
        mycoeffs(2) = hermitev(2);
        mycoeffs(3) = hermitev(1);
        mycoeffs(4) = hermitev(4) + hermitev(7) + hermitev(9);
    }
    else if(REP == HERMITE)
    {
        mycoeffs(0) = coeff(0);
        mycoeffs(1) = coeff(3);
        mycoeffs(2) = coeff(2);
        mycoeffs(3) = coeff(1);
        mycoeffs(4) = coeff(4) + coeff(7) + coeff(9);
    }
    double c = sqrt(sqrt( M_PI*M_PI*M_PI ) * ansatztemp * ansatztemp * ansatztemp);
    double density = mycoeffs(0) * c;
    double vx = ansatzv(0) + mycoeffs(1) * ansatztemp * c/(sqrt(2.0)*density);
    double vy = ansatzv(1) + mycoeffs(2) * ansatztemp * c/(sqrt(2.0)*density);
    double vz = ansatzv(2) + mycoeffs(3) * ansatztemp * c/(sqrt(2.0)*density);
    double Temp = 1.0/3.0 * ( 1.5 * ansatztemp + c/(density*sqrt(2.0)) * pow(ansatztemp,2.5) * mycoeffs(4) - 
                            ( pow(ansatzv(0)-vx,2) + pow(ansatzv(1)-vy,2) +pow(ansatzv(2)-vz,2)   ) );
    result(0) = density;
    result(1) = vx; result(2) = vy; result(3) = vz;
    result(4) = Temp;
}

template<>
template<Representation REP>
void Distribution<1>::Macroscopics(FlatVector<> coeff,FlatVector<> result,double aansatztemp,FlatVector<> aansatzv,LocalHeap & lh) const
{
    FlatVector<> mycoeffs(3,lh);
    if(REP == NODAL)
    {
        FlatVector<> coeff1(coeff.Size(),lh);
        coeff1 = coeff;
        FlatVector<> hermite(GetNDof<NODAL>(),lh);
        //N2H(coeff1,hermite,lh);
        if(use_hm_funcs == false)
            for(int i = 0;i<coeff1.Size();i++)
                coeff1(i)*=scale(i);
        hermite = Trans(nodal2hermite)*coeff1;
        mycoeffs(0) = hermite(0);
        mycoeffs(1) = hermite(1);
        mycoeffs(2) = hermite(2);
    }
    else if(REP == HERMITE)
    {
        mycoeffs(0) = coeff(0);
        mycoeffs(1) = coeff(1);
        mycoeffs(2) = coeff(2);
    }
    double rho0 = sqrt(sqrt( M_PI ) * aansatztemp ); //   \sqrt(T)*\int     H_0(v)*exp(-v^2)
    double v0 = rho0*sqrt(aansatztemp*0.5);          //          T*\int v   H_1(v)*exp(-v^2)
    double e0 = rho0*aansatztemp*sqrt(0.5);          // \sqrt(T)*T*\int v^2 H_2(v)*exp(-v^2)
    double m0 = mycoeffs(0) * rho0; // passt
    double m10 = mycoeffs(0)*aansatzv(0)*rho0 + mycoeffs(1)* v0;
    double m222= ( rho0*aansatzv(0)*aansatzv(0) + 0.5*aansatztemp*rho0 )*mycoeffs(0) + 2.0*aansatzv(0)*v0*mycoeffs(1) + e0*mycoeffs(2);
    result(0) = m0;
    result(1) = m10;
    result(2) = m222;
}

template<int dim>
double Distribution<dim>::BKWSolution(double s, IntegrationPoint & ip) const
{
    // BKW solution
    return pow(2*M_PI*s,-dim/2.0)*(1.0-(1.0-s)/(2.0*s)* (dim*s-ip(0)*ip(0)-ip(1)*ip(1)-ip(2)*ip(2))/s )*exp(- 1.0/(2.0*s)*( ip(0)*ip(0) + ip(1)*ip(1) + ip(2)*ip(2) ) );
}

template<int dim>
double Distribution<dim>::TwoMaxwelliansSolution(double rho1, double rho2, Vec<dim> V1, Vec<dim> V2, double T1, double T2, IntegrationPoint & ip) const
{
    // Sum of 2 Maxwellian peaks
    if(dim == 3)
    {
      Vec<3> vminusV1(ip(0)-V1(0),ip(1)-V1(1),ip(2)-V1(2));
      Vec<3> vminusV2(ip(0)-V2(0),ip(1)-V2(1),ip(2)-V2(2));
      return rho1*pow(2.0*M_PI*T1, -dim/2.0)*exp( - 1.0/(2.0*T1)*L2Norm2(vminusV1) ) + rho2*pow(2.0*M_PI*T2, -dim/2.0)*exp( - 1.0/(2.0*T2)*L2Norm2(vminusV2) );
    }
    if(dim == 2)
    {
      Vec<2> vminusV1(ip(0)-V1(0),ip(1)-V1(1));
      Vec<2> vminusV2(ip(0)-V2(0),ip(1)-V2(1));
      return rho1*pow(2.0*M_PI*T1, -dim/2.0)*exp( - 1.0/(2.0*T1)*L2Norm2(vminusV1) ) + rho2*pow(2.0*M_PI*T2, -dim/2.0)*exp( - 1.0/(2.0*T2)*L2Norm2(vminusV2) );
    }
    if(dim == 1)
    {
      Vec<1> vminusV1(ip(0)-V1(0));
      Vec<1> vminusV2(ip(0)-V2(0));
      return rho1*pow(2.0*M_PI*T1, -dim/2.0)*exp( - 1.0/(2.0*T1)*L2Norm2(vminusV1) ) + rho2*pow(2.0*M_PI*T2, -dim/2.0)*exp( - 1.0/(2.0*T2)*L2Norm2(vminusV2) );
    }    
}
/*
template<>
void Distribution<3>::SetHighOrder(int highorder,const std::function<double(IntegrationPoint &)> function, FlatVector<>coeffs, LocalHeap &lh)
{
    Distribution<3> disttemp(highorder,false);
    Vector<> coeffstempn(disttemp.GetNDof<NODAL>());
    disttemp.SetAnsatzTemp(GetAnsatzTemp());
    Vec<3> AnsatzV;
    GetAnsatzV(AnsatzV);
    disttemp.SetAnsatzV(AnsatzV);
    disttemp.Set( [&](IntegrationPoint & ip){return function(ip);},coeffstempn );

    Tensor<3> coeffstempht(highorder+1,highorder+1,highorder+1);
    FlatTensor<3> coeffstempnt(highorder+1,highorder+1,highorder+1,&coeffstempn(0));

    disttemp.N2H(coeffstempnt,coeffstempht,lh);
    Tensor<3> coeffstemphtsmall(order+1,order+1,order+1);
    for(int i=0;i<order+1;i++)
        coeffstemphtsmall(i,STAR,STAR) = coeffstempht(i,STAR,STAR).Cols(0,order+1).Rows(0,order+1);
    FlatTensor<3> coeffsnt(order+1,order+1,order+1,&coeffs(0));
    H2N(coeffsnt,coeffstempht,lh);

    Array<double> nodes,weights;
    disttemp.HermiteRule(nodes,weights);
    Vector<> error(pow(nodes.Size(),3));
    Array<IntegrationPoint * > pts;
    for(int i=0,ii=0;i<nodes.Size();i++)
        for(int j=0;j<nodes.Size();j++)
            for(int k=0;k<nodes.Size();k++)
            {
                pts.Append(new IntegrationPoint(nodes[i],nodes[j],nodes[k]));
                IntegrationPoint ip(nodes[i],nodes[j],nodes[k]);
                error(ii++) = function(ip);
            }
    Vector<> values(pts.Size());
    Evaluate<NODAL>(coeffs,pts,values);
    error -=values;
    double L2error = 0.0;
    for(int i=0,ii=0;i<nodes.Size();i++)
        for(int j=0;j<nodes.Size();j++)
            for(int k=0;k<nodes.Size();k++)
                L2error += sqr(error(ii++))*weights[i]*weights[j]*weights[k]*exp(L2Norm2(Vec<3>(nodes[i],nodes[j],nodes[k])));

    cout<<"-------------------- Maximum nodal error in Initial Distribution = "<<MaxNorm(error)<<" --------------------"<<endl;
    cout<<"-------------------- Maximum L2 error in Initial Distribution =    "<<L2error<<" --------------------"<<endl;
}
*/
template  class Distribution<3>;
template  class Distribution<2>;
template  class Distribution<1>;

template void Distribution<3>::Evaluate<NODAL>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<3>::Evaluate<HERMITE>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<3>::Evaluate<POLAR>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<3>::Evaluate<POLARHERMITE>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<3>::Evaluate<POLARHERMITE2>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template double Distribution<3>::Evaluate<HERMITE>(FlatVector<>, IntegrationPoint);
template double Distribution<3>::Evaluate<NODAL>(FlatVector<>, IntegrationPoint);

template void Distribution<2>::Evaluate<NODAL>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<2>::Evaluate<HERMITE>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<2>::Evaluate<POLAR>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<2>::Evaluate<POLARHERMITE>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<2>::Evaluate<POLARHERMITE2>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template double Distribution<2>::Evaluate<HERMITE>(FlatVector<>, IntegrationPoint);
template double Distribution<2>::Evaluate<NODAL>(FlatVector<>, IntegrationPoint);

template void Distribution<1>::Evaluate<NODAL>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<1>::Evaluate<HERMITE>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<1>::Evaluate<POLAR>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<1>::Evaluate<POLARHERMITE>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template void Distribution<1>::Evaluate<POLARHERMITE2>(FlatVector<>, IntegrationRule &, FlatVector<>, double, FlatVector<>) const;
template double Distribution<1>::Evaluate<HERMITE>(FlatVector<>, IntegrationPoint);
template double Distribution<1>::Evaluate<NODAL>(FlatVector<>, IntegrationPoint);

template void Distribution<3>::Macroscopics<NODAL>(FlatVector<> coeff,FlatVector<> result,double ansatztemp, FlatVector<> ansatzv,LocalHeap & lh) const;
template void Distribution<3>::Macroscopics<HERMITE>(FlatVector<> coeff,FlatVector<> result,double ansatztemp, FlatVector<> ansatzv,LocalHeap & lh) const;

//template void Distribution<2>::Macroscopics<NODAL>(FlatVector<> coeff,FlatVector<> result,LocalHeap & lh) const;
//template void Distribution<2>::Macroscopics<HERMITE>(FlatVector<> coeff,FlatVector<> result,LocalHeap & lh) const;
template void Distribution<1>::Macroscopics<NODAL>(FlatVector<> coeff,FlatVector<> result,double ansatztemp, FlatVector<> ansatzv,LocalHeap & lh) const;
template void Distribution<1>::Macroscopics<HERMITE>(FlatVector<> coeff,FlatVector<> result,double ansatztemp, FlatVector<> ansatzv,LocalHeap & lh) const;

//template void Distribution<3>::CalcShape1<POLAR>(const IntegrationPoint & ip, SliceVector<> shape);
//template void Distribution<2>::CalcShape1<HERMITE>(const IntegrationPoint & ip, SliceVector<> shape);
//template void Distribution<1>::CalcShape1<POLAR>(const IntegrationPoint & ip, SliceVector<> shape);
