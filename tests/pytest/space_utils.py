from ngsolve import *

xvec = [x,y,z]

def GetDiffOp(name, order, dim=1, dims=[], sym=False, dev=False, vb=VOL):
    if vb == BND and name != "id" and name.find("boundary") == -1:
        name = name +"boundary"
    #print("name = ", name, ", order = ", order, ", dim = ", dim, ", dims = ", dims, "sym = ", sym, ", vb =", vb)
    n_3d = specialcf.normal(3)
    t_3d = specialcf.tangential(3)
    Ptau_3d = Id(3) - OuterProduct(n_3d,n_3d)
    
    if name == "id":
        if len(dims) == 0 or (len(dims) == 1 and dims[0] == 1):
            if dim == 1:
                return CF( 3*x**order - x**(int(order/2)) )
            elif dim == 2:
                return CF( x**order + y**order + (x*y)**(int(order/2)) )
            elif dim == 3:
                return CF( x**order + y**order + z**order + (x*y)**(int(order/2)) + (x*z)**(int(order/2)) + (y*z)**(int(order/2)) + (x*y*z)**(int(order/3)) )
        elif len(dims) == 1:
            if dim == 2:
                return CF( (x**order+ 3*y**order-x**(int(order/2)),(x*y)**(int(order/2))) )
            elif dim == 3:
                return CF( (x**order-3*y**order,z**order+(x*y)**(int(order/2)), (x*z)**(int(order/2))+(y*z)**(int(order/2))+(x*y*z)**(int(order/3))) )
        elif len(dims) == 2:
            cf = CF( tuple( [xvec[i]**order - 3*xvec[j]**order + 5*(xvec[(i+1)%dim]*xvec[(j+2)%dim])**int(order/2) + (4 if i == j else 0) + (1/3*xvec[(i+1)%dim]**2 if i == 0 and j == 0 else 0) for i in range(dim) for j in range(dim)]), dims=(dim,dim) )
            if sym: cf = Sym(cf)
            if dev: cf = Deviator(cf)
            return cf
    elif name == "grad":
        cf = GetDiffOp("id", order, dim, dims, sym, dev)
        if len(dims) == 0 or (len(dims) == 1 and dims[0] == 1):
            return CF( tuple( [cf.Diff(xvec[i]) for i in range(dim)] ) )
        elif len(cf.dims) == 1:
            return CF( tuple( [cf[i].Diff(xvec[j]) for i in range(cf.dims[0]) for j in range(cf.dims[0])] ), dims=(cf.dims[0],cf.dims[0]))
        else:
            return CF( tuple( [cf[j,k].Diff(xvec[i]) for i in range(cf.dims[0]) for j in range(cf.dims[0]) for k in range(cf.dims[0])] ), dims=(cf.dims[0], cf.dims[0]**2) )
    elif name == "div":
        cf = GetDiffOp("id", order, dim, dims, sym, dev)
        if len(cf.dims) == 1:
            return sum(cf[i].Diff(xvec[i]) for i in range(cf.dims[0]))
        elif len(cf.dims) == 2:
            return CF( tuple( [sum(cf[i,j].Diff(xvec[j]) for j in range(cf.dims[0])) for i in range(cf.dims[0])] ) )
    elif name == "curl":
        cf = GetDiffOp("id", order, dim, dims, sym, dev)
        if len(cf.dims) == 1:
            if cf.dims[0] == 2:
                return cf[1].Diff(x) - cf[0].Diff(y)
            elif cf.dims[0] == 3:
                return CF( (cf[2].Diff(y)-cf[1].Diff(z), cf[0].Diff(z)-cf[2].Diff(x), cf[1].Diff(x)-cf[0].Diff(y)) )
        elif len(cf.dims)==2:
            if cf.dims[0] == 2:
                return CF( (cf[1,0].Diff(x)-cf[0,0].Diff(y), cf[1,1].Diff(x)-cf[0,1].Diff(y)) )
            elif cf.dims[0] == 3:
                return CF( (cf[0,2].Diff(y)-cf[0,1].Diff(z), cf[0,0].Diff(z)-cf[0,2].Diff(x), cf[0,1].Diff(x)-cf[0,0].Diff(y),
                            cf[1,2].Diff(y)-cf[1,1].Diff(z), cf[1,0].Diff(z)-cf[1,2].Diff(x), cf[1,1].Diff(x)-cf[1,0].Diff(y),
                            cf[2,2].Diff(y)-cf[2,1].Diff(z), cf[2,0].Diff(z)-cf[2,2].Diff(x), cf[2,1].Diff(x)-cf[2,0].Diff(y)), dims=(3,3) )
    elif name == "inc":
        cfcurl = GetDiffOp("curl", order, dim, dims, sym, dev)
        if dim == 2:
            return cfcurl[1].Diff(x) - cfcurl[0].Diff(y)
        elif dim == 3:
            cfcurl_T = cfcurl.trans
            return CF( (cfcurl_T[0,2].Diff(y)-cfcurl_T[0,1].Diff(z), cfcurl_T[0,0].Diff(z)-cfcurl_T[0,2].Diff(x), cfcurl_T[0,1].Diff(x)-cfcurl_T[0,0].Diff(y),
                         cfcurl_T[1,2].Diff(y)-cfcurl_T[1,1].Diff(z), cfcurl_T[1,0].Diff(z)-cfcurl_T[1,2].Diff(x), cfcurl_T[1,1].Diff(x)-cfcurl_T[1,0].Diff(y),
                         cfcurl_T[2,2].Diff(y)-cfcurl_T[2,1].Diff(z), cfcurl_T[2,0].Diff(z)-cfcurl_T[2,2].Diff(x), cfcurl_T[2,1].Diff(x)-cfcurl_T[2,0].Diff(y)), dims=(3,3) )
    elif name == "hesse":
        cf = GetDiffOp("id", order, dim, dims, sym, dev)
        if len(dims) == 0 or (len(dims) == 1 and dims[0] == 1):
            if dim == 1:
                return cf.Diff(x).Diff(x)
            elif dim == 2:
                return CF( (cf.Diff(x).Diff(x), cf.Diff(x).Diff(y), cf.Diff(y).Diff(x), cf.Diff(y).Diff(y)), dims=(2,2) )
            elif dim == 3:
                return CF( (cf.Diff(x).Diff(x), cf.Diff(x).Diff(y), cf.Diff(x).Diff(z),
                            cf.Diff(y).Diff(x), cf.Diff(y).Diff(y), cf.Diff(y).Diff(z),
                            cf.Diff(z).Diff(x), cf.Diff(z).Diff(y), cf.Diff(z).Diff(z)), dims=(3,3)  )
        else:
            return CF( tuple( [cf[i].Diff(xvec[j]).Diff(xvec[k]) for i in range(cf.dims[0]) for j in range(cf.dims[0]) for k in range(cf.dims[0])] ), dims=(dim,dim**2) )
    elif name == "gradboundary":
        cf = GetDiffOp("grad", order, dim, dims, sym, dev, vb=VOL)
        if len(dims) == 0 or (len(dims) == 1 and dims[0] == 1):
            return Ptau_3d*cf
        else:
            return cf*Ptau_3d
    elif name == "divboundary":
        return Trace(GetDiffOp("gradboundary", order, dim, dims, sym, dev, vb))
    elif name == "hesseboundary":
        cfgrad  = GetDiffOp("grad", order, dim, dims, sym, dev, vb=VOL)
        cfhesse = GetDiffOp("hesse", order, dim, dims, sym, dev, vb=VOL)
        n = specialcf.normal(dim)
        if len(dims) == 0 or (len(dims) == 1 and dims[0] == 1):
            return (-(cfgrad*n)*Grad(n) - OuterProduct(n,Grad(n)*cfgrad) + Ptau_3d*cfhesse)*Ptau_3d
        else:
            raise Exception("hesseboundary only for scalar implemented yet!")
    else:
        raise Exception("In GetDiffOp: Something went wrong: name =", name, ", order =", order, ", dim =", dim, ", dim =", dims, ", sym =", sym, ", dev =", dev, ", vb =", vb)

    # surface curl = surface div ( u times n) = (surface grad times u) cdot n
