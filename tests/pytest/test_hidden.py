from netgen.geom2d import unit_square
from ngsolve import *


def test_hidden():
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))

    elim_options = [False, True]
    bddc_options = [False, True]
    hidden_options = [False, True]
    compress_options = [False, True]

    solutions = {}
    ndof = {}
    nze = {}
    nze_total = {}

    i=0

    option_tuples = []

    for elim_internal in elim_options:
        for use_bddc in bddc_options:
            if not elim_internal and use_bddc:
                continue
            for use_hidden in hidden_options:
                for compress in compress_options:
                    current_tuple = (elim_internal,use_bddc,use_hidden,compress)
                    option_tuples.append(current_tuple)
    for elim_internal,use_bddc,use_hidden,compress in option_tuples:
        current_tuple = (elim_internal,use_bddc,use_hidden,compress)
        i +=1
        order = 4
        fes1 = HDiv(mesh, order=order, discontinuous=True, hide_all_dofs=use_hidden)
        fes2 = L2(mesh, order=order-1)
        fes3 = FacetFESpace(mesh, order=order, dirichlet=[1,2,3])
        fes = FESpace([fes1,fes2,fes3])
        if compress:
            fes = CompressCompound(fes)

        sigma,u,uhat = fes.TrialFunction()
        tau,v,vhat = fes.TestFunction()

        n = specialcf.normal(mesh.dim)

        a = BilinearForm(fes, symmetric=False,
                         eliminate_hidden = use_hidden,
                         eliminate_internal = elim_internal)
        a += SymbolicBFI(sigma*tau + div(sigma)*v + div(tau)*u)
        a += SymbolicBFI(sigma*n*vhat+tau*n*uhat, element_boundary=True)

        if use_bddc:
            c = Preconditioner(type="bddc", bf=a)
        else:
            c = Preconditioner(type="direct", bf=a)
            
        a.Assemble()

        f = LinearForm(fes)
        f += SymbolicLFI(-1*v)
        f.Assemble()

        gfu = GridFunction(fes)

        inv = CGSolver(a.mat, c.mat, maxsteps=1000)

        if elim_internal:
            f.vec.data += a.harmonic_extension_trans * f.vec

        gfu.vec.data = inv * f.vec

        if elim_internal:
            gfu.vec.data += a.harmonic_extension * gfu.vec
            gfu.vec.data += a.inner_solve * f.vec
        
        solutions[current_tuple] = gfu
        ndof[current_tuple] = fes.ndof
        nze[current_tuple] = a.mat.nze
        if (elim_internal):
            nze_total[current_tuple] = a.mat.nze + a.harmonic_extension.nze + a.harmonic_extension_trans.nze + a.inner_solve.nze
        else:
            nze_total[current_tuple] = a.mat.nze
    
    
    for elim_internal,use_bddc,use_hidden,compress in option_tuples:
        print("({:1},{:1},{:1},{:1}), nze A(Schur) {:7}".format(elim_internal,use_bddc,use_hidden,compress,nze[(elim_internal,use_bddc,use_hidden,compress)]))

    for elim_internal,use_bddc,use_hidden,compress in option_tuples:
        print("({:1},{:1},{:1},{:1}), nze A(Schur) {:7}".format(elim_internal,use_bddc,use_hidden,compress,nze[(elim_internal,use_bddc,use_hidden,compress)]))
        
    for elim_internal,use_bddc,use_hidden,compress in option_tuples:
        print("({:1},{:1},{:1}), nze A(total) {:7}".format(elim_internal,use_bddc,use_hidden,compress,nze_total[(elim_internal,use_bddc,use_hidden,compress)]))

    for elim_internal,use_bddc,use_hidden,compress in option_tuples:
        print("({:1},{:1},{:1}), nze A(total) {:7}".format(elim_internal,use_bddc,use_hidden,compress,ndof[(elim_internal,use_bddc,use_hidden,compress)]))


    for elim_internal,use_bddc,use_hidden,compress in option_tuples:
        my = solutions[(elim_internal,use_bddc,use_hidden,compress)].components[1].vec
        diff = my.CreateVector()
        for elim_internal2,use_bddc2,use_hidden2,compress2 in option_tuples:
            diff.data = my - solutions[(elim_internal2,use_bddc2,use_hidden2,compress2)].components[1].vec
            error = Norm(diff)
            print("comparing ({:1},{:1},{:1},{:1}) with ({:1},{:1},{:1},{:1}), difference is {}".format(elim_internal,use_bddc,use_hidden,compress,elim_internal2,use_bddc2,use_hidden2,compress2,error))
            assert error < 1e-9

if __name__ == "__main__":
    test_hidden()
