from netgen import *
from ngsolve import *
from sys import argv

def printMeasurement( name ):
    timers = Timers()
    timer = [ t for t in timers if t['name'] == name ]
    if(len(timer) == 1):
        timer = timer[0]
    else:
        return
    time =  timer['time']

    name = name.replace(' ', '')
    name = name.replace(':', '')
    name = name.replace('-', '')

    if(time == 0):
        return
    print('<DartMeasurement name="{}"'.format(name))
    print('type="numeric/double">{}</DartMeasurement>'.format(time))

    mflops = timer['flops']*1e-6 #NgProfiler::GetFlops(nr) / NgProfiler::GetTime(nr) * 1e-6;
    if(mflops == 0):
        return
    print('<DartMeasurement name="{}MFlops"'.format(name))
    print('type="numeric/double">{}</DartMeasurement>'.format(mflops))

if(__name__ == '__main__'):
    if(len(argv)) < 2:
        print('Usage" {} filename [niterations] '.format(argv[0]))
        raise RuntimeError("invalid arguments")
    niterations = 1
    if(len(argv)>2):
        niterations = int(argv[2])
    print("CTEST_FULL_OUTPUT")

    pde = PDE(argv[1])
    for i in range(niterations):
        pde.Solve();

    printMeasurement( "Matrix assembling" );
    printMeasurement( "Solver - Total" );
    printMeasurement( "CG solver" );
    printMeasurement( "SparseMatrixSymmetric::MultAdd" );
    printMeasurement( "SparseMatrixSymmetric::MultAdd1" );
    printMeasurement( "SparseMatrixSymmetric::MultAdd2" );
    printMeasurement( "BaseVector::Scale (taskhandler)" );
    printMeasurement( "BaseVector::InnerProduct (taskhandler)" );
    printMeasurement( "BlockJacobiPrecond::GSSmoothBack" );
    printMeasurement( "BaseVector::Add (taskhandler)" );
    printMeasurement( "BlockJacobiPrecond::GSSmooth" );
    printMeasurement( "BaseVector::SetScalar (taskhandler)" );
    printMeasurement( "SparseMatrix::MultAdd (taskhandler)" );
    printMeasurement( "BaseVector::Set (taskhandler)" );

    for t in Timers():
        if t['flops']:
            t['name'] = 'MFlops = {:.2} {}'.format(t['flops']*1e-6, t['name'])
        print('calls{:>8}, time {:.4f} sec {}'.format(t['counts'], t['time'], t['name']))
