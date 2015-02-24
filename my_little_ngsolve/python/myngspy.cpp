#include <solve.hpp>
using namespace ngsolve;

#include <python_ngstd.hpp>


// a global function ...
void Hello()
{
  cout << "Hello everybody" << endl;
}




class NumProcPyDemo : public NumProc
{

public:
    
  NumProcPyDemo (PDE & apde, const Flags & flags)
    : NumProc (apde) { ; }

  virtual void Do(LocalHeap & lh)
  {
    cout << "solving NumProcPyDemo" << endl;
  }
  
  virtual void PrintReport (ostream & ost)
  {
    ost << "I am NumProcPyDemo" << endl;
  }

  virtual void Hello (string name)
  {
    cout << "Hi " << name << endl;
  }

  virtual double Sum (double a, double b)
  {
    return a+b;
  }
};


static RegisterNumProc<NumProcPyDemo> npinit1("demopy");





BOOST_PYTHON_MODULE(libmyngspy) {
  
  bp::def("Hello", &Hello);

  bp::class_<NumProcPyDemo,bp::bases<NumProc>> 
    ("NumProcPyDemo", bp::no_init)

    .def ("Hello", &NumProcPyDemo::Hello)
    .def ("Sum", &NumProcPyDemo::Sum)
    ;
}





struct Init {
  Init() 
  { 
    cout << endl << endl 
         << "************ please execute in py-shell: *************" 
         << endl << endl << endl

         << "import myngspy\n"
         << "np = pde.numprocs['np1']\n"
         << "print (np)\n"
         << "np.Hello ('Matthias')\n"
         << "np.Sum (3,8)"

         << endl << endl << endl;

    PyImport_AppendInittab("myngspy", PyInit_libmyngspy); 
  }
};
static Init init;





