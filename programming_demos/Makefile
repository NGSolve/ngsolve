all: demo_std demo_bla  demo_fem  demo_comp demo_comp1d demo_solve

demo_std:  demo_std.cpp
	ngscxx demo_std.cpp -o demo_std -lngstd

demo_bla:  demo_bla.cpp
	ngscxx demo_bla.cpp -o demo_bla -lngstd -lngbla

demo_fem:  demo_fem.cpp
	ngscxx demo_fem.cpp -o demo_fem -lngfem -lngstd

demo_comp:  demo_comp.cpp
	ngscxx demo_comp.cpp -o demo_comp  -lngcomp -lngsolve -lngla -lngfem -lparallel -lngstd -lnglib -linterface

demo_comp1d:  demo_comp1d.cpp netgen_dummy.cpp
	ngscxx demo_comp1d.cpp netgen_dummy.cpp -o demo_comp1d  -lngcomp -lngsolve -lngla -lngfem -lparallel -lngstd -lnglib -linterface




demo_solve:  demo_solve.cpp
	ngscxx demo_solve.cpp -o demo_solve   -lngcomp -lngsolve -lngla -lngfem -lparallel -lngstd -lnglib -linterface



install:

clean:
	rm demo_std demo_bla demo_fem demo_comp demo_solve