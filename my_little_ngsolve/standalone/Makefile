NGSCXX = /opt/netgen/bin/ngscxx


all: demo_std demo_bla  demo_fem  demo_comp  demo_solve demo_parallel
# demo_comp1d

demo_std:  demo_std.cpp
	$(NGSCXX) demo_std.cpp -o demo_std -lngstd

demo_bla:  demo_bla.cpp
	$(NGSCXX) demo_bla.cpp -o demo_bla -lngstd -lngbla

demo_fem:  demo_fem.cpp
	$(NGSCXX) demo_fem.cpp -o demo_fem -lngfem -lngstd

demo_comp:  demo_comp.cpp
	$(NGSCXX) demo_comp.cpp -o demo_comp  -lngcomp -lngsolve -lngla -lngfem  -lngstd -lnglib -linterface

demo_comp1d:  demo_comp1d.cpp netgen_dummy.cpp
	$(NGSCXX) demo_comp1d.cpp netgen_dummy.cpp -o demo_comp1d  -lngcomp -lngsolve -lngla -lngfem  -lngstd -lnglib -linterface




demo_solve:  demo_solve.cpp
	$(NGSCXX) demo_solve.cpp -o demo_solve   -lngcomp -lngsolve -lngla -lngfem -lngstd -lnglib -linterface


demo_parallel:  demo_parallel.cpp
	$(NGSCXX) demo_parallel.cpp -o demo_parallel -lngstd



install:

clean:
	rm demo_std demo_bla demo_fem demo_comp demo_solve demo_parallel
