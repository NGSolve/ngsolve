objects = mixed_diffusion.o hybridmixed.o

%.o : %.cpp
	gcc  -O2 -fopenmp -fpic -DNETGEN_ELTRANS -DUSE_TIMEOFDAY -DLAPACK -I. -I$(NETGENDIR)/../include -c $? -o $@

libmixed.so : $(objects)
	gcc -shared -fopenmp -L/opt/netgen/lib -fpic $(objects) -lngcomp -lngsolve -o $@

clean:
	rm *.o libmixed.so

