
notwendig ?


//      LAPACK++ (V. 1.1)
//      (C) 1992-1996 All Rights Reserved.


// Linkage names between C, C++, and Fortran (platform dependent)

#ifndef _ARCH_H_
#define _ARCH_H_

// #include <generic.h> 
#define F77NAME(x) name2(x,_)
#endif

#if defined(SGI) && !defined(SGI_DEC)
#define SGI_DEC

extern "C" {
	void mkidxname() {}
	void mkdatname() {}
}
#endif

