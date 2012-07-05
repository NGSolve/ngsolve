das braucht keiner mehr



#ifndef FILE_PARALLELFUNC
#define FILE_PARALLELFUNC


void ParallelRun();

void LoadPDEParallel ( const char*  filename );


#ifdef NGSOLVE
// void NGS_ParallelRun ( const string & message);
#endif

#endif
