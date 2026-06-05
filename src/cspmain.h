#include "IProgressListener.h"
#include "cspexport.h"

// Different implementation for different solvers
// Assigned to different namespaces to be able to choose solver at runtime
namespace SCIPv
{
EXPORTFUNC void CSPSetFileNames(const char*);    
EXPORTFUNC void CSPFreeFileNames();
EXPORTFUNC void CSPSetDoubleConstant(const int, double);
EXPORTFUNC double CSPGetDoubleConstant(const int);
EXPORTFUNC void CSPSetIntegerConstant(const int, int);
EXPORTFUNC int CSPGetIntegerConstant(const int);
EXPORTFUNC int CSPloadprob(int,double*,int,double*,int*,char*,double*,double*,double*,double*,char**,int*,int*,signed char*);
EXPORTFUNC int CSPoptimize(IProgressListener*);
EXPORTFUNC int CSPfreeprob();
EXPORTFUNC int CSPsolution(int*,int*,char*);
EXPORTFUNC int CSPrelbounds(int,int*,double*,double*,char);
}  //namespace SCIPv end

// Functions to be used inside own code
int CSPtestprob(int,double*,int,double*,int*,char*,double*,double*,double*,double*,char**,int*,int*,signed char*);
int CSPabsbounds(int,int*,double*,double*,char);
int CSPpartialbounds();
int CSPwrite(char*);
