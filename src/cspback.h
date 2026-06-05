#include "cspexport.h"

int CSPstopcondition(void);
int CSPstoptime(void);
int CSPnewsolution(void);
int CSPexit(int);
int CSPdefinestop(int (*)(void));
int CSPdefineheur(int (*)(void));

EXPORTFUNC int CSPdefineexit(int (*)(int));
EXPORTFUNC int CSPdefinestoptime(int (*)(void));

#define EXIT_MEMO 1
#define EXIT_ERROR 2
#define EXIT_LPSOLVER 3
