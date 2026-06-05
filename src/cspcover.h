#include "cspdefns.h"

int separa_cover(int *, CONSTRAINT **);
int extend_cover(CONSTRAINT *, VARIABLE **, double *);
double violation_cover(CONSTRAINT *);
double get_coeficient_cover(VARIABLE *, CONSTRAINT *);
