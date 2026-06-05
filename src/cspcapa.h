#include "cspdefns.h"

int separa_capacity(int *, CONSTRAINT **);
double violation_capacity(CONSTRAINT *);
double get_coeficient_capacity(VARIABLE *, CONSTRAINT *);
CONSTRAINT *capacity_constraint(int, VARIABLE **, double *, VARIABLE *, double);
