void   load_network(double *,char);
void   unload_network(void);
void   bounding_1(void);
double protection_level(VARIABLE *,int,int*,VARIABLE**,double*,double*,char);
int    protected_flow(int,VARIABLE**);

/* for multi-capacity separation */

void   free_col(int,double*);
void   unfree_col(int,double*);
