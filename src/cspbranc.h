int    load_branch_tree(void);
int    unload_branch_tree(void);
int    read_prob(void);
int    var_branching(void);
int    con_branching(void);
double get_coeficient_branching(VARIABLE *,CONSTRAINT *);
double violation_branch(CONSTRAINT *);

extern struct BRANCH *tree;
