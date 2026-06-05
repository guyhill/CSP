#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../../src/cspmain.h"

// Return codes.
enum {
    SUCCESS = 0,
    SUCCESS_SUBOPTIMAL,
    ERROR_UNDEFINED = 1000,
    ERROR_INPUT_FILE_NOT_FOUND,
    ERROR_INPUT_FILE_INVALID,
    ERROR_CSP_LOAD,
    ERROR_CSP_OPTIMIZE,
    ERROR_CSP_SOLUTION,
    ERROR_COST_INF,
    ERROR_OUTPUT_FILE
};

// ID's for JJ variables.
#define JJZERO          101
#define JJZERO1         102
#define JJZERO2         103
#define JJINF           104
#define JJMAXCOLSLP     105
#define JJMAXROWSLP     106
#define JJMAXCUTSPOOL   107
#define JJMAXCUTSITER   108
#define JJMINVIOLA      109
#define JJMAXSLACK      110
#define JJFEASTOL       111
#define JJOPTTOL        112
#define JJMAXTIME       113

// CPS problem struct.
typedef struct {

    int INumVar;
    char **names;
    double *data;
    int *weight;
    char *states;
    double *lb;
    double *ub;
    double *lpl;
    double *upl;

    int INumCons;
    int *ncard;
    double *rhs;
    int* list;
    signed char *val;

} csp_problem_t;

// Free memory for CSP problem struct.
void destroy_csp_problem(csp_problem_t *p) {
    for(int i = 0; i < p->INumVar; i++) {
        free(p->names[i]);
    }
    free(p->names);
    free(p->data);
    free(p->lpl);
    free(p->upl);
    free(p->lb);
    free(p->ub);
    free(p->weight);
    free(p->states);

    free(p->ncard);
    free(p->rhs);
    free(p->list);
    free(p->val);
}

// Read CSP problem.
int read_problem(const char *input_filename, csp_problem_t *p) {
    // Open input file.
    FILE *input_file = fopen(input_filename, "r");
    if (input_file == NULL) {
        return ERROR_INPUT_FILE_NOT_FOUND;
    }

    // Input file should start with 0.
    int dum;
    fscanf(input_file, "%d\n", &dum);
    if (dum != 0) {
        fclose(input_file);
        return ERROR_INPUT_FILE_INVALID;
    }

    // Read rest of input file. Note that valid contents is assumed.
    int variable_cnt;
    fscanf(input_file, "%d\n", &variable_cnt);
    p->INumVar = variable_cnt;
    p->names = (char **)malloc(p->INumVar * sizeof(char *));
    p->data = (double *)malloc(p->INumVar * sizeof(double));
    p->weight = (int *)malloc(p->INumVar * sizeof(int));
    p->states = (char *)malloc(p->INumVar * sizeof(char));
    p->lb = (double *)malloc(p->INumVar * sizeof(double));
    p->ub = (double *)malloc(p->INumVar * sizeof(double));
    p->lpl = (double *)malloc(p->INumVar * sizeof(double));
    p->upl = (double *)malloc(p->INumVar * sizeof(double));

    for (int i = 0; i < p->INumVar; i++) {
        double d_, lb_, ub_, lpl_, upl_, spl_;
        int w_, x_;
        char s_;
        fscanf(input_file, "%d %lf %d %c %lf %lf %lf %lf %lf\n", &x_, &d_, &w_, &s_, &lb_, &ub_, &lpl_, &upl_, &spl_);
        char name[16];  // Large enough to fit "x<int>\0" for any integer <int>.
        int name_len = sprintf(name, "x%d", x_);
        p->names[i] = (char *)malloc((name_len + 1) * sizeof(char));
        strcpy(p->names[i], name);
        p->data[i] = d_;
        p->weight[i] = w_;
        p->states[i] = s_;
        p->lb[i] = lb_;
        p->ub[i] = ub_;
        p->lpl[i] = lpl_;
        p->upl[i] = upl_;
    }

    int constraint_cnt;
    fscanf(input_file, "%d\n", &constraint_cnt);

    p->INumCons = constraint_cnt;
    p->ncard = (int *)malloc(p->INumCons * sizeof(int));
    p->rhs = (double *)malloc(p->INumCons * sizeof(double));
    p->list = NULL;
    p->val  = NULL;

    int t = 0;
    int listsize = 0;
    for (int i = 0; i < p->INumCons; i++) {
        double d_;
        int terms;
        fscanf(input_file,"%lf %d :", &d_, &terms);

        p->ncard[i] = terms;
        p->rhs[i] = d_;
        listsize += terms;
        p->list = (int *)realloc(p->list, listsize * sizeof(int));
        p->val = (signed char *)realloc(p->val, listsize * sizeof(signed char));

        for (int j = 0; j < terms; j++) {
            int Icell, v_;
            fscanf(input_file," %d (%d)",&Icell,&v_);
            p->list[t] = Icell;
            p->val[t] = v_;
            t++;
        }
    }

    // Close input file.
    fclose(input_file);

    return SUCCESS;
}

// Solve CSP problem.
int solve_problem(csp_problem_t *p, int *lcost, int *ucost, long max_time_in_sec, const char *output_path) {
    int status = ERROR_UNDEFINED;

    // Configure CSP.
    SCIPv::CSPSetDoubleConstant(JJMAXTIME, max_time_in_sec);
    SCIPv::CSPSetDoubleConstant(JJZERO, 1E-7);
    SCIPv::CSPSetDoubleConstant(JJINF, 2.14E9);
    SCIPv::CSPSetIntegerConstant(JJMAXCOLSLP, 50000);
    SCIPv::CSPSetIntegerConstant(JJMAXROWSLP, 15000);
    SCIPv::CSPSetIntegerConstant(JJMAXCUTSITER, 50);
    SCIPv::CSPSetIntegerConstant(JJMAXCUTSPOOL, 500000);
    SCIPv::CSPSetDoubleConstant(JJMINVIOLA, 0.0001);
    SCIPv::CSPSetDoubleConstant(JJMAXSLACK, 0.0001);
    SCIPv::CSPSetDoubleConstant(JJFEASTOL, 1E-6);
    SCIPv::CSPSetDoubleConstant(JJOPTTOL, 1E-9);
    SCIPv::CSPSetFileNames(output_path);

    // Solve problem using CSP.
    do {
        int load_status = SCIPv::CSPloadprob(
            p->INumCons, p->rhs, p->INumVar, p->data, p->weight, p->states,
            p->lpl, p->upl, p->lb, p->ub, p->names, p->ncard, p->list, p->val
        );
        if (load_status != SUCCESS) {
            status = load_status;
            break;
        }

        do {
            int opt_status = SCIPv::CSPoptimize(NULL);
            if (opt_status != SUCCESS) {
                status = opt_status;
                break;
            }

            int sol_status = SCIPv::CSPsolution(lcost, ucost, p->states);
            if (sol_status != SUCCESS) {
                status = sol_status;
                break;
            }

            if (*lcost >= SCIPv::CSPGetDoubleConstant(JJINF) - 0.1) {
                status = ERROR_COST_INF;
                break;
            }

            status = SUCCESS;
        } while (0);

        SCIPv::CSPfreeprob();

    } while (0);

    SCIPv::CSPFreeFileNames();

    return status;
}

int write_solution(const char *output_filename, const csp_problem_t *p) {
    // Open output file.
    FILE *output_file = fopen(output_filename, "w");
    if (output_file == NULL) {
        return ERROR_OUTPUT_FILE;
    }

    // Write result to output file.
    for (int i = 0; i < p->INumVar; i++)
        if (p->states[i] == 'm')
            fprintf(output_file, "%5d\n", i);

    // Close output file.
    fclose(output_file);

    return SUCCESS;
}

int main(void) {

    const char input_filename[] = "in.jj";
    const char output_filename[] = "out.jj";
    const char output_path[] = "./";
    long max_time_in_sec = 60;

    printf("Reading problem...\n");
    csp_problem_t csp_problem;
    int read_status = read_problem(input_filename, &csp_problem);
    if (read_status != SUCCESS) {
        printf("ERROR: %d\n", read_status);
        return read_status;
    }

    int status = ERROR_UNDEFINED;
    do {
        printf("Solving problem...\n");
        int lcost, ucost;
        int solve_status = solve_problem(&csp_problem, &lcost, &ucost, max_time_in_sec, output_path);
        if (solve_status != SUCCESS) {
            printf("ERROR: %d\n", solve_status);
            status = solve_status;
            break;
        }

        printf("Writing result...\n");
        int write_status = write_solution(output_filename, &csp_problem);
        if (write_status != SUCCESS) {
            printf("ERROR: %d\n", write_status);
            status = write_status;
            break;
        }

        printf("lcost = %d, ucost = %d\n", lcost, ucost);
        status = SUCCESS;
    } while (0);

    destroy_csp_problem(&csp_problem);

    return status;
}
