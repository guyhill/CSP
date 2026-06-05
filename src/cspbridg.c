/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPBRIDG.c                                       */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/

#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "cspdefns.h"
#include "cspglob2.h"
#include "cspnet.h"
#include "cspsep.h"
#include "cspback.h"
#include "cspbridg.h"

static CONSTRAINT *bridge_constraint(int, VARIABLE **, int);
static double bridge_lhs(int, int, VARIABLE **, VARIABLE **, double *, double *,
                         int *, VARIABLE **, double *);

/* -----------------------  BRIDGE SEPARATION --------------------------*/

int separa_bridge(int *card, CONSTRAINT **stack)

{

    int k, l, nvar, nvar1, nvar2, num;
    double *primal1, *primal2;
    double *cvar1, *cvar2, *cvar;
    VARIABLE *var, *var0;
    VARIABLE **xvar1, **xvar2, **xvar;
    CONSTRAINT *con;
    double minvalue, maxvalue, value;
    char *sleep;

    if (*card == MAX_CUTS_ITER) {
        return (0);
    }

    cvar = (double *)malloc(Rncells * sizeof(double));
    cvar1 = (double *)malloc(Rncells * sizeof(double));
    cvar2 = (double *)malloc(Rncells * sizeof(double));
    xvar1 = (VARIABLE **)malloc(Rncells * sizeof(VARIABLE *));
    xvar2 = (VARIABLE **)malloc(Rncells * sizeof(VARIABLE *));
    xvar = (VARIABLE **)malloc(Rncells * sizeof(VARIABLE *));
    primal1 = (double *)malloc(Rncells * sizeof(double));
    primal2 = (double *)malloc(Rncells * sizeof(double));
    if (!cvar1 || !cvar2 || !xvar || !primal1 || !primal2) {
        std::cout << "There is not enough memory for CAPACITY" << std::endl;
        CSPexit(EXIT_MEMO); // exit(1);
    }
    for (k = 0; k < Rncells; k++) {
        primal2[k] = primal1[k] = columns[k].val;
    }
    load_network(primal1, 'C');

    num = *card;
    sleep = (char *)calloc(nsupport, sizeof(char));
    for (k = 0; k < nsupport; k++) {
        if (support[k]->sensitive) {
            sleep[k] = 1;
        }
    }

    l = nsupport;
    while (l--) {
        var = support[l];
        if (!sleep[l]) {
            bounding_1();
            maxvalue =
                protection_level(var, 1, &nvar1, xvar1, cvar1, primal1, 'C');
            minvalue =
                -protection_level(var, -1, &nvar2, xvar2, cvar2, primal2, 'C');
            value = bridge_lhs(nvar1, nvar2, xvar1, xvar2, cvar1, cvar2, &nvar,
                               xvar, cvar);
            if (maxvalue - minvalue < ZERO ||
                value < var->val + ZERO - MIN_VIOLA) {
                con = bridge_constraint(nvar, xvar, var->index);
                if (con && new_row(con)) {
                    stack[num++] = con;
                    if (num == MAX_CUTS_ITER) {
                        break;
                    }
                } else {
                    remove_row(con);
                }
            }
            k = l;
            while (k--) {
                var0 = support[k];
                if (fabs(primal2[var0->index] - primal1[var0->index]) >
                    var0->val - MIN_VIOLA + ZERO) {
                    sleep[k] = 1;
                }
            }
        }
    }
    free(sleep);
    sleep = NULL; /*PWOF*/
    num -= *card;
    *card += num;
    unload_network();

    free(primal1);
    primal1 = NULL; /*PWOF*/
    free(primal2);
    primal2 = NULL; /*PWOF*/
    free(xvar);
    xvar = NULL; /*PWOF*/
    free(xvar1);
    xvar1 = NULL; /*PWOF*/
    free(xvar2);
    xvar2 = NULL; /*PWOF*/
    free(cvar);
    cvar = NULL; /*PWOF*/
    free(cvar1);
    cvar1 = NULL; /*PWOF*/
    free(cvar2);
    cvar2 = NULL; /*PWOF*/
    cbrid += num;
    return (num);
}

static CONSTRAINT *bridge_constraint(int nvar, VARIABLE **xvar, int index)

{
    int k;
    CONSTRAINT *row;
    VARIABLE **x;
    double *c;

    row = (CONSTRAINT *)malloc(sizeof(CONSTRAINT));
    if (row == NULL) {
        std::cout << "There is not enough memory for ROW" << std::endl;
        CSPexit(EXIT_MEMO); // exit(1);
    }
    x = (VARIABLE **)malloc((nvar + 1) * sizeof(VARIABLE *));
    if (x == NULL) {
        std::cout << " ERROR: not enough memory for the X (bridge.c)"
                  << std::endl;
        CSPexit(EXIT_MEMO); // exit(1);
    }
    c = (double *)malloc((nvar + 1) * sizeof(double));
    if (c == NULL) {
        std::cout << " ERROR: not enough memory for the C" << std::endl;
        CSPexit(EXIT_MEMO); // exit(1);
    }
    for (k = 0; k < nvar; k++) {
        x[k] = xvar[k];
        c[k] = 1.0;
    }

    x[nvar] = columns + index;
    c[nvar] = -1.0;

    row->rhs = 0;
    row->index = index;
    row->card = nvar + 1;
    row->stack = x;
    row->coef = c;
    row->sense = 'G';
    row->type = BRIDGE;
    row->stat = LP_BA;
    row->lp = -1;
    row->con = NULL;
    return (row);
}

double violation_bridge(CONSTRAINT *con)

{
    int k;
    double viola;
    VARIABLE **xvar;
    double *cvar;

    if (con->type != BRIDGE) {
        std::cout << " ERROR: not bridge constraint" << std::endl;
        CSPexit(EXIT_ERROR); // exit(1);
    }
    xvar = con->stack;
    cvar = con->coef;
    k = con->card;
    viola = 0;
    while (k--) {
        viola -= cvar[k] * xvar[k]->val;
    }
    return (viola);
}

double get_coeficient_bridge(VARIABLE *col, CONSTRAINT *con)

{
    int k;
    VARIABLE **xvar;
    double *cvar;

    xvar = con->stack;
    cvar = con->coef;
    k = con->card;
    while (k--) {
        if (xvar[k] == col) {
            return (cvar[k]);
        }
    }
    return 0;
}

static double bridge_lhs(int nvar1, int nvar2, VARIABLE **xvar1,
                         VARIABLE **xvar2, double *cvar1, double *cvar2,
                         int *nvar, VARIABLE **xvar, double *cvar)

{
    int i, k;
    double value;

    for (i = 0; i < Rncells; i++) {
        cvar[i] = 0.0;
    }
    for (i = 0; i < nvar1; i++) {
        cvar[xvar1[i]->index] += cvar1[i];
    }
    for (i = 0; i < nvar2; i++) {
        cvar[xvar2[i]->index] += cvar2[i];
    }
    for (k = i = 0, value = 0; i < Rncells; i++) {
        if (cvar[i] > ZERO) {
            value += columns[i].val;
            xvar[k] = columns + i;
            k++;
        }
    }
    *nvar = k;
    return value;
}
