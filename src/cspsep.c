/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPSEP.c                                         */
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
#include "cspcapa.h"
#include "cspbridg.h"
#include "cspcover.h"
#include "cspgomo.h"
#include "cspbranc.h"
#include "cspback.h"
#include "cspsep.h"

/*  PROTOTYPES OF FUNCTIONS */

static int separa_pool(int *, CONSTRAINT **);
static int iden_rows(CONSTRAINT *, CONSTRAINT *);
static int pool_row(CONSTRAINT *);
static long int hash(CONSTRAINT *);
static int enumeration(int, double, VARIABLE **, double *, int, double,
                       VARIABLE **, double *);
static int integrability(void);
static int tailingoff(void);

/*  FUNCTIONS DEFINITIONS */

void separa(int *card, CONSTRAINT **stack)

{
    if (tailingoff()) {
        return;
    }
    separa_pool(card, stack);

    if (*card) {
        return;
    }

    separa_capacity(card, stack);

    if (*card == 0 && integrability()) {
        return;
    }

    separa_bridge(card, stack);
    separa_cover(card, stack);
}

double violated(CONSTRAINT *con)

{
    switch (con->type) {
    case CAPACITY:
        return (violation_capacity(con));
    case BRIDGE:
        return (violation_bridge(con));
    case COVER:
        return (violation_cover(con));
    case BRANCHCUT:
        return (violation_branch(con));
    case GOMORY:
        return (violation_gomory(con, nsupport, support));
    default:
        std::cout << "ERROR: un-considered constraint " << con->type
                  << std::endl;
        CSPexit(EXIT_ERROR); // exit(1);
    }
    return 0;
}

/* --------------------------  POOL SEPARATION ------------------- */

static int separa_pool(int *card, CONSTRAINT **stack)

{
    int k, num;
    CONSTRAINT *con;

    if (*card == MAX_CUTS_ITER) {
        return (0);
    }
    num = 0;
    for (k = 0; k < nrows; k++) {
        con = rows[k];
        if (con->type != BRANCHCUT && con->stat < 0 &&
            violated(con) > MIN_VIOLA) {
            num++;
            stack[(*card)++] = con;
            con->stat = LP_BA;
            if (*card == MAX_CUTS_ITER) {
                break;
            }
        }
    }
    cpool += num;
    return (num);
}

static long int hash(CONSTRAINT *con)

{
    int k;
    long int val = 0;
    double rhs;

    rhs = con->rhs;
    if (fabs(rhs) < ZERO) {
        rhs = 1;
    }

    switch (con->type) {
    case CAPACITY:
    case BRIDGE:
    case BRANCHCUT:
    case GOMORY:
        val = 0;
        k = con->card;
        while (k--) {
            srand(con->stack[k]->index);
            val +=
                rand() *
                ((long unsigned int)(con->stack[k]->index) * (con->coef[k])) /
                rhs;
        }
        break;
    case COVER:
        val = (long int)con->rhs;
        k = con->con->card;
        while (k--) {
            srand(con->con->stack[k]->index);
            val += (long int)rand() * 2 * con->con->stack[k]->index;
        }
        k = con->card;
        while (k--) {
            srand(con->stack[k]->index);
            val -= (long int)rand() * 2 * con->stack[k]->index;
        }
        break;
    default:
        std::cout << "ERROR: un-considered constraint " << con->type
                  << std::endl;
        CSPexit(EXIT_ERROR); // exit(1);
    }
    return (val);
}

int new_row(CONSTRAINT *con)

{
    con->hash = hash(con);

    if (pool_row(con) < 0) {
        rows[nrows++] = con;
        if (nrows == MAX_CUTS_POOL) {
            std::cout << "ERROR: Too much global rows" << std::endl;
            CSPexit(EXIT_ERROR); // exit(1);
        }
        return (1);
    }

    return (0);
}

static int pool_row(CONSTRAINT *con)

{
    int i;

    i = nrows;
    while (i--) {
        if (iden_rows(con, rows[i])) {
            return (i);
        }
    }
    return (-1);
}

static int iden_rows(CONSTRAINT *con1, CONSTRAINT *con2)

{
    int card1, card2, identical;
    double *coef1, *coef2;
    VARIABLE **stack1, **stack2;

    // PWOF added initialisation
    stack1 = NULL;
    stack2 = NULL;
    coef1 = NULL;
    coef2 = NULL;
    card1 = 0;
    card2 = 0;

    if (con1->hash != con2->hash) {
        return (0);
    }

    switch (con1->type) {
    case CAPACITY:
    case BRIDGE:
    case BRANCHCUT:
    case GOMORY:
        stack1 = con1->stack;
        coef1 = con1->coef;
        card1 = con1->card;
        break;
    case COVER:
        stack1 = (VARIABLE **)malloc(ncols * sizeof(VARIABLE *)); /*PWOF*/
        coef1 = (double *)malloc(ncols * sizeof(double));
        card1 = extend_cover(con1, stack1, coef1);
        break;
    default:
        std::cout << "ERROR: unknown type " << con1->type << " of constraint"
                  << std::endl;
        CSPexit(EXIT_ERROR); // exit(1);
    }

    switch (con2->type) {
    case CAPACITY:
    case BRIDGE:
    case BRANCHCUT:
    case GOMORY:
        stack2 = con2->stack;
        coef2 = con2->coef;
        card2 = con2->card;
        break;
    case COVER:
        stack2 = (VARIABLE **)malloc(ncols * sizeof(VARIABLE *)); /*PWOF*/
        coef2 = (double *)malloc(ncols * sizeof(double));
        card2 = extend_cover(con2, stack2, coef2);
        break;
    default:
        std::cout << "ERROR: unknown type " << con2->type << " of constraint"
                  << std::endl;
        CSPexit(EXIT_ERROR); // exit(1);
    }

    identical = enumeration(card1, con1->rhs, stack1, coef1, card2, con2->rhs,
                            stack2, coef2);

    switch (con1->type) {
    case CAPACITY:
    case BRIDGE:
    case BRANCHCUT:
    case GOMORY:
        break;
    case COVER:
        free(stack1);
        stack1 = NULL; /*PWOF*/
        free(coef1);
        coef1 = NULL; /*PWOF*/
        break;
    default:
        std::cout << "ERROR: unknown type " << con1->type << " of constraint"
                  << std::endl;
        CSPexit(EXIT_ERROR); // exit(1);
    }

    switch (con2->type) {
    case CAPACITY:
    case BRIDGE:
    case BRANCHCUT:
    case GOMORY:
        break;
    case COVER:
        free(stack2);
        stack2 = NULL; /*PWOF*/
        free(coef2);
        coef2 = NULL; /*PWOF*/
        break;
    default:
        std::cout << "ERROR: unknown type " << con2->type << " of constraint"
                  << std::endl;
        CSPexit(EXIT_ERROR); // exit(1);
    }
    return (identical);
}

static int enumeration(int card1, double rhs1, VARIABLE **stack1, double *coef1,
                       int card2, double rhs2, VARIABLE **stack2, double *coef2)

{
    int i, j;
    VARIABLE *ptr;

    if (card1 != card2) {
        return (0);
    }
    if (fabs(rhs1) < ZERO) {
        rhs1 = 1;
    }
    if (fabs(rhs2) < ZERO) {
        rhs2 = 1;
    }
    for (i = 0; i < card1; i++) {
        ptr = stack1[i];
        for (j = 0; j < card2; j++) {
            if (ptr == stack2[j]) {
                if (fabs(coef1[i] * rhs2 - coef2[j] * rhs1) > ZERO) {
                    return (0);
                }
                break;
            }
        }
        if (j == card2) {
            return (0);
        }
    }
    return (1);
}

void remove_row(CONSTRAINT *con)

{
    if (con == NULL) {
        return;
    }
    switch (con->type) {
    case CAPACITY:
    case BRIDGE:
    case BRANCHCUT:
    case GOMORY:
        free(con->stack);
        con->stack = NULL; /*PWOF*/
        free(con->coef);
        con->coef = NULL; /*PWOF*/
        break;
    case COVER:
        free(con->stack);
        con->stack = NULL; /*PWOF*/
        break;

    default:
        std::cout << " not known constraint in UNLOAD LP " << std::endl;
        CSPexit(EXIT_ERROR); // exit(1);
    }

    free((void *)con);
    con = NULL; /*PWOF*/
}

static int integrability() {
    int k = nsupport;
    while (k--) {
        if (fabs(support[k]->val - 0.5) < 0.5 - ZERO) {
            return (0);
        }
    }
    return (1);
}

#define NTAIL 10
#define ETAIL 0.01

static int tailingoff() {
    static double lbound[NTAIL];
    int i;

    i = NTAIL;
    while (--i) {
        lbound[i] = lbound[i - 1];
    }
    lbound[0] = lowerb;
    ntail++;
    if (ntail > NTAIL && lbound[0] < lbound[NTAIL - 1] + ETAIL) {
        if (integrability()) {
            return (0);
        }
        for (i = 1; i < NTAIL; i++) {
            if (lbound[i - 1] > lbound[i] + ETAIL) {
                return (0);
            }
        }
        return (1);
    }
    return (0);
}
