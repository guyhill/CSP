/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPPREP.c                                        */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/

#include <stdlib.h>
#include <iostream>

#include "cspdefns.h"
#include "cspglob2.h"
#include "cspnet.h"
#include "cspsep.h"
#include "cspcapa.h"
#include "cspback.h"
#include "cspprep.h"

int l0u0, l1u0, l0u1, l1u1;

static void clean_prot_level(void);

/*  FUNCTIONS DEFINITIONS */

static int sort_pl(const void *p, const void *q /*PROT_LEVEL *p,PROT_LEVEL *q*/)

{
  if (((PROT_LEVEL *)p)->level < ((PROT_LEVEL *)q)->level) {
    return (-1);
  }
  if (((PROT_LEVEL *)p)->level > ((PROT_LEVEL *)q)->level) {
    return (+1);
  }
  return (0);
}

int preprocessing() {
    int k, l, nvar;
    double lpl, upl;
    PROT_LEVEL *pro, *pro0;
    double *primal, *cvar;
    VARIABLE **xvar;
    CONSTRAINT *con;

    clean_prot_level();
    if (nprot_level == 0) {
        return 0;
    }

    qsort((void *)prot_level, nprot_level, sizeof(PROT_LEVEL), sort_pl);

    cvar = (double *)malloc(Rncells * sizeof(double));
    xvar = (VARIABLE **)malloc(Rncells * sizeof(VARIABLE *));
    primal = (double *)calloc(sizeof(double), Rncells);
    if (!cvar || !xvar || !primal) {
        std::cout << "There is not enough memory in PREPROCESSING" << std::endl;
        CSPexit(EXIT_MEMO); // exit(1);
    }
    for (k = 0; k < nsensitive; k++) {
        primal[sensitive[k].var->index] = sensitive[k].var->val;
    }
    load_network(primal, 'C');

    l = nprot_level;
    while (l--) {
        pro = prot_level + l;
        if (pro->level > ZERO && pro->study) {
            if (protection_level(pro->sen->var, pro->sense, &nvar, xvar, cvar,
                                 primal, 'C') <
                pro->sense * pro->sen->var->nominal + pro->level - ZERO) {

                if (nvar == 0) {
                    break;
                }

                pro->study = 0;

                con = capacity_constraint(nvar, xvar, cvar, pro->sen->var,
                                          pro->level);

                if (!new_row(con)) {
                    remove_row(con);
                } else {
                    rows[nrows - 1]->stat = FIX_LB;
                    rows[nrows - 1]->lp = -1;
                }

            } else {
                pro->level = 0;
            }

            k = l;
            while (k--) {
                pro0 = prot_level + k;
                if (pro0->level > ZERO &&
                    pro0->sense * primal[pro0->sen->var->index] >
                        pro0->sense * pro0->sen->var->nominal + pro0->level -
                            ZERO) {
                    pro0->level = 0;
                }
            }
        }
    }

    free(primal);
    primal = NULL; /*PWOF*/
    free(cvar);
    cvar = NULL; /*PWOF*/
    free(xvar);
    xvar = NULL; /*PWOF*/
    unload_network();

    if (l != -1) {
        return (-l - 1); /* not feasible sol. exits */
    }

    clean_prot_level();
    if (nprot_level == 0) {
        return 0; /* auto-protection */
    }

    l1u1 = l1u0 = l0u1 = l0u0 = 0;
    for (k = 0; k < nsensitive; k++) {
        lpl = sensitive[k].lpl;
        upl = sensitive[k].upl;
        if (lpl > ZERO && upl > ZERO) {
            l1u1++;
        } else if (lpl > ZERO) {
            l1u0++;
        } else if (upl > ZERO) {
            l0u1++;
        } else {
            l0u0++;
        }
    }
    return (l1u1 + l1u0 + l0u1);
}

static void clean_prot_level() {
    int k, l;

    l = nprot_level;
    while (l--) {
        if (prot_level[l].level < ZERO) {
            if (prot_level[l].sense > 0) {
                prot_level[l].sen->upl = 0;
            } else {
                prot_level[l].sen->lpl = 0;
            }
            for (k = l + 1; k < nprot_level; k++) {
                prot_level[k - 1].sen = prot_level[k].sen;
                prot_level[k - 1].sense = prot_level[k].sense;
                prot_level[k - 1].level = prot_level[k].level;
                prot_level[k - 1].study = prot_level[k].study;
            }
            nprot_level--;
        }
    }
}
