/*******************************************************************/
/*  Interface to LP-solvers:                                       */
/*       CPLEX 3.0, CPLEX 5.0, CPLEX 7.0  XPRESS 10.0  XPRESS 11.0 */
/*                                                                 */
/*  J.J.Salazar                       Rotterdam, August 10, 1997   */
/*                                                                 */
/*  Interface to LP-solver: SCIP 3.0.1    (VSCIP)                  */
/*  M-S Hernï¿½ndez-Garcï¿½a              Spain, 2013                  */
/*******************************************************************/

#include <iostream>
#include <scip/scip.h>

#include "jjsolver.h"

static unsigned int lpJJ = 0; /* for counting the loaded LP's   */

static bool dual = false;
static JJLPptr *Env = NULL;

JJLPptr JJloadprob(char *probname, int numcols, int numrows, int numrims,
                   int objsen, double *obj, double *rhs, char *sense,
                   int *matbeg, int *matcnt, int *matind, double *matval,
                   double *lb, double *ub, double *rngval, int *freerowind,
                   int *rimtype, int *rimbeg, int *rimcnt, int *rimind,
                   double *rimval, char *dataname, char *objname, char *rhsname,
                   char *rngname, char *bndname, char **colname,
                   char *colnamestore, char **rowname, char *rownamestore,
                   char **rimname, char *rimnamestore, int colspace,
                   int rowspace, int nzspace, int rimspace, int rimnzspace,
                   unsigned colnamespace, unsigned rownamespace,
                   unsigned rimnamespace) {
    if (Env == NULL) {
        Env = new JJLPptr[2];
        Env[0] = NULL;
        Env[1] = NULL;
        lpJJ = 0;
    }
    JJLPptr lp = NULL;

    SCIPlpiCreate(&lp, NULL, probname, (SCIP_OBJSEN)objsen);
    SCIPlpiSetIntpar(lp, SCIP_LPPAR_FASTMIP, TRUE); // Salomé added 31-01-2014

    Env[lpJJ] = lp;

    double *rowLower = new double[numrows];
    double *rowUpper = new double[numrows];
    for (int i = 0; i < numrows; i++) {
        if (sense[i] == 'L') {
            rowLower[i] = -SCIPlpiInfinity(lp);
            rowUpper[i] = rhs[i];
        } else if (sense[i] == 'E') {
            rowLower[i] = rhs[i];
            rowUpper[i] = rhs[i];
        } else if (sense[i] == 'G') {
            rowLower[i] = rhs[i];
            rowUpper[i] = SCIPlpiInfinity(lp);
        }
    }

    SCIPlpiLoadColLP(lp, (SCIP_OBJSEN)objsen, numcols, obj, lb, ub, colname,
                     numrows, rowLower, rowUpper, rowname,
                     matbeg[numcols - 1] + matcnt[numcols - 1], matbeg, matind,
                     matval);

    delete[] rowLower;
    delete[] rowUpper;

    if (lp) {
        lpJJ++;
    }
    return lp;
}

int JJcopyctype(JJLPptr lp, char *xctype) {
    std::cout << "\nJJcopyctype?" << std::endl;
    system("pause");
    return 1;
}

void JJfreeprob(JJLPptr *lp) {
    SCIPlpiFree(lp);
    lpJJ--;
    if (lpJJ == 0) {  // no lp left in Env[]
        free(*lp);    // PWOF added 27-08-2013
        *lp = NULL;   // already done by SCIPlpiFree(lp)???
        delete[] Env; // PWOF added 23-08-2013
        Env = NULL;
    } else {
        *lp = Env[lpJJ - 1];
    }
}

int JJoptimize(JJLPptr lp) {
    dual = false;
    if (SCIPlpiSolvePrimal(lp) == SCIP_OKAY) {
        return 0;
    }
    return 1;
}

int JJdualopt(JJLPptr lp) {
    dual = true;
    if (SCIPlpiSolveDual(lp) == SCIP_OKAY) {
        return 0;
    }
    return 1;
}

int JJhybnetopt(JJLPptr lp, char method) {
    std::cout << "\nJJhybnetopt?" << std::endl;
    system("pause");
    return 1;
}

int JJnetopt(JJLPptr lp, int *netstatus_p, int *numnodes_p, int *numarcs_p,
             int *itcnt_p) {
    *numnodes_p = 0;
    *numarcs_p = 0;
    *itcnt_p = 0;
    dual = true;
    if (SCIPlpiSolveDual(lp) == SCIP_OKAY) {
        *netstatus_p = 0;
    } else {
        *netstatus_p = 1;
    }
    return *netstatus_p;
}

int JJmipopt(JJLPptr lp) {
    dual = false;
    if (SCIPlpiSolvePrimal(lp) == SCIP_OKAY) {
        return 0;
    }
    return 1;
}

int JJsetscr_ind(JJLPptr lp, int scr_ind) {
    SCIPlpiSetIntpar(Env[lpJJ - 1], SCIP_LPPAR_LPINFO, scr_ind);
    return 0;
}

int JJgetmac(JJLPptr lp) {
    int ncols;
    SCIPlpiGetNCols(lp, &ncols);
    return ncols;
}

int JJgetmar(JJLPptr lp) {
    int nrows;
    SCIPlpiGetNRows(lp, &nrows);
    return nrows;
}

int JJgetmat(JJLPptr lp) {
    int num;
    SCIPlpiGetNNonz(lp, &num);
    return num;
}

int JJlpiterlimit(int val) {
    SCIPlpiSetIntpar(Env[lpJJ - 1], SCIP_LPPAR_LPITLIM, val);
    return 0;
}

int JJgetobjsen(JJLPptr lp) {
    std::cout << "\nJJgetobjsen?" << std::endl;
    system("pause");
    return 1;
}

int JJgetobj(JJLPptr lp, double *obj, int begin, int end) {
    SCIPlpiGetObj(lp, begin, end, obj);
    return 0;
}

int JJgetrhs(JJLPptr lp, double *rhs, int begin, int end) {
    SCIP_Real *lhss = new SCIP_Real[end - begin + 1];
    SCIP_Real *rhss = new SCIP_Real[end - begin + 1];
    SCIPlpiGetSides(lp, begin, end, lhss, rhss);
    for (int i = 0; i <= end - begin; i++) {
        if (SCIPlpiIsInfinity(lp, rhss[i])) {
            rhs[i] = lhss[i];
        } else {
            rhs[i] = rhss[i];
        }
    }
    delete[] lhss;
    delete[] rhss;
    return 0;
}

int JJgetsense(JJLPptr lp, char *sense, int begin, int end) {
    SCIP_Real *lhss = new SCIP_Real[end - begin + 1];
    SCIP_Real *rhss = new SCIP_Real[end - begin + 1];
    SCIPlpiGetSides(lp, begin, end, lhss, rhss);
    for (int i = 0; i <= end - begin; i++) {
        if (SCIPlpiIsInfinity(lp, -lhss[i])) {
            sense[i] = 'L';
        } else if (SCIPlpiIsInfinity(lp, rhss[i])) {
            sense[i] = 'G';
        } else {
            sense[i] = 'E';
        }
    }
    delete[] lhss;
    delete[] rhss;
    return 0;
}

int JJgetcols(JJLPptr lp, int *nzcnt, int *cmatbeg, int *cmatind,
              double *cmatval, int cmatspace, int *surplus, int begin,
              int end) {
    SCIPlpiGetCols(lp, begin, end, NULL, NULL, nzcnt, cmatbeg, cmatind,
                   cmatval);
    *surplus = cmatspace - *nzcnt;
    return 0;
}

int JJgetrows(JJLPptr lp, int *nzcnt, int *rmatbeg, int *rmatind,
              double *rmatval, int rmatspace, int *surplus, int begin,
              int end) {
    SCIPlpiGetRows(lp, begin, end, NULL, NULL, nzcnt, rmatbeg, rmatind,
                   rmatval);
    *surplus = rmatspace - *nzcnt;
    return 0;
}

int JJsolution(JJLPptr lp, int *lpstat_p, double *objval_p, double *x,
               double *pi, double *slack, double *dj) {
    *lpstat_p = JJgetstat(lp);
    SCIPlpiGetSol(lp, objval_p, x, pi, NULL /*slack*/, dj);
    JJgetslack(lp, slack, 0, JJgetmar(lp) - 1);
    return 0;
}

int JJgetstat(JJLPptr lp) {
    int m_stat = SCIPlpiGetInternalStatus(lp);
    switch (m_stat) // soplex
    {
    case 1:
        return 1; // optimal
    case 2:
        return 2; // unbounded
    case 3:
        return 3; // infeasible
    case -7:
        return 11; // timelimExc
    case -6:
        return 10; // iterlimExc
    case -5:
        return 12; // objlimExc
    default:
        // std::cout << "error " << m_stat << std::endl;
        return 1101;
    }
}

int JJgetobjval(JJLPptr lp, double *objval_p) {
    SCIPlpiGetObjval(lp, objval_p);
    return 0;
}

int JJgetx(JJLPptr lp, double *x, int begin, int end) {
    SCIP_Real *primsol = new SCIP_Real[JJgetmac(lp)];
    SCIPlpiGetSol(lp, NULL, primsol, NULL, NULL, NULL);

    for (int i = begin; i <= end; i++) {
        x[i - begin] = primsol[i];
    }
    delete[] primsol;
    return 0;
}

int JJgetpi(JJLPptr lp, double *pi, int begin, int end) {
    SCIP_Real *dualsol = new SCIP_Real[JJgetmar(lp)];
    SCIPlpiGetSol(lp, NULL, NULL, dualsol, NULL, NULL);

    for (int i = begin; i <= end; i++) {
        pi[i - begin] = dualsol[i];
    }
    delete[] dualsol;
    return 0;
}

int JJgetslack(JJLPptr lp, double *slack, int begin, int end) {
    SCIP_Real *rhs = new SCIP_Real[end - begin + 1];
    SCIP_Real *lhs = new SCIP_Real[end - begin + 1];
    SCIP_Real *activ = new SCIP_Real[JJgetmar(lp)];
    SCIPlpiGetSol(lp, NULL, NULL, NULL, activ, NULL);
    SCIPlpiGetSides(lp, begin, end, lhs, rhs);
    for (int i = begin; i <= end; i++) {
        if (SCIPlpiIsInfinity(lp, rhs[i - begin])) {
            slack[i - begin] = lhs[i - begin] - activ[i];
        } else {
            slack[i - begin] = rhs[i - begin] - activ[i];
        }
    }
    delete[] activ;
    delete[] rhs;
    delete[] lhs;
    return 0;
}

int JJgetdj(JJLPptr lp, double *dj, int begin, int end) {
    SCIP_Real *redcost = new SCIP_Real[JJgetmac(lp)];
    SCIPlpiGetSol(lp, NULL, NULL, NULL, NULL, redcost);

    for (int i = begin; i <= end; i++) {
        dj[i - begin] = redcost[i];
    }
    delete[] redcost;
    return 0;
}

int JJgetitc(JJLPptr lp) {
    int iterations;
    SCIPlpiGetIterations(lp, &iterations);
    return iterations;
}

int JJgetitci(JJLPptr lp) {
    int status = JJgetstat(lp);
    if ((status != 1) && (status != 2) && (status != 3)) {
        return 0;
    }
    int iterations;
    SCIPlpiGetIterations(lp, &iterations);
    return iterations;
}

int JJgetbase(JJLPptr lp, int *cstat, int *rstat) {
    SCIPlpiGetBase(lp, cstat, rstat);
    return 0;
}

int JJloadbase(JJLPptr lp, int *cstat, int *rstat) {
    SCIPlpiSetBase(lp, cstat, rstat);
    return 0;
}

int JJaddrows(JJLPptr lp, int ccnt, int rcnt, int nzcnt, double *rhs,
              char *sense, int *rmatbeg, int *rmatind, double *rmatval,
              char **colname, char **rowname) {
    if (ccnt != 0) // add cols
    {
        int *mybeg = new int[ccnt + 1];
        double *lb = new double[ccnt];
        double *ub = new double[ccnt];
        double *obj = new double[ccnt];
        for (int j = 0; j < ccnt; ++j) {
            mybeg[j] = 0;
            lb[j] = 0.0;
            ub[j] = SCIPlpiInfinity(lp);
            obj[j] = 0.0;
        }
        mybeg[ccnt] = 0;
        // add columns
        SCIPlpiAddCols(lp, ccnt, obj, lb, ub, colname, nzcnt, mybeg, 0, 0);

        delete[] mybeg;
        delete[] obj;
        delete[] lb;
        delete[] ub;
    }
    double *rowLower = new double[rcnt];
    double *rowUpper = new double[rcnt];
    for (int i = 0; i < rcnt; i++) {
        if (sense[i] == 'L') {
            rowLower[i] = -SCIPlpiInfinity(lp);
            rowUpper[i] = rhs[i];
        } else if (sense[i] == 'E') {
            rowLower[i] = rhs[i];
            rowUpper[i] = rhs[i];
        } else if (sense[i] == 'G') {
            rowLower[i] = rhs[i];
            rowUpper[i] = SCIPlpiInfinity(lp);
        }
    }

    SCIPlpiAddRows(lp, rcnt, rowLower, rowUpper, rowname, nzcnt, rmatbeg,
                   rmatind, rmatval);

    delete[] rowLower;
    delete[] rowUpper;
    return 0;
}

int JJdelrows(JJLPptr lp, int begin, int end) {
    SCIPlpiDelRows(lp, begin, end);
    return 0;
}

int JJdelsetrows(JJLPptr lp, int *delstat) {
    SCIPlpiDelRowset(lp, delstat);
    return 0;
}

int JJaddcols(JJLPptr lp, int ccnt, int nzcnt, double *obj, int *cmatbeg,
              int *cmatind, double *cmatval, double *lb, double *ub,
              char **colname) {
    SCIPlpiAddCols(lp, ccnt, obj, lb, ub, colname, nzcnt, cmatbeg, cmatind,
                   cmatval);
    return 0;
}

int JJdelcols(JJLPptr lp, int begin, int end) {
    SCIPlpiDelCols(lp, begin, end);
    return 0;
}

int JJdelsetcols(JJLPptr lp, int *delstat) {
    SCIPlpiDelColset(lp, delstat);
    return 0;
}

int JJchgcoef(JJLPptr lp, int i, int j, double newvalue) {
    int ind;

    if (i == -1) {
        ind = j;
        SCIP_Real opc;
        opc = newvalue;
        SCIPlpiChgObj(lp, 1, &ind, &opc);

    } else {
        if (j == -1) {
            ind = i;
            SCIP_Real lhs;
            SCIP_Real rhs;
            SCIP_Real rowlower;
            SCIP_Real rowupper;
            SCIPlpiGetSides(lp, i, i, &rowlower, &rowupper);
            if (SCIPlpiIsInfinity(lp, -rowlower)) {
                lhs = -SCIPlpiInfinity(lp);
                rhs = newvalue;
                SCIPlpiChgSides(lp, 1, &ind, &lhs, &rhs);
            } else if (SCIPlpiIsInfinity(lp, rowupper)) {
                lhs = newvalue;
                rhs = SCIPlpiInfinity(lp);
                SCIPlpiChgSides(lp, 1, &ind, &lhs, &rhs);
            } else { // if (sense[i] == 'E') {
                lhs = newvalue;
                rhs = newvalue;
                SCIPlpiChgSides(lp, 1, &ind, &lhs, &rhs);
            }
        } else {
            if (j == -2) {
                system("pause");
                return 1;
            } else {
                SCIPlpiChgCoef(lp, i, j, newvalue);
            }
        }
    }

    return 0;
}

int JJchgbds(JJLPptr lp, int cnt, int *index, char *lu, double *bd) {
    SCIP_Real lb;
    SCIP_Real ub;
    SCIP_Real lbAux;
    SCIP_Real ubAux;
    int indAux;
    for (int i = 0; i < cnt; i++) {
        SCIPlpiGetBounds(lp, index[i], index[i], &lbAux, &ubAux);
        indAux = index[i];
        if (lu[i] == 'L') {
            lb = bd[i];
            ub = ubAux;
            SCIPlpiChgBounds(lp, 1, &indAux, &lb, &ub);
        } else if (lu[i] == 'U') {
            lb = lbAux;
            ub = bd[i];
            SCIPlpiChgBounds(lp, 1, &indAux, &lb, &ub);
        } else {
            lb = bd[i];
            ub = bd[i];
            SCIPlpiChgBounds(lp, 1, &indAux, &lb, &ub);
        }
    }
    return 0;
}

int JJlpwrite(JJLPptr lp, std::string filename) {
    SCIPlpiWriteLP(lp, filename.c_str());
    return 0;
}

int JJmpswrite(JJLPptr lp, char *filename) {
    SCIPlpiWriteLP(lp, filename);
    return 0;
}

int JJgetbdl(JJLPptr lp, double *xlb, int begin, int end) {
    double *xub = new double[end - begin + 1];
    SCIPlpiGetBounds(lp, begin, end, xlb, xub);
    delete[] xub;
    return 0;
}

int JJgetbdu(JJLPptr lp, double *xub, int begin, int end) {
    double *xlb = new double[end - begin + 1];
    SCIPlpiGetBounds(lp, begin, end, xlb, xub);
    delete[] xlb;
    return 0;
}

int JJloadctype(JJLPptr lp, char *ctype) {
    std::cout << "\nJJloadctype?" << std::endl;
    system("pause");
    return 1;
}

int JJgetmx(JJLPptr lp, double *x, int begin, int end) {
    JJgetx(lp, x, begin, end);
    return 0;
}

int JJmipoptimize(JJLPptr lp) {
    JJoptimize(lp);
    return 0;
}

int JJbinvrow(JJLPptr lp, int i, double *y) {
    SCIPlpiGetBInvRow(lp, i, y, NULL, NULL);
    return 0;
}

int JJbinvarow(JJLPptr lp, int i, double *z) {
    double *binvrow = new SCIP_Real[JJgetmar(lp)];
    SCIPlpiGetBInvRow(lp, i, binvrow, NULL, NULL);
    SCIPlpiGetBInvARow(lp, i, binvrow, z, NULL, NULL);
    delete[] binvrow;
    return 0;
}

int JJgetbhead(JJLPptr lp, int *head, double *x) {
    std::cout << "\nJJgetbhead?" << std::endl;
    system("pause");
    return 1;
}

double JJinfinity(JJLPptr lp) {
  return SCIPlpiInfinity(lp);
}
