/*******************************************************************/
/*  Interface to LP-solvers:                                       */
/*       CPLEX 3.0, CPLEX 5.0, CPLEX 7.0  XPRESS 10.0  XPRESS 11.0 */
/*                                                                 */
/*  J.J.Salazar                       Rotterdam, August 10, 1997   */
/*                                                                 */
/*  Interface to LP-solver: SCIP 3.0.1    (VSCIP)                  */
/*  M-S Hernï¿½ndez-Garcï¿½a              Spain, 2013                  */
/*******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include "jjsolver.h"
extern char* fproblemlp;

extern double FEAS_TOL;
extern double OPT_TOL;

static unsigned int lpJJ = 0;           /* for counting the loaded LP's   */

#ifdef VSCIP
static bool dual = false;
static JJLPptr *Env = NULL;
#endif

/*#ifdef CPLEX3
#define  SYSWAT386
#include "c:\cplex\cpxdefs.inc"
typedef CPXLPptr JJLPptr;
#endif

#ifdef CPLEX7
#include <ilcplex/cplex.h>
#define CPLEX5
#endif

#ifdef CPLEX5
#define  SYSWAT386
#ifndef CPLEX7
#include "c:\cplex50\cplex.h"
#endif
typedef CPXLPptr JJLPptr;
static CPXENVptr Env = NULL;         // CPLEX enviroment
//static unsigned int lpJJ;            // for counting the loaded LP's
#endif


#ifdef XPRESS11
#define XPRESS
#endif


#ifdef XPRESS
//#define DLL

#ifdef XPRESS11
#include "xpresso.h"
//#include "c:\xpress11\xpresso.h"        Attentive with output(,)
//#define XOSLDIR "c:\\xpress11"
   #define XOSLDIR NULL            //set XPRESSMP in  "autoexec.bat"
   #define XOSLMEM 0
#else
   #include "c:\xpress10\w32_osl\xpresso.h"
   #define XOSLDIR "c:\\xpress10\\w32_osl"
   #define XOSLMEM 16000000
#endif
struct JJLP {
        int ref;
        int objsen;
};
typedef struct JJLP *JJLPptr;
static JJLPptr current_lp = NULL;
//static unsigned int lpJJ;            // for counting the loaded LP's
static void JJcheck(JJLPptr);

static void JJcheck(lp)
JJLPptr lp;
{
    if(lp==NULL){
        fprintf(stderr,"JJ ERROR: not lp\n");
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }                                         
    if(lp != current_lp){ 
        if( current_lp ){
            if( savmat( &(current_lp->ref) ) ){
                fprintf(stderr,"JJ ERROR: not savmat()\n");
                CSPexit(EXIT_LPSOLVER); //exit(1);
            }
        }                                         
        current_lp = lp;                         
        if( resmat( current_lp->ref ) ){ 
            fprintf(stderr,"JJ ERROR: not resmat()\n");
            CSPexit(EXIT_LPSOLVER); //exit(1);
        }
    }
}
#endif */


/*****
void JJfreedata ( probname, obj, rhs, sense, matbeg, matcnt,
                    matind, matval, lb, ub, rngval,
                    freerowind, rimtype, rimbeg, rimcnt, rimind, rimval,
                    dataname, objname, rhsname, rngname, bndname,
                    colname, colnamestore, rowname, rownamestore,
                    rimname, rimnamestore)
char *probname;
double *obj, *rhs;
char *sense;
int *matbeg, *matcnt, *matind;
double *matval, *lb, *ub, *rngval;
int *freerowind, *rimtype, *rimbeg, *rimcnt, *rimind;
double *rimval;
char *dataname, *objname, *rhsname, *rngname, *bndname;
char **colname;
char *colnamestore;
char **rowname;
char *rownamestore;
char **rimname;
char *rimnamestore;
{
#ifdef CPLEX3
	if(probname)   free(probname);
	if(obj)        free(obj);
	if(rhs)        free(rhs);
	if(sense)      free(sense);
	if(matbeg)     free(matbeg);
	if(matcnt)     free(matcnt);
	if(matind)     free(matind);
	if(matval)     free(matval);
	if(lb)         free(lb);
	if(ub)         free(ub);
	if(rngval)     free(rngval);
	if(freerowind) free(freerowind);
	if(rimtype)    free(rimtype);
	if(rimbeg)     free(rimbeg);
	if(rimcnt)     free(rimcnt);
	if(rimind)     free(rimind);
	if(rimval)     free(rimval);
	if(dataname)   free(dataname);
	if(objname)    free(objname);
	if(rhsname)    free(rhsname);
	if(rngname)    free(rngname);
	if(bndname)    free(bndname);
	if(colname)    free(colname);
	if(colnamestore) free(colnamestore);
	if(rowname)    free(rowname);
	if(rownamestore) free(rownamestore);
	if(rimname)    free(rimname);
	if(rimnamestore) free(rimnamestore);
#endif
}
*****/



JJLPptr JJloadprob (
char *probname,
int numcols, int numrows, int numrims, int objsen,
double *obj, double *rhs,
char *sense,
int *matbeg, int *matcnt, int *matind,
double *matval, double *lb, double *ub, double *rngval,
int *freerowind, int *rimtype, int *rimbeg, int *rimcnt, int *rimind,
double *rimval,
char *dataname, char *objname, char *rhsname, char *rngname, char *bndname,
char **colname,
char *colnamestore,
char **rowname,
char *rownamestore,
char **rimname,
char *rimnamestore,
int colspace, int rowspace, int nzspace,
int rimspace, int rimnzspace,
unsigned colnamespace, unsigned rownamespace, unsigned rimnamespace
)
{
#ifdef CPLEX3
    return loadprob (probname, numcols, numrows, numrims,
                    objsen, obj, rhs, sense, matbeg, matcnt,
                    matind, matval, lb, ub, rngval,
                    freerowind, rimtype, rimbeg, rimcnt, rimind, rimval,
                    dataname, objname, rhsname, rngname, bndname,
                    colname, colnamestore, rowname, rownamestore,
                    rimname, rimnamestore, colspace, rowspace, nzspace,
                    rimspace, rimnzspace, colnamespace, rownamespace,
                    rimnamespace);
#endif


#ifdef CPLEX5    
   int  status;
   char errmsg[1024];
   JJLPptr lp;

   if( CPLEXv::Env==NULL ){
       CPLEXv::Env = CPXopenCPLEX (&status);
       if ( CPLEXv::Env == NULL ) {
           CPXgeterrorstring (CPLEXv::Env, status, errmsg);
           fprintf (stderr, "JJ ERROR: not CPLEX environment.\n%s",errmsg);
           return(NULL);
       }
       lpJJ = 0;
   }
/****
   lp = CPXloadprob (Env, probname, numcols, numrows, numrims,
                    objsen, obj, rhs, sense, matbeg, matcnt,
                    matind, matval, lb, ub, rngval,
                    freerowind, rimtype, rimbeg, rimcnt, rimind, rimval,
                    dataname, objname, rhsname, rngname, bndname,
                    colname, colnamestore, rowname, rownamestore,
                    rimname, rimnamestore, colspace, rowspace, nzspace,
                    rimspace, rimnzspace, colnamespace, rownamespace,
                    rimnamespace);
*****/  // modified on December 2006 to link with CPLEX 10

    lp = CPXcreateprob (CPLEXv::Env, &status, probname);
    status = CPXcopylp (CPLEXv::Env, lp, numcols, numrows, objsen, obj, rhs,
                     sense, matbeg, matcnt, matind, matval, lb, ub, rngval);


/*****
    if(probname)   free(probname);
    if(obj)        free(obj);
	if(rhs)        free(rhs);
	if(sense)      free(sense);
	if(matbeg)     free(matbeg);
	if(matcnt)     free(matcnt);
	if(matind)     free(matind);
	if(matval)     free(matval);
	if(lb)         free(lb);
	if(ub)         free(ub);
	if(rngval)     free(rngval);
	if(freerowind) free(freerowind);
	if(rimtype)    free(rimtype);
	if(rimbeg)     free(rimbeg);
	if(rimcnt)     free(rimcnt);
	if(rimind)     free(rimind);
	if(rimval)     free(rimval);
	if(dataname)   free(dataname);
	if(objname)    free(objname);
	if(rhsname)    free(rhsname);
	if(rngname)    free(rngname);
	if(bndname)    free(bndname);
	if(colname)    free(colname);
	if(colnamestore) free(colnamestore);
	if(rowname)    free(rowname);
	if(rownamestore) free(rownamestore);
    if(rimname)    free(rimname);
    if(rimnamestore) free(rimnamestore);
*****/
   if( lp ) lpJJ++;
   return lp;
#endif


#ifdef XPRESS
    int k,total;
    if( lpJJ == 0 ){
        if( initlz(XOSLDIR,XOSLMEM) ){
            fprintf (stderr, "JJ ERROR: XPRESS not opened." << std::endl;
            return NULL;
        }
    }
    if( current_lp != NULL ){
        if( savmat( &(current_lp->ref) ) ){
            std::cout << "JJ WARNING: not savmat()" << std::endl;
            return NULL;
        }
    }
    current_lp = (JJLPptr) malloc(sizeof(struct JJLP));
    if( current_lp == NULL ){
        std::cout << "JJ WARNING: not memory" << std::endl;
        return NULL;
    }
    current_lp->ref    = 0;
    current_lp->objsen = objsen;

    seticv(N_NCXTRA,colspace-numcols);
    seticv(N_NRXTRA,rowspace-numrows);
    for(total=k=0;k<numcols;k++) total += matcnt[k];
    seticv(N_NMXTRA,nzspace-total);
    
    if( loadprob(probname,numcols,numrows,sense,rhs,NULL,
                  obj,matbeg,matcnt,matind,matval,lb,ub) ){
        free( current_lp );
        current_lp = NULL; /*PWOF*/
        return NULL;
    }
    lpJJ++;
/*****
    if(probname)   free(probname);
    if(obj)        free(obj);
	if(rhs)        free(rhs);
	if(sense)      free(sense);
	if(matbeg)     free(matbeg);
	if(matcnt)     free(matcnt);
	if(matind)     free(matind);
	if(matval)     free(matval);
	if(lb)         free(lb);
	if(ub)         free(ub);
	if(rngval)     free(rngval);
	if(freerowind) free(freerowind);
	if(rimtype)    free(rimtype);
	if(rimbeg)     free(rimbeg);
	if(rimcnt)     free(rimcnt);
	if(rimind)     free(rimind);
	if(rimval)     free(rimval);
	if(dataname)   free(dataname);
	if(objname)    free(objname);
	if(rhsname)    free(rhsname);
	if(rngname)    free(rngname);
	if(bndname)    free(bndname);
	if(colname)    free(colname);
	if(colnamestore) free(colnamestore);
	if(rowname)    free(rowname);
	if(rownamestore) free(rownamestore);
    if(rimname)    free(rimname);
    if(rimnamestore) free(rimnamestore);
*****/
    return current_lp;
#endif 
    
// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
    int k,total,status;
    JJLPptr lp=NULL;
    
    lp = (JJLPptr)malloc(sizeof(JJLP));
    if (lp == NULL)
    {
        std::cout << "JJ Warning: not memory" << std::endl;
        return NULL;
    }
    lp->prob = NULL;
    lp->objsen = objsen;
    
    if (XPRScreateprob(&(lp->prob)))
    {
        lp->prob = NULL;
        std::cout << "JJ ERROR: problem with XPRScreateprob" << std::endl;
        free(lp);
        lp = NULL;
        return NULL;
    }
    
    XPRSsetintcontrol(lp->prob,XPRS_EXTRACOLS,colspace-numcols);
    XPRSsetintcontrol(lp->prob,XPRS_EXTRAROWS,rowspace-numrows);
    for (total=k=0;k<numcols;k++) total += matcnt[k];
    XPRSsetintcontrol(lp->prob,XPRS_EXTRAELEMS,nzspace-total);
    // Needed ????
    //XPRSsetintcontrol(lp->prob,XPRS_SOLUTIONFILE,0);
    
    XPRSsetdblcontrol(lp->prob,XPRS_FEASTOL, FEAS_TOL);
    XPRSsetdblcontrol(lp->prob,XPRS_OPTIMALITYTOL, OPT_TOL);
    
    XPRSsetdblcontrol(lp->prob,XPRS_MARKOWITZTOL, 1.0E-4);
    XPRSsetdblcontrol(lp->prob,XPRS_PERTURB, 1.0E-8);
    XPRSsetdblcontrol(lp->prob,XPRS_MATRIXTOL, 1.0E-10);
    XPRSsetdblcontrol(lp->prob,XPRS_PIVOTTOL, 1.0E-10);
    
    XPRSsetintcontrol(lp->prob,XPRS_PRESOLVE,0);
    XPRSsetintcontrol(lp->prob,XPRS_SCALING,0);
    
    status = XPRSloadlp(lp->prob,probname,numcols,numrows,sense,rhs,NULL,obj,matbeg,matcnt,matind,matval,lb,ub);
    
    if (status)
    {
        std::cout << "XPRESS ERROR " << status << std::endl;
        XPRSdestroyprob(lp->prob);
        lp->prob = NULL;
        free(lp);
        lp=NULL;
        return NULL;
    }
    lpJJ++;
    JJlpwrite(lp,fproblemlp);
    return lp;
#endif
// End adding PWOF    

#ifdef VSCIP
    if (Env == NULL) {
        Env = new JJLPptr[2];
        Env[0] = NULL;
        Env[1] = NULL;
        lpJJ = 0;
    }
    //SCIP_LPI *lp = NULL;
    JJLPptr lp=NULL;
    
    SCIPlpiCreate(&lp,NULL,probname,(SCIP_OBJSEN)objsen);
    SCIPlpiSetIntpar(lp, SCIP_LPPAR_FASTMIP, TRUE);     // Salomé added 31-01-2014

    Env[lpJJ] = lp;

    double *rowLower = new double[numrows];
    double *rowUpper = new double[numrows];
    for (int i = 0; i <  numrows; i++)
    {
         if (sense[i] == 'L') {
             rowLower[i] = -SCIPlpiInfinity(lp);
             rowUpper[i] = rhs[i];
         }else
             if (sense[i] == 'E') {
                 rowLower[i] = rhs[i];
                 rowUpper[i] = rhs[i];
             }else
                 if (sense[i] == 'G') {
                     rowLower[i] = rhs[i];
                     rowUpper[i] = SCIPlpiInfinity(lp);
                 }
    }

    SCIPlpiLoadColLP(lp,(SCIP_OBJSEN)objsen,numcols,obj,lb,ub,colname,numrows,rowLower,rowUpper,rowname,matbeg[numcols-1]+matcnt[numcols-1],matbeg,matind,matval);

    delete[] rowLower;
    delete[] rowUpper;
    //SCIPlpiSetRealpar(lp,SCIP_LPPAR_FEASTOL,1.0E-5);   // feasibility tolerence for primal variables and slacks
    //SCIPlpiSetRealpar(lp,SCIP_LPPAR_DUALFEASTOL,1.0E-5); // feasibility tolerance for dual variables and reduced costs

    if( lp ) lpJJ++;
    return lp;
#endif
}


int JJcopyctype(
JJLPptr lp,
char    *xctype
)
{
#ifdef CPLEX5    
    return CPXcopyctype(CPLEXv::Env, lp, xctype);
#endif 
#ifdef XPRESS_13
    std::cout << "\nJJcopyctype?" << std::endl;
    system("pause");
    return 1;
#endif
#ifdef VSCIP
    std::cout << "\nJJcopyctype?" << std::endl;
    system("pause");
    return 1;
#endif
}


void JJfreeprob (
JJLPptr *lp)
{
#ifdef CPLEX3
    freeprob (lp);
#endif

    
#ifdef CPLEX5    
   char  errmsg[1024];
   int   status;

   if ( lp != NULL ) {
//    status = CPXunloadprob ( Env , lp);
//    Modified on December 2006 to link with CPLEX 10
      status = CPXfreeprob ( CPLEXv::Env , lp);
      if ( status ) {
        fprintf (stderr, "CPXunloadprob failed, error code %d.\n", status);
        return;
      }
      lp = NULL; /*PWOF*/
      lpJJ--;
      if ( lpJJ==0 ) {
        status = CPXcloseCPLEX (&CPLEXv::Env);
        if ( status ) {
          CPXgeterrorstring ( CPLEXv::Env , status, errmsg);
          fprintf (stderr, "Could not close CPLEX environment.\n%s",errmsg);
          return;
        }
        CPLEXv::Env = NULL;
      }
   }
#endif


#ifdef XPRESS
    if(*lp != current_lp){
        if( delmat( (*lp)->ref ) ){
            fprintf(stderr,"JJ ERROR: not savmat()\n");
            CSPexit(EXIT_LPSOLVER); //exit(1);
        }
        free( *lp );
        *lp = NULL;
    } else {
        free(current_lp);
        current_lp = NULL;
    }
    lpJJ--;
    if( lpJJ==0 ) {
        freexo();
    }
#endif 
    
// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
    int status;
    status = XPRSdestroyprob((*lp)->prob);
    if (status)
    {
        std::cout << "XPRSdestroyprob failed, error code " << status << std::endl;
        return;
    }
    free(*lp);
    *lp = NULL;
    lpJJ--;
#endif
// End adding PWOF

#ifdef VSCIP
    //std::cout << "\nJJfreeprob" << std::endl;
    SCIPlpiFree(lp);
    lpJJ--;
    if ( lpJJ == 0) { // no lp left in Env[]
        free(*lp); // PWOF added 27-08-2013
        *lp = NULL; // already done by SCIPlpiFree(lp)???
        delete[] Env;//PWOF added 23-08-2013
        Env = NULL;
    }else
        *lp = Env[lpJJ-1];
#endif
}

int JJoptimize (
JJLPptr lp)
{
#ifdef CPLEX3
   return optimize (lp);
#endif


#ifdef CPLEX5      
   return CPXoptimize (CPLEXv::Env, lp);
#endif


#ifdef XPRESS
   JJcheck(lp);
   seticv(N_IALG,3);
   if( lp->objsen == 1 ) return minim("l");
   return maxim("l");
#endif 
   
// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   if (lp->objsen == 1) status = XPRSminim(lp->prob,"l");
   else                 status = XPRSmaxim(lp->prob,"l");
   if (status)
       std::cout << "XPRSoptimize failed with error code " << status << std::endl;
   return status;
#endif
// End adding PWOF

#ifdef VSCIP
    //std::cout << "\nJJoptimize" << std::endl;
    dual = false;
    if (SCIPlpiSolvePrimal(lp) == SCIP_OKAY) return 0;
    return 1;
#endif
}

int JJdualopt  (
JJLPptr lp)
{    
#ifdef CPLEX3      
   return dualopt (lp);
#endif

#ifdef CPLEX5      
   return CPXdualopt (CPLEXv::Env, lp);
#endif


#ifdef XPRESS
   JJcheck(lp);
   seticv(N_IALG,1);
   if( lp->objsen == 1 ) return minim("l");
   return maxim("l");
#endif 

   // Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   if (lp->objsen == 1) status = XPRSminim(lp->prob,"dl");
   else                 status = XPRSmaxim(lp->prob,"dl");
   if (status)
       std::cout << "XPRSdualopt failed with error code " << status << std::endl;
   return status;
#endif
// End adding PWOF

#ifdef VSCIP
    //std::cout << "\nJJdualopt" << std::endl;
    dual = true;
    if (SCIPlpiSolveDual(lp) == SCIP_OKAY) return 0;
    return 1;
#endif
}

int JJhybnetopt (
JJLPptr lp,
char method)
{
#ifdef CPLEX3      
   return hybnetopt (lp, method);
#endif

#ifdef CPLEX5      
   return CPXhybnetopt (CPLEXv::Env, lp, method);
#endif


#ifdef XPRESS
   JJcheck(lp);
   seticv(N_IALG,1);
   if( lp->objsen == 1 ) return minim("l");
   return maxim("l");
#endif 

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   if (lp->objsen == 1) status = XPRSminim(lp->prob,"dl");
   else                 status = XPRSmaxim(lp->prob,"l");
   if (status)
       std::cout << "JJhypnetopt failed with error code " << status << std::endl;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    std::cout << "\nJJhybnetopt?" << std::endl;
    system("pause");
    return 1;
#endif
}

int JJnetopt (
JJLPptr lp,
int *netstatus_p,
int *numnodes_p,
int *numarcs_p,
int *itcnt_p
)
{
#ifdef CPLEX3      
   return netopt (lp,netstatus_p,numnodes_p,numarcs_p,itcnt_p);
#endif


#ifdef CPLEX5      
#ifndef CPLEX7
   return CPXnetopt (Env,lp,netstatus_p,numnodes_p,numarcs_p,itcnt_p);
#endif
#ifdef CPLEX7
   *numnodes_p  = 0;
   *numarcs_p   = 0;
   *itcnt_p     = 0;
   *netstatus_p = CPXdualopt(CPLEXv::Env,lp);
   return *netstatus_p;
   // return CPXhybnetopt(Env,lp,CPX_ALG_DUAL);
#endif
#endif


#ifdef XPRESS
   int status;

   JJcheck(lp);
   seticv(N_IALG,1);
   if( lp->objsen == 1 ) status = minim("l");
   else                  status = maxim("l");
   if( status ) {
	   *netstatus_p = 1101;     // JJ_INForUNB
	   return status;
   }
   return getipv(N_STATUS, netstatus_p);
#endif 

   // Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   if (lp->objsen == 1) status = XPRSminim(lp->prob,"l");
   else                 status = XPRSmaxim(lp->prob,"l");
   if (status)
       std::cout << "JJnetopt failed with error code " << status << std::endl;
   return status;
#endif
// End adding PWOF

#ifdef VSCIP
    //std::cout << "\nJJnetopt" << std::endl;
    *numnodes_p  = 0;
    *numarcs_p   = 0;
    *itcnt_p     = 0;
    dual = true;
    if (SCIPlpiSolveDual(lp) == SCIP_OKAY) *netstatus_p = 0;
    else *netstatus_p = 1;
    return *netstatus_p;
#endif
}

int JJmipopt  (
JJLPptr lp)
{    

#ifdef CPLEX5      
   return CPXmipopt (CPLEXv::Env, lp);
#endif 

#ifdef XPRESS_13
    std::cout << "\nJJmipopt?" << std::endl;
    system("pause");
    return 1;
#endif
   
#ifdef VSCIP
    //std::cout << "\nJJmipopt" << std::endl;
    dual = false;
    if (SCIPlpiSolvePrimal(lp) == SCIP_OKAY) return 0;
    return 1;
#endif
}

int JJsetscr_ind(
JJLPptr lp,
int scr_ind)
{
#ifdef CPLEX3
    return setscr_ind(scr_ind);
#endif


#ifdef CPLEX5      
    if( scr_ind )
        return CPXsetintparam( CPLEXv::Env , CPX_PARAM_SCRIND , CPX_ON );
    return CPXsetintparam( CPLEXv::Env , CPX_PARAM_SCRIND , CPX_OFF );
#endif


#ifdef XPRESS
   int status;
   status = setoptlog("xosl.log");
//   std::cout << "JJ WARNING: not seticv" << std::endl;
//   status = seticv(N_IFMSG,1);
   status = seticv(N_PRTMSG,scr_ind);
   return 0;
#endif 

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   status = XPRSsetlogfile(lp->prob,"xosl.log");
   if (status)
   {
       std::cout << "Error in XPRESS JJsetscr_ind() - 1" << std::endl;
       return status;
   }
   status = XPRSsetintcontrol(lp->prob, XPRS_OUTPUTLOG, scr_ind);
   if (status)
       std::cout << "Error in XPRESS JJsetscr_ind() - 2" << std::endl;
   return status;
#endif
// End adding PWOF

   
#ifdef VSCIP
    //std::cout << "\nJJsetscr_ind" << std::endl;
    SCIPlpiSetIntpar(Env[lpJJ-1],SCIP_LPPAR_LPINFO,scr_ind);
    return 0;
#endif
}

int JJgetmac(
JJLPptr lp)
{
#ifdef CPLEX3
   return getmac (lp);
#endif


#ifdef CPLEX5
   return CPXgetnumcols (CPLEXv::Env, lp);
#endif


#ifdef XPRESS
   int mac;
   JJcheck(lp);
   getipv(N_NCOL, &mac);
   return mac;
#endif 

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int mac, status;
   status = XPRSgetintattrib(lp->prob, XPRS_COLS, &mac);
   if (status)
       std::cout << "Error in XPRESS JJgetmac()" << std::endl;
   return mac;
#endif
// End adding PWOF

#ifdef VSCIP
    //std::cout << "\nJJgetmac" << std::endl;
    int ncols;
    SCIPlpiGetNCols(lp,&ncols);
    return ncols;
#endif
}

int JJgetmar(
JJLPptr lp)
{
#ifdef CPLEX3
   return getmar (lp);
#endif


#ifdef CPLEX5
   return CPXgetnumrows (CPLEXv::Env, lp);
#endif

#ifdef XPRESS
   int mar;
   JJcheck(lp);
   getipv(N_NROW, &mar);
   return mar;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int mar, status;
   status = XPRSgetintattrib(lp->prob, XPRS_ROWS, &mar);
   if (status)
       std::cout << "Error in XPRESS JJgetmar()" << std::endl;
   return mar;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJgetmar" << std::endl;
    int nrows;
    SCIPlpiGetNRows(lp,&nrows);
    return nrows;
#endif
}

int JJgetmat(
JJLPptr lp)
{
#ifdef CPLEX3
   return getmat (lp);
#endif


#ifdef CPLEX5
   return CPXgetnumnz (CPLEXv::Env, lp);
#endif


#ifdef XPRESS
   int mat;
   JJcheck(lp);
   getipv(N_NELEM, &mat);
   return mat;
#endif
   
// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int mat, status;
   status = XPRSgetintattrib(lp->prob, XPRS_ELEMS, &mat);
   if (status)
       std::cout << "Error in XPRESS JJgetmat()" << std::endl;
   return mat;
#endif
// End adding PWOF


#ifdef VSCIP
    //std::cout << "\nJJgetmat" << std::endl;
    int num;
    SCIPlpiGetNNonz(lp,&num);    
    return num;
#endif
}


int JJlpiterlimit(
/*JJLPptr lp,*/
int val)
{
#ifdef CPLEX3
   return 1;
#endif


#ifdef CPLEX5
   return CPXsetintparam (CPLEXv::Env, CPX_PARAM_ITLIM, val);
#endif


#ifdef XPRESS
   return 1;   //XPRSsetintcontrol(  ,  XPRS_LPITERLIMIT, val);
#endif

#ifdef XPRESS_13
   return 1;   //XPRSsetintcontrol(  ,  XPRS_LPITERLIMIT, val);
#endif

#ifdef VSCIP
    //std::cout << "\nJJlpiterlimit" << std::endl;
    SCIPlpiSetIntpar(Env[lpJJ-1],SCIP_LPPAR_LPITLIM,val);
    return 0;
#endif
}




int JJgetobjsen(
JJLPptr lp)
{
#ifdef CPLEX3
   return getobjsen (lp);
#endif


#ifdef CPLEX5
   return CPXgetobjsen (CPLEXv::Env,lp);
#endif


#ifdef XPRESS
   JJcheck(lp);
   return lp->objsen;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   return lp->objsen;
#endif
// End adding PWOF
   
#ifdef VSCIP
    std::cout << "\nJJgetobjsen?" << std::endl;
    system("pause");
    return 1;
#endif
}

int JJgetobj(
JJLPptr lp,
double *obj,
int    begin, int end
)
{
#ifdef CPLEX3
   return getobj (lp, obj, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetobj (CPLEXv::Env, lp, obj, begin, end);
#endif


#ifdef XPRESS
   JJcheck(lp);
   return getobj(obj,begin,end);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   status = XPRSgetobj(lp->prob,obj,begin,end);
   if (status) 
       std::cout << "Error XPRESS JJgetobj() " << status << std::endl;
   return status;
#endif
// End adding PWOF

#ifdef VSCIP
    //std::cout << "\nJJgetobj" << std::endl;
    SCIPlpiGetObj(lp,begin,end,obj);
    return 0;
#endif
}

int JJgetrhs (
JJLPptr lp,
double *rhs,
int    begin,int end
)
{
#ifdef CPLEX3
   return getrhs (lp, rhs, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetrhs (CPLEXv::Env, lp, rhs, begin, end);
#endif


#ifdef XPRESS
   JJcheck(lp);
   return getrhs (rhs, begin, end);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   status = XPRSgetrhs(lp->prob,rhs,begin,end);
   if (status) 
       std::cout << "Error XPRESS JJgetrhs() " << status << std::endl;
   return status;
#endif
// End adding PWOF

#ifdef VSCIP
    //std::cout << "\nJJgetrhs" << std::endl;
    SCIP_Real *lhss = new SCIP_Real[end-begin+1];
    SCIP_Real *rhss = new SCIP_Real[end-begin+1];
    SCIPlpiGetSides(lp,begin,end,lhss,rhss);
    for (int i = 0; i <= end-begin; i++)
    {
        if (SCIPlpiIsInfinity(lp,rhss[i]))
            rhs[i] = lhss[i];
        else
            rhs[i] = rhss[i];
    }
    //PWOF change 23-08-2013
    //delete lhss;
    //delete rhss;
    delete[] lhss;
    delete[] rhss;
    return 0;
#endif
}

int JJgetsense (
JJLPptr lp,
char   *sense,
int    begin,int end
)
{
#ifdef CPLEX3
   return getsense (lp, sense, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetsense (CPLEXv::Env, lp, sense, begin, end);
#endif


#ifdef XPRESS
   JJcheck(lp);
   return getrowtype (sense, begin, end);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   status = XPRSgetrowtype(lp->prob,sense,begin,end);
   if (status) 
       std::cout << "Error XPRESS JJgetsense() " << status << std::endl;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJgetsense" << std::endl;
    SCIP_Real *lhss = new SCIP_Real[end-begin+1];
    SCIP_Real *rhss = new SCIP_Real[end-begin+1];
    SCIPlpiGetSides(lp,begin,end,lhss,rhss);
    for (int i = 0; i <= end-begin; i++)
    {
        if (SCIPlpiIsInfinity(lp,-lhss[i]))
            sense[i] = 'L';
        else
            if (SCIPlpiIsInfinity(lp,rhss[i]))
                sense[i] = 'G';
            else
                sense[i] = 'E';
    }
    //PWOF change 23-08-2013
    //delete lhss;
    //delete rhss;
    delete[] lhss;
    delete[] rhss;
    return 0;
#endif
}


int JJgetcols (
JJLPptr lp,
int *nzcnt, int *cmatbeg, int *cmatind, 
double *cmatval,
int cmatspace, 
int *surplus,
int begin,  int end
)
{
#ifdef CPLEX3
   return getcols (lp, nzcnt, cmatbeg, cmatind, cmatval,
                       cmatspace, surplus, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetcols (CPLEXv::Env, lp, nzcnt, cmatbeg, cmatind, cmatval,
                       cmatspace, surplus, begin, end);
#endif


#ifdef XPRESS
   int status;

   JJcheck(lp);
   status = getcols( cmatbeg, cmatind, cmatval, cmatspace, nzcnt, begin, end);
   *surplus = cmatspace - *nzcnt;
   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   status = XPRSgetcols(lp->prob, cmatbeg, cmatind, cmatval, cmatspace, nzcnt, begin, end);
   if (status) 
       std::cout << "Error XPRESS JJgetcols() " << status << std::endl;
   *surplus = cmatspace - *nzcnt;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJgetcols" << std::endl;
    SCIPlpiGetCols(lp,begin,end,NULL,NULL,nzcnt,cmatbeg,cmatind,cmatval);
    *surplus = cmatspace - *nzcnt;
    return 0;
#endif
}

int JJgetrows(
JJLPptr lp,
int *nzcnt, int *rmatbeg, int *rmatind, double *rmatval,
int rmatspace, 
int *surplus,
int begin,  int end
)
{
#ifdef CPLEX3
   return getrows (lp, nzcnt, rmatbeg, rmatind, rmatval,
                       rmatspace, surplus, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetrows (CPLEXv::Env, lp, nzcnt, rmatbeg, rmatind, rmatval,
                       rmatspace, surplus, begin, end);
#endif


#ifdef XPRESS
   int status;

   JJcheck(lp);
   status = getrows( rmatbeg, rmatind, rmatval, rmatspace, nzcnt, begin, end);
   *surplus = rmatspace - *nzcnt;
   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   status = XPRSgetrows(lp->prob, rmatbeg, rmatind, rmatval, rmatspace, nzcnt, begin, end);
   if (status) 
       std::cout << "Error XPRESS JJgetrows() " << status << std::endl;
   *surplus = rmatspace - *nzcnt;
   return status;
#endif
// End adding PWOF

   
#ifdef VSCIP
    //std::cout << "\nJJgetrows" << std::endl;
    SCIPlpiGetRows(lp,begin,end,NULL,NULL,nzcnt,rmatbeg,rmatind,rmatval);
    *surplus = rmatspace - *nzcnt;
    return 0;
#endif
}

int JJsolution (
JJLPptr lp,
int *lpstat_p,
double *objval_p, double *x, double *pi, double *slack, double *dj
)
{
#ifdef CPLEX3
   return solution (lp, lpstat_p, objval_p, x, pi, slack, dj);
#endif


#ifdef CPLEX5
   return CPXsolution (CPLEXv::Env, lp, lpstat_p, objval_p, x, pi, slack, dj);
#endif


#ifdef XPRESS
#ifndef XPRESS11
   int k;
#endif
   int status;

   JJcheck(lp);
   if( lpstat_p ){
       status = getipv(N_STATUS, lpstat_p);
       if( status ) return status;
   }
   if( objval_p ){
       status = getdpv(N_DOBJVL, objval_p);
       if( status ) return status;
   }
   status = solution (x, slack, pi, dj);
#ifndef XPRESS11
   if( pi && lp->objsen== -1 ){
//     fprintf(stderr,"INFO: changing dual variable signs\n");
       getipv(N_NROW,&k);
       while(k--) pi[k] = -pi[k];
   }
   if( dj && lp->objsen== -1 ){
       fprintf(stderr,"INFO: changing reduced cost signs\n");
       getipv(N_NCOL,&k);
       while(k--) dj[k] = -dj[k];
   }
#endif
   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   if (lpstat_p)
   {
        status = XPRSgetintattrib(lp->prob, XPRS_LPSTATUS, lpstat_p);
        if (status) 
                std::cout << "Error XPRESS JJsolution() " << status << std::endl;
        return status;
   }
   if (objval_p)
   {
        status = XPRSgetdblattrib(lp->prob, XPRS_LPOBJVAL, objval_p);
        if (status) 
                std::cout << "Error XPRESS JJsolution() " << status << std::endl;
        return status;
   }
   status = XPRSgetsol(lp->prob, x, slack, pi, dj);
   if (status) 
        std::cout << "Error XPRESS JJsolution() " << status << std::endl;
   return status;
#endif
// End adding PWOF

   
#ifdef VSCIP
    //std::cout << "\nJJsolution" << std::endl;
    *lpstat_p = JJgetstat(lp);
    SCIPlpiGetSol(lp,objval_p,x,pi,NULL/*slack*/,dj);
    JJgetslack(lp,slack, 0,JJgetmar(lp)-1);
    return 0;
#endif
}

int JJgetstat (
JJLPptr lp)
{
#ifdef CPLEX3
   return getstat (lp);
#endif


#ifdef CPLEX5
   return CPXgetstat (CPLEXv::Env, lp);
#endif


#ifdef XPRESS
   int lpstat;
   
   JJcheck(lp);
   getipv(N_STATUS, &lpstat);
   return lpstat;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status, lpstat;
   status = XPRSgetintattrib(lp->prob, XPRS_LPSTATUS, &lpstat);
   if (status) 
           std::cout << "Error XPRESS JJsolution() " << status << std::endl;
   if (lpstat == 5) lpstat = 3;   /* JJ_UNBOUNDED */
   if (lpstat == 3) lpstat = -1000;
   return lpstat;
#endif
// End adding PWOF   
   
#ifdef VSCIP    // Salomé modified 31/01/2014
    //std::cout << "\nJJgetstat" << std::endl;
    int m_stat = SCIPlpiGetInternalStatus(lp);
	#ifdef soplex    
    switch( m_stat )  //soplex
    {
        case 1: return 1;     // optimal
        case 2: return 2;     // unbounded
        case 3: return 3;     // infeasible
        case -7: return 11;   // timelimExc
        case -6: return 10;   // iterlimExc
        case -5: return 12;   // objlimExc
        default:
           //std::cout << "error " << m_stat << std::endl;
           return 1101;
    }
	#endif
    #ifdef clp
    /* We first check if status is ok, i.e., is one of the following:
    * 0 - optimal
    * 1 - primal infeasible
    * 2 - dual infeasible
    * 3 - stopped on iterations or time
    * 4 - stopped due to errors
    * 5 - stopped by event handler (virtual int ClpEventHandler::event())
    */
    if (m_stat == 0) return 1;   // optimal
    if (m_stat == 1) return 3;   // infeasible
    if (m_stat == 2) return 2;   // unbounded
    if (SCIPlpiIsTimelimExc(lp)) return 11;  // timelimExc
    if (SCIPlpiIsIterlimExc(lp)) return 10;  // iterlimExc
    if (SCIPlpiIsObjlimExc(lp)) return 12;   // objlimExc
    return 1101;
	#endif
#endif
}

int JJgetobjval (
JJLPptr lp,
double *objval_p)
{
#ifdef CPLEX3
   return getobjval (lp, objval_p);
#endif


#ifdef CPLEX5
   return CPXgetobjval (CPLEXv::Env, lp, objval_p);
#endif


#ifdef XPRESS
   JJcheck(lp);
   return getdpv(N_DOBJVL, objval_p);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   status = XPRSgetdblattrib(lp->prob,XPRS_LPOBJVAL,objval_p);
   if (status) 
       std::cout << "Error XPRESS JJgetobjval() " << status << std::endl;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJgetobjval" << std::endl;
    SCIPlpiGetObjval(lp,objval_p);    
    return 0;
#endif
}

int JJgetx ( 
JJLPptr lp,
double *x,
int    begin,int end
)
{
#ifdef CPLEX3
   return getx (lp, x, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetx (CPLEXv::Env, lp, x, begin, end);
#endif


#ifdef XPRESS
   int k,status;
   double *ptr;

   k = JJgetmac(lp); 
   ptr = (double *)malloc( k * sizeof( double ) );
   if( ptr==NULL ) return 1;
   status = JJsolution (lp, NULL, NULL, ptr, NULL, NULL, NULL);
   if( status==0 )
       for( k=begin ; k<=end ; k++ )
           x[k-begin] = ptr[k];
   free( ptr );
   ptr = NULL; /*PWOF*/
   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int k, status;
   double *ptr=NULL;
   
   status = XPRSgetintattrib(lp->prob,XPRS_COLS,&k);
   if (status) 
   {
       std::cout << "Error XPRESS JJgetx() " << status << std::endl;
       return status;
   }
   ptr = (double *)malloc(k*sizeof(double));
   if (ptr==NULL) return 1;
   
   status = XPRSgetsol(lp->prob,ptr, NULL,NULL,NULL);
   if (status) std::cout << "Error XPRESS JJgetx() " << status << std::endl;
   if (status==0)
       for (k=begin;k<=end;k++)
           x[k-begin] = ptr[k];
   free(ptr);
   ptr=NULL;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJgetx" << std::endl;
    SCIP_Real *primsol = new SCIP_Real[JJgetmac(lp)];
    SCIPlpiGetSol(lp,NULL,primsol,NULL,NULL,NULL);

    for (int i = begin; i <= end; i++) {
        x[i-begin] = primsol[i];        
    }
    //PWOF change 23-08-2013
    //delete primsol;
    delete[] primsol;
    return 0;
#endif
}

int JJgetpi ( 
JJLPptr lp,
double *pi,
int    begin,int end
)
{
#ifdef CPLEX3
   return getpi (lp, pi, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetpi (CPLEXv::Env, lp, pi, begin, end);
#endif


#ifdef XPRESS
   int k,status;
   double *ptr;

   k = JJgetmar(lp); 
   ptr = (double *)malloc( k * sizeof( double ) );
   if( ptr==NULL ) return 1;
   status = JJsolution (lp, NULL, NULL, NULL, ptr, NULL, NULL);
   if( status==0 )
       for( k=begin ; k<=end ; k++ )
           pi[k-begin] = ptr[k];
   free( ptr );
   ptr = NULL; /*PWOF*/
   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int k, status;
   double *ptr=NULL;
   
   status = XPRSgetintattrib(lp->prob,XPRS_ROWS,&k);
   if (status) 
   {
       std::cout << "Error XPRESS JJgetpi() " << status << std::endl;
       return status;
   }
   ptr = (double *)malloc(k*sizeof(double));
   if (ptr==NULL) return 1;
   
   status = XPRSgetsol(lp->prob,NULL,NULL,ptr,NULL);
   if (status) std::cout << "Error XPRESS JJgetpi() " << status << std::endl;
   if (status==0)
       for (k=begin;k<=end;k++)
           pi[k-begin] = ptr[k];
   free(ptr);
   ptr=NULL;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJgetpi" << std::endl;

    SCIP_Real *dualsol = new SCIP_Real[JJgetmar(lp)];
    SCIPlpiGetSol(lp,NULL,NULL,dualsol,NULL,NULL);

    for (int i = begin; i <= end; i++) {
        pi[i-begin] = dualsol[i];        
    }
    //PWOF change 23-08-2013
    //delete dualsol;    
    delete[] dualsol;    
    return 0;
#endif
}

int JJgetslack ( 
JJLPptr lp,
double *slack,
int    begin,int end
)
{
#ifdef CPLEX3
   return getslack (lp, slack, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetslack (CPLEXv::Env, lp, slack, begin, end);
#endif


#ifdef XPRESS
   int k,status;
   double *ptr;

   k = JJgetmar(lp); 
   ptr = (double *)malloc( k * sizeof( double ) );
   if( ptr==NULL ) return 1;
   status = JJsolution (lp, NULL, NULL, NULL, NULL, ptr, NULL);
   if( status==0 )
       for( k=begin ; k<=end ; k++ )
           slack[k-begin] = ptr[k];
   free( ptr );
   ptr = NULL; /*PWOF*/
   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int k, status;
   double *ptr=NULL;
   
   status = XPRSgetintattrib(lp->prob,XPRS_ROWS,&k);
   if (status) 
   {
       std::cout << "Error XPRESS JJgetslack() " << status << std::endl;
       return status;
   }
   ptr = (double *)malloc(k*sizeof(double));
   if (ptr==NULL) return 1;
   
   status = XPRSgetsol(lp->prob,NULL,ptr,NULL,NULL);
   if (status) std::cout << "Error XPRESS JJgetslack() " << status << std::endl;
   if (status==0)
       for (k=begin;k<=end;k++)
           slack[k-begin] = ptr[k];
   free(ptr);
   ptr=NULL;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJgetslack" << std::endl;//  Ax <= b => slack = b - Ax*
    SCIP_Real *rhs = new SCIP_Real[end-begin+1];
    SCIP_Real *lhs = new SCIP_Real[end-begin+1];
    SCIP_Real *activ = new SCIP_Real[JJgetmar(lp)];
    SCIPlpiGetSol(lp,NULL,NULL,NULL,activ,NULL);
    SCIPlpiGetSides(lp,begin,end,lhs,rhs);
    for (int i = begin; i <= end; i++) {
        if (SCIPlpiIsInfinity(lp,rhs[i-begin]))
            slack[i-begin] = lhs[i-begin] - activ[i];
        else
            slack[i-begin] = rhs[i-begin] - activ[i];        
    }
    //PWOF change 23-08-2013
    //delete activ;    
    //delete rhs;
    //delete lhs;
    delete[] activ;    
    delete[] rhs;
    delete[] lhs;
    return 0;
#endif
}

int JJgetdj ( 
JJLPptr lp,
double *dj,
int    begin,int end
)
{
#ifdef CPLEX3
   return getdj (lp, dj, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetdj (CPLEXv::Env, lp, dj, begin, end);
#endif


#ifdef XPRESS
   int k,status;
   double *ptr;

   k = JJgetmac(lp); 
   ptr = (double *)malloc( k * sizeof( double ) );
   if( ptr==NULL ) return 1;
   status = JJsolution (lp, NULL, NULL, NULL, NULL, NULL, ptr);
   if( status==0 )
       for( k=begin ; k<=end ; k++ )
           dj[k-begin] = ptr[k];
   free( ptr );
   ptr = NULL; /*PWOF*/
   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int k, status;
   double *ptr=NULL;
   
   status = XPRSgetintattrib(lp->prob,XPRS_COLS,&k);
   if (status) 
   {
       std::cout << "Error XPRESS JJgetdj() " << status << std::endl;
       return status;
   }
   ptr = (double *)malloc(k*sizeof(double));
   if (ptr==NULL) return 1;
   
   status = XPRSgetsol(lp->prob,NULL,NULL,NULL,ptr);
   if (status) std::cout << "Error XPRESS JJgetdj() " << status << std::endl;
   if (status==0)
       for (k=begin;k<=end;k++)
           dj[k-begin] = ptr[k];
   free(ptr);
   ptr=NULL;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJgetdj" << std::endl;
    SCIP_Real *redcost = new SCIP_Real[JJgetmac(lp)];
    SCIPlpiGetSol(lp,NULL,NULL,NULL,NULL,redcost);

    for (int i = begin; i <= end; i++) {
        dj[i-begin] = redcost[i];
    }
    //PWOF change 23-08-2013
    //delete redcost;
    delete[] redcost;
    return 0;
#endif
}



int JJgetitc ( 
JJLPptr lp)
{
#ifdef CPLEX3
   return getitc (lp);
#endif


#ifdef CPLEX5
   return CPXgetitcnt (CPLEXv::Env, lp);
#endif


#ifdef XPRESS
   int itcnt;
   JJcheck(lp);
   getipv(N_ITCNT, &itcnt);
   return itcnt;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int itcnt, status;
      
   status = XPRSgetintattrib(lp->prob,XPRS_SIMPLEXITER,&itcnt);
   if (status) 
       std::cout << "Error XPRESS JJgetitc() " << status << std::endl;
   return itcnt;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJgetitc" << std::endl;
    int iterations;
    SCIPlpiGetIterations(lp,&iterations);    
    return iterations;
#endif
}

int JJgetitci ( 
JJLPptr lp)
{
#ifdef CPLEX3
   return getitci (lp);
#endif


#ifdef CPLEX5
   return CPXgetphase1cnt (CPLEXv::Env, lp);
#endif


#ifdef XPRESS
   JJcheck(lp);
   return 0;
#endif
   
// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   return 0;
#endif
// End adding PWOF

#ifdef VSCIP
    //std::cout << "\nJJgetitci" << std::endl;
    int status = JJgetstat(lp);
    if (( status != 1 ) && (status != 2) && (status != 3))
        return 0;
    int iterations;
    SCIPlpiGetIterations(lp,&iterations);
    return iterations;
#endif
}

int JJgetbase ( 
JJLPptr lp,
int *cstat, int *rstat
)
{
#ifdef CPLEX3
   return getbase (lp, cstat, rstat);
#endif


#ifdef CPLEX5
   return CPXgetbase (CPLEXv::Env, lp, cstat, rstat);
#endif


#ifdef XPRESS
   JJcheck(lp);
   return getbasis(rstat, cstat);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   status = XPRSgetbasis(lp->prob,rstat,cstat);
   if (status)
   {
       std::cout << "Error XPRESS JJgetbase() " << status << std::endl;
       return status;
   }
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJgetbase" << std::endl;
    SCIPlpiGetBase(lp,cstat,rstat);
    return 0;
#endif
}

int JJloadbase( 
JJLPptr lp,
int *cstat, int *rstat
)
{
#ifdef CPLEX3
   return loadbase (lp, cstat, rstat);
#endif


#ifdef CPLEX5
   return CPXcopybase (CPLEXv::Env, lp, cstat, rstat);
#endif


#ifdef XPRESS
   JJcheck(lp);
   return loadbasis (rstat, cstat);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   status = XPRSloadbasis(lp->prob,rstat,cstat);
   if (status)
   {
       std::cout << "Error XPRESS JJloadbase() " << status << std::endl;
       return status;
   }
   return status;
#endif
// End adding PWOF
   
   
#ifdef VSCIP
    //std::cout << "\nJJloadbase" << std::endl;
    SCIPlpiSetBase(lp,cstat,rstat);
    return 0;
#endif
}


int JJaddrows (
JJLPptr lp,
int ccnt, int rcnt, int nzcnt,
double *rhs,
char *sense,
int *rmatbeg, int *rmatind,
double *rmatval,
char **colname, char **rowname
)
{
#ifdef CPLEX3
   return addrows (lp, ccnt, rcnt, nzcnt, rhs, sense, rmatbeg,
                  rmatind, rmatval, colname, rowname);
#endif


#ifdef CPLEX5
   return CPXaddrows (CPLEXv::Env, lp, ccnt, rcnt, nzcnt, rhs, sense, rmatbeg,
                  rmatind, rmatval, colname, rowname);
#endif


#ifdef XPRESS
   JJcheck(lp);
   if ( ccnt ) fprintf(stderr,"JJ WARNING: not new columns for new rows!!!\n");
   if ( colname ) fprintf(stderr,"JJ WARNING: not name to new columns!!!\n");
   if ( rowname ) fprintf(stderr,"JJ WARNING: not name to new rows!!!\n");
   return addrows (rcnt, nzcnt, sense, rhs, NULL, rmatbeg, rmatind, rmatval);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   if (ccnt) std::cout << "JJ WARNING: not new columns for new rows!!!" << std::endl;
   if (colname) std::cout << "JJ WARNING: not name to new columns!!!" << std::endl;
   if (rowname) std::cout << "JJ WARNING: not name to new rows!!!" << std::endl;
   status = XPRSaddrows(lp->prob, rcnt, nzcnt, sense, rhs, NULL, rmatbeg, rmatind, rmatval);
   if (status)
       std::cout << "Error XPRESS JJaddrows() " << status << std::endl;
   return status;
#endif
// End adding PWOF

   
#ifdef VSCIP
    //std::cout << "\nJJaddrows" << std::endl;
    if (ccnt != 0)   // add cols
    {
        int* mybeg = new int[ccnt + 1];
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
        SCIPlpiAddCols(lp,ccnt,obj,lb,ub,colname,nzcnt,mybeg,0,0);

        delete[] mybeg;
        delete[] obj;
        delete[] lb;
        delete[] ub;
    }
    double *rowLower = new double[rcnt];
    double *rowUpper = new double[rcnt];
    for (int i = 0; i <  rcnt; i++)
    {
        if (sense[i] == 'L') {
            rowLower[i] = -SCIPlpiInfinity(lp);
            rowUpper[i] = rhs[i];
        }else
            if (sense[i] == 'E') {
                rowLower[i] = rhs[i];
                rowUpper[i] = rhs[i];
            }else
                if (sense[i] == 'G') {
                    rowLower[i] = rhs[i];
                    rowUpper[i] = SCIPlpiInfinity(lp);
                }
    }

    SCIPlpiAddRows(lp,rcnt,rowLower,rowUpper,rowname,nzcnt,rmatbeg,rmatind,rmatval);

    //PWOF change 23-08-2013
    //delete rowLower;
    //delete rowUpper;
    delete[] rowLower;
    delete[] rowUpper;
    return 0;
#endif
}

int JJdelrows( 
JJLPptr lp,
int    begin,int end
)
{
#ifdef CPLEX3
   return delrows (lp, begin, end);
#endif


#ifdef CPLEX5
   return CPXdelrows (CPLEXv::Env, lp, begin, end);
#endif


#ifdef XPRESS
   int nrows,status,i;
   int *mindex;

   JJcheck(lp);
   nrows = end-begin+1;
   mindex = (int *)malloc( nrows*sizeof(int) );
   if( mindex==NULL ){
       fprintf(stderr,"JJ WARNING: not space in del-rows!!!\n");
       return 1;
   }
//   for(i=begin;i<=end;i++) mindex[i] = i;   //  it was way before!!!
   for(i=begin;i<=end;i++) mindex[i-begin] = i;    // Rafael marz 2001
   status = delrows(nrows, mindex);
   free( mindex );
   mindex = NULL; /*PWOF*/
   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int nrows_, status, i;
   int *mindex = NULL;
   
   nrows_ = end - begin + 1;
   mindex = (int *)malloc( nrows_*sizeof(int) );
   if (mindex==NULL)
   {
       std::cout << "JJ WARNING: not space in del-rows!!!" << std::endl;
       return 1;
   }
   for (i=begin;i<=end;i++) mindex[i-begin] = i;
   status = XPRSdelrows(lp->prob, nrows_, mindex);
   if (status)
       std::cout << "Error XPRESS JJdelrows() " << status << std::endl;
   free(mindex);
   mindex = NULL; /*PWOF*/
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJdelrows" << std::endl;
    SCIPlpiDelRows(lp,begin,end);
    return 0;
#endif
}

int JJdelsetrows( 
JJLPptr lp,
int    *delstat
)
{
#ifdef CPLEX3
   return delsetrows (lp, delstat);
#endif

#ifdef CPLEX5
   return CPXdelsetrows (CPLEXv::Env, lp, delstat);
#endif

#ifdef XPRESS
   int nrows,status,mar;
   int *mindex;

   JJcheck(lp);
   nrows = 0;
   getipv (N_NROW, &mar);
   while(mar--)
       if( delstat[mar] ) nrows++;
   if( nrows==0 ) return 0;
   mindex = (int *)malloc( nrows*sizeof(int) );
   if( mindex==NULL ){
       fprintf(stderr,"JJ WARNING: not space for del-set-rows!!!\n");
       return 1;
   }
   nrows = 0;
   getipv (N_NROW, &mar);
   while(mar--)
       if( delstat[mar] ) mindex[nrows++] = mar;
   status = delrows(nrows, mindex);
   free( mindex );
   mindex = NULL; /*PWOF*/
   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int nrows_, status, mar;
   int *mindex = NULL;
   
   nrows_ = 0;
   status = XPRSgetintattrib(lp->prob, XPRS_ROWS, &mar);
   if (status) 
   {
       std::cout << "Error XPRESS JJdelsetrows() " << status << std::endl;
       return status;
   }
   
   while (mar--)
       if (delstat[mar]) nrows_++;
   if (nrows_==0) return 0;
   
   mindex = (int *)malloc( nrows_*sizeof(int) );
   if (mindex==NULL)
   {
       std::cout << "JJ WARNING: not space for del-set-rows!!!" << std::endl;
       return 1;
   }
   nrows_=0;
   status = XPRSgetintattrib(lp->prob, XPRS_ROWS, &mar);
   if (status) 
   {
       std::cout << "Error XPRESS JJdelsetrows() " << status << std::endl;
       free(mindex);
       mindex=NULL;
       return status;
   }
   while (mar--)
       if (delstat[mar]) mindex[nrows_++] = mar;
   
   status = XPRSdelrows(lp->prob, nrows_, mindex);
   
   if (status)
       std::cout << "Error XPRESS JJdelsetrows() " << status << std::endl;
   free(mindex);
   mindex=NULL;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJdelsetrows" << std::endl;
    SCIPlpiDelRowset(lp,delstat);
    return 0;
#endif
}

int JJaddcols (
JJLPptr lp,
int ccnt, int nzcnt,
double *obj,
int *cmatbeg, int *cmatind,
double *cmatval,double *lb,double *ub,
char **colname
)
{
#ifdef CPLEX3
   return addcols (lp, ccnt, nzcnt, obj, cmatbeg,
                  cmatind, cmatval, lb, ub, colname );
#endif


#ifdef CPLEX5
   return CPXaddcols (CPLEXv::Env, lp, ccnt, nzcnt, obj, cmatbeg,
                  cmatind, cmatval, lb, ub, colname );
#endif


#ifdef XPRESS
   JJcheck(lp);
   if ( colname ) fprintf(stderr,"JJ WARNING: not name to new columns!!!\n");
   return addcols (ccnt, nzcnt, obj, cmatbeg, cmatind, cmatval, lb, ub);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   if (colname) std::cout << "JJ WARNING: not name to new columns!!!" << std::endl;
   status = XPRSaddcols(lp->prob, ccnt, nzcnt, obj, cmatbeg, cmatind, cmatval, lb, ub);
   if (status)
       std::cout << "Error XPRESS JJaddcols() " << status << std::endl;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJaddcols" << std::endl;
    SCIPlpiAddCols(lp,ccnt,obj,lb,ub,colname,nzcnt,cmatbeg,cmatind,cmatval);
    return 0;
#endif
}

int JJdelcols( 
JJLPptr lp,
int    begin,int end
)
{
#ifdef CPLEX3
   return delcols (lp, begin, end);
#endif


#ifdef CPLEX5
   return CPXdelcols (CPLEXv::Env, lp, begin, end);
#endif


#ifdef XPRESS
   int ncols,status,i;
   int *mindex;

   JJcheck(lp);
   ncols = end-begin+1;
   mindex = (int *)malloc( ncols*sizeof(int) );
   if( mindex==NULL ){
       fprintf(stderr,"JJ WARNING: not space for del-columns!!!\n");
       return 1;
   }
//   for(i=begin;i<=end;i++) mindex[i] = i;  //   before!!!
   for(i=begin;i<=end;i++) mindex[i-begin] = i;    // Rafael marzo 2001
   status = delcols(ncols, mindex);
   free( mindex );
   mindex = NULL; /*PWOF*/
   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int ncols_, status, i;
   int *mindex = NULL;
   
   ncols_ = end - begin + 1;
   mindex = (int *)malloc( ncols_*sizeof(int) );
   if (mindex==NULL)
   {
       std::cout << "JJ WARNING: not space in del-rows!!!" << std::endl;
       return 1;
   }
   for (i=begin;i<=end;i++) mindex[i-begin] = i;
   status = XPRSdelcols(lp->prob, ncols_, mindex);
   if (status)
       std::cout << "Error XPRESS JJdelcols() " << status << std::endl;
   free(mindex);
   mindex=NULL;
   return status;
#endif
// End adding PWOF

   
#ifdef VSCIP
    //std::cout << "\nJJdelcols" << std::endl;
    SCIPlpiDelCols(lp,begin,end);
    return 0;
#endif
}

int JJdelsetcols( 
JJLPptr lp,
int    *delstat
)
{
#ifdef CPLEX3
   return delsetcols (lp, delstat);
#endif


#ifdef CPLEX5
   return CPXdelsetcols (CPLEXv::Env, lp, delstat);
#endif


#ifdef XPRESS
   int ncols,status,mac;
   int *mindex;

   JJcheck(lp);
   ncols = 0;
   getipv (N_NCOL, &mac);
   while(mac--)
       if( delstat[mac] ) ncols++;
   if( ncols==0 ) return 0;
   mindex = (int *)malloc( ncols*sizeof(int) );
   if( mindex==NULL ){
       fprintf(stderr,"JJ WARNING: not space for del-columns!!!\n");
       return 1;
   }
   ncols = 0;
   getipv (N_NCOL, &mac);
   while(mac--)
       if( delstat[mac] ) mindex[ncols++] = mac;
   status = delcols(ncols, mindex);
   free( mindex );
   mindex = NULL; /*PWOF*/
   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int ncols_, status, mac;
   int *mindex = NULL;
   
   ncols_ = 0;
   status = XPRSgetintattrib(lp->prob, XPRS_COLS, &mac);
   if (status) 
   {
       std::cout << "Error XPRESS JJdelsetcols() " << status << std::endl;
       return status;
   }
   
   while (mac--)
       if (delstat[mac]) ncols_++;
   if (ncols_==0) return 0;
   
   mindex = (int *)malloc( ncols_*sizeof(int) );
   if (mindex==NULL)
   {
       std::cout << "JJ WARNING: not space for del-set-columns!!!" << std::endl;
       return 1;
   }
   ncols_=0;
   status = XPRSgetintattrib(lp->prob, XPRS_COLS, &mac);
   if (status) 
   {
       std::cout << "Error XPRESS JJdelsetcols() " << status << std::endl;
       free(mindex);
       mindex=NULL;
       return status;
   }
   while (mac--)
       if (delstat[mac]) mindex[ncols_++] = mac;
   
   status = XPRSdelcols(lp->prob, ncols_, mindex);
   
   if (status)
       std::cout << "Error XPRESS JJdelsetcols() " << status << std::endl;
   free(mindex);
   mindex=NULL;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJdelsetcols" << std::endl;
    SCIPlpiDelColset(lp,delstat);
    return 0;
#endif
}

int JJchgcoef (
JJLPptr lp,
int i,int j,
double newvalue
)
{
#ifdef CPLEX3
   return chgcoef (lp, i, j, newvalue);
#endif


#ifdef CPLEX5
   return CPXchgcoef (CPLEXv::Env, lp, i, j, newvalue);
#endif


#ifdef XPRESS
   JJcheck(lp);
   if( i == -1 )
       return chgobj(1, &j, &newvalue);
   if( j == -1 )
       return chgrhs(1, &i, &newvalue);
   if( j == -2 )
       return chgrng(1, &i, &newvalue);
   return chgcof(i, j, newvalue);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   
   if (i == -1)
       status = XPRSchgobj(lp->prob,1,&j,&newvalue);
   else if (j == -1)
       status = XPRSchgrhs(lp->prob, 1, &i, &newvalue);
   else if (j == -2)
       status = XPRSchgrhsrange(lp->prob, 1, &i, &newvalue);
   else
       status = XPRSchgcoef(lp->prob, i, j, newvalue);
   
   if (status) 
       std::cout << "Error XPRESS JJchgcoef() " << status << std::endl;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJchgcoef" << std::endl;
    int ind;

    if( i == -1 ) {        
        ind = j;
        SCIP_Real opc;
        opc = newvalue;
        SCIPlpiChgObj(lp,1,&ind,&opc);

    }else {
      if( j == -1 ){        
        ind = i;
        SCIP_Real lhs;
        SCIP_Real rhs;
        SCIP_Real rowlower;
        SCIP_Real rowupper;
        SCIPlpiGetSides(lp,i,i,&rowlower,&rowupper);
        if (SCIPlpiIsInfinity(lp,-rowlower)){
            lhs = -SCIPlpiInfinity(lp);
            rhs = newvalue;
            SCIPlpiChgSides(lp,1,&ind,&lhs,&rhs);
        }else
            if (SCIPlpiIsInfinity(lp,rowupper)){
                lhs = newvalue;
                rhs = SCIPlpiInfinity(lp);
                SCIPlpiChgSides(lp,1,&ind,&lhs,&rhs);
            }else {//if (sense[i] == 'E') {                
                lhs = newvalue;
                rhs = newvalue;
                SCIPlpiChgSides(lp,1,&ind,&lhs,&rhs);
            }
      }else {
        if( j == -2 ) {
            system("pause"); return 1;
        }else
        {            
            SCIPlpiChgCoef(lp,i,j,newvalue);
        }
      }
    }

    return 0;
#endif
}

int JJchgbds (
JJLPptr lp,
int cnt,
int *index,
char *lu,
double *bd
)
{
#ifdef CPLEX3
   return chgbds (lp, cnt, index, lu, bd);
#endif


#ifdef CPLEX5
   return CPXchgbds (CPLEXv::Env, lp, cnt, index, lu, bd);
#endif


#ifdef XPRESS
   JJcheck(lp);
   return chgbds (cnt, index, lu, bd);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   
   status = XPRSchgbounds(lp->prob, cnt, index, lu, bd);
   
   if (status) 
       std::cout << "Error XPRESS JJchgbds() " << status << std::endl;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJchgbds" << std::endl;
    SCIP_Real lb;
    SCIP_Real ub;
    SCIP_Real lbAux;
    SCIP_Real ubAux;    
    int indAux;
    for (int i = 0; i < cnt; i++)
    {
        SCIPlpiGetBounds(lp,index[i],index[i],&lbAux,&ubAux);
        indAux = index[i];
        if (lu[i] == 'L') {
            lb = bd[i];
            ub = ubAux;            
            SCIPlpiChgBounds(lp,1,&indAux,&lb,&ub);
        }else
            if (lu[i] == 'U') {
                lb = lbAux;
                ub = bd[i];                
                SCIPlpiChgBounds(lp,1,&indAux,&lb,&ub);
            }else {
                lb = bd[i];
                ub = bd[i];                
                SCIPlpiChgBounds(lp,1,&indAux,&lb,&ub);
            }        
    }
    return 0;
#endif
}

int JJlpwrite( 
JJLPptr lp,
//char *filename
std::string filename
)
{
#ifdef CPLEX3
   return lpwrite (lp, filename.c_str());
#endif

#ifdef CPLEX5
   char *filenameCplex = new char[filename.length()+1];
   std::strcpy(filenameCplex,filename.c_str());
   int ret = CPXlpwrite (CPLEXv::Env, lp, filenameCplex);
   delete[] filenameCplex;
   return ret;
#endif


#ifdef XPRESS
   JJcheck(lp);
   fprintf(stderr,"JJ WARNING: not implemented!!!" << std::endl;
   return 1;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   
   status = XPRSwriteprob(lp->prob, filename.c_str(), "lops");
   
   if (status) 
       std::cout << "Error XPRESS JJlpwrite() " << status << std::endl;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJlpwrite\n");
    SCIPlpiWriteLP(lp,filename.c_str());
    return 0;
#endif
}

int JJmpswrite( 
JJLPptr lp,
char *filename
)
{
#ifdef CPLEX3
   return mpswrite (lp, filename);
#endif


#ifdef CPLEX5
   return CPXmpswrite (CPLEXv::Env, lp, filename);
#endif


#ifdef XPRESS
   int status;

#ifdef XPRESS11
   char flist = 'o';
   const char *flaglist = &flist;
   JJcheck(lp);
   status = output ( filename, flaglist );
#else
   JJcheck(lp);
   status = output ( filename);
#endif

   return status;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   
   status = XPRSwriteprob(lp->prob, filename, "ops");
   
   if (status) 
       std::cout << "Error XPRESS JJmpswrite() " << status << std::endl;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJmpswrite" << std::endl;
    SCIPlpiWriteLP(lp,filename);
    return 0;
#endif
}


int JJgetbdl (
JJLPptr  lp,
double    *xlb,
int       begin,
int       end
)
{
#ifdef CPLEX3
   return getbdl (lp, xlb, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetlb (CPLEXv::Env, lp, xlb, begin, end);
#endif


#ifdef XPRESS
   JJcheck(lp);
   return getbdl (xlb, begin, end);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   
   status = XPRSgetlb(lp->prob, xlb, begin, end);
   
   if (status) 
       std::cout << "Error XPRESS JJgetbdl() " << status << std::endl;
   return status;
#endif
// End adding PWOF

   
#ifdef VSCIP
    //std::cout << "\nJJgetbdl" << std::endl;
    double *xub = new double[end-begin+1];
    SCIPlpiGetBounds(lp,begin,end,xlb,xub);
    //PWOF change 23-08-2013
    //delete xub;
    delete[] xub;
    return 0;
#endif
}


int JJgetbdu (
JJLPptr  lp,
double    *xub,
int       begin,
int       end
)
{
#ifdef CPLEX3
   return getbdu (lp, xub, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetub (CPLEXv::Env, lp, xub, begin, end);
#endif


#ifdef XPRESS
   JJcheck(lp);
   return getbdu (xub, begin, end);
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int status;
   
   status = XPRSgetub(lp->prob, xub, begin, end);
   
   if (status) 
       std::cout << "Error XPRESS JJgetbdu() " << status << std::endl;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    //std::cout << "\nJJgetbdu\n");
    double *xlb = new double[end-begin+1];
    SCIPlpiGetBounds(lp,begin,end,xlb,xub);
    //PWOF change 23-08-2013
    //delete xlb;
    delete[] xlb;
    return 0;
#endif
}

int JJloadctype(
JJLPptr  lp,
char     *ctype
)
{
#ifdef CPLEX3
   return loadctype (lp, ctype);
#endif

#ifdef CPLEX5
//   return CPXloadctype (Env, lp, ctype);
//   Modified on December 2006 to link with CPLEX 10
   return CPXcopyctype (CPLEXv::Env, lp, ctype);

#endif

#ifdef XPRESS
   int  i,mac;
   int *mindex;

   JJcheck(lp);
   mac = JJgetmac(lp);
   mindex = (int *)malloc( mac*sizeof(int) );
   if( mindex==NULL ){
       fprintf(stderr,"JJ ERROR: not enough memory in JJloadctype\n");
       return 1;
   }
   for(i=0;i<mac;i++) mindex[i] = i;
   i = chgcoltype(mac, mindex, ctype);
   free( mindex );
   mindex = NULL; /*PWOF*/
   return i;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   int i, mac, status;
   int *mindex = NULL;
   
   status = XPRSgetintattrib(lp->prob, XPRS_COLS, &mac);
   
   if (status) 
   {
       std::cout << "Error XPRESS JJloadctype() " << status << std::endl;
       return status;
   }
   mindex = (int*) malloc(mac*sizeof(int));
   if (mindex == NULL)
   {
       std::cout << "JJ ERROR: not enough memory in JJloadctype" << std::endl;
       return 1;
   }   
   for (i=0;i<mac;i++) mindex[i] = i;
   status = XPRSchgcoltype(lp->prob, mac, mindex, ctype);
   if (status) 
       std::cout << "Error XPRESS JJloadctype() " << status << std::endl;
   free(mindex);
   mindex = NULL;
   return status;
#endif
// End adding PWOF
   
#ifdef VSCIP
    std::cout << "\nJJloadctype?" << std::endl;
    system("pause");
    return 1;
#endif
}

int JJgetmx ( 
JJLPptr lp,
double *x,
int    begin,int end
)
{
#ifdef CPLEX3
   return getmx (lp, x, begin, end);
#endif


#ifdef CPLEX5
   return CPXgetmipx (CPLEXv::Env, lp, x, begin, end);
#endif


#ifdef XPRESS
   int k,status;
   double *ptr;

   k = JJgetmac(lp); 
   ptr = (double *)malloc( k * sizeof( double ) );
   if( ptr==NULL ) return 1;
   status = JJsolution (lp, NULL, NULL, ptr, NULL, NULL, NULL);
   if( status==0 )
       for( k=begin ; k<=end ; k++ )
           x[k-begin] = ptr[k];
   free( ptr );
   ptr = NULL; /*PWOF*/
   return status;
#endif

#ifdef XPRESS_13
   int k,status;
   double *ptr;

   k = JJgetmac(lp); 
   ptr = (double *)malloc( k * sizeof( double ) );
   if( ptr==NULL ) return 1;
   status = JJsolution (lp, NULL, NULL, ptr, NULL, NULL, NULL);
   if( status==0 )
       for( k=begin ; k<=end ; k++ )
           x[k-begin] = ptr[k];
   free( ptr );
   ptr = NULL; /*PWOF*/
   return status;
#endif
   
#ifdef VSCIP
    //std::cout << "\nJJgetmx" << std::endl;
    JJgetx(lp,x,begin,end);
    return 0;
#endif
}

int JJmipoptimize (
JJLPptr lp)
{
#ifdef CPLEX3
   return mipoptimize (lp);
#endif


#ifdef CPLEX5      
   return CPXmipopt (CPLEXv::Env, lp);
#endif

#ifdef XPRESS
   JJcheck(lp);
   return global();
#endif

#ifdef XPRESS_13
   JJoptimize(lp);
   return 0;
#endif
   
#ifdef VSCIP
    //std::cout << "\nJJmipoptimize" << std::endl;
    JJoptimize(lp);
    return 0;
#endif
}

/* Not used
JJLPptr JJloadmprob (
char *probname,
int numcols, int numrows, int numrims, int objsen,
double *obj, double *rhs,
char *sense,
int *matbeg, int *matcnt, int *matind,
double *matval, double *lb, double *ub, double *rngval,
int *freerowind, int *rimtype, int *rimbeg, int *rimcnt, int *rimind,
double *rimval,
char *dataname, char *objname, char *rhsname, char *rngname, char *bndname,
char **colname,
char *colnamestore,
char **rowname,
char *rownamestore,
char **rimname,
char *rimnamestore,
int colspace, int rowspace, int nzspace,
int rimspace, int rimnzspace,
unsigned colnamespace, unsigned rownamespace, unsigned rimnamespace,
char *ctype
)
{
#ifdef CPLEX3
    return loadmprob (probname, numcols, numrows, numrims,
                    objsen, obj, rhs, sense, matbeg, matcnt,
                    matind, matval, lb, ub, rngval,
                    freerowind, rimtype, rimbeg, rimcnt, rimind, rimval,
                    dataname, objname, rhsname, rngname, bndname,
                    colname, colnamestore, rowname, rownamestore,
                    rimname, rimnamestore, colspace, rowspace, nzspace,
                    rimspace, rimnzspace, colnamespace, rownamespace,
                    rimnamespace,ctype);
#endif


#ifdef CPLEX5    
   int  status;
   char errmsg[1024];
   JJLPptr lp;

   if( CPLEXv::Env==NULL ){
       CPLEXv::Env = CPXopenCPLEX (&status);
       if ( CPLEXv::Env == NULL ) {
           CPXgeterrorstring (CPLEXv::Env, status, errmsg);
           fprintf (stderr, "JJ ERROR: not CPLEX environment.\n%s",errmsg);
           return(NULL);
       }
       lpJJ = 0;
   }

  // lp = CPXloadmipprob (Env, probname, numcols, numrows, numrims,
  //                 objsen, obj, rhs, sense, matbeg, matcnt,
  //                  matind, matval, lb, ub, rngval,
  //                  freerowind, rimtype, rimbeg, rimcnt, rimind, rimval,
  //                  dataname, objname, rhsname, rngname, bndname,
  //                  colname, colnamestore, rowname, rownamestore,
  //                  rimname, rimnamestore, colspace, rowspace, nzspace,
  //                  rimspace, rimnzspace, colnamespace, rownamespace,
  //                  rimnamespace,ctype);

//  Modified on December 2006 to link with CPLEX 10

    lp = CPXcreateprob (CPLEXv::Env, &status, probname);
    status = CPXcopylp (CPLEXv::Env, lp, numcols, numrows, objsen, obj, rhs,
                     sense, matbeg, matcnt, matind, matval, lb, ub, rngval);
    status = CPXcopyctype (CPLEXv::Env, lp, ctype);
    
  
  
   if( lp ) lpJJ++;
   return lp;
#endif

#ifdef XPRESS
    int k,total,ngents;
    int *mgcols;

    if( lpJJ == 0 ){
        if( initlz(XOSLDIR,XOSLMEM) ){
            fprintf (stderr, "JJ ERROR: XPRESS not opened.\n");
            return NULL;
        }
    }
    if( current_lp != NULL ){
        if( savmat( &(current_lp->ref) ) ){
            std::cout << "JJ WARNING: not savmat()" << std::endl;
            return NULL;
        }
    }
    current_lp = (JJLPptr) malloc(sizeof(struct JJLP));
    if( current_lp == NULL ){
        std::cout << "JJ WARNING: not memory" << std::endl;
        return NULL;
    }
    current_lp->ref    = 0;
    current_lp->objsen = objsen;

    seticv(N_NCXTRA,colspace-numcols);
    seticv(N_NRXTRA,rowspace-numrows);
    for(total=k=0;k<numcols;k++) total += matcnt[k];
    seticv(N_NMXTRA,nzspace-total);

    mgcols = (int *)malloc( numcols * sizeof(int) );
    if( mgcols==NULL ){
        free( current_lp );
        return NULL;
    }
    for(ngents=k=0;k<numcols;k++)
        if( ctype[k]=='B' || ctype[k]=='I')
            mgcols[ngents++] = k;

    //seticv(N_NGXTRA,ngents);

    if( loadglobal(probname,numcols,numrows,sense,rhs,NULL,
                  obj,matbeg,matcnt,matind,matval,lb,ub,
                  ngents,0,ctype,mgcols,NULL,NULL,NULL,NULL,NULL) ){
        free( current_lp );
        current_lp = NULL;
    } else
        lpJJ++;
    free( mgcols );
    return current_lp;
#endif

#ifdef XPRESS_13
    int k,total,ngents;
    int *mgcols;

    if( lpJJ == 0 ){
        if( initlz(XOSLDIR,XOSLMEM) ){
            fprintf (stderr, "JJ ERROR: XPRESS not opened.\n");
            return NULL;
        }
    }
    if( current_lp != NULL ){
        if( savmat( &(current_lp->ref) ) ){
            std::cout << "JJ WARNING: not savmat()" << std::endl;
            return NULL;
        }
    }
    current_lp = (JJLPptr) malloc(sizeof(struct JJLP));
    if( current_lp == NULL ){
        std::cout << "JJ WARNING: not memory" << std::endl;
        return NULL;
    }
    current_lp->ref    = 0;
    current_lp->objsen = objsen;

    seticv(N_NCXTRA,colspace-numcols);
    seticv(N_NRXTRA,rowspace-numrows);
    for(total=k=0;k<numcols;k++) total += matcnt[k];
    seticv(N_NMXTRA,nzspace-total);

    mgcols = (int *)malloc( numcols * sizeof(int) );
    if( mgcols==NULL ){
        free( current_lp );
        return NULL;
    }
    for(ngents=k=0;k<numcols;k++)
        if( ctype[k]=='B' || ctype[k]=='I')
            mgcols[ngents++] = k;

    //seticv(N_NGXTRA,ngents);

    if( loadglobal(probname,numcols,numrows,sense,rhs,NULL,
                  obj,matbeg,matcnt,matind,matval,lb,ub,
                  ngents,0,ctype,mgcols,NULL,NULL,NULL,NULL,NULL) ){
        free( current_lp );
        current_lp = NULL;
    } else
        lpJJ++;
    free( mgcols );
    return current_lp;
#endif

    
#ifdef VSCIP
    //std::cout << "\nJJloadmprob" << std::endl;
    return JJloadprob (probname,numcols,numrows,numrims,objsen,obj,rhs,sense,matbeg,matcnt,matind,
          matval,lb,ub,rngval,freerowind,rimtype,rimbeg,rimcnt,rimind,rimval,dataname,objname,rhsname,
          rngname,bndname,colname,colnamestore,rowname,rownamestore,rimname,rimnamestore,colspace,rowspace,
          nzspace,rimspace,rimnzspace,colnamespace,rownamespace,rimnamespace);
#endif
}
*/

int JJbinvrow (
JJLPptr   lp,
int       i,
double    *y
)
{
#ifdef CPLEX3
    return binvrow (lp, i, y);
#endif

#ifdef CPLEX5
    return CPXbinvrow (CPLEXv::Env, lp, i, y);
#endif

#ifdef XPRESS
    std::cout << "JJ ERROR: binvrow not implemented" << std::endl;
    return 1;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   std::cout << "JJ ERROR: binvrow not implemented" << std::endl;
   return 1;
#endif
// End adding PWOF
    
#ifdef VSCIP
    //std::cout << "\nJJbinvrow" << std::endl;
    SCIPlpiGetBInvRow(lp,i,y, NULL, NULL);
    return 0;
#endif
}

int JJbinvarow (
JJLPptr   lp,
int       i,
double    *z
)
{
#ifdef CPLEX3
    return binvarow (lp, i, z);
#endif

#ifdef CPLEX5
    return CPXbinvarow (CPLEXv::Env, lp, i, z);
#endif

#ifdef XPRESS
    std::cout << "JJ ERROR: binvarow not implemented" << std::endl;
    return 1;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   std::cout << "JJ ERROR: binvarow not implemented" << std::endl;
   return 1;
#endif
// End adding PWOF

    
#ifdef VSCIP
    //std::cout << "\nJJbinvarow" << std::endl;
    double *binvrow = new SCIP_Real[JJgetmar(lp)];
    SCIPlpiGetBInvRow(lp,i,binvrow, NULL, NULL);
    SCIPlpiGetBInvARow(lp,i,binvrow,z, NULL, NULL);
    //PWOF change 23-08-2013
    //delete binvrow;
    delete[] binvrow;
    return 0;
#endif
}

int JJgetbhead (
JJLPptr   lp,
int       *head,
double    *x
)
{
#ifdef CPLEX3
    return getbhead (lp, head, x);
#endif

#ifdef CPLEX5
    return CPXgetbhead (CPLEXv::Env, lp, head, x);
#endif

#ifdef XPRESS
    std::cout << "JJ ERROR: getbhead not implemented" << std::endl;
    return 1;
#endif

// Start adding PWOF, taken from RAFA's code
#ifdef XPRESS_13
   std::cout << "JJ ERROR: getbhead not implemented" << std::endl;
   return 1;
#endif
// End adding PWOF

    
#ifdef VSCIP
    std::cout << "\nJJgetbhead?" << std::endl;
    system("pause");
    return 1;
#endif
}
