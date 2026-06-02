/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPGOMO.c                                        */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cspdefns.h"
#include "cspglob2.h"
#include "cspgomo.h"
#include "cspcover.h"
#include "cspsep.h"
#include "jjsolver.h"
#include "cspback.h"



#ifdef GOMORY


/* -------------------------- GOMORY CUT SEPARATION ---------------------- */


#define    RE_GOMORY          1
#define    MAX_CUTS_GOMO     10 
#define    MAX_DEN_GOMORY  3000
#define    MAX_COEF_GOMORY 3000

extern JJLPptr lp;

static int  chvatal(int,double *,CONSTRAINT **,int *,CONSTRAINT **,int,int,int*,double,int);
static int  different_gomory(int,VARIABLE **,int *,int,CONSTRAINT **);
static int  adding(CONSTRAINT *,double,double *,double *);

/*******************************************************************/

int        separa_gomory(int        *card,CONSTRAINT **stack)

{
    int        i,j,k,cont;
    double     viola;
    double     *rb,*ry;
    int        *head,*cstat;
    CONSTRAINT **ri;
    static int set_gomory = 1;

    if( *card ==MAX_CUTS_ITER) return(0);
    if( set_gomory%RE_GOMORY ) {
        set_gomory++;
        return(0);
    }

    JJdualopt(lp);

    if(mar != JJgetmar(lp)){        
        std::cout << "ERROR: mar=" << mar << "  getmar(lp)=" << JJgetmar(lp) << std::endl;
        CSPexit(EXIT_LPSOLVER); //exit(1);
    }
    rb=(double *)malloc( mar*sizeof(double) );
    if(rb==NULL){        
        std::cout << "There is not enough memory for RB" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    ri=(CONSTRAINT **)malloc( mar*sizeof(CONSTRAINT *) );
    if(ri==NULL){        
        std::cout << "There is not enough memory for RI" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    ry=(double *)malloc( mar*sizeof(double) );
    if(ry==NULL){        
        std::cout << "There is not enough memory for RY" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    head=(int *)malloc( mar*sizeof(int) );
    if(head==NULL){        
        std::cout << "There is not enough memory for HEAD" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    cstat=(int *)malloc( mac*sizeof(int) );
    if(cstat==NULL){        
        std::cout << "There is not enough memory for CSTAT" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    if ( JJgetbase(lp,cstat,NULL) ){
        CSPexit(EXIT_LPSOLVER);
        //exit(1);
    }
    
    cont=0;

    JJgetbhead( lp, head, ry);
    for( i=1 ; i<mar ; i++){
        viola = ceil(ry[i]-ZERO) - ry[i];
        if( viola > MIN_VIOLA+ZERO ){
             cont++;
             JJbinvrow( lp, i, rb);
             k=0;
             for( j=1 ; j<mar ; j++ ) {
                  if( rb[j] - floor( rb[j]+ZERO ) > ZERO ){
                      ri[k] = rind[j];
                      rb[k] = rb[j];
                      k++;
                  }
             }
             if( k && chvatal(k,rb,ri,card,stack,GOMORY,MAX_CUTS_GOMO,cstat,viola,i) < 0 )
                  break;
        }
    }

    if( i<mar && (*card)<MAX_CUTS_GOMO ) set_gomory=1;

    free((void *)cstat);
    cstat = NULL; /*PWOF*/
    free((void *)head);
    head = NULL; /*PWOF*/
    free((void *)ry);
    ry = NULL; /*PWOF*/
    free((void *)rb);
    rb = NULL; /*PWOF*/
    free((void *)ri);
    ri = NULL; /*PWOF*/
    cgomo += *card;
    return( *card );
}


static int chvatal(int k,double *rb,CONSTRAINT **con,int *card,CONSTRAINT **stack,int type,int max_local_cuts,int *cstat,double viola,int number)

{
    int        j,v,coef,maxcoef;
    double     *ra;
    int        *rc;
    double     rhs;
    VARIABLE   **ri;
    CONSTRAINT *row;
    double     vio;

    ra=(double *)calloc( (size_t)ncols,sizeof(double) );
    if(ra==NULL){        
        std::cout << "There is not enough memory for RA" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    rhs = 0.0;
    for( j=0 ; j<k ; j++ )
        if( adding( con[j] , rb[j] , ra , &rhs ) ){
            free( (void *)ra );
            ra = NULL; /*PWOF*/
            return(0);
        }

    maxcoef = 0;
    v = 0;
    for( j=0 ; j<ncols ; j++ ){
        switch( columns[j].stat ){
            case FIX_UB:
                coef = (int)floor(ra[j]+ZERO);
                break;
            case FIX_LB:
                coef = (int)ceil(ra[j]-ZERO);
                break;
            default:
                if( cstat[columns[j].lp] == 2 )
                    coef = (int)floor(ra[j]+ZERO);
                else
                    coef = (int)ceil(ra[j]-ZERO);
        }
        if( coef ) v++;
        if( coef>maxcoef ) maxcoef = coef;
    }
    if( v==0 ){
        free( (void *)ra );
        ra = NULL; /*PWOF*/
        return(0);
    }
    if( v>MAX_DEN_GOMORY || maxcoef>MAX_COEF_GOMORY ){
        return(-2);
    }


    rc=(int *)calloc( (size_t)v,sizeof(int) );
    if(rc==NULL){        
        std::cout << "There is not enough memory for RC" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    ri=(VARIABLE **)malloc( (size_t)v*sizeof(VARIABLE *) );
    if(ri==NULL){
        std::cout << "There is not enough memory for RI" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    v = 0;
    for( j=0 ; j<ncols ; j++ )
        switch( columns[j].stat ){
            case FIX_UB:
                coef = (int)floor(ra[j]+ZERO);
                if( coef ){
                    ri[v] = &columns[j];
                    rc[v] = coef;
                    v++;
                }
                rhs -= ra[j]-coef;
                break;
            case FIX_LB:
                coef = (int)ceil(ra[j]-ZERO);
                if( coef ){
                    ri[v] = &columns[j];
                    rc[v] = coef;
                    v++;
                }
                break;
            default:
                if( cstat[columns[j].lp] == 2 ){
                    coef = (int)floor(ra[j]+ZERO);
                    if( coef ){
                        ri[v] = &columns[j];
                        rc[v] = coef;
                        v++;
                    }
                    rhs -= ra[j]-coef;
                }else{
                    coef = (int)ceil(ra[j]-ZERO);
                    if( coef ){
                        ri[v] = &columns[j];
                        rc[v] = coef;
                        v++;
                    }
                }
        }
    rhs = ceil( rhs-ZERO );
    free( (void *)ra );
    ra = NULL; /*PWOF*/


    if( type==GOMORY && different_gomory(v,ri,rc,*card,stack) ){
        free( (void *)rc );
        rc = NULL; /*PWOF*/
        free( (void *)ri );
        ri = NULL; /*PWOF*/
        return(0);
    }
    row = (CONSTRAINT *)malloc( sizeof(CONSTRAINT) );
    if(row==NULL){        
        std::cout << "There is not enough memory for ROW" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    row->stack = (VARIABLE **)malloc( v * sizeof(VARIABLE *) );
    if(row->stack==NULL){        
        std::cout << "There is not enough memory for ROW->stack " << v << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    row->coef = (double *)malloc( v * sizeof(double) );
    if(row->coef==NULL){        
        std::cout << "There is not enough memory for ROW->coef" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    row->type  = type;
    row->stat  = LP_BA;
    row->lp    = -1;
    row->card  = v;
    row->con   = NULL;
    row->rhs   = rhs;
    row->sense = 'G';
    for( j=0 ; j<v ; j++){
            row->stack[j] = ri[j];
            row->coef[j]  = rc[j];
    }
    free( (void *)rc );
    rc = NULL; /*PWOF*/
    free( (void *)ri );
    ri = NULL; /*PWOF*/

    vio = violated( row );
    if( fabs( vio - viola ) > ZERO ){
        free( (void *)(row->stack) );
        row->stack = NULL; /*PWOF*/
        free( (void *)(row->coef) );
        row->coef = NULL; /*PWOF*/
        free( (void *)row );
        row = NULL; /*PWOF*/
        return(0);
    }
    if( vio<MIN_VIOLA ){
        free( (void *)(row->stack) );
        row->stack = NULL; /*PWOF*/
        free( (void *)(row->coef) );
        row->coef = NULL; /*PWOF*/
        free( (void *)row );
        row = NULL; /*PWOF*/
        return(0);
    }

    if( new_row( row ) ){
        stack[ (*card)++ ] = row;
        if( (*card) == max_local_cuts || (*card)==MAX_CUTS_ITER ){
            std::cout << "WARNING: Too much local rows (" << max_local_cuts << ")" << std::endl;
            return(-1);
        }
    }
    return(v);
}



static int different_gomory(int v,VARIABLE    **ri,int         *rc,int card,CONSTRAINT  **stack)

{
    int      l;
    double   *ptr_coe;
    VARIABLE **ptr_var;

    while(card--){
        if( (l=stack[card]->card) == v ){
            ptr_var = stack[card]->stack;
            ptr_coe = stack[card]->coef;
            while(l-- && ri[l]==ptr_var[l] && fabs(rc[l]-ptr_coe[l])<ZERO );
            if(l<0)return(1);
        }
    }
    return(0);
}



static      int adding(CONSTRAINT  *con,double rb,double      *ra,double      *rhs)

{
    int       card;
    VARIABLE  **stack;
    double    *coef;
    

    if( con->sense != 'G' ){
         return(1);
    }

    switch( con->type ){
    case COVER:
        stack = (VARIABLE **)malloc( ncols*sizeof(VARIABLE *) ); /*PWOF*/
        coef  = (double *)malloc( ncols*sizeof(double) );
        card  = extend_cover(con,stack,coef);
        break;
    case CAPACITY:
    case BRIDGE:
    case GOMORY:
        stack = con->stack;
        coef  = con->coef;
        card  = con->card;
        break;
    default:
        return(1);
    }


    rb = rb - floor(rb+ZERO);

    if( rb>ZERO ){
        while(card--) ra[ stack[card]->index ] += coef[card] * rb;
        *rhs += con->rhs * rb;
    }


    switch( con->type ){
    case COVER:
        free(stack);
        stack = NULL; /*PWOF*/
        free(coef);
        coef = NULL; /*PWOF*/
        break;
    case BRIDGE:
    case CAPACITY:
    case GOMORY:
        break;
    default:
        return(1);
    }

    return 0;
}

double     violation_gomory(CONSTRAINT *con,int n,VARIABLE   **var)

{
    int       card,i;
    double    viola;
    VARIABLE  **stack;
    double    *coef;
    VARIABLE  *ptr;

    viola = con->rhs;
    coef  = con->coef;
    stack = con->stack;
    card  = con->card;

    while(n--){
        ptr = var[n];
        for(i=0;i<card;i++)
            if( stack[i]==ptr ){
                viola -= coef[i] * ptr->val;
                break;
            }
    }
    return( viola );
}

    
void control_gomory(int  i ,double * ra , double rhs )

{
    int j;
    double *z;
    z = (double *)malloc( mac * sizeof(double) );

    JJbinvarow(lp,i,z);
    for(j=0;j<mac;j++)
        if( fabs( z[j] - ra[j] )>ZERO )
            std::cout << "ERROR in Gomory: row=" << i << " col=" << j << " z=" << z[j] << " ra=" << ra[j] << std::endl;

    free(z); 
    z = NULL; /*PWOF*/
}

double     get_coeficient_gomory(VARIABLE   *col,CONSTRAINT *con)

{
    int i;
    for(i=0;i<con->card;i++)
        if( col==con->stack[i] )
            return( con->coef[i] );
    return(0);
}

#endif
