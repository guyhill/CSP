/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPCAPA.c                                        */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "cspdefns.h"
#include "cspglob2.h"
#include "cspcapa.h"
#include "cspsep.h"
#include "cspnet.h"
#include "cspback.h"



/* -----------------------  CAPACITY SEPARATION --------------------------*/


int        separa_capacity(int        *card,CONSTRAINT **stack)

{
    
    int       k,l,nvar,num,index;
    double    *primal;
    double    *cvar;
    PROT_LEVEL*pro,*pro0;
    VARIABLE  **xvar;
    CONSTRAINT*con;
    char      *sleep;
    double    bd[2];
    double    maxc;

    if( *card ==MAX_CUTS_ITER) return(0);

    cvar = (double *)malloc( Rncells * sizeof(double) );
    xvar = (VARIABLE **)malloc( Rncells * sizeof(VARIABLE *) );
    primal = (double *)malloc( Rncells * sizeof(double) );
    if( !cvar || !xvar || !primal ){        
        std::cout << "There is not enough memory for CAPACITY" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    for(k=0;k<Rncells;k++) primal[k] = columns[k].val;
    load_network(primal,'C');

    num = *card;
    sleep = (char *)calloc( nprot_level , sizeof(char) );
    l = nprot_level;
    while(l--){
        pro = prot_level+l;
        if( !sleep[l] ){
			index = pro->sen->var->index;
            if( protection_level(pro->sen->var,pro->sense,&nvar,xvar,cvar,primal,'C') < pro->sense * pro->sen->var->nominal + pro->level-ZERO ){
                con = capacity_constraint(nvar,xvar,cvar,pro->sen->var,pro->level);
                if( con && new_row(con) ){
                    stack[ num++] = con;
                    if( num==MAX_CUTS_ITER) break;
                }else
                    remove_row( con );
            }
            k = l;
            while(k--){
                pro0 = prot_level+k;
                if( pro0->sense * primal[ pro0->sen->var->index ] > pro0->sense * pro0->sen->var->nominal + pro0->level - ZERO )
                    sleep[k] = 1;
            }
//                    for introducing an alternative Capacity Cut
//                         a cell must be choosen to be super-suppressed

            index = -1;
            maxc  = 0.0;
            for(k=0;k<nvar;k++)
                if( xvar[k]->val>ZERO && cvar[k]>maxc ){
                    index = xvar[k]->index;
                    maxc  = cvar[k];
                }
            if( index != -1 ){
                free_col(index,bd);
                if( protection_level(pro->sen->var,pro->sense,&nvar,xvar,cvar,NULL,'C') < pro->sense * pro->sen->var->nominal + pro->level-ZERO ){
                    con = capacity_constraint(nvar,xvar,cvar,pro->sen->var,pro->level);
                    if( con && new_row(con) ){
                        stack[ num++] = con;
                        if( num==MAX_CUTS_ITER) break;
                    }else
                        remove_row( con );
                }
                unfree_col(index,bd);
            }
        }
    }
    free( sleep );
    sleep = NULL; /*PWOF*/
    num -= *card;
    *card += num;
    unload_network();

    free(primal);
    primal = NULL; /*PWOF*/
    free(cvar);
    cvar = NULL; /*PWOF*/
    free(xvar);
    xvar = NULL; /*PWOF*/
    ccapa += num;
    return( num );
}




CONSTRAINT *capacity_constraint(int nvar,VARIABLE  **xvar,double    *cvar,VARIABLE  *var,double rhs)

{
        int        k,l;
        CONSTRAINT *row;
        VARIABLE   **x;
        double     *c;
        double     maxlhs;

        for(l=k=0;k<nvar;k++)
            if( cvar[k] > rhs-ZERO ) l++;
        if(l==nvar && rhs>ZERO){
            rhs = 1;
            for(k=0;k<nvar;k++) cvar[k]=1;
        }


        row = (CONSTRAINT *)malloc( sizeof(CONSTRAINT) );
        if(row==NULL){            
            std::cout << "There is not enough memory for ROW" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }

        x = (VARIABLE **)malloc( nvar*sizeof(VARIABLE *) );
        if( x==NULL ){            
            std::cout << " ERROR: not enough memory for the X (cspcapa.c)" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }
        c = (double *)malloc( nvar*sizeof(double) );
        if( c==NULL ){            
            std::cout << " ERROR: not enough memory for the C" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }
        maxlhs=rhs;
        for(k=0;k<nvar;k++)
            if( xvar[k]->stat==FIX_UB && xvar[k]->lp==0 )
                maxlhs -= cvar[k];

        for(k=0;k<nvar;k++){
            x[k] = xvar[k];
            c[k] = cvar[k];

            if( xvar[k]->stat==FIX_UB && xvar[k]->lp==0 )
                c[k] = ( cvar[k]>rhs ? rhs : cvar[k] );
            else   
                c[k] = ( cvar[k]>maxlhs ? maxlhs : cvar[k] );

        }
   
        row->rhs    = rhs;
        row->index  = var->index;
        row->card   = nvar;
        row->stack  = x;
        row->coef   = c;
        row->sense  = 'G';
        row->type   = CAPACITY;
        row->stat   = LP_BA;
        row->lp     = -1;
        row->con    = NULL;
        return(row);
}




double     violation_capacity(CONSTRAINT *con)

{
    int       k;
    double    viola;
    VARIABLE  **xvar;
    double    *cvar;

    if( con->type != CAPACITY ){        
        std::cout << " ERROR: not capacity constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    xvar  = con->stack;
    cvar  = con->coef;
    k     = con->card;
    viola = con->rhs;
    while(k--)
        viola -= cvar[k] * xvar[k]->val;
    return( viola );
}


double     get_coeficient_capacity(VARIABLE   *col,CONSTRAINT *con)

{
    int       k;
    VARIABLE  **xvar;
    double    *cvar;

    xvar = con->stack;
    cvar = con->coef;
    k = con->card;
    while(k--)
        if( xvar[k] == col )
            return( cvar[k] );
    return 0;
}
