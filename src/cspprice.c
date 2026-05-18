/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPPRICE.c                                       */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "cspdefns.h"
#include "cspglob2.h"
#include "cspprice.h"
#include "cspsolve.h"
#include "cspcapa.h"
#include "cspbridg.h"
#include "cspcover.h"
#include "cspgomo.h"
#include "cspbranc.h"
#include "cspback.h"


#ifdef STAMP
static void control_pricing(void);
#endif

/*  PRIVATE GLOBAL VARIABLES */

float   tpricing = 0.0;
int     npricing = 0;     /* number of WAITING variables */

#define MAX_ADD  50
#define MAX_DEL  100000




int pricing(int cleaning)

{
    int      i,card,card2;
    VARIABLE *col;
    VARIABLE **stack,**stack2;
    double   redcost,lb;
    struct   PRICE *p,*anterior,*posterior;
    double   *u;
    float    t1;

    anterior = (struct PRICE*) malloc(sizeof(struct PRICE)); // PWOF added
    
    if( list_pricing==NULL ) return(0);
/*
    if( mac+npricing < MAX_COLS_LP ){
        add_cols_all();
        return(1);
    }
*/
    t1 = seconds();
#ifdef STAMP
    for( p=list_pricing,card=0 ; p ; p=p->next,card++ )
        if( p->col->stat != WAITING ){
            std::cout << "ERROR: not waiting columns in LIST PRIcing"  << std::endl;
            CSPexit(EXIT_ERROR); //exit(1);
        }
    std::cout << "  << pricing " << card << "(=" << npricing << ") variables>>" << std::endl;
    if( card!=npricing ){
        std::cout << "ERROR: different numbers of pricing variables"  << std::endl ;
/*
        CSPexit(EXIT_ERROR); //exit(1);
*/
    }
    control_pricing();
#endif
    u=(double *)malloc(sizeof(double)*mar);
    if(u==NULL){	
        std::cout << "Not enough memory for U" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    stack=(VARIABLE**)malloc(sizeof(VARIABLE*)*MAX_ADD);
    if(stack==NULL){	
        std::cout << "Not enough memory for STACK" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    stack2=(VARIABLE**)malloc(sizeof(VARIABLE*)*MAX_DEL);
    if(stack2==NULL){	
        std::cout << "Not enough memory for STACK2" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }

    lb = dualcost(u,NULL);

    card = 0;
    card2 = 0;
    p = list_pricing;
    while( p ){
        posterior = p->next;
        col       = p->col;
        redcost   = col->weight;
        for(i=1;i<mar;i++)
            if( u[i]>ZERO )
                redcost -= get_coeficient(col,rind[i]) * u[i];
        if( redcost < -ZERO ){
            if( card == MAX_ADD ) break;
            stack[card++] = col;
            if( p==list_pricing ) list_pricing = posterior;
            else anterior->next = posterior;
            free( p );
            p = NULL; /*PWOF*/
        } else {
            if( ceil( lb+redcost-ZERO ) > upperb && card2<MAX_DEL )
                stack2[card2++] = col;
            anterior = p;
        }
        p = posterior;
    }
    free((char *)u);
    u = NULL; /*PWOF*/

    if(card){
        npricing -= card;
        add_cols(card,stack);
    } else if( cleaning ){
        npricing -= card2;
        rem_cols(card2,stack2);
    }
    free((char *)stack);
    stack = NULL; /*PWOF*/
    free((char *)stack2);
    stack2 = NULL; /*PWOF*/
    tpricing += seconds()-t1;
    return(card);
}



#ifdef STAMP
static void control_pricing()
{
    int      i,k;
    double   *u,*dj;
    double   redcost;
    VARIABLE *col;

    control();
    dj=(double *)malloc(sizeof(double)*mac);
    if(dj==NULL){	
        std::cout << "Not enough memory for DJ" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    u=(double *)malloc(sizeof(double)*mar);
    if(u==NULL){	
        std::cout << "Not enough memory for U" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    for(k=0;k<mar;k++)
        if( u[k]< -ZERO ){	    
            std::cout << "ERROR: negative dual variable in row " << k << std::endl;
            CSPexit(EXIT_ERROR); //exit(1);
        }
    dualcost(u,dj);
    for( k=1;k<mac;k++ ){
        col = cind[k];
        redcost = col->weight;
        for(i=1;i<mar;i++)
            if( u[i]>ZERO )
                redcost -= get_coeficient(col,rind[i]) * u[i];
        if( fabs(dj[col->lp]-redcost) > ZERO ) {	    
            std::cout << "\n ERROR pricing the variable k=" << col->index << " v=" << col->lp << std::endl;
            CSPexit(EXIT_ERROR); //exit(1);
        }
    }
    free((char *)dj);
    dj = NULL; /*PWOF*/
    free((char *)u);
    u = NULL; /*PWOF*/
}
#endif


double     get_coeficient(VARIABLE   *col,CONSTRAINT *con)

{
    switch( con->type ){
    case CAPACITY:
         return( get_coeficient_capacity(col,con) );
    case BRIDGE:
         return( get_coeficient_bridge(col,con) );
    case COVER:
         return( get_coeficient_cover(col,con) );
    case BRANCHCUT:
         return( get_coeficient_branching(col,con) );
    case GOMORY:
         return( get_coeficient_gomory(col,con) );
    default:	
        std::cout << "ERROR (price): unknown type " << con->type << " of constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    return(0);
}


void insert_list_pricing(VARIABLE *col)

{
    struct PRICE *p;

    p = (struct PRICE *)malloc( sizeof( struct PRICE ) );
    if( p==NULL ){	
        std::cout << "ERROR: not memory for LIST PRICING" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    p->col       = col;
    p->next      = list_pricing;
    list_pricing = p;
    npricing ++;
}

        
