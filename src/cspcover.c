/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPCOVER.c                                       */
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
#include "cspcover.h"
#include "mt1rc.h"
#include "cspsep.h"
#include "cspback.h"
#include "cspcapa.h"


static  CONSTRAINT *cover_constraint(CONSTRAINT *);
static  int        cover_extension(int*,VARIABLE**,CONSTRAINT*,double*,double*);
static  double     kp(int,double *,double *,double,int *);

typedef struct{
        VARIABLE *var;
        double   coef;
        double   val;
}ITEM;

/* -----------------------  COVER SEPARATION --------------------------*/

int        separa_cover(int *card,CONSTRAINT **stack)

{
    int        i;
    int        num=0;
    CONSTRAINT *con;
    
    if( *card ==MAX_CUTS_ITER) return(0);
    i = nrows;
    while(i--)
        if(rows[i]->type == CAPACITY && rows[i]->rhs > 1+ZERO){
            con = cover_constraint(rows[i]);
            if( con && new_row(con) ){
                stack[ (*card)++] = con;
                num++;
                if( (*card)==MAX_CUTS_ITER) break;
            }else
                remove_row( con );
        }
    ccove += num;
    return(num);
}

static CONSTRAINT *cover_constraint(CONSTRAINT *con)

{
    int        j,k,nitem,card=0,new_card,ext;
    double     profit,profit0;
    double     rhs,viola,viola0,val;
    double     *coef = NULL;
    int        *x = NULL;
    double     *p = NULL,*w = NULL;
    VARIABLE   **stack = NULL,**new_stack = NULL;
    VARIABLE   *var = NULL;
    CONSTRAINT *inequality = NULL;
    ITEM       *item = NULL;
    int        sort_item(const void*,const void*);

    item       = (ITEM *)malloc( ncols*sizeof(ITEM) );
    //new_stack  = (VARIABLE **)malloc( ncols*sizeof(ITEM) );
    new_stack  = (VARIABLE **)malloc( ncols*sizeof(VARIABLE *) ); /*PWOF*/

    switch( con->type ){
    case CAPACITY:
        k = con->card;
        coef  = (double *)malloc( k*sizeof(double) );
        stack = (VARIABLE **)malloc( k*sizeof(VARIABLE *) );
        card = 0;
        while(k--)
            if( con->stack[k]->val > ZERO ){
                stack[card] = con->stack[k];
                coef[card]  = con->coef[k];
                card++;
            }
        break;
    default:        
        std::cout << "ERROR cover: unknown type " << con->type << " of constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    profit  = 0;
    profit0 = con->rhs;
    viola   = 0;
    viola0  = -1 + MIN_VIOLA + ZERO;
    rhs     = con->rhs -1;
    new_card=nitem=0;
    for(j=0;j<card;j++){
        var = stack[j];
        val = var->val;
        if( val > 1.0-ZERO ){
             viola0 += val;
             rhs    -= coef[j];
             new_stack[ new_card++ ] = var;
             profit += coef[j];
             viola  += val;
        } else if( val > ZERO ){
             viola0 += val;
             item[ nitem ].var  = var;
             item[ nitem ].coef = coef[j];
             item[ nitem ].val  = val;
             nitem++;
        } else {             
             std::cout << " WARNING: very rare in KP" << std::endl;
             CSPexit(EXIT_ERROR); //exit(1);
        }
    }
    switch( con->type ){
    case CAPACITY:
        free(coef);
        coef = NULL; /*PWOF*/
        free(stack);
        stack = NULL; /*PWOF*/
        break;
    default:        
        std::cout << "ERROR cover: unknown type " << con->type << " of constraint\n" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }


    for(j=0;j<nitem;j++)
        if( item[j].coef > rhs+ZERO ){
             nitem--;
             item[j].var  = item[nitem].var;
             item[j].coef = item[nitem].coef;
             item[j].val  = item[nitem].val;
             j--;
        }
                 
    if(rhs<ZERO && nitem){        
        std::cout << "ERROR: not empty knapsack in SDCCOVER!" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    if(nitem>1){
        //qsort( (char *)item , nitem , sizeof(ITEM) , sort_item );
        qsort( (void *)item , nitem , sizeof(ITEM) , sort_item );

        x =(int *)malloc( (nitem+2)*sizeof(int) );
        if (x==NULL){            
            std::cout << "Not enough memory for EXACT COVER" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }
        p =(double *)malloc( (nitem+2)*sizeof(double) );
        if (p==NULL){            
            std::cout << "Not enough memory for EXACT COVER" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }
        w =(double *)malloc( (nitem+2)*sizeof(double) );
        if (w==NULL){            
            std::cout << "Not enough memory for EXACT COVER" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }

        for(j=0;j<nitem;j++){
            p[j] = item[j].val;
            w[j] = item[j].coef;
            x[j] = 0;
            if(w[j]>rhs+ZERO || w[j]<-ZERO){                
                std::cout << "WARNING: item " << j << "(p=" << p[j] << ",w=" << w[j] << ") for KP with rhs=" << rhs << std::endl;
                CSPexit(EXIT_ERROR); //exit(1);
            }
        }

        kp(nitem,p,w,rhs,x);
        
        for(j=0;j<nitem;j++)
            if( x[j] ){
                profit += item[j].coef;
                viola  += item[j].val;
                new_stack[ new_card++ ] = item[j].var;
            }

        //free((void *)p);
        free(p);
        p = NULL; /*PWOF*/
        //free((void *)w);
        free(w);
        w = NULL; /*PWOF*/
        //free((void *)x);
        free(x);
        x = NULL; /*PWOF*/
    } else if(nitem==1){
        profit += item[0].coef;
        viola  += item[0].val;
        new_stack[ new_card++ ] = item[0].var;
    }

    free(item);
    item = NULL; /*PWOF*/

    if ( profit < profit0-ZERO ){
        ext = cover_extension( &new_card , new_stack , con , &viola , &viola0 );
        if( viola > viola0 ){
            stack = (VARIABLE **)malloc( new_card * sizeof(VARIABLE *) ); /*PWOF*/
            for(j=0;j<new_card;j++)
                stack[j] = new_stack[j];
            inequality=(CONSTRAINT *)malloc( sizeof(CONSTRAINT) );
            inequality->rhs    = 1+ext;
            inequality->index  = -1;
            inequality->card   = new_card;
            inequality->stack  = stack;
            inequality->coef   = NULL;
            inequality->sense  = 'G';
            inequality->type   = COVER;
            inequality->stat   = LP_BA;
            inequality->lp     = -1;
            inequality->con    = con;
        }
    }
    free(new_stack);
    new_stack = NULL; /*PWOF*/

    return(inequality);
}


static double kp(int N,double *P,double *W,double C,int    *X)

{
    int   j,test;
    int   dim;
    int   *INT;
    double z;
    double *FLOAT;

    dim=N+2;
    INT=(int *)malloc(2*dim*sizeof(int));
    if (INT==NULL){        
        std::cout << "Not enough memory for KP" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    FLOAT=(double *)malloc(5*dim*sizeof(double));
    if (FLOAT==NULL){        
        std::cout << "Not enough memory for KP" << std::endl;
        CSPexit(EXIT_MEMO); //exit(1);
    }
    for(j=0;j<N;j++){
        if( W[j]>C ) W[j]=C;
        else if( W[j]<0.0 ) W[j]=0.0;
    }  
    test = 0;
    MT1RC(N,P-1,W-1,C,0.0001,&z,X-1,dim,test,
        INT,FLOAT,FLOAT+dim,FLOAT+(2*dim),INT+dim,FLOAT+(3*dim),FLOAT+(4*dim));
    free((void *)INT);
    INT = NULL; /*PWOF*/
    free((void *)FLOAT);
    FLOAT = NULL; /*PWOF*/
    return(z);
}



int sort_item(const void * i, const void *j/*ITEM *i,ITEM *j*/)

{
  double vi,vj;

  vi = ((ITEM *)i)->val / ((ITEM *)i)->coef;
  vj = ((ITEM *)j)->val/ ((ITEM *)j)->coef;

  if(vi<vj) return(1);
  if(vi>vj) return(-1);
  return(0);
}



static int  cover_extension(int * new_card ,VARIABLE    ** new_stack ,CONSTRAINT  * con ,double      * viola ,double      * viola0 )

{
    int      k,l,card=0,ext;
    double   maxi;
    VARIABLE **stack=NULL; // Initialized to NULL-pointer
    VARIABLE *var;
    double   *coef=NULL,*new_coef=NULL; // Initialized to NULL-pointer
    
    new_coef  = (double *)malloc( (*new_card)*sizeof(double) );

    switch( con->type ){
    case CAPACITY:
        card  = con->card;
        coef  = (double *)malloc( card*sizeof(double) );
        stack = (VARIABLE **)malloc( card*sizeof(VARIABLE *) );
        for(k=0;k<card;k++) stack[k] = con->stack[k];
        for(k=0;k<card;k++) coef[k]  = con->coef[k];
        break;
    default:        
        std::cout << "ERROR cover: unknown type " << con->type << " of constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }

    l  = *new_card;
    while(l--){
        var = new_stack[l];
        k   = card;
        while( k-- )
            if( var == stack[k] ){
                new_coef[l] = coef[k];
                --card;
                stack[k] = stack[card];
                coef[k]  = coef[card];
                break;
            }
    }

    maxi = 0;
    l = card;
    while(l--)
        if( coef[l] > maxi ) maxi=coef[l];

    ext = 0;
    l = *new_card;
    while(l--)
        if( new_coef[l] >= maxi && new_stack[l]->stat!=FIX_UB ){
            ext++;
            *viola -= new_stack[l]->val;
            (*viola0) -= 1;
            --(*new_card);
            new_stack[l] = new_stack[ *new_card ];
            new_coef[l]  = new_coef[ *new_card ];
        }
    free( new_coef );
    new_coef = NULL; /*PWOF*/
    switch( con->type ){
    case CAPACITY:
        free(coef);
        coef = NULL; /*PWOF*/
        free(stack);
        stack = NULL; /*PWOF*/
        break;
    default:        
        std::cout << "ERROR cover: unknown type " << con->type << " of constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    return( ext );
}



/************************/

int        extend_cover(CONSTRAINT *con,VARIABLE   **stack,double     *coef)

{
   int      k,card,num;
   VARIABLE **ptr;
   VARIABLE *var;

   card  = con->con->card;
   ptr   = con->con->stack;
   k = card;
   while(k--)
       stack[k] = ptr[k];

   num  = con->card;
   ptr  = con->stack;
   while(num--){
       var = ptr[num];
       k   = card;
       while( k-- )
           if( var == stack[k] ){
               stack[k] = stack[--card];
               break;
           }
   }
   k = card;
   while(k--)
       coef[k] = 1;
   return(card);
}

double     violation_cover(CONSTRAINT *con)

{
    double    viola;
    VARIABLE  **stack;
    int       card;

    if( con->type != COVER ){        
        std::cout << " ERROR: not cover constraint" << std::endl;
        CSPexit(EXIT_ERROR); //exit(1);
    }
    viola = con->rhs;
    stack = con->con->stack;
    card  = con->con->card;
    while(card--)
        viola -= stack[card]->val;
    stack = con->stack;
    card  = con->card;
    while(card--)
        viola += stack[card]->val;
    return( viola );
}


double     get_coeficient_cover(VARIABLE   *col,CONSTRAINT *con)

{
    int i;

    if( get_coeficient_capacity(col,con->con)<ZERO ) return(0);
    for(i=0;i<con->card;i++)
        if( col==con->stack[i] ) return(0);
    return(1);
}
