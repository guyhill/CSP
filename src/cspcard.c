/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPCARD.c                                        */
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
#include "cspcard.h"
#include "cspback.h"



void    insert_cutcard()
{
        int        rhs,card;
        int        k;
        CONSTRAINT *row;
        VARIABLE   **x;
        double     *c;
        char       *inbetter;
        int        NFREE = nbetter;

        std::cout << " >>>>>>>>>>> neighbour constraint:  NFREE =" << NFREE << std::endl;

        row = (CONSTRAINT *)malloc( sizeof(CONSTRAINT) );
        if(row==NULL){            
            std::cout << "There is not enough memory for ROW" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }
        x = (VARIABLE **)malloc( ncols*sizeof(VARIABLE *) );
        if( x==NULL ){            
            std::cout << " ERROR: not enough memory for the SORT" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }
        c = (double *)malloc( ncols*sizeof(double) );
        if( c==NULL ){            
            std::cout << " ERROR: not enough memory for the SORT" << std::endl;
            CSPexit(EXIT_MEMO); //exit(1);
        }
        if( NFREE == -1 ){
            inbetter = (char *)calloc( ncols , sizeof(char) );
            if(inbetter==NULL){                
                std::cout << "There is not enough memory for ROW" << std::endl;
                CSPexit(EXIT_MEMO); //exit(1);
            }
            for(k=0;k<nbetter;k++) inbetter[ better[k]->index ] = 1;
            card = 0;
            for(k=0;k<ncols;k++)
                if( inbetter[k]==0 ){
                    x[card] = columns+k;
                    c[card] = -1.0;
                    card++;
                }
            free(inbetter);
            inbetter = NULL; /*PWOF*/
            rhs = 0;
        } else {            
            for(k=0;k<nbetter;k++){
                x[k] = better[k];
                c[k] = 1.0;
            }
            rhs  = nbetter - NFREE;
            card = nbetter;
        }
   
        row->rhs    = rhs;
        row->card   = card;
        row->stack  = x;
        row->coef   = c;
        row->sense  = 'G';
        row->type   = BRANCHCUT;
        row->stat   = FIX_LB;
        row->lp     = -1;
        row->con    = NULL;
        rows[nrows++] = row;
}





