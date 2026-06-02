/****************************************************/
/* The Cell Suppression Problem for Linked Tables   */
/* in Statistical Disclosure Control (SDC)          */
/*                                                  */
/* Funded by IST-2000-25069-CASC                    */
/* to be used only inside tau-ARGUS                 */
/*                                                  */
/* CSPBACK.c                                       */
/* Version 1.0.0                                    */
/*                                                  */
/* M.D. Montes de Oca & J.J. SALAZAR                */
/*                                                  */
/* Last modified September 10, 2001                 */
/****************************************************/


#include <stdlib.h>
#include "cspback.h"

static int (*external_stop)(void) = {NULL};
static int (*external_heur)(void) = {NULL};
static int (*external_exit)(int) = {NULL};
static int (*external_time)(void) = {NULL};

int CSPstoptime()
{
    if ( external_time==NULL )
        return 0;
    else
        return (*external_time)();
}

int CSPstopcondition()
{
    if( external_stop==NULL )
        return 0;
    else
        return (*external_stop)();
}

int CSPnewsolution()
{
    if( external_heur==NULL )
        return 0;
    else
        return (*external_heur)();
}

int CSPexit(int flag)

{
    if( external_exit==NULL )
        return 0;
    else
        return (*external_exit)(flag);
}

int CSPdefinestop(int (*stop)(void))

{
    external_stop = stop;
    return 0;
}

int CSPdefineheur(int (*heur)(void))

{
    external_heur = heur;
    return 0;
}

int CSPdefineexit(int (*fexit)(int))

{
    external_exit = fexit;
    return 0;
}

int CSPdefinestoptime(int (*fstoptime)(void))
{
    external_time = fstoptime;
    return 0;
}
