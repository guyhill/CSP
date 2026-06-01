/***********************************************************
 *                                                         *
 *   cputime.c                          J.J.Salazar        *
 *                                                         *
 *   Different functions for getting CPU time since        *
 *   the program start for different computers with        *
 *   different C compilers:                                *
 *                                                         *
 *      BORLAND    for PC with DOS and Borland C           *
 *      WATCOM     for PC with DOS and Watcom C            *
 *      HPUX       for HP with HP-UX                       *
 *      SUN        for SUN with UNIX                       *
 *      ULTRIX     for DEC Station with ULTRIX             *
 *      VAX        for DEC Station with VMS                *
 *                                                         *
 *   float seconds(void)     return the CPU time           *
 *                                                         *
 ***********************************************************/

#include <time.h>
#include "my_time.h"


float seconds(void) {
    return( (float)clock()/CLOCKS_PER_SEC );
}
