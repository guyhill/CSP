#include "mt1rc.h"

static void   CHMT1RC(int, double *, double *, double, double *, int);
//static void MT1RC(int,double *,double *,double,double,double *,int *,
  //     int,int,int *,double *,double *,double *,int *,double *,double *);
/***********************************************************/

void MT1RC(int N, double *P, double *W, double C, double EPS, double *Z, int *X, int JDIM, int JCK,
           int *XX, double *MIN, double *PSIGN, double *WSIGN, int *ZSIGN, double *CRC, double *CRP)

{
/***************************************************************
 * Questa subroutine risolve "0-1 KNAPSACK PROBLEM" con        *
 * parametri reali:                                            *
 *                                                             *
 *   maximize Z = P[1]*X[1]+...+P[N]*X[N]                      *
 *   s.t.                                                      *
 *                W[1]*X[1]+...+W[N]*X[N] <= C,                *
 *                X[J] = 0 or 1 per J = 1,..,N.                *
 *                                                             *
 * Il programma � una versione modificata della subroutine     *
 * MT1C.                                                       *
 *                                                             *
 * I dati di ingresso del problema devono soddisfare le        *
 * seguenti condizioni:                                        *
 *   1) 2 <= N <= JDIM-2;                                      *
 *   2) P[J],W[J],C numeri reali positivi;                     *
 *   3) MAX(W[J]) <= C;                                        *
 *   4) W[1]+...+W[N] > C;                                     *
 *   5) P[J]/W[J] >= P[J+1]/W[J+1] per J = 1,...,N-1.          *
 *                                                             *
 * MT1RC chiama 1 subroutine:CHMT1RC.                          *
 *                                                             *
 * La comunicazione con il programma � solamente possibile     *
 * attraverso la lista dei parametri.                          *
 * Il programma � stato scritto utilizzando il linguaggio di   *
 * programmazione C,non sono state usate costanti dipendenti   *
 * dalla macchina che si intende utilizzare.                   *
 * Il programma � stato collaudato su DIGITAL VAX_STATION 2000.*
 *                                                             *
 * MT1RC ha bisogno di 10 vettori:                             *
 *   P , W , X , XX , MIN , PSIGN , WSIGN , ZSIGN , CRC, CRP   *
 * la cui dimensione deve essere almeno N+2.                   *
 *                                                             *
 * I principali parametri di ingresso sono:                    *
 *   N    = numero di oggetti;                                 *
 *   P[J] = profitto dell'oggetto J (J = 1,...,N);             *
 *   W[J] = peso dell'oggetto J (J = 1,...N);                  *
 *   C    = capacit�;                                          *
 *   EPS  = tolleranza(due valori positivi Q e R sono          *
 *          considerati uguali if ABS(Q-R)/MAX(Q,R) <= EPS;    *
 *   JDIM = dimensione dei 10  arrays;                         *
 *   JCK  = 1 se � richiesto il controllo dei dati di ingresso,*
 *          0 altrimenti.                                      *
 *                                                             *
 * I principali parametri di uscita sono:                      *
 *   Z =    valore della soluzione ottima se Z > 0,            *
 *          errore nei dati di ingresso (quando JCK = 1)       *
 *          se Z < 0,la condizione -Z � stata viol�ta;         *
 *   X[J] = 1 se l'oggetto J � nella soluzione ottima,         *
 *          0 altrimenti.                                      *
 *                                                             *
 * I vettori XX,MIN,PSIGN,WSIGN,ZSIGN,CRC,CRP sono ausiliari.  *
 *                                                             *
 * I parametri N,X,JDIM,JCK,XX,ZSIGN sono interi,i parametri   *
 * P,W,C,Z,MIN,PSIGN,WSIGN,CRC,CRP,EPS sono reali.AL ritorno   *
 * al chiamante tutti i parametri di ingresso sono invariati.  *
 **************************************************************/

int LL,KK,NM2,JJ,LOLD,II,JJ1,JP1,NEL,J1,IN,J,NN,N1;
double LIM,LIM1,IP,MINK,IU,CH,CHS,PROFIT,R,DIFF,T,A,B,EPSP;
//void CHMT1RC();

*Z = 0.0;
if (JCK == 1) CHMT1RC(N,P,W,C,Z,JDIM);

if (*Z == -4) {
  *Z = 0;
  for (J=1; J<=N; J++) { 
       if (W[J] > 0) {
              X[J] = 1; *Z += W[J];
       }else{ X[J]=0; }
  } return;
}

if (*Z < 0.0) return;

// Initialize

IP=0.0;
CHS=C+EPS*C;
for (LL=1;LL<=N;LL++)
  {
  if (W[LL]<=CHS)
    {
    IP+=P[LL];
    CHS-=W[LL];
    }
  else
    break;
  }
LL--;
if (CHS > 0.0)
  {
  P[N+1]= -4.0*IP*C;
  W[N+1]=2.0*C;
  LIM=IP+CHS*P[LL+2]/W[LL+2];
  A=W[LL+1]-CHS;
  B=IP+P[LL+1];
  LIM1=B-A*P[LL]/W[LL];
  if (LIM1 > LIM)
    LIM=LIM1;
  EPSP=EPS*IP;
  LIM-=EPSP;
  MINK=2.0*C;
  MIN[N]=MINK;
  for (J=N;J>1;J--)
    {
    if (W[J] < MINK)
      MINK=W[J];
    MIN[J-1]=MINK;
    }
  for (J=1;J<=N;J++)
    XX[J]=0;
  CH=C+EPS*C;
  PROFIT=0.0;
  LOLD=N;
  JJ=1;
  NM2=N-2;
  goto L180;
  }
*Z=IP;
for (J=1;J<=LL;J++)
  X[J]=1;
NN=LL+1;
for (J=NN;J<=N;J++)
  X[J]=0;
return;

L80:

// Try to insert the II_th item into the current solution

JJ=II;
CH=CRC[II];
PROFIT=CRP[II];
while (W[JJ] > CH)
  {
  JJ1=JJ+1;
  if ((*Z+EPSP) >= (CH*P[JJ1]/W[JJ1] + PROFIT))
    goto L290;
  JJ=JJ1;
  }

// Build a new current solution

IP=PSIGN[JJ];
CHS=CH-WSIGN[JJ];
IN=ZSIGN[JJ];
for (LL=IN;LL<=N;LL++)
  {
  if (W[LL] > CHS)
    goto L170;
  IP+=P[LL];
  CHS-=W[LL];
  }
LL=N;

L120:

if (*Z >= (PROFIT+IP))
  goto L290;
*Z=PROFIT+IP;
NN=JJ-1;
for (J=1;J<=NN;J++)
  X[J]=XX[J];
for (J=JJ;J<=LL;J++)
  X[J]=1;
if (LL != N)
  {
  NN=LL+1;
  for (J=NN;J<=N;J++)
    X[J]=0;
  }
if (*Z <= LIM)
  goto L290;
else
  return;

L170:

IU=CHS*P[LL]/W[LL];
LL--;
if (IU <= EPSP)
  goto L120;
if ((*Z+EPSP) >= (PROFIT+IP+IU))
  goto L290;

L180:

// Save the current solution

II=JJ;
CRC[II]=CH;
CRP[II]=PROFIT;
CRC[II+1]=CRC[II]-W[II];
CRP[II+1]=CRP[II]+P[II];
NN=LL-1;
J1=LL+1;
if (NN >= II)
  {
  for (J=II;J<=NN;J++)
    {
    JP1=J+1;
    CRC[J+2]=CRC[JP1]-W[JP1];
    CRP[J+2]=CRP[JP1]+P[JP1];
    }
  }
PROFIT=CRP[LL+1];
for (J=J1;J<=LOLD;J++)
  {
  WSIGN[J]=0.0;
  PSIGN[J]=0.0;
  ZSIGN[J]=J;
  }
LOLD=LL;
NEL=LL-II+1;
for (JJ=1;JJ<=NEL;JJ++)
  {
  J=J1-JJ;
  WSIGN[J]=WSIGN[J+1]+W[J];
  PSIGN[J]=PSIGN[J+1]+P[J];
  ZSIGN[J]=J1;
  XX[J]=1;
  }
if (LL > NM2)
  {
  II=N;
  goto L260;
  }
CRC[LL+2]=CRC[LL+1];
CRP[LL+2]=CRP[LL+1];
if (LL >= NM2)
  {
  if (CRC[N] >= W[N])
    {
    PROFIT+=P[N];
    XX[N]=1;
    }
  II=N-1;
  goto L260;
  }
II=LL+2;
if (CRC[LL+2] >= MIN[II-1])
  goto L80;

L260:

// Save the current optimal solution

if (*Z < PROFIT)
  {
  *Z=PROFIT;
  for (J=1;J<=N;J++)
    X[J]=XX[J];
  if (*Z >= LIM)
    return;
  }
XX[N]=0;

L290:

// Backtrack

NN=II-1;
if (NN == 0)
  return;
for (J=1;J<=NN;J++)
  {
  KK=II-J;
  if (XX[KK] == 1)
    goto L310;
  }
return;

L310:

R=CRC[KK+1];
CRC[KK+1]=CRC[KK];
CRP[KK+1]=CRP[KK];
CH=CRC[KK];
PROFIT=CRP[KK];
XX[KK]=0;
if (R >= MIN[KK])
  {
  II=KK+1;
  goto L80;
  }
NN=KK+1;
II=KK;

L330:

// Try to substitute the NN_th item for the KK_th

if ((*Z+EPSP) >= (PROFIT+CH*P[NN]/W[NN]))
  goto L290;
DIFF=W[NN]-W[KK];
if (DIFF == 0.0)
  {

  L340:

  NN++;
  goto L330;
  }
if (DIFF > 0.0)
  {
  if (DIFF > R)
    goto L340;
  if ((*Z+EPSP) >= (PROFIT+P[NN]))
    goto L340;
  *Z=PROFIT+P[NN];
  for (J=1;J<=KK;J++)
    X[J]=XX[J];
  JJ=KK+1;
  for (J=JJ;J<=N;J++)
    X[J]=0;
  X[NN]=1;
  if (*Z >= LIM)
    return;
  else
    {
    R=CH-W[NN];
    KK=NN;
    NN++;
    goto L330;
    }
  }
if (DIFF < 0.0)
  {
  T=CH-W[NN];
  if (T < MIN[NN])
    goto L340;
  if ((*Z+EPSP) >= (PROFIT+P[NN]+T*P[NN+1]/W[NN+1]))
    goto L290;
  CRC[NN]=CH;
  CRP[NN]=PROFIT;
  CRC[NN+1]=CH-W[NN];
  CRP[NN+1]=PROFIT+P[NN];
  XX[NN]=1;
  II=NN+1;
  WSIGN[NN]=W[NN];
  PSIGN[NN]=P[NN];
  ZSIGN[NN]=II;
  N1=NN+1;
  for (J=N1;J<=LOLD;J++)
    {
    WSIGN[J]=0.0;
    PSIGN[J]=0.0;
    ZSIGN[J]=J;
    }
  LOLD=NN;
  goto L80;
  }
return;
}

static void CHMT1RC(int N,double *P,double *W,double C,double *Z,int JDIM)

{
// check the input data

int J;
double JSW,RR,R;
if (N < 2)
  {
  *Z= -1;
  return;
  }
if (N > (JDIM-2))
  {
  *Z= -1;
  return;
  }
if (C <= 0)
  {

  L20:

  *Z= -2;
  return;
  }
JSW=0.0;
RR=P[1]/W[1];
for (J=1;J<=N;J++)
  {
  R=RR;
  if (P[J] <= 0.0)
    goto L20;
  if (W[J] <= 0.0)
    goto L20;
  JSW+=W[J];
  if (W[J] > C)
    {
    *Z= -3;
    return;
    }
  RR=P[J]/W[J];
  if (RR > R)
    {
    *Z= -5;
    return;
    }
  }
if (JSW > C)
  return;
   *Z= -4; 
return;
}
